#!/usr/bin/env python3
"""Integrate research designs into ResearchRuns and link the local FrameWork."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = SCRIPT_DIR / "workspace_links.yaml"
RUN_ENVIRONMENT_DIR = SCRIPT_DIR.parent / "RunEnvironment"
sys.path.insert(0, str(RUN_ENVIRONMENT_DIR))

from yaml_support import YamlFormatError, load_yaml  # noqa: E402


class WorkspaceLinkError(RuntimeError):
    """Raised when workspace integration cannot be completed safely."""


@dataclass(frozen=True)
class DesignIntegration:
    legacy_path: Path
    integrated_path: Path
    keep_legacy_link: bool


@dataclass(frozen=True)
class LinkSpec:
    key: str
    name: str
    target: Path


@dataclass(frozen=True)
class LinkConfig:
    source: Path
    execution_root: Path
    link_type: str
    designs: DesignIntegration
    links: tuple[LinkSpec, ...]


@dataclass(frozen=True)
class LinkState:
    spec: LinkSpec
    link: Path
    status: str
    actual_target: Path | None = None


_ENV_PATTERN = re.compile(r"\$\{([^}]+)\}|\$([A-Za-z_][A-Za-z0-9_]*)|%([^%]+)%")


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise WorkspaceLinkError(f"{label} must be a mapping")
    return value


def _expand_environment(text: str) -> str:
    def replace(match: re.Match[str]) -> str:
        name = next(group for group in match.groups() if group is not None)
        value = os.environ.get(name)
        if value is None:
            raise WorkspaceLinkError(f"environment variable is not set: {name}")
        return value

    return _ENV_PATTERN.sub(replace, text)


def _path(value: Any, base: Path, label: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise WorkspaceLinkError(f"{label} must be a non-empty path")
    expanded = Path(_expand_environment(value.strip())).expanduser()
    if not expanded.is_absolute():
        expanded = base / expanded
    return Path(os.path.abspath(expanded))


def _safe_link_name(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise WorkspaceLinkError(f"{label} must be a non-empty name")
    name = value.strip()
    if name in {".", ".."} or Path(name).name != name or "/" in name or "\\" in name:
        raise WorkspaceLinkError(f"{label} must be a single directory name")
    return name


def _normalized(path: Path) -> str:
    return os.path.normcase(os.path.abspath(path))


def _is_within(path: Path, parent: Path) -> bool:
    try:
        return os.path.commonpath([_normalized(path), _normalized(parent)]) == _normalized(
            parent
        )
    except ValueError:
        return False


def _load_config(path: Path) -> LinkConfig:
    source = path.resolve()
    data = _mapping(load_yaml(source), "workspace link design")
    if data.get("schema_version") != 1:
        raise WorkspaceLinkError("workspace link design.schema_version must be 1")

    execution_root = _path(
        data.get("execution_root"), source.parent, "execution_root"
    )
    link_type = str(data.get("link_type", "auto")).strip().lower()
    if link_type not in {"auto", "junction", "symlink"}:
        raise WorkspaceLinkError("link_type must be auto, junction, or symlink")
    if link_type == "junction" and os.name != "nt":
        raise WorkspaceLinkError("junction links are available only on Windows")

    design_data = _mapping(data.get("design_environment"), "design_environment")
    designs = DesignIntegration(
        legacy_path=_path(
            design_data.get("legacy_path"),
            source.parent,
            "design_environment.legacy_path",
        ),
        integrated_path=_path(
            design_data.get("integrated_path"),
            source.parent,
            "design_environment.integrated_path",
        ),
        keep_legacy_link=bool(design_data.get("keep_legacy_link", True)),
    )
    if not _is_within(designs.integrated_path, execution_root):
        raise WorkspaceLinkError(
            "design_environment.integrated_path must be inside execution_root"
        )
    if _normalized(designs.integrated_path) == _normalized(execution_root):
        raise WorkspaceLinkError(
            "design_environment.integrated_path must not equal execution_root"
        )
    if _normalized(designs.legacy_path) == _normalized(designs.integrated_path):
        raise WorkspaceLinkError(
            "design_environment legacy_path and integrated_path must differ"
        )

    raw_links = _mapping(data.get("links"), "links")
    if not raw_links:
        raise WorkspaceLinkError("links must contain at least one entry")

    specs: list[LinkSpec] = []
    names: set[str] = set()
    for key, raw in raw_links.items():
        entry = _mapping(raw, f"links.{key}")
        name = _safe_link_name(entry.get("name"), f"links.{key}.name")
        folded = os.path.normcase(name)
        if folded in names:
            raise WorkspaceLinkError(f"duplicate link name: {name}")
        names.add(folded)
        target = _path(entry.get("target"), source.parent, f"links.{key}.target")
        if _is_within(target, execution_root) or _is_within(execution_root, target):
            raise WorkspaceLinkError(
                f"links.{key}.target must be outside execution_root "
                "to prevent recursive links"
            )
        specs.append(LinkSpec(key=str(key), name=name, target=target))

    return LinkConfig(
        source=source,
        execution_root=execution_root,
        link_type=link_type,
        designs=designs,
        links=tuple(specs),
    )


def _lexists(path: Path) -> bool:
    return os.path.lexists(path)


def _is_junction(path: Path) -> bool:
    checker = getattr(path, "is_junction", None)
    return bool(checker and checker())


def _is_managed_link(path: Path) -> bool:
    return path.is_symlink() or _is_junction(path)


def _resolved(path: Path) -> Path | None:
    try:
        return path.resolve(strict=True)
    except (FileNotFoundError, OSError, RuntimeError):
        return None


def _same_path(left: Path, right: Path) -> bool:
    try:
        left_text = _normalized(left.resolve(strict=True))
        right_text = _normalized(right.resolve(strict=True))
    except (FileNotFoundError, OSError, RuntimeError):
        return False
    return left_text == right_text


def _effective_link_type(config: LinkConfig) -> str:
    if config.link_type == "auto":
        return "junction" if os.name == "nt" else "symlink"
    return config.link_type


def _create_link(link: Path, target: Path, link_type: str) -> None:
    link.parent.mkdir(parents=True, exist_ok=True)
    if link_type == "junction":
        result = subprocess.run(
            ["cmd.exe", "/d", "/c", "mklink", "/J", str(link), str(target)],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            detail = (result.stderr or result.stdout).strip()
            raise WorkspaceLinkError(f"failed to create junction {link}: {detail}")
        return
    try:
        link.symlink_to(target, target_is_directory=True)
    except OSError as exc:
        raise WorkspaceLinkError(f"failed to create symlink {link}: {exc}") from exc


def _remove_link(path: Path) -> None:
    if path.is_symlink():
        path.unlink()
    elif _is_junction(path):
        os.rmdir(path)
    else:
        raise WorkspaceLinkError(f"refusing to remove a real file or directory: {path}")


def _inspect_links(config: LinkConfig) -> list[LinkState]:
    states: list[LinkState] = []
    for spec in config.links:
        link = config.execution_root / spec.name
        if not spec.target.is_dir():
            states.append(LinkState(spec, link, "target_missing"))
        elif not _lexists(link):
            states.append(LinkState(spec, link, "link_missing"))
        elif not _is_managed_link(link):
            states.append(LinkState(spec, link, "occupied"))
        else:
            actual = _resolved(link)
            status = "ok" if _same_path(link, spec.target) else "wrong_target"
            states.append(LinkState(spec, link, status, actual))
    return states


def _design_status(config: LinkConfig) -> str:
    legacy = config.designs.legacy_path
    integrated = config.designs.integrated_path
    if not integrated.is_dir() or _is_managed_link(integrated):
        if _lexists(integrated):
            return "invalid_integrated"
        if _lexists(legacy) and not _is_managed_link(legacy) and legacy.is_dir():
            return "legacy_only"
        return "designs_missing"
    if not config.designs.keep_legacy_link:
        return "ready"
    if not _lexists(legacy):
        return "legacy_link_missing"
    if not _is_managed_link(legacy):
        return "merge_required" if legacy.is_dir() else "legacy_occupied"
    return "ready" if _same_path(legacy, integrated) else "wrong_legacy_link"


def _merge_conflicts(source: Path, destination: Path) -> list[Path]:
    if not source.is_dir() or _is_managed_link(source):
        return []
    return [
        child
        for child in source.iterdir()
        if _lexists(destination / child.name)
    ]


def _integrate_designs(
    config: LinkConfig, replace_links: bool, dry_run: bool
) -> None:
    legacy = config.designs.legacy_path
    integrated = config.designs.integrated_path
    link_type = _effective_link_type(config)
    legacy_will_be_available = not _lexists(legacy)

    if _lexists(integrated) and (
        _is_managed_link(integrated) or not integrated.is_dir()
    ):
        raise WorkspaceLinkError(
            f"integrated design path must be a real directory: {integrated}"
        )

    if _lexists(legacy) and _is_managed_link(legacy):
        if _same_path(legacy, integrated):
            if not integrated.is_dir():
                raise WorkspaceLinkError(
                    f"legacy design link target is missing: {integrated}"
                )
            print(f"[OK] design compatibility link: {legacy} -> {integrated}")
            return
        if not replace_links:
            raise WorkspaceLinkError(
                f"legacy design link points elsewhere: {legacy} -> {_resolved(legacy)}; "
                "use --replace-links to replace links only"
            )
        if dry_run:
            print(f"[WOULD REPLACE DESIGN LINK] {legacy} -> {integrated}")
            legacy_will_be_available = True
        else:
            _remove_link(legacy)
            print(f"[UNLINKED] old design link: {legacy}")
            legacy_will_be_available = True

    if _lexists(legacy) and not _is_managed_link(legacy):
        if not legacy.is_dir():
            raise WorkspaceLinkError(
                f"legacy design path is not a directory: {legacy}"
            )
        if integrated.is_dir():
            conflicts = _merge_conflicts(legacy, integrated)
            if conflicts:
                names = ", ".join(path.name for path in conflicts[:5])
                raise WorkspaceLinkError(
                    f"cannot merge design directories because names already exist "
                    f"in {integrated}: {names}"
                )
            if dry_run:
                print(f"[WOULD MERGE DESIGNS] {legacy} -> {integrated}")
                legacy_will_be_available = True
            else:
                for child in legacy.iterdir():
                    shutil.move(str(child), str(integrated / child.name))
                legacy.rmdir()
                print(f"[MERGED] designs: {legacy} -> {integrated}")
                legacy_will_be_available = True
        elif dry_run:
            print(f"[WOULD MOVE DESIGNS] {legacy} -> {integrated}")
            legacy_will_be_available = True
        else:
            integrated.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(legacy), str(integrated))
            print(f"[MOVED] designs: {legacy} -> {integrated}")
            legacy_will_be_available = True
    elif not integrated.is_dir():
        if dry_run:
            print(f"[WOULD CREATE DESIGNS] {integrated}")
        else:
            integrated.mkdir(parents=True)
            print(f"[CREATED] designs: {integrated}")

    if config.designs.keep_legacy_link and legacy_will_be_available:
        if dry_run:
            print(f"[WOULD LINK LEGACY PATH] {legacy} -> {integrated}")
        else:
            _create_link(legacy, integrated, link_type)
            print(f"[LINKED] legacy design path: {legacy} -> {integrated}")


def _apply(config: LinkConfig, replace_links: bool, dry_run: bool) -> None:
    link_type = _effective_link_type(config)
    if dry_run:
        print("[DRY-RUN] No directories, files, or links will be changed.")
    elif not config.execution_root.exists():
        config.execution_root.mkdir(parents=True)
        print(f"[CREATED] execution root: {config.execution_root}")
    elif not config.execution_root.is_dir():
        raise WorkspaceLinkError(
            f"execution_root is not a directory: {config.execution_root}"
        )

    _integrate_designs(config, replace_links, dry_run)

    for spec in config.links:
        link = config.execution_root / spec.name
        if not spec.target.is_dir():
            raise WorkspaceLinkError(f"link target does not exist: {spec.target}")
        if not _lexists(link):
            if dry_run:
                print(f"[WOULD LINK] {link} -> {spec.target}")
            else:
                _create_link(link, spec.target, link_type)
                print(f"[LINKED] {link} -> {spec.target}")
            continue
        if not _is_managed_link(link):
            raise WorkspaceLinkError(
                f"link path is occupied by a real file or directory: {link}"
            )
        if _same_path(link, spec.target):
            print(f"[OK] {link} -> {spec.target}")
            continue
        if not replace_links:
            raise WorkspaceLinkError(
                f"link points elsewhere: {link} -> {_resolved(link)}; "
                "use --replace-links to replace links only"
            )
        if dry_run:
            print(f"[WOULD REPLACE LINK] {link} -> {spec.target}")
        else:
            _remove_link(link)
            _create_link(link, spec.target, link_type)
            print(f"[RELINKED] {link} -> {spec.target}")


def _remove(config: LinkConfig, dry_run: bool) -> None:
    if dry_run:
        print("[DRY-RUN] Design files, source files, and run results will not be removed.")

    candidates = [config.execution_root / spec.name for spec in config.links]
    if config.designs.keep_legacy_link:
        candidates.append(config.designs.legacy_path)

    for link in candidates:
        if not _lexists(link):
            print(f"[MISSING] {link}")
        elif not _is_managed_link(link):
            raise WorkspaceLinkError(
                f"refusing to remove a real file or directory: {link}"
            )
        elif dry_run:
            print(f"[WOULD UNLINK] {link}")
        else:
            _remove_link(link)
            print(f"[UNLINKED] {link}")

    manifest = config.execution_root / ".workspace_links.json"
    if manifest.is_file() and not dry_run:
        manifest.unlink()


def _print_status(config: LinkConfig) -> bool:
    design_status = _design_status(config)
    print(f"execution_root:    {config.execution_root}")
    print(f"integrated_design: {config.designs.integrated_path}")
    print(f"legacy_design:     {config.designs.legacy_path}")
    print(f"link_type:         {_effective_link_type(config)}")
    print(f"[{design_status.upper():19}] design environment")

    states = _inspect_links(config)
    for state in states:
        detail = f" -> {state.spec.target}"
        if state.actual_target is not None and state.status == "wrong_target":
            detail += f" (actual: {state.actual_target})"
        print(f"[{state.status.upper():19}] {state.link}{detail}")
    return design_status == "ready" and all(state.status == "ok" for state in states)


def _write_manifest(config: LinkConfig) -> None:
    manifest = {
        "schema_version": 1,
        "config": str(config.source),
        "execution_root": str(config.execution_root),
        "link_type": _effective_link_type(config),
        "design_environment": {
            "legacy_path": str(config.designs.legacy_path),
            "integrated_path": str(config.designs.integrated_path),
            "keep_legacy_link": config.designs.keep_legacy_link,
        },
        "links": {
            spec.key: {
                "name": spec.name,
                "target": str(spec.target),
            }
            for spec in config.links
        },
    }
    path = config.execution_root / ".workspace_links.json"
    path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Move C:\\ResearchDesigns into ResearchRuns\\Designs and link the "
            "local FrameWork from the execution root."
        )
    )
    parser.add_argument(
        "config",
        nargs="?",
        default=str(DEFAULT_CONFIG),
        help="Workspace link YAML; default: workspace_links.yaml",
    )
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--apply", action="store_true", help="Integrate and link")
    action.add_argument("--remove", action="store_true", help="Remove managed links")
    action.add_argument("--status", action="store_true", help="Show current status")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show changes without modifying directories or links",
    )
    parser.add_argument(
        "--replace-links",
        action="store_true",
        help="Replace only links that point to a different target",
    )
    args = parser.parse_args()

    try:
        config = _load_config(Path(args.config))
        if args.remove:
            _remove(config, args.dry_run)
        elif args.apply:
            _apply(config, args.replace_links, args.dry_run)
            if not args.dry_run:
                _write_manifest(config)
                if not _print_status(config):
                    raise WorkspaceLinkError("workspace verification failed")
                print("[OK] ResearchRuns integration and links are ready.")
        elif not _print_status(config):
            return 1
    except (OSError, ValueError, YamlFormatError, WorkspaceLinkError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
