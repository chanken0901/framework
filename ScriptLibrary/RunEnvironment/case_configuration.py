#!/usr/bin/env python3
"""Resolve a case design and its explicitly referenced extension files."""

from __future__ import annotations

import copy
import json
import os
import uuid
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Any

from yaml_support import YamlFormatError, load_yaml


class CaseConfigurationError(ValueError):
    """Raised when a modular case configuration is invalid or unsafe."""


# Register extension ownership explicitly. Adding an extension must not silently
# change how an existing YAML key is interpreted.
EXTENSION_TARGETS: dict[str, tuple[str, ...]] = {
    "multicomponent": ("physics", "multicomponent"),
    "thermodynamics": ("thermodynamics",),
    "transport": ("transport",),
    "chemistry": ("chemistry",),
}


@dataclass(frozen=True)
class ResolvedCaseConfiguration:
    """A standalone case document and every source file used to build it."""

    document: dict[str, Any]
    source_paths: tuple[Path, ...]
    extension_paths: dict[str, Path]


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise CaseConfigurationError(f"{label} must be a YAML mapping")
    return value


def _extension_path(case_path: Path, name: str, reference: Any) -> Path:
    if not isinstance(reference, str) or not reference.strip():
        raise CaseConfigurationError(
            f"extensions.{name} must be a non-empty relative YAML path"
        )
    relative = Path(reference)
    if relative.is_absolute():
        raise CaseConfigurationError(
            f"extensions.{name} must be relative to case.yaml: {reference}"
        )
    case_directory = case_path.parent.resolve()
    lexical = (case_directory / relative).absolute()
    candidate = lexical.resolve()
    try:
        candidate.relative_to(case_directory)
    except ValueError as exc:
        raise CaseConfigurationError(
            f"extensions.{name} escapes the case directory: {reference}"
        ) from exc
    if candidate != lexical:
        raise CaseConfigurationError(
            f"extensions.{name} must not use '..', symbolic links, or junctions: "
            f"{reference}"
        )
    if candidate == case_path.resolve():
        raise CaseConfigurationError(f"extensions.{name} cannot reference case.yaml")
    if candidate.suffix.lower() not in {".yaml", ".yml"}:
        raise CaseConfigurationError(
            f"extensions.{name} must reference a .yaml or .yml file: {reference}"
        )
    if not candidate.is_file():
        raise CaseConfigurationError(
            f"extensions.{name} file was not found: {candidate}"
        )
    return candidate


def _install_extension(
    document: dict[str, Any], name: str, target: tuple[str, ...], config: dict[str, Any]
) -> None:
    parent = document
    for key in target[:-1]:
        existing = parent.get(key)
        if existing is None:
            existing = {}
            parent[key] = existing
        if not isinstance(existing, dict):
            dotted = ".".join(target[:-1])
            raise CaseConfigurationError(
                f"cannot install extensions.{name}; {dotted} must be a YAML mapping"
            )
        parent = existing
    leaf = target[-1]
    if leaf in parent:
        dotted = ".".join(target)
        raise CaseConfigurationError(
            f"duplicate configuration for {dotted}; define it in case.yaml or "
            f"extensions.{name}, not both"
        )
    parent[leaf] = copy.deepcopy(config)


def resolve_case_configuration(case_path: str | Path) -> ResolvedCaseConfiguration:
    """Load a legacy case or compose a schema-v2 case with explicit extensions.

    Legacy single-file cases need no schema version and are returned unchanged.
    A modular case declares ``schema_version: 2`` and an ``extensions`` mapping.
    Each sidecar declares schema version 1, a matching extension name, and a
    ``config`` mapping. Sidecars must remain inside the case directory.
    """

    source = Path(case_path).resolve()
    if not source.is_file():
        raise CaseConfigurationError(f"case YAML was not found: {source}")
    try:
        loaded = load_yaml(source)
    except (OSError, YamlFormatError) as exc:
        raise CaseConfigurationError(f"failed to load case YAML {source}: {exc}") from exc
    document = copy.deepcopy(_mapping(loaded, f"case YAML {source}"))
    if "extensions" not in document:
        return ResolvedCaseConfiguration(document, (source,), {})
    raw_extensions = document.pop("extensions", None)
    if document.get("schema_version") != 2:
        raise CaseConfigurationError(
            "case.yaml must declare schema_version: 2 when extensions are used"
        )
    extensions = _mapping(raw_extensions, "extensions")
    if not extensions:
        raise CaseConfigurationError(
            "extensions must not be empty; omit it for a single-file case"
        )

    extension_paths: dict[str, Path] = {}
    used_paths: set[Path] = set()
    for raw_name, reference in extensions.items():
        name = str(raw_name).strip().lower()
        if name != raw_name or name not in EXTENSION_TARGETS:
            choices = ", ".join(sorted(EXTENSION_TARGETS))
            raise CaseConfigurationError(
                f"unsupported extension {raw_name!r}; supported extensions: {choices}"
            )
        path = _extension_path(source, name, reference)
        if path in used_paths:
            raise CaseConfigurationError(
                f"extension file is referenced more than once: {path}"
            )
        used_paths.add(path)
        try:
            sidecar = _mapping(load_yaml(path), f"extension YAML {path}")
        except (OSError, YamlFormatError) as exc:
            raise CaseConfigurationError(
                f"failed to load extensions.{name} from {path}: {exc}"
            ) from exc
        unexpected = set(sidecar) - {"schema_version", "extension", "config"}
        if unexpected:
            keys = ", ".join(sorted(str(key) for key in unexpected))
            raise CaseConfigurationError(
                f"extensions.{name} contains unsupported top-level keys: {keys}"
            )
        if sidecar.get("schema_version") != 1:
            raise CaseConfigurationError(
                f"extensions.{name} must declare schema_version: 1"
            )
        if sidecar.get("extension") != name:
            raise CaseConfigurationError(
                f"extensions.{name} file must declare extension: {name}"
            )
        config = _mapping(sidecar.get("config"), f"extensions.{name}.config")
        _install_extension(document, name, EXTENSION_TARGETS[name], config)
        extension_paths[name] = path

    return ResolvedCaseConfiguration(
        document,
        (source, *extension_paths.values()),
        extension_paths,
    )


def _json_yaml_scalar(value: Any) -> str:
    """Preserve YAML timestamps as ISO text in the JSON-compatible snapshot."""
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def write_resolved_case(
    configuration: ResolvedCaseConfiguration,
    output_path: str | Path | None = None,
) -> Path:
    """Atomically write the standalone merged document as JSON-compatible YAML."""

    destination = (
        Path(output_path).resolve()
        if output_path is not None
        else configuration.source_paths[0].with_name("resolved_case.yaml")
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        f".{destination.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp"
    )
    try:
        temporary.write_text(
            json.dumps(
                configuration.document, indent=2, ensure_ascii=False,
                default=_json_yaml_scalar,
            ) + "\n",
            encoding="utf-8",
        )
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination
