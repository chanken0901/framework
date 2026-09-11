#!/usr/bin/env python3
"""Safely manage the ScriptLibrary and SolverLibrary Git repositories.

The command is configuration driven and intentionally excludes automatic merges,
rebases, resets, and force pushes. Mutating operations require ``--apply``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
YAML_SUPPORT_DIR = SCRIPT_DIR.parents[1] / "SetupCase"
if str(YAML_SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(YAML_SUPPORT_DIR))

try:
    from yaml_support import YamlFormatError, load_yaml
except ImportError as exc:  # pragma: no cover - deployment integrity check
    raise SystemExit(
        f"yaml_support.py was not found under {YAML_SUPPORT_DIR}. "
        "Deploy the complete ScriptLibrary integration first."
    ) from exc


class GitAutomationError(RuntimeError):
    """A safe workflow precondition was not satisfied."""


@dataclass(frozen=True)
class RepositoryConfig:
    name: str
    path: Path
    branch: str
    remote: str
    remote_url: str | None
    enabled: bool


@dataclass(frozen=True)
class Settings:
    git_command: str
    max_changed_file_mb: float
    user_name: str | None
    user_email: str | None


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise GitAutomationError(f"{label} must be a YAML mapping")
    return value


def _expand(text: Any, framework_root: str) -> str:
    value = "" if text is None else str(text)
    value = os.path.expandvars(os.path.expanduser(value))
    return value.replace("{framework_root}", framework_root)


def load_configuration(
    path: Path, framework_override: str | None
) -> tuple[Settings, list[RepositoryConfig]]:
    try:
        raw = _mapping(load_yaml(path), "configuration root")
    except (OSError, YamlFormatError, ValueError) as exc:
        raise GitAutomationError(f"failed to read configuration {path}: {exc}") from exc

    configured_root = raw.get("framework_root", "")
    framework_text = framework_override or str(configured_root or "")
    if not framework_text:
        raise GitAutomationError("framework_root is empty; set it in YAML or use --framework-root")
    framework_text = os.path.expandvars(os.path.expanduser(framework_text))

    settings_raw = _mapping(raw.get("settings", {}), "settings")
    settings = Settings(
        git_command=str(settings_raw.get("git_command", "git")),
        max_changed_file_mb=float(settings_raw.get("max_changed_file_mb", 50.0)),
        user_name=(str(settings_raw["user_name"]) if settings_raw.get("user_name") else None),
        user_email=(str(settings_raw["user_email"]) if settings_raw.get("user_email") else None),
    )
    if settings.max_changed_file_mb <= 0:
        raise GitAutomationError("settings.max_changed_file_mb must be positive")

    default_branch = str(settings_raw.get("default_branch", "main"))
    default_remote = str(settings_raw.get("default_remote", "origin"))
    repositories_raw = _mapping(raw.get("repositories", {}), "repositories")
    repositories: list[RepositoryConfig] = []
    for name, value in repositories_raw.items():
        repo = _mapping(value, f"repositories.{name}")
        path_text = _expand(repo.get("path", ""), framework_text)
        if not path_text:
            raise GitAutomationError(f"repositories.{name}.path is empty")
        remote_url = _expand(repo.get("remote_url", ""), framework_text).strip() or None
        repositories.append(
            RepositoryConfig(
                name=str(name),
                path=Path(path_text),
                branch=str(repo.get("branch", default_branch)),
                remote=str(repo.get("remote", default_remote)),
                remote_url=remote_url,
                enabled=bool(repo.get("enabled", True)),
            )
        )
    if not repositories:
        raise GitAutomationError("no repositories are defined")
    return settings, repositories


class GitRunner:
    def __init__(self, command: str) -> None:
        self.executable = command

    def executable_exists(self) -> bool:
        candidate = Path(self.executable)
        return (
            candidate.is_file()
            if candidate.parent != Path(".")
            else shutil.which(self.executable) is not None
        )

    def _build_command(self, repo: RepositoryConfig, args: list[str]) -> list[str]:
        # Native Windows processes cannot reliably use a UNC path as cwd.
        # Git's -C option handles UNC repositories without changing process cwd.
        return [
            self.executable,
            "-c",
            "core.quotepath=false",
            "-C",
            str(repo.path),
            *args,
        ]

    def display(self, repo: RepositoryConfig, args: list[str]) -> str:
        command = self._build_command(repo, args)
        return subprocess.list2cmdline(command) if os.name == "nt" else shlex.join(command)

    def run(
        self,
        repo: RepositoryConfig,
        args: list[str],
        *,
        check: bool = True,
        show_output: bool = False,
        mutate: bool = False,
        apply: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        command = self._build_command(repo, args)
        prefix = "[DRY-RUN]" if mutate and not apply else "[CMD]"
        print(f"{prefix} {repo.name}: {self.display(repo, args)}")
        if mutate and not apply:
            return subprocess.CompletedProcess(command, 0, "", "")

        result = subprocess.run(
            command,
            cwd=os.environ.get("TEMP") or str(Path.home()),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if show_output and result.stdout.strip():
            print(result.stdout.rstrip())
        if show_output and result.stderr.strip():
            print(result.stderr.rstrip(), file=sys.stderr)
        if check and result.returncode != 0:
            detail = result.stderr.strip() or result.stdout.strip() or "unknown Git error"
            if "dubious ownership" in detail.lower():
                detail += (
                    "\nThis NAS path is not trusted by Git. Review it, then run:"
                    f"\n  git config --global --add safe.directory \"{repo.path}\""
                )
            raise GitAutomationError(
                f"{repo.name}: Git command failed ({self.display(repo, args)}):\n{detail}"
            )
        return result


def is_git_repository(runner: GitRunner, repo: RepositoryConfig) -> bool:
    if not repo.path.is_dir():
        return False
    result = runner.run(repo, ["rev-parse", "--is-inside-work-tree"], check=False)
    if result.returncode != 0 and (repo.path / ".git").exists():
        detail = result.stderr.strip() or result.stdout.strip() or "Git metadata is not readable"
        if "dubious ownership" in detail.lower():
            detail += (
                "\nReview the path and mark it safe with:"
                f"\n  git config --global --add safe.directory \"{repo.path}\""
            )
        raise GitAutomationError(f"{repo.name}: .git exists but Git cannot use it:\n{detail}")
    return result.returncode == 0 and result.stdout.strip() == "true"


def require_repository(runner: GitRunner, repo: RepositoryConfig) -> None:
    if not repo.path.is_dir():
        raise GitAutomationError(f"{repo.name}: directory does not exist: {repo.path}")
    if not is_git_repository(runner, repo):
        raise GitAutomationError(
            f"{repo.name}: not a Git repository: {repo.path}. "
            "If GitHub already has commits, run adopt --apply. "
            "Use connect --apply only for an empty/new GitHub repository."
        )


def git_value(runner: GitRunner, repo: RepositoryConfig, args: list[str]) -> str | None:
    result = runner.run(repo, args, check=False)
    return result.stdout.strip() if result.returncode == 0 and result.stdout.strip() else None


def current_branch(runner: GitRunner, repo: RepositoryConfig) -> str | None:
    return git_value(runner, repo, ["branch", "--show-current"])


def has_head(runner: GitRunner, repo: RepositoryConfig) -> bool:
    return runner.run(repo, ["rev-parse", "--verify", "HEAD"], check=False).returncode == 0


def remote_url(runner: GitRunner, repo: RepositoryConfig) -> str | None:
    return git_value(runner, repo, ["remote", "get-url", repo.remote])


def probe_remote(runner: GitRunner, repo: RepositoryConfig) -> dict[str, Any]:
    """Check authentication/reachability without changing local Git references."""
    result = runner.run(
        repo,
        ["ls-remote", "--exit-code", repo.remote, f"refs/heads/{repo.branch}"],
        check=False,
    )
    reachable = result.returncode in (0, 2)
    branch_exists = result.returncode == 0
    remote_head = None
    if branch_exists and result.stdout.strip():
        remote_head = result.stdout.split()[0]
    return {
        "reachable": reachable,
        "branch_exists": branch_exists,
        "head": remote_head,
        "error": None if reachable else (result.stderr.strip() or result.stdout.strip()),
    }


def is_clean(runner: GitRunner, repo: RepositoryConfig) -> bool:
    result = runner.run(repo, ["status", "--porcelain"])
    return not result.stdout.strip()


def remote_ref_exists(runner: GitRunner, repo: RepositoryConfig) -> bool:
    ref = f"refs/remotes/{repo.remote}/{repo.branch}"
    return runner.run(repo, ["show-ref", "--verify", "--quiet", ref], check=False).returncode == 0


def divergence(runner: GitRunner, repo: RepositoryConfig) -> tuple[int, int] | None:
    if not has_head(runner, repo) or not remote_ref_exists(runner, repo):
        return None
    remote_ref = f"refs/remotes/{repo.remote}/{repo.branch}"
    result = runner.run(repo, ["rev-list", "--left-right", "--count", f"HEAD...{remote_ref}"])
    values = result.stdout.strip().split()
    if len(values) != 2:
        raise GitAutomationError(f"{repo.name}: could not parse ahead/behind counts")
    return int(values[0]), int(values[1])


def ensure_expected_branch(runner: GitRunner, repo: RepositoryConfig) -> None:
    branch = current_branch(runner, repo)
    if branch and branch != repo.branch:
        raise GitAutomationError(
            f"{repo.name}: current branch is {branch!r}, expected {repo.branch!r}. "
            "Switch branches explicitly before continuing."
        )
    if branch is None and has_head(runner, repo):
        raise GitAutomationError(f"{repo.name}: detached HEAD is not supported")


def configured_identity(runner: GitRunner, repo: RepositoryConfig) -> tuple[str | None, str | None]:
    name = git_value(runner, repo, ["config", "--get", "user.name"])
    email = git_value(runner, repo, ["config", "--get", "user.email"])
    return name, email


def doctor_repository(runner: GitRunner, repo: RepositoryConfig) -> dict[str, Any]:
    result: dict[str, Any] = {
        "name": repo.name,
        "path": str(repo.path),
        "exists": repo.path.is_dir(),
        "git_repository": False,
        "expected_branch": repo.branch,
        "remote": repo.remote,
    }
    if not repo.path.is_dir():
        print(f"[WARN] {repo.name}: directory not found: {repo.path}")
        return result
    if not is_git_repository(runner, repo):
        print(f"[WARN] {repo.name}: not a Git repository")
        return result

    result["git_repository"] = True
    result["branch"] = current_branch(runner, repo)
    result["remote_url"] = remote_url(runner, repo)
    result["clean"] = is_clean(runner, repo)
    name, email = configured_identity(runner, repo)
    result["user_name"] = name
    result["user_email"] = email
    counts = divergence(runner, repo)
    result["ahead"] = counts[0] if counts else None
    result["behind"] = counts[1] if counts else None
    if result["remote_url"]:
        probe = probe_remote(runner, repo)
        result["remote_reachable"] = probe["reachable"]
        result["remote_branch_exists"] = probe["branch_exists"]
        result["remote_head"] = probe["head"]
        result["remote_error"] = probe["error"]
    else:
        result["remote_reachable"] = False
        result["remote_branch_exists"] = False
        result["remote_head"] = None
        result["remote_error"] = "remote is not configured"

    print(f"[{repo.name}]")
    print(f"  path:       {repo.path}")
    print(f"  branch:     {result['branch'] or '(no commits)'}")
    print(f"  remote:     {result['remote_url'] or '(not configured)'}")
    print(f"  reachable:  {result['remote_reachable']}")
    print(f"  remote main:{' present' if result['remote_branch_exists'] else ' absent/unknown'}")
    print(f"  clean:      {result['clean']}")
    print(f"  ahead:      {result['ahead'] if result['ahead'] is not None else 'unknown'}")
    print(f"  behind:     {result['behind'] if result['behind'] is not None else 'unknown'}")
    print(f"  identity:   {name or '(unset)'} <{email or 'unset'}>")
    if result["remote_error"]:
        print(f"  remote error: {result['remote_error']}")
    return result


def connect_repository(
    runner: GitRunner,
    repo: RepositoryConfig,
    settings: Settings,
    *,
    apply: bool,
    update_remote: bool,
) -> dict[str, Any]:
    if not repo.path.is_dir():
        raise GitAutomationError(f"{repo.name}: directory does not exist: {repo.path}")

    already_repository = is_git_repository(runner, repo)
    if not already_repository:
        runner.run(repo, ["init", "-b", repo.branch], mutate=True, apply=apply)
        if not apply:
            if repo.remote_url:
                runner.run(
                    repo,
                    ["remote", "add", repo.remote, repo.remote_url],
                    mutate=True,
                    apply=False,
                )
            return {"initialized": True, "remote": repo.remote_url, "dry_run": True}

    require_repository(runner, repo)
    ensure_expected_branch(runner, repo)

    existing_url = remote_url(runner, repo)
    if existing_url is None:
        if not repo.remote_url:
            raise GitAutomationError(
                f"{repo.name}: remote {repo.remote!r} is missing and remote_url is empty in YAML"
            )
        runner.run(
            repo,
            ["remote", "add", repo.remote, repo.remote_url],
            mutate=True,
            apply=apply,
        )
        existing_url = repo.remote_url
    elif repo.remote_url and existing_url != repo.remote_url:
        if not update_remote:
            raise GitAutomationError(
                f"{repo.name}: configured URL differs from existing remote. "
                "Review both values and rerun with --update-remote --apply.\n"
                f"  existing:   {existing_url}\n  configured: {repo.remote_url}"
            )
        runner.run(
            repo,
            ["remote", "set-url", repo.remote, repo.remote_url],
            mutate=True,
            apply=apply,
        )
        existing_url = repo.remote_url

    if settings.user_name:
        runner.run(
            repo,
            ["config", "user.name", settings.user_name],
            mutate=True,
            apply=apply,
        )
    if settings.user_email:
        runner.run(
            repo,
            ["config", "user.email", settings.user_email],
            mutate=True,
            apply=apply,
        )

    return {
        "initialized": not already_repository,
        "remote": existing_url,
        "branch": current_branch(runner, repo) or repo.branch,
        "dry_run": not apply,
    }


def worktree_hashes(root: Path) -> dict[str, str]:
    """Hash existing files without reading or following anything under .git."""
    hashes: dict[str, str] = {}
    for path in root.rglob("*"):
        try:
            relative = path.relative_to(root)
        except ValueError:
            continue
        if not relative.parts or relative.parts[0].casefold() == ".git":
            continue
        if not path.is_file() or path.is_symlink():
            continue
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        hashes[relative.as_posix()] = digest.hexdigest()
    return hashes


def _chunks(items: list[str], size: int = 100) -> list[list[str]]:
    return [items[index : index + size] for index in range(0, len(items), size)]


def adopt_existing_remote_history(
    runner: GitRunner,
    repo: RepositoryConfig,
    settings: Settings,
    *,
    apply: bool,
    update_remote: bool,
) -> dict[str, Any]:
    """Attach a non-repository folder to an existing remote without overwriting files."""
    if not repo.path.is_dir():
        raise GitAutomationError(f"{repo.name}: directory does not exist: {repo.path}")
    if is_git_repository(runner, repo) and has_head(runner, repo):
        raise GitAutomationError(
            f"{repo.name}: local commit history already exists. Use doctor/sync; adopt is only "
            "for a folder with no local commits."
        )

    before = worktree_hashes(repo.path)
    existing_remote = remote_url(runner, repo) if is_git_repository(runner, repo) else None
    target = repo.remote_url or existing_remote
    if not target:
        raise GitAutomationError(
            f"{repo.name}: remote_url is required in YAML because no existing {repo.remote!r} exists"
        )

    if not apply:
        probe = runner.run(
            repo,
            ["ls-remote", "--exit-code", target, f"refs/heads/{repo.branch}"],
            check=False,
        )
        if probe.returncode != 0:
            detail = probe.stderr.strip() or "configured branch was not found"
            raise GitAutomationError(
                f"{repo.name}: GitHub branch {repo.branch!r} is not reachable at {target}: {detail}"
            )
        if not is_git_repository(runner, repo):
            runner.run(repo, ["init", "-b", repo.branch], mutate=True, apply=False)
            runner.run(
                repo,
                ["remote", "add", repo.remote, target],
                mutate=True,
                apply=False,
            )
        runner.run(repo, ["fetch", "--prune", repo.remote], mutate=True, apply=False)
        remote_ref = f"refs/remotes/{repo.remote}/{repo.branch}"
        runner.run(repo, ["reset", "--mixed", remote_ref], mutate=True, apply=False)
        runner.run(
            repo,
            ["branch", "--set-upstream-to", f"{repo.remote}/{repo.branch}", repo.branch],
            mutate=True,
            apply=False,
        )
        return {
            "dry_run": True,
            "existing_local_files": len(before),
            "remote_url": target,
        }

    connect_repository(
        runner,
        repo,
        settings,
        apply=True,
        update_remote=update_remote,
    )
    fetch_repository(runner, repo, True)
    remote_ref = f"refs/remotes/{repo.remote}/{repo.branch}"
    if not remote_ref_exists(runner, repo):
        raise GitAutomationError(
            f"{repo.name}: remote branch {repo.remote}/{repo.branch} does not exist"
        )

    tree = runner.run(repo, ["ls-tree", "-r", "--name-only", "-z", remote_ref]).stdout
    remote_files = sorted(path for path in tree.split("\0") if path)
    root_resolved = repo.path.resolve()
    missing: list[str] = []
    for relative in remote_files:
        candidate = (repo.path / relative).resolve()
        try:
            candidate.relative_to(root_resolved)
        except ValueError as exc:
            raise GitAutomationError(
                f"{repo.name}: remote path escapes repository root: {relative}"
            ) from exc
        if not candidate.exists() and not candidate.is_symlink():
            missing.append(relative)

    # --mixed adopts GitHub's commit/index while deliberately leaving the worktree untouched.
    runner.run(repo, ["reset", "--mixed", remote_ref], mutate=True, apply=True)
    for group in _chunks(missing):
        runner.run(
            repo,
            ["restore", "--source", remote_ref, "--worktree", "--", *group],
            mutate=True,
            apply=True,
        )
    runner.run(
        repo,
        ["branch", "--set-upstream-to", f"{repo.remote}/{repo.branch}", repo.branch],
        mutate=True,
        apply=True,
    )

    after = worktree_hashes(repo.path)
    changed_existing = [path for path, digest in before.items() if after.get(path) != digest]
    if changed_existing:
        raise GitAutomationError(
            f"{repo.name}: existing files changed unexpectedly during adopt: "
            + ", ".join(changed_existing)
        )
    status = runner.run(
        repo,
        ["status", "--short", "--branch"],
        show_output=True,
    ).stdout.rstrip()
    return {
        "dry_run": False,
        "remote_url": target,
        "existing_local_files_verified": len(before),
        "remote_only_files_restored": len(missing),
        "status": status,
    }


def fetch_repository(runner: GitRunner, repo: RepositoryConfig, apply: bool) -> None:
    if not remote_url(runner, repo):
        raise GitAutomationError(f"{repo.name}: remote {repo.remote!r} is not configured")
    runner.run(
        repo,
        ["fetch", "--prune", repo.remote],
        mutate=True,
        apply=apply,
        show_output=True,
    )


def push_repository(runner: GitRunner, repo: RepositoryConfig, apply: bool) -> None:
    runner.run(
        repo,
        ["push", "-u", repo.remote, repo.branch],
        mutate=True,
        apply=apply,
        show_output=True,
    )


def sync_repository(
    runner: GitRunner, repo: RepositoryConfig, *, apply: bool, no_push: bool
) -> dict[str, Any]:
    require_repository(runner, repo)
    ensure_expected_branch(runner, repo)
    if not is_clean(runner, repo):
        raise GitAutomationError(
            f"{repo.name}: working tree is not clean. Commit or stash changes before sync."
        )
    if not has_head(runner, repo):
        raise GitAutomationError(f"{repo.name}: repository has no local commits")

    fetch_repository(runner, repo, apply)
    counts = divergence(runner, repo)
    if counts is None:
        if not no_push:
            push_repository(runner, repo, apply)
        return {"ahead": None, "behind": None, "initial_push": not no_push}

    ahead, behind = counts
    if ahead and behind:
        raise GitAutomationError(
            f"{repo.name}: local and remote histories diverged (ahead={ahead}, behind={behind}). "
            "Resolve this repository manually; no automatic merge or rebase was attempted."
        )
    if behind:
        runner.run(
            repo,
            ["pull", "--ff-only", repo.remote, repo.branch],
            mutate=True,
            apply=apply,
            show_output=True,
        )
    if ahead and not no_push:
        push_repository(runner, repo, apply)
    return {"ahead": ahead, "behind": behind, "pushed": bool(ahead and not no_push)}


def changed_paths(runner: GitRunner, repo: RepositoryConfig) -> list[str]:
    commands = [
        ["ls-files", "--others", "--exclude-standard", "-z"],
        ["diff", "--name-only", "-z"],
        ["diff", "--cached", "--name-only", "-z"],
    ]
    paths: set[str] = set()
    for command in commands:
        output = runner.run(repo, command).stdout
        paths.update(item for item in output.split("\0") if item)
    return sorted(paths)


def check_large_changed_files(
    runner: GitRunner,
    repo: RepositoryConfig,
    max_changed_file_mb: float,
) -> list[dict[str, Any]]:
    limit = int(max_changed_file_mb * 1024 * 1024)
    root = repo.path.resolve()
    blocked: list[dict[str, Any]] = []
    for relative in changed_paths(runner, repo):
        candidate = (repo.path / relative).resolve()
        try:
            candidate.relative_to(root)
        except ValueError:
            raise GitAutomationError(f"{repo.name}: changed path escapes repository: {relative}")
        if not candidate.is_file() or candidate.is_symlink():
            continue
        size = candidate.stat().st_size
        if size > limit:
            blocked.append({"path": relative, "size": size})
    if blocked:
        lines = "\n".join(
            f"  {item['path']}: {item['size'] / (1024 * 1024):.2f} MiB" for item in blocked
        )
        raise GitAutomationError(
            f"{repo.name}: changed files exceed the {max_changed_file_mb:g} MiB safety limit:\n"
            f"{lines}\nAdd appropriate .gitignore rules or use external data storage."
        )
    return blocked


def index_has_changes(runner: GitRunner, repo: RepositoryConfig) -> bool:
    result = runner.run(repo, ["diff", "--cached", "--quiet"], check=False)
    if result.returncode not in (0, 1):
        raise GitAutomationError(f"{repo.name}: failed to inspect staged changes")
    return result.returncode == 1


def snapshot_repository(
    runner: GitRunner,
    repo: RepositoryConfig,
    settings: Settings,
    *,
    message: str,
    apply: bool,
    no_push: bool,
) -> dict[str, Any]:
    require_repository(runner, repo)
    ensure_expected_branch(runner, repo)
    fetch_repository(runner, repo, apply)

    local_has_head = has_head(runner, repo)
    remote_has_branch = remote_ref_exists(runner, repo)
    if not local_has_head and remote_has_branch:
        raise GitAutomationError(
            f"{repo.name}: GitHub already has branch {repo.branch!r}, but this local folder "
            "has no commit history. Do not create a second initial commit. Clone the existing "
            "GitHub repository to a separate folder, compare the files, and migrate changes "
            "into that clone."
        )

    counts = divergence(runner, repo)
    dirty = not is_clean(runner, repo)
    if counts is not None:
        ahead, behind = counts
        if ahead and behind:
            raise GitAutomationError(
                f"{repo.name}: histories diverged (ahead={ahead}, behind={behind}); resolve manually"
            )
        if behind:
            if dirty:
                raise GitAutomationError(
                    f"{repo.name}: remote is ahead and local changes exist. "
                    "Commit/stash deliberately before pulling."
                )
            runner.run(
                repo,
                ["pull", "--ff-only", repo.remote, repo.branch],
                mutate=True,
                apply=apply,
                show_output=True,
            )

    candidates = changed_paths(runner, repo)
    check_large_changed_files(runner, repo, settings.max_changed_file_mb)
    if candidates:
        print(f"[{repo.name}] changed paths:")
        for path in candidates:
            print(f"  {path}")
    else:
        print(f"[{repo.name}] no working-tree changes")

    runner.run(repo, ["add", "--all"], mutate=True, apply=apply)
    committed = False
    if apply:
        if index_has_changes(runner, repo):
            runner.run(
                repo,
                ["commit", "-m", message],
                mutate=True,
                apply=True,
                show_output=True,
            )
            committed = True
    elif candidates:
        runner.run(repo, ["commit", "-m", message], mutate=True, apply=False)

    pushed = False
    if not no_push:
        if apply:
            post_counts = divergence(runner, repo)
            should_push = committed or post_counts is None or bool(post_counts[0])
            if should_push:
                push_repository(runner, repo, True)
                pushed = True
        else:
            push_repository(runner, repo, False)
    return {
        "changed_paths": candidates,
        "committed": committed,
        "pushed": pushed,
        "dry_run": not apply,
    }


def status_repository(runner: GitRunner, repo: RepositoryConfig) -> dict[str, Any]:
    require_repository(runner, repo)
    ensure_expected_branch(runner, repo)
    result = runner.run(
        repo,
        ["status", "--short", "--branch"],
        show_output=True,
    )
    return {"status": result.stdout.rstrip()}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Manage ScriptLibrary and SolverLibrary Git repositories safely."
    )
    parser.add_argument(
        "action",
        choices=("doctor", "status", "connect", "adopt", "sync", "snapshot"),
    )
    parser.add_argument(
        "--config",
        default=str(SCRIPT_DIR / "library_repositories.yaml"),
        help="Repository configuration YAML.",
    )
    parser.add_argument("--framework-root", help="Override framework_root from YAML.")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--repository",
        action="append",
        help="Operate only on this repository name; may be repeated.",
    )
    selection.add_argument(
        "--all",
        action="store_true",
        help="Operate on every enabled repository. This must be requested explicitly.",
    )
    parser.add_argument("--message", help="Commit message used by snapshot.")
    parser.add_argument("--apply", action="store_true", help="Execute mutating commands.")
    parser.add_argument("--no-push", action="store_true", help="Do not push during sync/snapshot.")
    parser.add_argument(
        "--update-remote",
        action="store_true",
        help="Allow connect to replace a differing remote URL.",
    )
    parser.add_argument("--report", help="Optional JSON report path.")
    return parser.parse_args(argv)


def select_repositories(
    repositories: list[RepositoryConfig], requested: list[str] | None
) -> list[RepositoryConfig]:
    enabled = [repo for repo in repositories if repo.enabled]
    if not requested:
        return enabled
    wanted = {name.casefold() for name in requested}
    selected = [repo for repo in enabled if repo.name.casefold() in wanted]
    missing = wanted - {repo.name.casefold() for repo in selected}
    if missing:
        raise GitAutomationError(f"unknown or disabled repositories: {', '.join(sorted(missing))}")
    return selected


def write_report(path: Path, report: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] Report: {path}")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        settings, configured = load_configuration(Path(args.config), args.framework_root)
        if not args.repository and not args.all:
            raise GitAutomationError(
                "select one repository with --repository, or use a dedicated wrapper script"
            )
        repositories = select_repositories(configured, args.repository)
        runner = GitRunner(settings.git_command)
        if not runner.executable_exists():
            raise GitAutomationError(f"Git executable was not found: {settings.git_command}")
        if args.action in {"connect", "adopt", "sync", "snapshot"} and not args.apply:
            print("[DRY-RUN] Mutating Git commands will only be displayed. Add --apply to execute them.")

        message = args.message or f"Library update {datetime.now().astimezone():%Y-%m-%d %H:%M}"
        report: dict[str, Any] = {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "action": args.action,
            "apply": bool(args.apply),
            "repositories": [],
        }
        failures = 0
        for repo in repositories:
            try:
                if args.action == "doctor":
                    detail = doctor_repository(runner, repo)
                elif args.action == "status":
                    detail = status_repository(runner, repo)
                elif args.action == "connect":
                    detail = connect_repository(
                        runner,
                        repo,
                        settings,
                        apply=args.apply,
                        update_remote=args.update_remote,
                    )
                elif args.action == "adopt":
                    detail = adopt_existing_remote_history(
                        runner,
                        repo,
                        settings,
                        apply=args.apply,
                        update_remote=args.update_remote,
                    )
                elif args.action == "sync":
                    detail = sync_repository(
                        runner,
                        repo,
                        apply=args.apply,
                        no_push=args.no_push,
                    )
                else:
                    detail = snapshot_repository(
                        runner,
                        repo,
                        settings,
                        message=message,
                        apply=args.apply,
                        no_push=args.no_push,
                    )
                report["repositories"].append(
                    {"name": repo.name, "path": str(repo.path), "result": "ok", **detail}
                )
            except (GitAutomationError, OSError) as exc:
                failures += 1
                print(f"[ERROR] {exc}", file=sys.stderr)
                report["repositories"].append(
                    {
                        "name": repo.name,
                        "path": str(repo.path),
                        "result": "error",
                        "error": str(exc),
                    }
                )

        if args.report:
            write_report(Path(args.report), report)
        if failures:
            print(f"[FAILED] {failures} repository operation(s) failed.", file=sys.stderr)
            return 2
        print(f"[OK] {args.action} completed for {len(repositories)} repositories.")
        return 0
    except (GitAutomationError, OSError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
