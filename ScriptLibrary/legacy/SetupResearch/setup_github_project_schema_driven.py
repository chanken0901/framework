#!/usr/bin/env python3
"""
setup_github_project_schema_driven.py

STEP2: init_research_project_schema_driven.py で生成済みの研究プロジェクトを
Git / GitHub に登録するためのスクリプト。

目的
----
- 既存の研究プロジェクトを Git リポジトリとして初期化する。
- 初回コミットを作成する。
- GitHub CLI (gh) が使える場合は GitHub リポジトリを作成する。
- remote を設定し，初回 push を行う。
- Windows の日本語パス・Git出力で起きやすい UnicodeDecodeError を避ける。

想定ワークフロー
----------------
STEP1:
    python init_research_project_schema_driven.py --schema project_schema.yaml --project tgv

STEP2:
    python setup_github_project_schema_driven.py --root tgv --private

または，プロジェクト内で実行:
    cd tgv
    python ..\\setup_github_project_schema_driven.py --private

GitHub作成なしでGit初期化だけ:
    python setup_github_project_schema_driven.py --root tgv --no-github

事前確認だけ:
    python setup_github_project_schema_driven.py --root tgv --private --dry-run

注意
----
- Git が必要。
- GitHub リポジトリ作成には GitHub CLI (gh) と認証が必要。
- VTK, HDF5, restart, checkpoint, build, solver.exe などの生成物は Git 管理しない想定。
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


DEFAULT_INITIAL_COMMIT_MESSAGE = "Initial schema-driven project structure"


# -----------------------------------------------------------------------------
# Command utilities
# -----------------------------------------------------------------------------


def run_command(
    cmd: List[str],
    cwd: Path,
    dry_run: bool = False,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """
    Run a command with robust Windows-safe decoding.

    Windows環境では，gitやghの出力をcp932として読もうとして
    UnicodeDecodeErrorが出ることがある。
    そのため encoding='utf-8', errors='replace' を明示する。
    """
    cmd_text = " ".join(cmd)
    print(f"[CMD] {cmd_text}")

    if dry_run:
        return subprocess.CompletedProcess(cmd, 0, "", "")

    result = subprocess.run(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )

    if result.stdout and result.stdout.strip():
        print(result.stdout.strip())

    if result.stderr and result.stderr.strip():
        print(result.stderr.strip(), file=sys.stderr)

    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd_text}")

    return result


def command_exists(command: str) -> bool:
    return shutil.which(command) is not None


# -----------------------------------------------------------------------------
# Project checks
# -----------------------------------------------------------------------------


def validate_project_root(root: Path, require_schema: bool) -> None:
    if not root.exists():
        raise RuntimeError(f"Project root does not exist: {root}")

    if not root.is_dir():
        raise RuntimeError(f"Project root is not a directory: {root}")

    if require_schema and not (root / "project_schema.yaml").exists():
        raise RuntimeError(
            "project_schema.yaml was not found. "
            "This does not look like a project created by "
            "init_research_project_schema_driven.py. "
            "Use --no-require-schema to skip this check."
        )


def ensure_gitignore_has_research_rules(root: Path, dry_run: bool = False) -> None:
    """
    Ensure .gitignore contains common research/HPC/generated-output rules.

    init_research_project_schema_driven.py already creates .gitignore,
    but this function makes the Git setup script safer even if .gitignore is old.
    """
    gitignore = root / ".gitignore"

    rules = [
        "",
        "# Added by setup_github_project_schema_driven.py",
        "# Build products",
        "build/",
        "*.o",
        "*.mod",
        "*.smod",
        "*.exe",
        "*.out",
        "",
        "# Large simulation data",
        "*.h5",
        "*.hdf5",
        "*.vtk",
        "*.vtu",
        "*.plt",
        "*.raw",
        "*.bin",
        "",
        "# Case runtime outputs",
        "run/",
        "cases/*/run_env/",
        "cases/*/output/",
        "cases/*/logs/",
        "cases/*/restart/",
        "cases/*/checkpoint/",
        "",
    ]

    existing = gitignore.read_text(encoding="utf-8") if gitignore.exists() else ""

    missing = [rule for rule in rules if rule and rule not in existing]

    if not missing:
        print("[INFO] .gitignore already contains common research rules.")
        return

    print("[INFO] Updating .gitignore with common research rules.")

    if dry_run:
        print("[DRY-RUN] .gitignore would be updated.")
        return

    with gitignore.open("a", encoding="utf-8") as f:
        f.write("\n")
        f.write("\n".join(rules))
        f.write("\n")


# -----------------------------------------------------------------------------
# Git helpers
# -----------------------------------------------------------------------------


def is_git_repository(root: Path) -> bool:
    return (root / ".git").exists()


def git_has_commits(root: Path) -> bool:
    result = subprocess.run(
        ["git", "rev-parse", "--verify", "HEAD"],
        cwd=str(root),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    return result.returncode == 0


def git_is_clean(root: Path) -> bool:
    result = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=str(root),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    return (result.stdout or "").strip() == ""


def current_branch(root: Path) -> Optional[str]:
    result = subprocess.run(
        ["git", "branch", "--show-current"],
        cwd=str(root),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if result.returncode != 0:
        return None
    branch = (result.stdout or "").strip()
    return branch or None


def remote_exists(root: Path, remote_name: str) -> bool:
    result = subprocess.run(
        ["git", "remote", "get-url", remote_name],
        cwd=str(root),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    return result.returncode == 0


def configure_git_user(
    root: Path,
    name: Optional[str],
    email: Optional[str],
    dry_run: bool,
) -> None:
    if name:
        run_command(["git", "config", "user.name", name], cwd=root, dry_run=dry_run)
    if email:
        run_command(["git", "config", "user.email", email], cwd=root, dry_run=dry_run)


def ensure_git_repository(root: Path, default_branch: str, dry_run: bool) -> None:
    if is_git_repository(root):
        print("[INFO] Git repository already exists.")
        return

    # git init -b is supported by modern Git.
    # If it fails, user likely has an old Git; the error message will be shown.
    run_command(["git", "init", "-b", default_branch], cwd=root, dry_run=dry_run)


def initial_commit(
    root: Path,
    message: str,
    dry_run: bool,
    allow_empty: bool,
) -> None:
    if git_has_commits(root):
        print("[INFO] Repository already has commits. Initial commit is skipped.")
        return

    run_command(["git", "add", "."], cwd=root, dry_run=dry_run)

    if dry_run:
        cmd = ["git", "commit", "-m", message]
        if allow_empty:
            cmd.insert(2, "--allow-empty")
        run_command(cmd, cwd=root, dry_run=True)
        return

    if git_is_clean(root) and not allow_empty:
        print("[INFO] No files to commit. Initial commit is skipped.")
        return

    cmd = ["git", "commit", "-m", message]
    if allow_empty:
        cmd.insert(2, "--allow-empty")

    run_command(cmd, cwd=root, dry_run=False)


# -----------------------------------------------------------------------------
# GitHub helpers
# -----------------------------------------------------------------------------


def ensure_github_auth(root: Path, dry_run: bool) -> bool:
    if not command_exists("gh"):
        print("[WARN] GitHub CLI 'gh' was not found. GitHub repository creation will be skipped.")
        print("[HINT] Install GitHub CLI or run with --no-github.")
        return False

    if dry_run:
        run_command(["gh", "auth", "status"], cwd=root, dry_run=True, check=False)
        return True

    result = run_command(["gh", "auth", "status"], cwd=root, dry_run=False, check=False)

    if result.returncode != 0:
        print("[WARN] GitHub CLI is installed, but authentication is not complete.")
        print("[HINT] Run: gh auth login")
        return False

    return True


def create_github_repo(
    root: Path,
    repo_name: str,
    visibility: str,
    remote_name: str,
    description: Optional[str],
    dry_run: bool,
    skip_if_remote_exists: bool,
) -> None:
    if remote_exists(root, remote_name):
        print(f"[INFO] Remote '{remote_name}' already exists.")
        if skip_if_remote_exists:
            print("[INFO] GitHub repository creation is skipped because remote already exists.")
            return
        raise RuntimeError(
            f"Remote '{remote_name}' already exists. "
            "Use --skip-if-remote-exists or remove it manually."
        )

    cmd = [
        "gh",
        "repo",
        "create",
        repo_name,
        f"--{visibility}",
        "--source=.",
        f"--remote={remote_name}",
    ]

    if description:
        cmd.extend(["--description", description])

    run_command(cmd, cwd=root, dry_run=dry_run)


def push_initial(root: Path, remote_name: str, branch: str, dry_run: bool) -> None:
    branch_to_push = current_branch(root) or branch
    run_command(["git", "push", "-u", remote_name, branch_to_push], cwd=root, dry_run=dry_run)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def print_final_message(root: Path, remote_name: str, branch: str, github_enabled: bool) -> None:
    print("\n[OK] Git/GitHub setup completed.")
    print(f"  Project root: {root.resolve()}")
    print(f"  Remote name:  {remote_name}")
    print(f"  Branch:       {current_branch(root) or branch}")

    if not github_enabled:
        print("\n[INFO] GitHub repository creation/push was skipped or unavailable.")

    print("\nDaily commands:")
    print("  git status")
    print("  git add .")
    print('  git commit -m "your message"')
    print("  git push")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Initialize Git and register a schema-driven research project on GitHub."
    )

    parser.add_argument(
        "--root",
        default=".",
        help="Project root directory generated by init_research_project_schema_driven.py. Default: current directory.",
    )
    parser.add_argument(
        "--repo-name",
        default=None,
        help="GitHub repository name. Default: project directory name.",
    )

    visibility = parser.add_mutually_exclusive_group()
    visibility.add_argument("--private", action="store_true", help="Create a private GitHub repository. Default.")
    visibility.add_argument("--public", action="store_true", help="Create a public GitHub repository.")

    parser.add_argument("--description", default=None, help="GitHub repository description.")
    parser.add_argument("--branch", default="main", help="Default branch name. Default: main.")
    parser.add_argument("--remote", default="origin", help="Git remote name. Default: origin.")
    parser.add_argument(
        "--message",
        default=DEFAULT_INITIAL_COMMIT_MESSAGE,
        help=f"Initial commit message. Default: {DEFAULT_INITIAL_COMMIT_MESSAGE!r}.",
    )

    parser.add_argument("--git-user-name", default=None, help="Set local git user.name for this repository.")
    parser.add_argument("--git-user-email", default=None, help="Set local git user.email for this repository.")

    parser.add_argument("--no-github", action="store_true", help="Only initialize Git and initial commit.")
    parser.add_argument("--no-push", action="store_true", help="Create GitHub repository but do not push.")
    parser.add_argument("--skip-if-remote-exists", action="store_true", help="Skip GitHub creation if remote exists.")
    parser.add_argument("--allow-empty", action="store_true", help="Allow an empty initial commit.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them.")

    parser.add_argument(
        "--no-require-schema",
        action="store_true",
        help="Do not require project_schema.yaml in project root.",
    )
    parser.add_argument(
        "--no-update-gitignore",
        action="store_true",
        help="Do not append common research rules to .gitignore.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()

    try:
        validate_project_root(root, require_schema=not args.no_require_schema)

        if not command_exists("git"):
            print("[ERROR] Git was not found. Please install Git first.")
            return 1

        repo_name = args.repo_name or root.name
        visibility = "public" if args.public else "private"

        if not args.no_update_gitignore:
            ensure_gitignore_has_research_rules(root, dry_run=args.dry_run)

        ensure_git_repository(root, default_branch=args.branch, dry_run=args.dry_run)
        configure_git_user(root, args.git_user_name, args.git_user_email, args.dry_run)
        initial_commit(
            root,
            message=args.message,
            dry_run=args.dry_run,
            allow_empty=args.allow_empty,
        )

        github_attempted = False

        if not args.no_github:
            if ensure_github_auth(root, dry_run=args.dry_run):
                github_attempted = True
                create_github_repo(
                    root=root,
                    repo_name=repo_name,
                    visibility=visibility,
                    remote_name=args.remote,
                    description=args.description,
                    dry_run=args.dry_run,
                    skip_if_remote_exists=args.skip_if_remote_exists,
                )
                if not args.no_push:
                    push_initial(root, remote_name=args.remote, branch=args.branch, dry_run=args.dry_run)
            else:
                print("[WARN] GitHub setup was skipped. Git initialization is complete.")

        print_final_message(root, args.remote, args.branch, github_enabled=github_attempted and not args.no_push)
        return 0

    except RuntimeError as exc:
        print(f"[ERROR] {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
