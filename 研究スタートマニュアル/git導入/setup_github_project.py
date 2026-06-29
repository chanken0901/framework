#!/usr/bin/env python3
"""
setup_github_project.py

Step 2 script for research project framework.

Purpose
-------
This script is intended to be executed AFTER the project directory has already
been created by init_research_project.py.

It initializes Git, creates the first commit, creates a GitHub repository using
GitHub CLI (gh), adds the GitHub remote, and pushes the initial project to GitHub.

Typical workflow
----------------
Step 1:
    python init_research_project.py --project ProjectName

Step 2:
    cd ProjectName
    python setup_github_project.py --repo-name ProjectName --private

Dry run:
    python setup_github_project.py --repo-name ProjectName --private --dry-run

Notes
-----
- This script requires Git to be installed.
- GitHub repository creation requires GitHub CLI (gh) and authentication.
- If gh is not available, the script can still initialize Git and make the first commit.
- This script does NOT handle large simulation data. Large data should remain outside Git.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


DEFAULT_INITIAL_COMMIT_MESSAGE = "Initial project structure"


def run_command(cmd: List[str], cwd: Path, dry_run: bool = False, check: bool = True) -> subprocess.CompletedProcess:
    """Run a command with clear logging."""
    cmd_text = " ".join(cmd)
    print(f"[CMD] {cmd_text}")

    if dry_run:
        return subprocess.CompletedProcess(cmd, 0, "", "")

    result = subprocess.run(
        cmd,
        cwd=str(cwd),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    if result.stdout.strip():
        print(result.stdout.strip())
    if result.stderr.strip():
        print(result.stderr.strip(), file=sys.stderr)

    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd_text}")

    return result


def command_exists(command: str) -> bool:
    """Return True if command exists in PATH."""
    return shutil.which(command) is not None


def is_git_repository(root: Path) -> bool:
    """Check whether root is already a Git repository."""
    return (root / ".git").exists()


def git_has_commits(root: Path) -> bool:
    """Check whether the repository already has at least one commit."""
    result = subprocess.run(
        ["git", "rev-parse", "--verify", "HEAD"],
        cwd=str(root),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.returncode == 0


def git_is_clean(root: Path) -> bool:
    """Return True if there are no uncommitted changes."""
    result = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=str(root),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stdout.strip() == ""


def remote_exists(root: Path, remote_name: str) -> bool:
    """Return True if a Git remote exists."""
    result = subprocess.run(
        ["git", "remote", "get-url", remote_name],
        cwd=str(root),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.returncode == 0


def configure_git_user(root: Path, name: Optional[str], email: Optional[str], dry_run: bool) -> None:
    """Optionally set local Git user.name and user.email for this repository."""
    if name:
        run_command(["git", "config", "user.name", name], cwd=root, dry_run=dry_run)
    if email:
        run_command(["git", "config", "user.email", email], cwd=root, dry_run=dry_run)


def ensure_git_repository(root: Path, default_branch: str, dry_run: bool) -> None:
    """Initialize Git repository if needed."""
    if is_git_repository(root):
        print("[INFO] Git repository already exists.")
        return

    run_command(["git", "init", "-b", default_branch], cwd=root, dry_run=dry_run)


def initial_commit(root: Path, message: str, dry_run: bool, allow_empty: bool) -> None:
    """Create the first commit if the repository has no commits."""
    if git_has_commits(root):
        print("[INFO] Repository already has commits. Initial commit is skipped.")
        return

    run_command(["git", "add", "."], cwd=root, dry_run=dry_run)

    if dry_run:
        run_command(["git", "commit", "-m", message], cwd=root, dry_run=True)
        return

    if git_is_clean(root) and not allow_empty:
        print("[INFO] No files to commit. Initial commit is skipped.")
        return

    cmd = ["git", "commit", "-m", message]
    if allow_empty:
        cmd.insert(2, "--allow-empty")
    run_command(cmd, cwd=root, dry_run=False)


def ensure_github_auth(root: Path, dry_run: bool) -> bool:
    """Check GitHub CLI authentication."""
    if not command_exists("gh"):
        print("[WARN] GitHub CLI 'gh' was not found. GitHub repository creation will be skipped.")
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
    """Create GitHub repository using gh repo create."""
    if remote_exists(root, remote_name):
        print(f"[INFO] Remote '{remote_name}' already exists.")
        if skip_if_remote_exists:
            print("[INFO] GitHub repository creation is skipped because remote already exists.")
            return
        raise RuntimeError(f"Remote '{remote_name}' already exists. Use --skip-if-remote-exists or remove it manually.")

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
    """Push current branch to GitHub."""
    run_command(["git", "push", "-u", remote_name, branch], cwd=root, dry_run=dry_run)


def print_final_message(root: Path, remote_name: str, branch: str) -> None:
    """Print final usage hints."""
    print("\n[OK] Step 2 Git/GitHub setup completed.")
    print(f"  Project root: {root.resolve()}")
    print(f"  Remote name:  {remote_name}")
    print(f"  Branch:       {branch}")
    print("\nNext daily commands:")
    print("  git status")
    print("  git add .")
    print('  git commit -m "your message"')
    print("  git push")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Initialize Git and register an existing research project on GitHub."
    )

    parser.add_argument(
        "--root",
        default=".",
        help="Project root directory. Default: current directory.",
    )
    parser.add_argument(
        "--repo-name",
        default=None,
        help="GitHub repository name. Default: project directory name.",
    )

    visibility = parser.add_mutually_exclusive_group()
    visibility.add_argument("--private", action="store_true", help="Create a private GitHub repository. Default.")
    visibility.add_argument("--public", action="store_true", help="Create a public GitHub repository.")

    parser.add_argument(
        "--description",
        default=None,
        help="GitHub repository description.",
    )
    parser.add_argument(
        "--branch",
        default="main",
        help="Default branch name. Default: main.",
    )
    parser.add_argument(
        "--remote",
        default="origin",
        help="Git remote name. Default: origin.",
    )
    parser.add_argument(
        "--message",
        default=DEFAULT_INITIAL_COMMIT_MESSAGE,
        help=f"Initial commit message. Default: {DEFAULT_INITIAL_COMMIT_MESSAGE!r}.",
    )
    parser.add_argument(
        "--git-user-name",
        default=None,
        help="Set local git user.name for this repository.",
    )
    parser.add_argument(
        "--git-user-email",
        default=None,
        help="Set local git user.email for this repository.",
    )
    parser.add_argument(
        "--no-github",
        action="store_true",
        help="Only initialize Git and initial commit. Do not create GitHub repository.",
    )
    parser.add_argument(
        "--no-push",
        action="store_true",
        help="Create GitHub repository but do not push.",
    )
    parser.add_argument(
        "--skip-if-remote-exists",
        action="store_true",
        help="Skip GitHub repository creation if the remote already exists.",
    )
    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="Allow an empty initial commit.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()

    if not root.exists():
        print(f"[ERROR] Project root does not exist: {root}")
        return 1

    if not command_exists("git"):
        print("[ERROR] Git was not found. Please install Git first.")
        return 1

    repo_name = args.repo_name or root.name
    visibility = "public" if args.public else "private"

    try:
        ensure_git_repository(root, default_branch=args.branch, dry_run=args.dry_run)
        configure_git_user(root, args.git_user_name, args.git_user_email, args.dry_run)
        initial_commit(root, message=args.message, dry_run=args.dry_run, allow_empty=args.allow_empty)

        if not args.no_github:
            if ensure_github_auth(root, dry_run=args.dry_run):
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

        print_final_message(root, args.remote, args.branch)
        return 0

    except RuntimeError as exc:
        print(f"[ERROR] {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
