#!/usr/bin/env python3
"""
setup_solverlibrary_github.py

SolverLibrary を Git / GitHub に登録するためのスクリプト。

使い方:
  python setup_solverlibrary_github.py --root SolverLibrary --private

SolverLibrary の中で実行する場合:
  python ../setup_solverlibrary_github.py --private

GitHubなしでGit初期化だけ:
  python setup_solverlibrary_github.py --root SolverLibrary --no-github

注意:
  - Git が必要です。
  - GitHub 登録には GitHub CLI (gh) と認証が必要です。
  - Windows の UnicodeDecodeError 対策として subprocess の encoding を明示しています。
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional


DEFAULT_INITIAL_COMMIT_MESSAGE = "Initial SolverLibrary structure"


def run_command(
    cmd: List[str],
    cwd: Path,
    dry_run: bool = False,
    check: bool = True,
) -> subprocess.CompletedProcess:
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


def validate_root(root: Path, create_root: bool) -> None:
    if root.exists() and not root.is_dir():
        raise RuntimeError(f"Root exists but is not a directory: {root}")

    if not root.exists():
        if create_root:
            root.mkdir(parents=True, exist_ok=True)
            print(f"[OK] Created SolverLibrary root: {root}")
        else:
            raise RuntimeError(f"SolverLibrary root does not exist: {root}")


def detect_solvers(root: Path) -> list[str]:
    solvers: list[str] = []
    for child in root.iterdir():
        if not child.is_dir():
            continue
        if child.name.startswith("."):
            continue
        if child.name.lower() in {"build", "run", "logs", "tmp", "temp", "__pycache__"}:
            continue
        if (child / "src").exists():
            solvers.append(child.name)
    return sorted(solvers)


def ensure_gitignore(root: Path, dry_run: bool = False) -> None:
    path = root / ".gitignore"
    rules = """# Build products
build/
**/build/
*.o
*.obj
*.mod
*.smod
*.a
*.lib
*.so
*.dll
*.dylib
*.exe
*.out

# Runtime outputs
run/
**/run/
logs/
**/logs/
output/
**/output/
restart/
**/restart/
checkpoint/
**/checkpoint/

# Large data
*.h5
*.hdf5
*.vtk
*.vtu
*.plt
*.raw
*.bin

# Python
__pycache__/
*.pyc
.venv/
venv/

# Editor / OS
.vscode/
.idea/
.DS_Store
Thumbs.db

# Temporary
tmp/
temp/
*.tmp
*.bak
"""
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if "Added by setup_solverlibrary_github.py" in existing:
        print("[INFO] .gitignore already updated.")
        return

    print("[INFO] Updating .gitignore for SolverLibrary.")
    if dry_run:
        print("[DRY-RUN] .gitignore would be updated.")
        return

    with path.open("a", encoding="utf-8") as f:
        if existing and not existing.endswith("\n"):
            f.write("\n")
        f.write("\n# Added by setup_solverlibrary_github.py\n")
        f.write(rules)


def ensure_readme(root: Path, repo_name: str, dry_run: bool = False) -> None:
    path = root / "README.md"
    if path.exists():
        print("[INFO] README.md already exists.")
        return

    solvers = detect_solvers(root)
    solver_lines = "\n".join(f"- {name}" for name in solvers) if solvers else "- NSE\n- GPE\n- Shared"

    text = f"""# {repo_name}

SolverLibrary is a shared solver development repository.

## Solvers

{solver_lines}

## Recommended structure

```text
SolverLibrary/
├── NSE/
│   └── src/
├── GPE/
│   └── src/
└── Shared/
```

## Role

This repository stores reusable solver source code.

Research projects import solver code from this library into their own `src/` directory, then generate a project-specific Makefile.
"""

    print("[INFO] Creating README.md.")
    if dry_run:
        print("[DRY-RUN] README.md would be created.")
        return

    path.write_text(text, encoding="utf-8")


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
    return (result.stdout or "").strip() or None


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


def get_git_config(root: Path, key: str) -> Optional[str]:
    result = subprocess.run(
        ["git", "config", "--get", key],
        cwd=str(root),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if result.returncode != 0:
        return None
    return (result.stdout or "").strip() or None


def configure_git_user(root: Path, name: Optional[str], email: Optional[str], dry_run: bool) -> None:
    if name:
        run_command(["git", "config", "user.name", name], cwd=root, dry_run=dry_run)
    if email:
        run_command(["git", "config", "user.email", email], cwd=root, dry_run=dry_run)


def ensure_git_identity(root: Path) -> None:
    name = get_git_config(root, "user.name")
    email = get_git_config(root, "user.email")

    if name and email:
        print(f"[INFO] Git user.name : {name}")
        print(f"[INFO] Git user.email: {email}")
        return

    raise RuntimeError(
        "Git user identity is not configured.\n\n"
        "Run:\n"
        "  git config --global user.name \"Kento Tanaka\"\n"
        "  git config --global user.email \"k.tanaka.fluid@gmail.com\"\n\n"
        "or pass:\n"
        "  --git-user-name \"Kento Tanaka\" --git-user-email \"k.tanaka.fluid@gmail.com\""
    )


def ensure_git_repository(root: Path, branch: str, dry_run: bool) -> None:
    if is_git_repository(root):
        print("[INFO] Git repository already exists.")
        return
    run_command(["git", "init", "-b", branch], cwd=root, dry_run=dry_run)


def initial_commit(root: Path, message: str, dry_run: bool, allow_empty: bool) -> None:
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
    run_command(cmd, cwd=root)


def ensure_github_auth(root: Path, dry_run: bool) -> bool:
    if not command_exists("gh"):
        print("[WARN] GitHub CLI 'gh' was not found. GitHub repository creation will be skipped.")
        return False

    if dry_run:
        run_command(["gh", "auth", "status"], cwd=root, dry_run=True, check=False)
        return True

    result = run_command(["gh", "auth", "status"], cwd=root, check=False)
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
        raise RuntimeError(f"Remote '{remote_name}' already exists.")

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Initialize Git and register SolverLibrary on GitHub.")

    parser.add_argument("--root", default=".", help="SolverLibrary root directory. Default: current directory.")
    parser.add_argument("--repo-name", default=None, help="GitHub repository name. Default: root directory name.")

    visibility = parser.add_mutually_exclusive_group()
    visibility.add_argument("--private", action="store_true", help="Create a private GitHub repository. Default.")
    visibility.add_argument("--public", action="store_true", help="Create a public GitHub repository.")

    parser.add_argument("--description", default="Reusable solver library for research projects.")
    parser.add_argument("--branch", default="main")
    parser.add_argument("--remote", default="origin")
    parser.add_argument("--message", default=DEFAULT_INITIAL_COMMIT_MESSAGE)

    parser.add_argument("--git-user-name", default=None)
    parser.add_argument("--git-user-email", default=None)

    parser.add_argument("--no-github", action="store_true", help="Only initialize Git and initial commit.")
    parser.add_argument("--no-push", action="store_true", help="Create GitHub repository but do not push.")
    parser.add_argument("--skip-if-remote-exists", action="store_true")
    parser.add_argument("--allow-empty", action="store_true")
    parser.add_argument("--dry-run", action="store_true")

    parser.add_argument("--create-root", action="store_true", help="Create SolverLibrary root if needed.")
    parser.add_argument("--no-readme", action="store_true", help="Do not create README.md.")
    parser.add_argument("--no-update-gitignore", action="store_true", help="Do not create/update .gitignore.")
    parser.add_argument("--no-require-identity", action="store_true", help="Do not check git user.name/email.")

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.root).resolve()

    try:
        validate_root(root, create_root=args.create_root)

        if not command_exists("git"):
            print("[ERROR] Git was not found. Please install Git first.")
            return 1

        repo_name = args.repo_name or root.name
        visibility = "public" if args.public else "private"

        if not args.no_update_gitignore:
            ensure_gitignore(root, dry_run=args.dry_run)

        if not args.no_readme:
            ensure_readme(root, repo_name=repo_name, dry_run=args.dry_run)

        ensure_git_repository(root, branch=args.branch, dry_run=args.dry_run)
        configure_git_user(root, args.git_user_name, args.git_user_email, args.dry_run)

        if not args.no_require_identity and not args.dry_run:
            ensure_git_identity(root)

        initial_commit(root, message=args.message, dry_run=args.dry_run, allow_empty=args.allow_empty)

        github_done = False
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
                    github_done = True
            else:
                print("[WARN] GitHub setup was skipped. Git initialization is complete.")

        print("\n[OK] SolverLibrary Git/GitHub setup completed.")
        print(f"  Root:   {root}")
        print(f"  Branch: {current_branch(root) or args.branch}")
        solvers = detect_solvers(root)
        if solvers:
            print("  Solvers:")
            for solver in solvers:
                print(f"    - {solver}")
        if not github_done:
            print("\n[INFO] GitHub repository creation/push was skipped or unavailable.")

        return 0

    except RuntimeError as exc:
        print(f"[ERROR] {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
