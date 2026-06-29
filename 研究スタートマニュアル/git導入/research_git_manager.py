#!/usr/bin/env python3
"""
research_git_manager.py

Research-oriented Git helper written in pure Python.

Purpose
-------
This script provides a safer Git workflow for simulation-based research projects.
It is designed to help prevent:

  - forgetting to commit changes before running simulations,
  - committing code that does not compile,
  - forgetting to push important commits to GitHub,
  - losing the connection between simulation results and source code versions.

Main commands
-------------
  status       Show Git status and basic repository information.
  check-clean  Stop if there are uncommitted changes.
  verify       Run verification commands such as compilation or tests.
  save         Run verification, then git add, git commit, and optionally git push.
  sync         Pull the latest changes from GitHub using --ff-only.
  info         Print current branch, commit hash, and remote information.

Recommended workflow
--------------------
  python scripts/research_git_manager.py status
  python scripts/research_git_manager.py save -m "Update solver" --cmd "make"
  python scripts/research_git_manager.py check-clean

Examples
--------
Commit only if compilation succeeds:

  python scripts/research_git_manager.py save -m "Fix boundary condition" --cmd "make"

Commit only if multiple checks succeed:

  python scripts/research_git_manager.py save -m "Update case framework" \
      --cmd "python scripts/generate_case_docs.py --case cases/case0001" \
      --cmd "make"

Commit but do not push:

  python scripts/research_git_manager.py save -m "Work in progress" --no-push

Check that the repository is clean before build/run:

  python scripts/research_git_manager.py check-clean

Notes
-----
- This script uses only Python standard library modules.
- Git itself must be installed and available from PATH.
- Verification commands are executed before git add/commit.
- If any verification command fails, the commit is cancelled.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Optional, Sequence


CONFIG_FILE = ".research_git.json"


class CommandError(RuntimeError):
    """Raised when an external command fails."""


class GitManager:
    def __init__(self, repo: Path, dry_run: bool = False, verbose: bool = True) -> None:
        self.repo = repo.resolve()
        self.dry_run = dry_run
        self.verbose = verbose

    def run(
        self,
        args: Sequence[str],
        *,
        check: bool = True,
        capture: bool = False,
        shell: bool = False,
    ) -> subprocess.CompletedProcess:
        """Run a command in the repository root."""
        if self.verbose:
            if shell:
                print(f"[cmd] {args if isinstance(args, str) else ' '.join(args)}")
            else:
                print("[cmd] " + " ".join(shlex.quote(str(a)) for a in args))

        if self.dry_run:
            return subprocess.CompletedProcess(args=args, returncode=0, stdout="", stderr="")

        result = subprocess.run(
            args,
            cwd=self.repo,
            text=True,
            shell=shell,
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.PIPE if capture else None,
        )

        if check and result.returncode != 0:
            if capture:
                if result.stdout:
                    print(result.stdout)
                if result.stderr:
                    print(result.stderr, file=sys.stderr)
            raise CommandError(f"Command failed with exit code {result.returncode}: {args}")

        return result

    def git(self, *args: str, check: bool = True, capture: bool = False) -> subprocess.CompletedProcess:
        return self.run(["git", *args], check=check, capture=capture)

    def ensure_git_repo(self) -> None:
        result = self.git("rev-parse", "--is-inside-work-tree", check=False, capture=True)
        if result.returncode != 0 or result.stdout.strip() != "true":
            raise CommandError(
                "This directory is not inside a Git repository.\n"
                "Run 'git init' first, or move to the project root."
            )

    def branch(self) -> str:
        result = self.git("branch", "--show-current", capture=True)
        return result.stdout.strip()

    def head_hash(self, short: bool = False) -> str:
        args = ["rev-parse"]
        if short:
            args.append("--short")
        args.append("HEAD")
        result = self.git(*args, capture=True)
        return result.stdout.strip()

    def remote_url(self) -> str:
        result = self.git("remote", "get-url", "origin", check=False, capture=True)
        return result.stdout.strip() if result.returncode == 0 else ""

    def status_porcelain(self) -> str:
        result = self.git("status", "--porcelain", capture=True)
        return result.stdout

    def has_uncommitted_changes(self) -> bool:
        return bool(self.status_porcelain().strip())

    def has_staged_changes(self) -> bool:
        result = self.git("diff", "--cached", "--quiet", check=False)
        return result.returncode != 0

    def print_status(self) -> None:
        self.ensure_git_repo()
        print("\n[Repository]")
        print(f"  path:   {self.repo}")
        print(f"  branch: {self.branch()}")
        print(f"  commit: {self.head_hash(short=True)}")
        remote = self.remote_url()
        print(f"  remote: {remote if remote else '(no origin remote)'}")
        print("\n[Git status]")
        self.git("status", check=True)

    def print_info(self) -> None:
        self.ensure_git_repo()
        info = {
            "repo": str(self.repo),
            "branch": self.branch(),
            "commit_hash": self.head_hash(short=False),
            "commit_hash_short": self.head_hash(short=True),
            "remote_origin": self.remote_url(),
            "generated_at": datetime.now().isoformat(timespec="seconds"),
        }
        print(json.dumps(info, indent=2, ensure_ascii=False))

    def check_clean(self) -> None:
        self.ensure_git_repo()
        if self.has_uncommitted_changes():
            print("\n[ERROR] Uncommitted changes detected.\n")
            self.git("status", "--short", check=False)
            print("\nPlease save the current state before build/run, for example:")
            print('  python scripts/research_git_manager.py save -m "your message" --cmd "make"')
            raise CommandError("Repository is not clean.")
        print("[OK] Repository is clean.")

    def verify(self, commands: Sequence[str], allow_empty: bool = False) -> None:
        self.ensure_git_repo()
        if not commands:
            if allow_empty:
                print("[WARN] No verification commands were specified. Skipping verification.")
                return
            raise CommandError(
                "No verification command was specified.\n"
                "Use --cmd, for example: --cmd \"make\""
            )

        print("\n[verify] Starting verification...")
        for cmd in commands:
            print(f"\n[verify] {cmd}")
            self.run(cmd, shell=True, check=True)
        print("\n[OK] Verification succeeded.")

    def save(
        self,
        message: str,
        commands: Sequence[str],
        paths: Sequence[str],
        push: bool = True,
        allow_empty_verify: bool = False,
        allow_empty_commit: bool = False,
    ) -> None:
        self.ensure_git_repo()

        # Verification is intentionally done before staging/commit.
        # If compilation or tests fail, no commit is created.
        self.verify(commands, allow_empty=allow_empty_verify)

        print("\n[save] Staging files...")
        if paths:
            self.git("add", *paths)
        else:
            self.git("add", ".")

        if not self.has_staged_changes() and not allow_empty_commit:
            print("[INFO] No staged changes. Nothing to commit.")
            return

        print("\n[save] Creating commit...")
        commit_args = ["commit", "-m", message]
        if allow_empty_commit:
            commit_args.append("--allow-empty")
        self.git(*commit_args)

        print("\n[save] Commit created.")
        print(f"  commit: {self.head_hash(short=True)}")

        if push:
            print("\n[save] Pushing to remote...")
            self.git("push")
            print("[OK] Push completed.")
        else:
            print("[INFO] Push skipped because --no-push was specified.")

    def sync(self) -> None:
        self.ensure_git_repo()
        self.check_clean()
        print("\n[sync] Pulling latest changes with --ff-only...")
        self.git("pull", "--ff-only")
        print("[OK] Sync completed.")


def load_config(repo: Path) -> dict:
    path = repo / CONFIG_FILE
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def commands_from_args_or_config(args: argparse.Namespace, repo: Path) -> List[str]:
    if args.cmd:
        return list(args.cmd)
    config = load_config(repo)
    return list(config.get("verify_commands", []))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Research-oriented Git helper for safer commit/build/run workflows."
    )
    parser.add_argument(
        "--repo",
        default=".",
        help="Repository root. Default: current directory.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Reduce command printing.",
    )

    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("status", help="Show git status.")
    sub.add_parser("info", help="Show repository information as JSON.")
    sub.add_parser("check-clean", help="Fail if the repository has uncommitted changes.")

    p_verify = sub.add_parser("verify", help="Run verification commands.")
    p_verify.add_argument(
        "--cmd",
        action="append",
        help="Verification command. Can be specified multiple times.",
    )
    p_verify.add_argument(
        "--allow-empty",
        action="store_true",
        help="Allow no verification commands.",
    )

    p_save = sub.add_parser("save", help="Verify, then git add/commit/push.")
    p_save.add_argument(
        "-m",
        "--message",
        required=True,
        help="Commit message.",
    )
    p_save.add_argument(
        "--cmd",
        action="append",
        help="Verification command. Can be specified multiple times.",
    )
    p_save.add_argument(
        "--path",
        action="append",
        help="Path to stage. Can be specified multiple times. Default: git add .",
    )
    p_save.add_argument(
        "--no-push",
        action="store_true",
        help="Create commit but do not push.",
    )
    p_save.add_argument(
        "--allow-empty-verify",
        action="store_true",
        help="Allow save without verification commands.",
    )
    p_save.add_argument(
        "--allow-empty-commit",
        action="store_true",
        help="Allow an empty commit.",
    )

    sub.add_parser("sync", help="Pull latest changes using --ff-only after clean check.")

    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo = Path(args.repo)
    manager = GitManager(repo, dry_run=args.dry_run, verbose=not args.quiet)

    try:
        if args.command == "status":
            manager.print_status()
        elif args.command == "info":
            manager.print_info()
        elif args.command == "check-clean":
            manager.check_clean()
        elif args.command == "verify":
            cmds = commands_from_args_or_config(args, manager.repo)
            manager.verify(cmds, allow_empty=args.allow_empty)
        elif args.command == "save":
            cmds = commands_from_args_or_config(args, manager.repo)
            paths = args.path or []
            manager.save(
                message=args.message,
                commands=cmds,
                paths=paths,
                push=not args.no_push,
                allow_empty_verify=args.allow_empty_verify,
                allow_empty_commit=args.allow_empty_commit,
            )
        elif args.command == "sync":
            manager.sync()
        else:
            raise CommandError(f"Unknown command: {args.command}")
    except CommandError as exc:
        print(f"\n[FAILED] {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\n[FAILED] Interrupted by user.", file=sys.stderr)
        return 130

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
