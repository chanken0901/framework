#!/usr/bin/env python3
"""Run one Git workflow for ScriptLibrary only."""

from __future__ import annotations

import sys

from manage_library_repositories import main


if __name__ == "__main__":
    raise SystemExit(main([*sys.argv[1:], "--repository", "ScriptLibrary"]))
