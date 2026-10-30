#!/usr/bin/env python3
"""Commit already staged research changes, optionally pushing main (no force)."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def run(*arguments: str, capture: bool = False) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *arguments], cwd=ROOT, text=True, capture_output=capture, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--message", required=True, help="Specific research checkpoint or WIP message.")
    parser.add_argument("--push", action="store_true", help="Push main to origin/main after the commit; also retries an already committed checkpoint.")
    args = parser.parse_args()
    try:
        branch = run("branch", "--show-current", capture=True).stdout.strip()
        if branch != "main":
            raise ValueError("Use the integration checkout on main. This script will not switch branches or discard work.")
        if args.push:
            remotes = run("remote", "get-url", "--push", "--all", "origin", capture=True).stdout.splitlines()
            allowed = {"https://github.com/uedakazushi/poisson-p4.git", "https://github.com/uedakazushi/poisson-p4", "git@github.com:uedakazushi/poisson-p4.git", "ssh://git@github.com/uedakazushi/poisson-p4.git"}
            if len(remotes) != 1 or remotes[0] not in allowed:
                raise ValueError("origin must have exactly one push URL for the intended poisson-p4 repository; inspect the configuration.")
        run("diff", "--cached", "--check")
        staged = run("diff", "--cached", "--name-only", capture=True).stdout.strip()
        if staged:
            print("Committing previously staged files:\n" + staged, flush=True)
            run("commit", "-m", args.message)
        else:
            print("No staged changes; no commit created.")
        if args.push:
            # A normal push rejects conflicting remote updates; no forced reset,
            # automatic rebase, or blind inclusion of unstaged files occurs.
            run("push", "origin", "main")
        return 0
    except (subprocess.CalledProcessError, ValueError) as exc:
        print(f"Checkpoint incomplete: {exc}. Keep local work; inspect and resolve the blocker.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
