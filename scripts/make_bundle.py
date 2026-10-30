#!/usr/bin/env python3
"""Build a versioned, self-contained continuation ZIP from intended repo files."""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import zipfile
from zoneinfo import ZoneInfo

from bootstrap_sources import load_manifest, verify_archive

ROOT = Path(__file__).resolve().parents[1]


def git(*arguments: str) -> str:
    return subprocess.check_output(["git", *arguments], cwd=ROOT, text=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", default="v26-s61", help="Version label (letters, digits, dot, underscore, hyphen).")
    parser.add_argument("--output", type=Path, help="ZIP path. Default: ignored bundles/ directory.")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9._-]+", args.version):
        parser.error("Invalid version label")
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    stem = f"poisson_p4_continuation_{now:%Y%m%d}_{args.version}"
    output = args.output or ROOT / "bundles" / f"{stem}.zip"
    if not output.is_absolute():
        output = ROOT / output
    if output.exists():
        parser.error(f"Refusing to overwrite an existing bundle: {output}")
    # Git's exclude rules keep dependencies, extraction, credentials, and local
    # build products out. Include intended new files so pre-commit bundles work.
    candidates = sorted(set(git("ls-files", "--cached", "--others", "--exclude-standard", "-z").split("\0")) - {""})
    required = {"AGENTS.md", "RESEARCH_POLICY.md", "research/STATE.md", "research/NEXT_TASKS.md", "research/CLAIMS.json", "research/COMPONENTS.json", "papers/README.md", "sources/manifest.json", "scripts/setup.sh", "scripts/bootstrap_sources.py", "scripts/doctor.py", "scripts/verify.py", "requirements.txt"}
    files = []
    for name in candidates:
        # The manifest is generated below for this exact snapshot; importing a
        # previous bundle's manifest would create a duplicate ZIP entry.
        if name == "BUNDLE_MANIFEST.json":
            continue
        path = ROOT / name
        if path.resolve() == output.resolve() or not path.is_file():
            continue
        if path.is_symlink():
            parser.error(f"Unexpected tracked symlink: {name}; materialized sources must stay in ignored legacy/.")
        files.append((name, path))
    try:
        source_manifest = load_manifest()
        for item in source_manifest["archives"]:
            verify_archive(item)
            required.add("sources/" + item["path"])
    except (OSError, ValueError, KeyError) as exc:
        parser.error(f"Source archive verification failed: {exc}")
    included = {name for name, _ in files}
    if not required.issubset(included):
        parser.error(f"Missing actual continuation files: {sorted(required - included)}")
    try:
        head = git("rev-parse", "HEAD").strip()
    except subprocess.CalledProcessError:
        head = None
    status = git("status", "--porcelain").splitlines()
    manifest = {
        "schema_version": 1, "date": now.isoformat(), "version": args.version,
        "repository": "https://github.com/uedakazushi/poisson-p4", "base_commit": head,
        "worktree_clean": not status,
        "scope": "Repository snapshot including all versioned research source archives. Run scripts/setup.sh after extraction. Not a new claim of mathematical completeness.",
        "library_save": "Not performed by this script; save with an available authenticated Library capability and record its receipt.",
        "files": [],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        with zipfile.ZipFile(output, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
            for name, path in files:
                payload = path.read_bytes()
                archive.writestr(f"{stem}/{name}", payload)
                manifest["files"].append({"path": name, "size_bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()})
            archive.writestr(f"{stem}/BUNDLE_MANIFEST.json", json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
        with zipfile.ZipFile(output) as archive:
            damaged = archive.testzip()
            if damaged:
                raise ValueError(f"CRC failure: {damaged}")
        print(json.dumps({"path": str(output), "files": len(files), "bytes": output.stat().st_size, "sha256": hashlib.sha256(output.read_bytes()).hexdigest(), "base_commit": head, "worktree_clean": not status}, ensure_ascii=False, indent=2))
        return 0
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        print(f"Bundle failed: {exc}. An incomplete output may remain at {output}.", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
