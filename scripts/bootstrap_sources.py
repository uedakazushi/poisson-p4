#!/usr/bin/env python3
"""Verify and materialize immutable v02--v24 source packages (stdlib only)."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "sources"
LEGACY = ROOT / "legacy"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_relative(name: str) -> Path:
    path = PurePosixPath(name)
    if not name or path.is_absolute() or any(p in ("..", ".", "") for p in path.parts):
        raise ValueError(f"Unsafe archive member: {name!r}")
    if "\\" in name or ":" in name:
        raise ValueError(f"Unsupported archive member: {name!r}")
    return Path(*path.parts)


def load_manifest() -> dict:
    manifest = json.loads((SOURCES / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported source manifest schema")
    versions = [a["version"] for a in manifest["archives"]]
    if versions != list(range(2, 25)):
        raise ValueError("Expected exactly the supplied versions v02 through v24")
    if sum(a["file_count"] for a in manifest["archives"]) != manifest["source_bundle"]["member_count"]:
        raise ValueError("Inconsistent source member count")
    return manifest


def verify_archive(item: dict) -> None:
    path = SOURCES / safe_relative(item["path"])
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"Missing regular archive: {path}")
    if path.stat().st_size != item["size_bytes"] or sha256(path) != item["sha256"]:
        raise ValueError(f"Archive checksum mismatch: {path}")


def verify_version(destination: Path, item: dict) -> None:
    if not destination.is_dir() or destination.is_symlink():
        raise ValueError(f"Missing regular version directory: {destination}")
    expected = set()
    for member in item["members"]:
        relative = safe_relative(member["path"])
        expected.add(relative.as_posix())
        path = destination / relative
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"Missing regular source file: {path}")
        if path.stat().st_size != member["size_bytes"] or sha256(path) != member["sha256"]:
            raise ValueError(f"Source changed; refusing to overwrite: {path}")
    # rglob deliberately does not follow the version-history symlinks.
    actual = {p.relative_to(destination).as_posix() for p in destination.rglob("*") if p.is_file()}
    if actual != expected:
        raise ValueError(f"Unexpected or missing files in {destination}: {sorted(actual ^ expected)[:8]}")


def extract_version(item: dict) -> None:
    destination = LEGACY / f"v{item['version']:02d}"
    if destination.exists() or destination.is_symlink():
        verify_version(destination, item)
        return
    with tempfile.TemporaryDirectory(prefix=".source-stage-", dir=LEGACY) as stage:
        materialized = Path(stage) / destination.name
        materialized.mkdir()
        expected = {member["path"]: member for member in item["members"]}
        if len(expected) != item["file_count"]:
            raise ValueError(f"Duplicate manifest members in {item['path']}")
        with zipfile.ZipFile(SOURCES / item["path"]) as archive:
            names = archive.namelist()
            if len(names) != len(set(names)) or set(names) != set(expected):
                raise ValueError(f"Archive membership mismatch: {item['path']}")
            for name in names:
                member = expected[name]
                relative = safe_relative(name)
                payload = archive.read(name)
                if len(payload) != member["size_bytes"] or hashlib.sha256(payload).hexdigest() != member["sha256"]:
                    raise ValueError(f"Archive member checksum mismatch: {item['path']}:{name}")
                path = materialized / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(payload)
        materialized.rename(destination)


def history_links(manifest: dict, check: bool) -> None:
    for item in manifest["archives"]:
        previous = item.get("previous_version")
        if previous is None:
            continue
        destination = LEGACY / f"v{item['version']:02d}"
        link = destination / f"previous_v{previous}"
        target = LEGACY / f"v{previous:02d}"
        if link.is_symlink():
            if link.resolve() != target.resolve():
                raise ValueError(f"Incorrect previous-version link: {link}")
        elif link.exists() or check:
            raise ValueError(f"Missing or conflicting previous-version link: {link}")
        else:
            link.symlink_to(Path("..") / target.name, target_is_directory=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Only verify archives, every extracted file, and history links.")
    args = parser.parse_args()
    try:
        manifest = load_manifest()
        for item in manifest["archives"]:
            verify_archive(item)
        if LEGACY.is_symlink():
            raise ValueError("legacy/ must not itself be a symlink")
        if args.check:
            for item in manifest["archives"]:
                verify_version(LEGACY / f"v{item['version']:02d}", item)
        else:
            LEGACY.mkdir(exist_ok=True)
            for item in manifest["archives"]:
                extract_version(item)
        history_links(manifest, args.check)
        print(json.dumps({"all_passed": True, "mode": "check" if args.check else "materialize", "versions": len(manifest["archives"]), "files": manifest["source_bundle"]["member_count"], "source_bundle_sha256": manifest["source_bundle"]["sha256"]}, ensure_ascii=False, indent=2))
        return 0
    except (OSError, ValueError, KeyError, zipfile.BadZipFile) as exc:
        print(json.dumps({"all_passed": False, "error": str(exc)}, ensure_ascii=False, indent=2), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
