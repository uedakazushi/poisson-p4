"""Verify v21 preservation and the v22 SHA256 manifest (Python stdlib)."""
import argparse
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path)
    ap.add_argument('--preservation-only', action='store_true')
    args = ap.parse_args()
    old = json.loads((BASE/'data/previous_v21_preservation.json').read_text())
    previous = BASE/'previous_v21'
    files = {str(p.relative_to(previous)) for p in previous.rglob('*') if p.is_file()}
    assert files == set(old['files']), 'Previous-edition file inventory changed'
    for name, expected in old['files'].items():
        assert digest(previous/name) == expected, f'Previous file changed: {name}'
    assert digest(BASE/'previous_v21.pdf') == old['pdf_sha256'], 'Previous PDF changed'
    count = 0
    if not args.preservation_only:
        manifest = BASE/'SHA256SUMS.txt'
        recorded = {}
        for line in manifest.read_text().splitlines():
            sha, name = line.split('  ', 1)
            path = BASE/name
            assert path.is_file() and digest(path) == sha, f'Manifest mismatch: {name}'
            recorded[name] = sha
        actual = {str(p.relative_to(BASE)) for p in BASE.rglob('*')
                  if p.is_file() and p != manifest
                  and not (p.relative_to(BASE).parts[0] == 'code'
                           and '__pycache__' in p.relative_to(BASE).parts)}
        # New Python caches are excluded; any caches in the preserved v21 tree
        # remain ordinary immutable files and are included in its own inventory.
        assert set(recorded) == actual, 'Manifest inventory mismatch'
        count = len(recorded)
    report = {'all_passed': True, 'previous_files_unchanged': len(files),
              'previous_pdf_unchanged': True, 'manifest_files_checked': count,
              'scope': 'Byte preservation and file integrity, not a mathematical proof'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
