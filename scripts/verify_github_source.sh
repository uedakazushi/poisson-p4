#!/usr/bin/env bash
# No GitHub Actions, no Library I/O, no network download; manual reproducibility.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
MODE="${1:-incidence}"
case "$MODE" in
  incidence)
    bash scripts/reproduce_subvariety_incidence.sh
    ;;
  full)
    mkdir -p work
    python3 code/audit_required_dependencies.py --output work/github-source-dependencies.json
    python3 code/check_component_catalogue.py --output work/github-component-catalogue.json
    python3 code/check_global_structure.py --output work/github-global-structure.json
    bash scripts/build_structural_manuscripts.sh
    bash scripts/reproduce_subvariety_incidence.sh
    ;;
  *)
    echo 'Usage: bash scripts/verify_github_source.sh [incidence|full]' >&2
    exit 2
    ;;
esac
