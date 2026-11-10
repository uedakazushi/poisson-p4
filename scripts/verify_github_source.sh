#!/usr/bin/env bash
# No GitHub Actions, Library I/O or network downloads; manual reproducibility.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
MODE="${1:-incidence}"
case "$MODE" in
  tangent-cone)
    bash scripts/verify_isotropic_tangent_cone.sh
    ;;
  incidence)
    bash scripts/reproduce_subvariety_incidence.sh
    bash scripts/verify_isotropic_tangent_cone.sh
    ;;
  full)
    mkdir -p work
    python3 code/audit_required_dependencies.py --output work/github-source-dependencies.json
    python3 code/check_component_catalogue.py --output work/github-component-catalogue.json
    python3 code/check_global_structure.py --output work/github-global-structure.json
    bash scripts/build_structural_manuscripts.sh poisson_divisor_symmetries_ja poisson_p4_classification_ja normal_actual_divisor_recognition_ja s61_first_jet_reduction_ja
    bash scripts/reproduce_subvariety_incidence.sh
    bash scripts/verify_isotropic_tangent_cone.sh
    ;;
  *)
    echo 'Usage: bash scripts/verify_github_source.sh [tangent-cone|incidence|full]' >&2
    exit 2
    ;;
esac
