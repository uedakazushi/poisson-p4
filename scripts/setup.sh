#!/usr/bin/env bash
# Install only into this checkout; no sudo, OS changes, or automatic full audit.
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
run_verification=0

usage() {
  cat <<'USAGE'
Usage: bash scripts/setup.sh [--verify]

Create/reuse .venv, install the pinned Python dependencies, and materialize
the immutable research sources. --verify also runs the v24 quick suite.
PYTHON may select a Python 3.11+ interpreter; it defaults to python3.
LaTeX is optional and is not installed by this script.
USAGE
}

for argument in "$@"; do
  case "$argument" in
    --verify) run_verification=1 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown option: %s\n' "$argument" >&2; usage >&2; exit 2 ;;
  esac
done

setup_python="${PYTHON:-python3}"
"$setup_python" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 11) else "Python 3.11 or newer is required (3.12 is recommended).")'

if [[ ! -x "$repo_root/.venv/bin/python" ]]; then
  if [[ -e "$repo_root/.venv" ]]; then
    printf '%s\n' 'An incomplete .venv already exists. Inspect or remove it before retrying setup.' >&2
    exit 1
  fi
  "$setup_python" -m venv "$repo_root/.venv"
fi

venv_python="$repo_root/.venv/bin/python"
"$venv_python" -m pip install --requirement "$repo_root/requirements.txt"
"$venv_python" "$repo_root/scripts/bootstrap_sources.py"
"$venv_python" "$repo_root/scripts/doctor.py" --strict

if [[ "$run_verification" -eq 1 ]]; then
  "$venv_python" "$repo_root/scripts/verify.py" --suite quick
fi

printf '\nReady. Activate with: source %s/.venv/bin/activate\n' "$repo_root"
printf 'Run the v24 audit with: %s scripts/verify.py --suite quick\n' "$venv_python"
