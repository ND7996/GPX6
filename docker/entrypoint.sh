#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "check" ]]; then
  echo "Python: $(python --version 2>&1)"
  python - <<'PY'
import importlib

modules = [
    "numpy",
    "matplotlib",
    "Bio",
    "pandas",
    "scipy",
    "seaborn",
    "networkx",
    "PIL",
    "adjustText",
]

for module in modules:
    importlib.import_module(module)
    print(f"ok: {module}")

try:
    import Qpyl
    print("ok: Qpyl")
except Exception as exc:
    print(f"warn: Qpyl import failed: {exc}")
PY

  for exe in pymol q_mapper.py q_analysefeps.py qprep5 qdyn5_r8 qfep5 latexmk pdflatex sbatch mpirun; do
    if command -v "$exe" >/dev/null 2>&1; then
      echo "ok: $exe -> $(command -v "$exe")"
    else
      echo "missing: $exe"
    fi
  done

  echo
  echo "Q6 binaries are expected to be supplied by mounting them into /opt/q6/bin or adding them to PATH."
  exit 0
fi

exec "$@"
