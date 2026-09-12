#!/usr/bin/env bash
# Regenerate everything in experiments/results/.
#
# Some experiments report wall-clock times (exp06, exp29, exp30, exp32, exp34,
# exp35, exp41, exp43, exp44, exp45, exp46, exp47). Run those on an otherwise
# idle machine: under load their tables measure the load, not the methods.
set -e
cd "$(dirname "$0")/experiments"
python3 render_figure.py > /dev/null && echo "=== render_figure.py ==="
python3 render_gaskets.py > /dev/null && echo "=== render_gaskets.py ==="
for f in exp*.py; do
  echo "=== $f ==="
  python3 "$f"
done
