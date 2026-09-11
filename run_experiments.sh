#!/usr/bin/env bash
# Regenerate everything in experiments/results/.
set -e
cd "$(dirname "$0")/experiments"
python3 render_figure.py > /dev/null && echo "=== render_figure.py ==="
python3 render_gaskets.py > /dev/null && echo "=== render_gaskets.py ==="
for f in exp*.py; do
  echo "=== $f ==="
  python3 "$f"
done
