#!/usr/bin/env python3
"""Zero-dependency test runner: ``python3 run_tests.py`` (pytest also works)."""

from __future__ import annotations

import importlib.util
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))


def load(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    files = sorted((ROOT / "tests").glob("test_*.py"))
    passed = failed = 0
    for path in files:
        module = load(path)
        for name in sorted(vars(module)):
            if not name.startswith("test_"):
                continue
            fn = getattr(module, name)
            if not callable(fn):
                continue
            start = time.time()
            try:
                fn()
            except Exception:
                failed += 1
                print(f"FAIL  {path.name}::{name}", flush=True)
                traceback.print_exc()
            else:
                passed += 1
                print(f"ok    {path.name}::{name}  ({time.time() - start:.2f}s)", flush=True)
    print(f"\n{passed} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
