"""Regenerate the row-support figure used in the README."""

from _common import ROOT, RESULTS  # noqa: F401  (sets sys.path)

from aksfactor.arith import is_prime, spf_trial
from aksfactor.pascal import row_exact

UPTO = 45

lines = [
    "  n │ interior of row n mod n   (· = 0,  █ = non-zero,  ▲ = first non-zero)",
    "────┼" + "─" * 56,
]
for n in range(4, UPTO):
    row = row_exact(n, n)
    cells, first = [], None
    for k in range(1, n):
        if row[k] and first is None:
            first = k
        cells.append("·" if not row[k] else ("▲" if k == first else "█"))
    tag = "prime" if is_prime(n) else f"spf={spf_trial(n)}  n/spf={n // spf_trial(n)}"
    lines.append(f"{n:3d} │ {''.join(cells):<44s} {tag}")

out = ROOT / "docs" / "row_map.txt"
out.write_text("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"\n-> wrote {out.relative_to(ROOT)}")
