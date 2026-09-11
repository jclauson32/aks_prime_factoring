"""Regenerate docs/figures/*.png: Pascal's and an elliptic triangle mod 77.

Each entry is coloured by which of 7 and 11 divide it: red for 7 only, blue for
11 only, dark for both, light for neither.  The coloured entries are exactly
those whose gcd with 77 is a proper factor.
"""

from math import comb

from _common import ROOT

from aksfactor.divseq import curve_eds, eds_mod, nomial, rank_of_apparition
from aksfactor.figures import classify, write_png, zero_pattern

P, Q, ROWS, SCALE = 7, 11, 192, 3
out = ROOT / "docs" / "figures"
out.mkdir(parents=True, exist_ok=True)

naturals = list(range(ROWS + 1))
zp, zq = zero_pattern(naturals, P, 5, ROWS), zero_pattern(naturals, Q, 4, ROWS)
assert all((comb(n, k) % P == 0) == zp[n][k] for n in range(80) for k in range(n + 1))
write_png(out / "pascal_mod_77.png", classify(zp, zq, ROWS, SCALE))

curve = (0, 17, 2, 5)                    # y^2 = x^3 + 17, P = (2, 5)
w = curve_eds(*curve, 4)
ezp = zero_pattern(eds_mod(*w[1:5], ROWS + 1, P ** 6), P, 6, ROWS)
ezq = zero_pattern(eds_mod(*w[1:5], ROWS + 1, Q ** 6), Q, 6, ROWS)
exact = curve_eds(*curve, 50)
assert all((nomial(exact, n, k) % P == 0) == ezp[n][k] for n in range(40) for k in range(n + 1))
assert all((nomial(exact, n, k) % Q == 0) == ezq[n][k] for n in range(40) for k in range(n + 1))
write_png(out / "elliptic_mod_77.png", classify(ezp, ezq, ROWS, SCALE))
print(f"cells: Pascal {P}, {Q}; elliptic {rank_of_apparition(exact, P)}, "
      f"{rank_of_apparition(exact, Q)}")
print(f"-> wrote {out.relative_to(ROOT)}/pascal_mod_77.png, elliptic_mod_77.png")
