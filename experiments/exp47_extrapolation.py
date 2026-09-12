"""Round 41: how far away the wall actually is.

Forty rounds have said "sub-exponential, not polynomial". This one puts a number
on it. The quadratic sieve and the toy number field sieve in this repository are
timed across the sizes they can reach, their `L[1/2]` and `L[1/3]` constants are
fitted to those times, and the fits are extrapolated to the sizes that matter.

The fits turn out to be worth nothing outside their range -- the toy sieve's
fitted constant is half the textbook one, and extrapolating it "predicts"
RSA-768 in about two core-years -- so the round throws them away and calibrates
on the published records instead. That extrapolation says something sharper
than "far away": 1024 bits is close enough that constants still decide it, and
2048 bits is nine orders of magnitude past that. One doubling, nine orders of
magnitude, and no exponent moved.
"""

import random
import statistics
import time
from math import exp, log, sqrt

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.qs import quadratic_sieve
from aksfactor.nfs import number_field_sieve

report = Report(
    "exp47_extrapolation",
    "How far the wall is, in years",
    "Round 41. Fitting this repository's sieves and extrapolating them, "
    "calibrated against the published records.",
)

rng = random.Random(47)
SECONDS_PER_YEAR = 365.25 * 24 * 3600


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


def l_notation(bits, alpha, c):
    """`exp(c (ln N)^alpha (ln ln N)^(1-alpha))`, the cost shape of the sieves."""
    ln_n = bits * log(2)
    return exp(c * ln_n ** alpha * log(ln_n) ** (1 - alpha))


def shape_of(bits, alpha):
    ln_n = bits * log(2)
    return ln_n ** alpha * log(ln_n) ** (1 - alpha)


def fit_line(points, alpha):
    """Least squares for `log t = intercept + c * shape`.

    Forcing the line through the origin instead -- one parameter, `log t = c *
    shape` -- makes `c` absorb the implementation's constant factor, which is
    exactly the quantity that must be kept separate from the exponent.
    """
    xs = [shape_of(bits, alpha) for bits, _ in points]
    ys = [log(seconds) for _, seconds in points]
    mx, my = statistics.mean(xs), statistics.mean(ys)
    denom = sum((x - mx) ** 2 for x in xs)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denom
    return slope, my - slope * mx


def predict(bits, alpha, slope, intercept):
    return exp(intercept + slope * shape_of(bits, alpha))


# ---------------------------------------------------------------- section 1
report.p("## What this repository's sieves cost")
report.p()
report.p("Balanced semiprimes, median of three, on one core of whatever machine "
         "this ran on. The quadratic sieve is the one that reaches furthest here; "
         "the number field sieve is a toy at degree 3 and loses to it at every "
         "size it can manage, which is a fact about the implementation and not "
         "about the algorithms.")
report.p()

qs_points, nfs_points = [], []
rows = []
for bits in (60, 70, 80, 90, 100):
    times = []
    for _ in range(3):
        p, q = randprime(bits // 2), randprime(bits - bits // 2)
        n = p * q
        t0 = time.perf_counter()
        got = quadratic_sieve(n)
        dt = time.perf_counter() - t0
        assert got in (p, q), f"quadratic sieve failed at {bits} bits"
        times.append(dt)
    med = statistics.median(times)
    qs_points.append((bits, med))
    rows.append([bits, f"{med:.2f} s"])
report.table(["bits of N", "quadratic sieve (median of 3)"], rows)

rows = []
for bits in (40, 50, 60, 70):
    times, ok = [], 0
    for _ in range(2):
        p, q = randprime(bits // 2), randprime(bits - bits // 2)
        n = p * q
        t0 = time.perf_counter()
        got = number_field_sieve(n, 3, 1500, 1500, 24, 8000, 3000)
        dt = time.perf_counter() - t0
        ok += got in (p, q)
        times.append(dt)
    med = statistics.median(times)
    if ok == 2:
        nfs_points.append((bits, med))
    rows.append([bits, f"{med:.2f} s", f"{ok}/2"])
report.table(["bits of N", "number field sieve (median of 2)", "factored"], rows)
report.p()

# ---------------------------------------------------------------- section 2
report.p("## Fitting the shape")
report.p()
report.p("The model is `t = A exp(c (ln N)^alpha (ln ln N)^(1-alpha))` with "
         "`alpha` held at the theoretical value, fitted as a straight line in "
         "the log. Both parameters are needed: forced through the origin, `c` "
         "silently absorbs the implementation's constant factor `A`, which is "
         "the one quantity this round has to keep separate.")
report.p()

qs_c, qs_a = fit_line(qs_points, 0.5)
nfs_fit = fit_line(nfs_points, 1 / 3) if len(nfs_points) >= 3 else None
rows = [["quadratic sieve", "1/2", f"{qs_c:.3f}", f"{exp(qs_a):.2e}", "1"]]
if nfs_fit:
    rows.append(["number field sieve (toy, degree 3)", "1/3", f"{nfs_fit[0]:.3f}",
                 f"{exp(nfs_fit[1]):.2e}", "1.923"])
report.table(["method", "alpha", "fitted c", "fitted A (seconds)",
              "textbook c"], rows)

rows = []
for bits, seconds in qs_points:
    got = predict(bits, 0.5, qs_c, qs_a)
    rows.append([bits, f"{seconds:.2f} s", f"{got:.2f} s", f"{got / seconds:.2f}"])
report.table(["bits of N", "measured", "fitted", "ratio"], rows)
resid = max(abs(log(predict(b, 0.5, qs_c, qs_a) / t)) for b, t in qs_points)
report.p(f"The fit's worst residual inside the measured range is a factor of "
         f"{exp(resid):.2f}. A curve that is already off by that much where it "
         f"was fitted cannot be read to two decimal places outside it, so what "
         f"follows uses it only for orders of magnitude -- and then stops using "
         f"it altogether.")
report.p()

# ---------------------------------------------------------------- section 3
report.p("## Extrapolating this implementation")
report.p()
rows = []
for bits in (128, 256, 512, 1024):
    qs_years = predict(bits, 0.5, qs_c, qs_a) / SECONDS_PER_YEAR
    nfs_years = (predict(bits, 1 / 3, *nfs_fit) / SECONDS_PER_YEAR
                 if nfs_fit else None)
    rows.append([bits, f"{qs_years:.2e}",
                 f"{nfs_years:.2e}" if nfs_years else "--"])
report.table(["bits of N", "this repository's QS (years, one core)",
              "this repository's toy NFS (years, one core)"], rows)
report.p(f"Those numbers are worthless and it is worth saying why, because the "
         f"reason is the whole difficulty with extrapolation. The toy sieve's "
         f"fitted `c` is {nfs_fit[0]:.3f} against a textbook 1.923"
         + (f", so the curve bends far too gently and 'predicts' RSA-768 in "
            f"{predict(768, 1 / 3, *nfs_fit) / SECONDS_PER_YEAR:.1f} core-years"
            if nfs_fit else "")
         + ". A constant fitted over 30 bits of range does not survive 700 bits "
           "of extrapolation. The next section discards these fits entirely.")
report.p()

# ---------------------------------------------------------------- section 4
report.p("## Calibrating on the records instead")
report.p()
report.p("The published factorisations are the honest data for an extrapolation "
         "to real key sizes: a serious number field sieve, on hardware that "
         "existed, with the effort reported. Fitting `A` to each record with `c` "
         "held at the textbook 1.923 gives a constant per record, and the useful "
         "check is whether the two records agree on it -- they are 61 bits and "
         "ten years apart, so if the `L[1/3]` shape were wrong between them, "
         "they would not.")
report.p()

records = [("RSA-768 (2009)", 768, 2000.0), ("RSA-250 (2020)", 829, 2700.0)]
constants = []
rows = []
for name, bits, core_years in records:
    a = log(core_years * SECONDS_PER_YEAR) - 1.923 * shape_of(bits, 1 / 3)
    constants.append(a)
    rows.append([name, bits, f"{core_years:,.0f}", f"{exp(a):.2e}"])
report.table(["record", "bits", "published core-years",
              "implied A (seconds) at c = 1.923"], rows)
spread = exp(max(constants) - min(constants))
report.p(f"The two implied constants differ by a factor of {spread:.2f}. Ten "
         f"years of hardware and software progress, and 61 bits, and the "
         f"`L[1/3]` shape absorbs the difference to within that -- which is the "
         f"evidence that extrapolating along this curve is a reasonable thing "
         f"to do at all.")
report.p()

calibrated = statistics.mean(constants)
rows = []
for bits in (1024, 2048, 4096):
    years = exp(calibrated + 1.923 * shape_of(bits, 1 / 3)) / SECONDS_PER_YEAR
    rows.append([bits, f"{years:.2e}", f"{years / 1e6:.2e}",
                 f"{years / 1e6 / 1e6:.2e}"])
report.table(["bits of N", "core-years at the records' quality",
              "years on a million cores",
              "years on a million cores, after a million-fold speedup"], rows)

report.p("So one doubling of the key, from 1024 to 2048 bits, costs nine orders "
         "of magnitude. That is the shape of `L[1/3]`, and it is the whole "
         "security argument for RSA -- not that factoring is hard, but that it "
         "gets hard faster than keys get long.")
report.p()
report.p("It also says where constants still matter. At 1024 bits the total is "
         "within reach of a large enough machine for long enough, which is why "
         "1024-bit RSA is deprecated rather than merely discouraged; a "
         "million-fold improvement in constants would finish it outright. The "
         "same million-fold leaves 2048 bits at a billion core-years. Near the "
         "edge, constants decide; past it, only the exponent does.")
report.p()
report.p("Which is why the search in this project was always for a change of "
         "exponent. Rounds 33 to 36 found four mechanisms and one exponent each. "
         "Rounds 37 to 40 found three kinds of extra information -- bits shared "
         "between moduli, batches of instances, a predicate that splits the "
         "range -- and priced every one of them above the wall. The table above "
         "is the size of the prize for moving it.")
report.write()
