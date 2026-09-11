"""Round 25: can Shor's interference be done classically?

Shor's algorithm reads the order r from a heavy Fourier coefficient.  Classical
algorithms can also find heavy Fourier coefficients -- Goldreich-Levin and
Kushilevitz-Mansour do it with about 1/tau queries for weight tau -- *if* some
efficiently computable function of a^x mod N puts heavy weight on a frequency
that reveals r.  This round measures the spectra of every such function this
repository knows how to compute.
"""

import cmath
import math
import random
import statistics
from math import gcd

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.central import jacobi
from aksfactor.shor import order_spectrum

report = Report(
    "exp33_dequantize",
    "Can Shor's interference be done classically?",
    "Round 25. The Fourier spectra of computable functions of a^x mod N.",
)

rng = random.Random(33)


def randprime(lo, hi):
    while True:
        v = rng.randrange(lo, hi) | 1
        if is_prime(v):
            return v


def order(a, n):
    r, v = 1, a % n
    while v != 1:
        v = v * a % n
        r += 1
    return r


def tests(n, c):
    return {
        "additive character e(cv/N)": lambda v: cmath.exp(2j * math.pi * c * v / n),
        "e(c v^2 / N)": lambda v: cmath.exp(2j * math.pi * c * v * v / n),
        "v mod 7 a square mod 7": lambda v: 1 if (v % 7) in (0, 1, 2, 4) else -1,
        "lowest bit of v": lambda v: 1 - 2 * (v & 1),
        "v < N/2": lambda v: 1 if v < n // 2 else -1,
        "middle bit of v": lambda v: 1 - 2 * ((v >> (n.bit_length() // 2)) & 1),
        "Jacobi symbol (v/N)": lambda v: jacobi(v, n),
    }


report.p("For `g(s) = psi(a^s mod N)` over one period `s < r`, the weight "
         "`|g^(j)|^2` at frequency `j` is what Shor's first register would show at "
         "`y ~ jQ/r`. A frequency `j` coprime to `r` reveals `r` by continued "
         "fractions. The table reports `r x` the largest such weight: `1` means "
         "perfectly flat, and a random spectrum gives about `ln r`.")
report.p()
bands = [(200, 700), (1000, 3000), (4000, 12000), (15000, 40000)]
rows = []
summary = {}
all_orders = []
for lo, hi in bands:
    got = {}
    rs = []
    for _ in range(3):
        while True:
            p = randprime(40, 400)
            q = randprime(40, 400)
            if p == q:
                continue
            n = p * q
            a = rng.randrange(2, n)
            if gcd(a, n) > 1:
                continue
            r = order(a, n)
            if lo <= r <= hi:
                break
        rs.append(r)
        all_orders.append(r)
        c = rng.randrange(1, n)
        for name, psi in tests(n, c).items():
            r, w = order_spectrum(n, a, psi, r)
            best = max((w[j] for j in range(1, r) if gcd(j, r) == 1), default=0.0)
            got.setdefault(name, []).append(best * r)
    med_r = statistics.median(rs)
    summary[(lo, hi)] = got
    rows.append([f"{min(rs)}-{max(rs)}", f"{math.log(med_r):.1f}"] +
                [f"{statistics.median(v):.2f}" for v in got.values()])
names = list(tests(15, 1).keys())
report.table(["r", "ln r"] + names, rows)

meds = [statistics.median(v) for band in summary.values() for k, v in band.items()
        if not k.startswith("Jacobi") and not k.startswith("e(c v^2")]
worst = max(max(v) for band in summary.values() for k, v in band.items()
            if not k.startswith("Jacobi"))
report.p(f"Every function except the Jacobi symbol is flat: the heaviest coefficient "
         f"that reveals `r` carries between {min(meds):.1f} and {max(meds):.1f} times "
         f"the uniform weight `1/r` (median over three orders per band; the single "
         f"largest seen is {worst:.1f}/r), the size of the maximum of a random "
         f"spectrum. A classical heavy-coefficient search at weight `tau ~ c/r` needs "
         f"on the order of `r` queries -- worse than the `sqrt r` of baby-step "
         f"giant-step on the order, and `r` is about `N`. The quadratic character's "
         f"medians are zero because {sum(r % 2 == 0 for r in all_orders)} of the "
         f"{len(all_orders)} orders drawn were even, and then `a^(2s)` has period "
         f"`r/2`: it can only reveal a divisor.")
report.p()
report.p("The Jacobi symbol is the one computable function with a heavy "
         "coefficient -- `(a^s / N) = (a/N)^s`, so all its weight sits at `j = 0` or "
         "`j = r/2` -- and that coefficient says only whether `r` is even. It is also, "
         "as far as anyone knows, the only multiplicative character of `(Z/N)*` "
         "computable without the factorisation; a character of order `m` would "
         "reveal `r mod m` the same way, and no more.")
report.p()
coset = []
for _ in range(3):
    while True:
        p, q = randprime(40, 400), randprime(40, 400)
        n = p * q
        a = rng.randrange(2, n)
        if p != q and gcd(a, n) == 1 and 1000 <= order(a, n) <= 3000:
            break
    r = order(a, n)
    x0 = rng.randrange(r)
    target = pow(a, x0, n)
    # the collapsed second register: the indicator of {x : a^x = a^x0}, over x < 8r
    ind = [1.0 if pow(a, x, n) == target else 0.0 for x in range(8 * r)]
    from aksfactor.shor import dft_any

    spec = dft_any(ind)
    energy = sum(abs(z) ** 2 for z in spec)
    peak = max(abs(spec[j]) ** 2 for j in range(1, len(spec)) if j % 8 == 0 and gcd(j // 8, r) == 1)
    coset.append((r, peak / energy * r, sum(ind) / len(ind) * r))
report.p("Shor's own register is no exception. After the second register collapses, "
         "the first holds the indicator of one coset `{x : a^x = c}` -- a function "
         "that *is* classically computable. Its spectrum:")
report.p()
report.table(["r", "r x largest revealing weight / total", "r x density of the indicator"],
             [[r, f"{w:.2f}", f"{d:.2f}"] for r, w, d in coset])
report.p("Perfectly flat over its `r` peaks, each carrying `1/r` of the energy, and "
         "supported on a `1/r` fraction of the inputs. What the quantum measurement "
         "supplies is not a heavy coefficient -- there is none, in this function or "
         "any other tested -- but a *sample* from a distribution spread over `r` "
         "frequencies, any one of which reveals `r`. Classical heavy-coefficient "
         "search cannot draw that sample; it can only find coefficients of weight "
         "`1/poly`, and none of the computable functions of `a^x mod N` tested here "
         "has one.")
report.write()
