"""Round 23: Pascal's triangle over divisibility sequences, and the one that is ECM.

The first rounds of this project looked at Pascal's triangle mod n and found a
Sierpinski gasket whose cells are exactly p wide -- which is why reading p off
it costs p.  Replace n by the n-th term of a divisibility sequence and the
triangle survives, integers and all; only the cell width changes.  One choice
of sequence makes the cell width redrawable, and that choice is Lenstra's
elliptic curve method.
"""

import random
from math import gcd

from _common import Report

from aksfactor.arith import factorize, is_prime, sieve
from aksfactor.divseq import (
    curve_eds,
    elliptic_triangle_factor,
    factor_with_point_count,
    fibonacci,
    gasket,
    mixed_radix_carries,
    naturals,
    nomial,
    point_order,
    q_integers,
    rank_of_apparition,
)

report = Report(
    "exp31_elliptic_pascal",
    "Pascal's triangle over divisibility sequences, and the one that is ECM",
    "Round 23. Four integer triangles, one Sierpinski rule, and a cell width you can redraw.",
)

rng = random.Random(31)
CURVE = (0, 17, 2, 5)          # y^2 = x^3 + 17, P = (2, 5), a point of infinite order
families = {
    "Pascal, a_n = n": naturals(90),
    "Gaussian, a_n = (2^n - 1)": q_integers(2, 90),
    "Fibonomial, a_n = F_n": fibonacci(90),
    "elliptic, a_n = psi_n(P) on y^2 = x^3 + 17": curve_eds(*CURVE, 90),
}

report.p("## 1. Four integer triangles")
report.p()
ell = families["elliptic, a_n = psi_n(P) on y^2 = x^3 + 17"]
report.p("For a sequence with `gcd(a_m, a_n) = a_gcd(m,n)`, the generalised binomial "
         "`[n, k]_a = a_n a_(n-1) .. a_(n-k+1) / (a_1 .. a_k)` is an integer. The "
         "fourth sequence is an *elliptic divisibility sequence*: the division "
         "polynomials of a point `P` of infinite order, "
         f"`{', '.join(str(v) for v in ell[1:5])}, ...`.")
report.p()
rows = []
for name, seq in families.items():
    total = ok = 0
    for n in range(40):
        for k in range(n + 1):
            total += 1
            try:
                nomial(seq, n, k)
                ok += 1
            except ArithmeticError:
                pass
    rows.append([name, f"{ok} of {total}"])
report.table(["triangle", "integral entries, rows 0-39"], rows)

report.p("## 2. One Sierpinski rule")
report.p()
report.p("Modulo a prime `p`, `[n, k]_a` vanishes exactly when adding `k` and `n - k` "
         "carries in the mixed radix `(r, p, p, ...)`, where `r = r(p)` is the rank "
         "of apparition -- the first `n` with `p | a_n` (Kummer's theorem, as "
         "generalised by Knuth and Wilf). Each triangle mod `p` is a gasket whose "
         "first-level cell is `r` wide.")
report.p()
rows = []
primes = [p for p in sieve(50) if p >= 5]
for name, seq in families.items():
    total = ok = used = 0
    for p in primes:
        if seq[2] % p == 0:
            continue
        r = rank_of_apparition(seq, p)
        if r is None or r > 40:
            continue
        used += 1
        for n in range(45):
            for k in range(n + 1):
                total += 1
                ok += (nomial(seq, n, k) % p == 0) == (mixed_radix_carries(n, k, r, p) > 0)
    rows.append([name, used, f"{ok} of {total}"])
report.table(["triangle", "primes 5..47 used", "entries where the carry rule is right"], rows)

report.p("The same prime, `p = 7`, four cells:")
report.p()
for name, seq in families.items():
    r = rank_of_apparition(seq, 7)
    report.p(f"**{name}** mod 7, cell `r = {r}`:")
    report.p()
    report.p("```")
    for line in gasket(seq, 7, 30):
        report.p(line)
    report.p("```")
    report.p()

report.p("## 3. The cell is an order")
report.p()
rows = []
checked = {"Pascal": 0, "Gaussian": 0, "Fibonomial": 0, "elliptic": 0}
total = 0
for p in [p for p in sieve(400) if p >= 7 and (27 * 17 * 17) % p]:
    total += 1
    checked["Pascal"] += rank_of_apparition(naturals(p + 1), p) == p
    r2 = rank_of_apparition(q_integers(2, p), p)
    checked["Gaussian"] += r2 == next(d for d in range(1, p) if pow(2, d, p) == 1)
    rf = rank_of_apparition(fibonacci(2 * p + 2), p)
    checked["Fibonomial"] += (p - (1 if pow(5, (p - 1) // 2, p) == 1 else -1)) % rf == 0
    re = rank_of_apparition(curve_eds(*CURVE, 2 * p + 4), p) if p != 2 else None
    checked["elliptic"] += re == point_order(*CURVE, p)
report.table(["triangle", "cell r(p)", "checked on primes 7..397"], [
    ["Pascal", "`p` itself", f"{checked['Pascal']} of {total}"],
    ["Gaussian (q = 2)", "the order of 2 mod p, a divisor of `p - 1`", f"{checked['Gaussian']} of {total}"],
    ["Fibonomial", "a divisor of `p - (5/p)`", f"{checked['Fibonomial']} of {total}"],
    ["elliptic", "the order of `P` in `E(F_p)`, anywhere in `p + 1 +- 2 sqrt p`",
     f"{checked['elliptic']} of {total}"],
])
report.p("The Pascal triangle's cell *is* `p`: to see the cell is to have factored. "
         "The Gaussian and Fibonomial cells are orders, but orders fixed by `p` -- "
         "divisors of `p - 1` or `p + 1`, the one ticket of Proposition 14. The "
         "elliptic cell is the order of a point, and a different curve is a "
         "different triangle.")
report.p()

report.p("## 4. Redrawing the triangle")
report.p()
p = 10007
cells = []
for _ in range(300):
    while True:
        a, x, y = rng.randrange(p), rng.randrange(p), rng.randrange(1, p)
        b = (y * y - x ** 3 - a * x) % p
        if (4 * a ** 3 + 27 * b * b) % p:
            break
    cells.append(point_order(a, b, x, y, p))
B = 50


def smooth(m):
    return max(factorize(m)) <= B if m > 1 else True


def show(m):
    return " x ".join(f"{a}^{e}" if e > 1 else str(a) for a, e in sorted(factorize(m).items()))


eps = 1 if pow(5, (p - 1) // 2, p) == 1 else -1
fib_m = p - eps
report.table(["triangle family", "possible cells at p = 10007", f"a {B}-smooth cell available?"], [
    ["Gaussian, any q", f"divisors of p - 1 = {show(p - 1)}", "yes" if smooth(p - 1) else "no"],
    ["Fibonomial", f"divisors of p - (5/p) = {show(fib_m)}", "yes" if smooth(fib_m) else
     f"only small divisors ({max(factorize(fib_m))} blocks the rest)"],
    ["elliptic, 300 random curves", f"{len(set(cells))} distinct orders in [{min(cells)}, {max(cells)}]",
     f"{sum(smooth(c) for c in cells)} of the 300 curves"],
])
big = max(factorize(p - 1))
report.p(f"`p - 1` has the prime factor {big}, so every Gaussian triangle mod {p} "
         f"has a cell that is a multiple of {big} or a divisor of "
         f"{(p - 1) // big}; the Fibonomial cells are stuck with the divisors of "
         f"{fib_m}. The elliptic triangles take {len(set(cells))} different cell "
         f"widths across the Hasse interval, and {sum(smooth(c) for c in cells)} of "
         f"them are `{B}`-smooth.")
report.p()

report.p("How many different cells can each family offer at one prime? Sampling "
         "many members of each family at `p = 10007`:")
report.p()


def curve_order(a, b, p):
    """#E(F_p) for y^2 = x^3 + a x + b by counting."""
    total = 1
    for x in range(p):
        r = (x ** 3 + a * x + b) % p
        total += 1 if r == 0 else (2 if pow(r, (p - 1) // 2, p) == 1 else 0)
    return total


pc = 10009 if 10009 % 3 == 1 else 10039          # a prime = 1 (mod 3) for the CM family
gauss_cells = {next(d for d in range(1, p) if pow(qq, d, p) == 1) for qq in range(2, 400)}
cm_orders = {curve_order(0, b, pc) for b in range(1, 200) if (27 * b * b) % pc}
generic_orders = set()
for _ in range(60):
    while True:
        a, b = rng.randrange(p), rng.randrange(p)
        if (4 * a ** 3 + 27 * b * b) % p:
            break
    generic_orders.add(curve_order(a, b, p))
report.table(["family", "members sampled", "distinct group orders (cells divide these)"], [
    ["Pascal", "--", "1: the cell is `p`"],
    ["Gaussian, q = 2..399", 398, f"{len(gauss_cells)} cells, all dividing `p - 1`: one order"],
    ["Fibonomial / Lucas", "--", "2: `p - 1` or `p + 1`"],
    [f"elliptic with CM, `y^2 = x^3 + b` at p = {pc}", 199, f"{len(cm_orders)}"],
    ["elliptic, generic", 60, f"{len(generic_orders)}"],
])
report.p(f"The family with complex multiplication offers only {len(cm_orders)} orders "
         f"-- the six twists of a curve with `j = 0` -- so it is barely more "
         f"redrawable than `p +- 1`. Generic curves offer a new order almost every "
         f"time: the Hasse interval holds about `4 sqrt p` of them. That count of "
         f"tickets per prime is the whole difference between the rigid order methods "
         f"and ECM.")
report.p()

report.p("## 5. Jumping to a row: ECM")
report.p()
report.p("Row `M` of an elliptic triangle mod `N` begins with `a_M = W_M`, and `p | W_M` "
         "exactly when the cell `r_p` divides `M`. Take `M = lcm(1..B)` -- reached in "
         "`O(log M)` steps by the sequence's own double-and-add -- and `gcd(W_M, N)` "
         "splits `N` whenever `r_p` is `B`-smooth. If not, draw another curve.")
report.p()
rows = []
for pbits in (20, 24, 28, 32):
    curves = []
    for _ in range(3):
        while True:
            pp = rng.randrange(1 << (pbits - 1), 1 << pbits) | 1
            if is_prime(pp):
                break
        while True:
            qq = rng.randrange(1 << 63, 1 << 64) | 1
            if is_prime(qq):
                break
        st = {}
        g = elliptic_triangle_factor(pp * qq, 2000, 2000, rng, st)
        assert g in (pp, qq)
        curves.append(st["curves"])
    rows.append([pbits, 64 + pbits, " / ".join(map(str, curves))])
report.table(["bits of p", "bits of N", "curves needed (three semiprimes)"], rows)
report.p("That is stage 1 of Lenstra's elliptic curve method, read as the original "
         "observation of this project: *look along a row of the triangle mod `N` for "
         "entries that share a factor with `N`.* In Pascal's triangle the row that "
         "works is `N` itself and the entries that work sit at multiples of `p`, so "
         "finding them is finding `p`. In the elliptic triangle the row that works "
         "is any multiple of the cell width, and the cell width can be redrawn until "
         "it is smooth. The gasket the project started from was the right picture; "
         "it needed a different sequence under it.")
report.p()

report.p("## 6. The cell widths, counted: equivalent to factoring")
report.p()
report.p("The elliptic cells `r_p`, `r_q` divide `#E(F_p)` and `#E(F_q)`, and "
         "`#E(Z/N) = #E(F_p) #E(F_q)`. An oracle for that one number factors `N`: "
         "`[#E(Z/N)] P = O` modulo both primes, and removing a small prime `l` that "
         "divides only one of the two orders leaves a multiple that vanishes modulo "
         "one prime and not the other.")
report.p()
ok = total = 0
used = []
for _ in range(40):
    while True:
        pp = rng.randrange(2000, 8000) | 1
        qq = rng.randrange(2000, 8000) | 1
        if pp != qq and is_prime(pp) and is_prime(qq):
            break
    nn = pp * qq
    while True:
        a, x, y = rng.randrange(nn), rng.randrange(nn), rng.randrange(nn)
        b = (y * y - x ** 3 - a * x) % nn
        if gcd((4 * a ** 3 + 27 * b * b) % nn, nn) == 1:
            break
    count = curve_order(a % pp, b % pp, pp) * curve_order(a % qq, b % qq, qq)   # the oracle
    g, ell = factor_with_point_count(nn, a, (x, y), count)
    total += 1
    ok += g in (pp, qq)
    if ell:
        used.append(ell)
report.p(f"With the count supplied by brute force as the oracle, one count per "
         f"number factored **{ok} of {total}** semiprimes; the prime removed was "
         f"at most {max(used)}. The reverse direction -- factoring gives the count -- "
         f"is Schoof's algorithm modulo each prime, so counting points on `E(Z/N)` is "
         f"equivalent to factoring `N` (Kunihiro and Koyama, 1998). It joins "
         f"`phi(N)`, `tau(N)` and `r_4(N)` on the list of quantities that carry the "
         f"factorisation openly and are polynomial-time equivalent to it.")
report.write()
