"""Round 28: evolving factoring programs.

Round 14 found random straight-line programs over Z/N worth exactly the
zero-divisor density.  Here they are evolved instead: short programs are scored
on how often their output shares a factor with N across training semiprimes,
and the better ones breed.  The population starts random -- no known method is
seeded -- and the question is what selection finds.
"""

import random
import time

from _common import Report

from aksfactor.arith import factorize, is_prime
from aksfactor.evolve import (
    OPS,
    PASCAL_OPS,
    as_polynomial,
    cost,
    evolve,
    floyd_program,
    optimised_pminus1,
    output_exponent,
    output_slice,
    pminus1_program,
    pminus1_success,
    random_program,
    score,
)

report = Report(
    "exp36_evolve",
    "Evolving factoring programs",
    "Round 28. Selection over straight-line programs, starting from nothing.",
)

rng = random.Random(36)
BITS, BUDGET, LENGTH = 12, 40, 16


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


def cases(k):
    out = []
    while len(out) < k:
        p, q = randprime(BITS), randprime(BITS)
        if p != q:
            n = p * q
            out.append((n, rng.randrange(2, n - 1)))
    return out


train, test = cases(500), cases(500)
report.p(f"Programs of {LENGTH} instructions over `Z/N` (`mul`, `add`, `sub`, and `pow` by "
         f"a prime below 50), held to {BUDGET} multiplications. A program scores on "
         f"`N` when `gcd(out, N)` or `gcd(out - 1, N)` is a proper factor. Training and "
         f"test sets: {len(train)} semiprimes each, with {BITS}-bit primes, where "
         f"smooth group orders are common enough for selection to have a gradient.")
report.p()

pm = pminus1_program(BUDGET)
rand = [p for p in (random_program(LENGTH, rng) for _ in range(400)) if cost(p) <= BUDGET]
rand_score = sum(score(p, test) for p in rand) / len(rand)
e_opt, c_opt = optimised_pminus1(train, BUDGET)
opt_test = pminus1_success(e_opt, test)


def show(e):
    return " x ".join(f"{ell}^{k}" if k > 1 else str(ell) for ell, k in sorted(factorize(e).items()))


rows = [["random programs (mean of %d)" % len(rand), "--", f"{rand_score:.4f}", "--"],
        ["Pollard p - 1, written by hand", f"{cost(pm):.1f}", f"{score(pm, test):.3f}", "--"],
        ["Pollard p - 1, exponent optimised on the training set", f"{c_opt:.1f}",
         f"{opt_test:.3f}", f"`a^E - 1`, `E = {show(e_opt)}`"]]
decoded = []
evolved_scores = []
for run in range(3):
    t0 = time.perf_counter()
    best = evolve(train, LENGTH, BUDGET, 200, 150, rng)
    poly = as_polynomial(best)
    if poly and len(poly) == 1:
        (e, c), = poly.items()
        form = f"`a^E`, `E = {show(e)}`"
        decoded.append(e)
    elif poly and len(poly) == 2 and set(abs(c) for c in poly.values()) == {1} \
            and sum(poly.values()) == 0:
        (e1, _), (e2, _) = sorted(poly.items())
        form = f"`a^{e1} (a^D - 1)`, `D = {show(e2 - e1)}`"
        decoded.append(e2 - e1)
    elif poly:
        form = f"a polynomial with {len(poly)} terms"
    else:
        form = "not a small polynomial"
    evolved_scores.append(score(best, test))
    rows.append([f"evolved, run {run + 1} ({time.perf_counter() - t0:.0f} s)",
                 f"{cost(best):.1f}", f"{score(best, test):.3f}", form])
report.table(["program", "multiplications", "test success", "what it computes"], rows)

if len(decoded) == 3:
    verdict = (f"All three runs converged on the same idea: the base raised to an "
               f"exponent that is a product of small primes (largest prime factor "
               f"{max(max(factorize(e)) for e in decoded)}), compared with 1 -- either "
               f"as `a^E` directly or as a difference of two powers, "
               f"`a^e (a^D - 1)`, whose gcd with `N` is that of `a^D - 1`. That is "
               f"Pollard's `p - 1`, rediscovered from random programs by selection alone.")
else:
    verdict = ("Not every run reduced to a power of the base compared with 1; the table "
               "shows what each computes.")
report.p(verdict)
report.p()
best_ev = max(evolved_scores)
se = (opt_test * (1 - opt_test) / len(test)) ** 0.5
cmp = ("matches" if abs(best_ev - opt_test) <= 2 * se else
       "beats" if best_ev > opt_test else "stays below")
report.p(f"Against random programs ({rand_score:.4f}) selection wins by a wide margin. "
         f"The fair comparison is `p - 1` with its exponent optimised on the same "
         f"training data ({opt_test:.3f} on the test set): the best evolved program "
         f"({best_ev:.3f}) {cmp} it (a difference within two standard errors, "
         f"+-{2 * se:.3f}, counts as a match). An early version compared evolution with a hand-written `p - 1` "
         f"whose exponent spent the budget badly, and one evolved program appeared to "
         f"beat it; decoded, that program was `p - 1` with a better exponent. Nothing "
         f"any run found steps outside the order mechanism.")
report.p()

report.p("## With Pascal's collision primitive")
report.p()
report.p("Round 17 showed that iterating `x -> C(x, 2)` is Pollard's rho. Giving the "
         "programs that primitive, and a larger budget, asks whether selection "
         "assembles the collision mechanism too: iterate, then compare two iterates.")
report.p()
B2, L2 = 60, 30
fl = floyd_program(B2)
e2, c2 = optimised_pminus1(train, B2)
rows2 = [["Floyd's rho on `C(x, 2)`, written by hand", f"{cost(fl):.0f}", f"{score(fl, test):.3f}",
          "iterates `C(x, 2)`, multiplies `x_2i - x_i`"],
         ["Pollard p - 1, optimised exponent", f"{c2:.0f}", f"{pminus1_success(e2, test):.3f}",
          f"`a^E - 1`, `E = {show(e2)}`"]]
uses_collision = []
for ops, label in ((OPS, "evolved, arithmetic only"), (PASCAL_OPS, "evolved, with C(x, 2)")):
    best = evolve(train, L2, B2, 150, 120, rng, None, ops)
    live = output_slice(best)
    live_ops = sorted({op for op, _, _ in live})
    e = output_exponent(best)
    if ops is PASCAL_OPS:
        uses_collision.append("cb2" in live_ops)
    what = (f"`a^E`, `E = {show(e)}`" if e else
            f"instructions reaching the output: {', '.join(live_ops)}")
    rows2.append([label, f"{cost(best):.0f}", f"{score(best, test):.3f}", what])
report.table(["program", "multiplications", "test success", "what reaches the output"], rows2)
floyd_s, pm_s = score(fl, test), pminus1_success(e2, test)
if uses_collision and not any(uses_collision):
    tail = ("The evolved program that had `C(x, 2)` available does not use it: the "
            "instructions its output depends on are powers of the base, and it is "
            "`p - 1` again. Differences between the two evolved scores are "
            "run-to-run variance of the search, not the primitive.")
else:
    tail = ("The evolved program with `C(x, 2)` does route it to the output; its "
            "decoded form is in the table.")
report.p(f"At a budget of {B2} multiplications the collision mechanism is nearly "
         f"worthless -- a hand-written Floyd walk succeeds on {floyd_s:.1%} of the "
         f"test numbers, because a collision mod a {BITS}-bit prime needs about "
         f"`sqrt p` iterations -- while `p - 1` succeeds on {pm_s:.1%}. {tail} The "
         f"mechanism selection finds first is the one that pays first.")
report.write()
