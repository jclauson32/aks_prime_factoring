"""Round 13: taking the Coppersmith door -- what the window actually needs.

Round 12 ended with an equivalence: polynomial-time factoring of a balanced
semiprime is the same problem as producing `p` to within `N^(1/4)`.  The
remaining lever was Coppersmith's `beta^2/d` bound itself.  This round attacks it
not by trying to widen the window -- the bound is tight for this lattice family
-- but by asking what *shape* the `N^(1/4)` budget can take, and whether any
shape is cheaper to fill.

The answer reduces factoring to something startlingly small: one residue modulo
a tiny prime.
"""

from math import gcd, isqrt, log

from _common import Report

from aksfactor.arith import is_prime, sieve
from aksfactor.lattice import (
    congruence_modulus_needed,
    crt_assembly_cost,
    factor_with_congruence,
    residue_candidates,
)

report = Report(
    "exp21_coppersmith_door",
    "The Coppersmith door: the budget is one residue per small prime",
    "Round 13. Not widening the window -- changing the shape of what fills it.",
)

import random  # noqa: E402

rng = random.Random(9)


def randprime(bits):
    while True:
        x = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if is_prime(x):
            return x


report.p("## The budget can be spent as a congruence")
report.p()
report.p("Coppersmith is usually stated with the *high bits* of `p`. It works just "
         "as well with a congruence: given `p = r + M x`, the unknown `x` is "
         "bounded by `sqrt(N)/M`, so the method succeeds once `M >= N^(1/4)`. "
         "(The polynomial `r + M x` is not monic, so it is multiplied by "
         "`M^(-1) mod N` first, which does not move the roots mod `p`.)")
report.p()
rows = []
for pb in (13, 15, 17):
    p = randprime(pb)
    q = randprime(pb)
    if p == q:
        continue
    n = p * q
    need = congruence_modulus_needed(n)
    for mult in (0.5, 1.0, 2.0, 4.0):
        m_mod = max(2, int(need * mult))
        while gcd(m_mod, n) != 1:
            m_mod += 1
        got = factor_with_congruence(n, p % m_mod, m_mod, m=4)
        rows.append([n.bit_length(), f"{m_mod:,}", f"{need:,}", f"{m_mod / need:.2f}",
                     "yes" if got and got[0] * got[1] == n else "no"])
report.table(["bits of N", "M", "N^(1/4)", "M / N^(1/4)", "factored"], rows)
report.p("At finite lattice dimension it needs `M` a small multiple of `N^(1/4)`, "
         "approaching `N^(1/4)` as the lattice grows -- the same budget, a "
         "different shape.")
report.p()

report.p("## Which matters, because M can be built from small primes")
report.p()
report.p("A congruence modulus is composite-friendly: `p mod M` follows by CRT "
         "from `p mod ell` for primes `ell` with `prod ell >= N^(1/4)`, i.e. all "
         "primes up to about `(1/4) ln N`. So the question becomes: **how much "
         "does `N` tell us about `p mod ell` for a single tiny prime?**")
report.p()
report.p("Over `F_ell` the residues of `p` and `q` are the roots of")
report.p()
report.p("```")
report.p("z^2 - s z + N,      s = (p + q) mod ell.")
report.p("```")
report.p()
report.p("`N mod ell` is free. `s` is the whole unknown -- and since "
         "`p + q = N + 1 - phi(N)`, knowing `s` is knowing `phi(N) mod ell`.")
report.p()
rows = []
for ell in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
    total = count = 0
    for _ in range(60):
        n = rng.randrange(2, 10**6)
        if n % ell == 0:
            continue
        total += residue_candidates(n, ell)["candidates"]
        count += 1
    mean = total / count
    rows.append([ell, f"{mean:.2f}", (ell - 1) // 2, f"{mean / ((ell - 1) / 2):.3f}"])
report.table(["ell", "mean candidates for p mod ell", "(ell-1)/2", "ratio"], rows)
report.p("The map `a -> a + N/a` is two-to-one, since `a` and `N/a` give the same "
         "sum. So **`N` gives away exactly the `p <-> q` symmetry -- a factor of "
         "two -- and nothing else.**")
report.p()

report.p("## Is any of it free?")
report.p()
rows = []
for ell in (3, 5, 7, 11, 13):
    free = []
    for residue in range(1, ell):
        info = residue_candidates(residue, ell)
        if info["free"]:
            free.append((residue, info["sums"][0]))
    rows.append([ell, str(free) if free else "none"])
report.table(["ell", "N mod ell values that pin (p+q) mod ell outright"], rows)
report.p("One case in the whole range: `N = 2 (mod 3)` forces `p + q = 0 (mod 3)`, "
         "because the only unordered pair with product `2` is `{1, 2}`. Beyond "
         "that, nothing is free.")
report.p()

report.p("## The cost of assembling the budget")
report.p()
rows = []
for bits in (64, 128, 256, 512, 1024, 2048):
    c = crt_assembly_cost(bits)
    rows.append([bits, c["primes"], c["largest_prime"],
                 f"2^{c['log2_search']:.1f}", f"2^{c['log2_target']:.0f}",
                 f"2^{c['saving_bits']:.0f}"])
report.table(["bits of N", "primes needed", "largest prime y",
              "search prod (ell-1)/2", "N^(1/4)", "symmetry saving"], rows)
report.p("`prod (ell-1)/2 = N^(1/4) / 2^pi(y)`, and `pi(y) = O(log N / log log N)`, "
         "so the saving is `N^o(1)`. At 2048 bits it is a factor of `2^73` -- "
         "enormous in absolute terms, invisible in the exponent.")
report.p()

report.p("## What the door leads to")
report.p()
report.p("Following the shape of the budget rather than its size compresses the "
         "problem further than anything else in this project:")
report.p()
report.p("> **Factoring a balanced semiprime in polynomial time is equivalent to "
         "computing `phi(N) mod ell` for every prime `ell` up to about "
         "`(1/4) ln N`.**")
report.p()
report.p("Each of those is a *single element of a field with `O(log N)` elements*. "
         "`N mod ell` is handed to you; `phi(N) mod ell` is the entire content of "
         "the problem. That is the smallest this project has been able to make "
         "factoring, and it is worth stating plainly how little room is left in "
         "it: for a 2048-bit modulus the whole difficulty is 76 residues, modulo "
         "primes no larger than 383.")
report.p()
report.p("It is also, of course, not easier. Knowing `phi(N) mod ell` for enough "
         "primes gives `phi(N)` outright, and `phi(N)` factors `N` immediately -- "
         "so the reduction is again an equivalence, the third one this project "
         "has run into. The pattern is consistent enough to be worth naming: "
         "**every reformulation of factoring that this repository found "
         "computable turned out to be symmetric in `p` and `q`, and every "
         "asymmetric one turned out to be equivalent to factoring.** The `p <-> q` "
         "symmetry that Theorem 8 first identified in the AKS ring is the same "
         "obstruction that survives here, in a setting with no binomial "
         "coefficients anywhere in it.")
report.write()
