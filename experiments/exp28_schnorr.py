"""Round 20: testing Schnorr's lattice factoring at small sizes.

In 2021 Schnorr claimed that lattice reduction finds the smooth relations of a
sieve fast enough to break RSA.  The claim did not survive scrutiny, but its
core is a clean, testable question, and this repository has the tools to ask
it: *are the residues r = u - vN produced by close vectors in the prime-number
lattice smooth more often than integers of the same size?*
"""

import math
import random
import statistics
import time

from _common import Report

from aksfactor.arith import is_prime, sieve
from aksfactor.schnorr import (
    prime_lattice_relations,
    residue_baselines,
    smooth_exponents,
    split_from_relations,
    support_coprime,
)

report = Report(
    "exp28_schnorr",
    "Schnorr's prime-number lattice, tested",
    "Round 20. Does lattice reduction find smooth relations better than chance?",
)

rng = random.Random(28)


def smooth(m, primes):
    return smooth_exponents(m, primes) is not None


def parity_baseline(r, u, v, primes, k=60):
    """The second attempt: within 10% of |r|, coprime to the used primes, r's parity."""
    used = [p for p in primes if u % p == 0 or v % p == 0]
    a = abs(r)
    lo, hi = max(2, int(a * 0.9)), int(a * 1.1) + 2
    hits = got = 0
    while got < k:
        x = rng.randrange(lo, hi)
        if (x - a) % 2 or any(x % p == 0 for p in used):
            continue
        got += 1
        hits += smooth(x, primes)
    return hits / k


report.p("## Setup")
report.p()
report.p("For primes `p_1 .. p_n`, the lattice spanned by "
         "`b_i = (0, .., f_i, .., 0, round(C ln p_i))` with `f` a random permutation "
         "of `ceil(i/2)` and `C = N`, target `t = (0, .., 0, round(C ln N))`. Kannan's "
         "embedding and an exact integer LLL give close vectors "
         "`sum e_i b_i`; each yields `u = prod p_i^(e_i > 0)`, `v = prod p_i^(-e_i < 0)` "
         "and `r = u - vN`. A smooth `r` is a relation `u = r (mod N)`; `n + 2` of "
         "them factor `N` by linear algebra mod 2.")
report.p()

report.p("## A structural handicap")
report.p()
checked = ok = 0
for _ in range(20):
    for u, v, r in prime_lattice_relations(1000003 * 999983, sieve(100)[:20], 1.0, rng):
        checked += 1
        ok += support_coprime(u, v, r, sieve(100)[:20])
report.p(f"`u` and `v` use disjoint primes, and `r` is divisible by none of them: "
         f"`p | u` forces `r = -vN != 0 (mod p)`, and `p | v` forces `r = u`. "
         f"Checked on {ok} of {checked} candidates. A close vector that uses half "
         f"the factor base leaves `r` only the other half to be smooth over.")
report.p()
report.p("There is a partial compensation at the primes the vector did *not* use. "
         "There `u` and `v` are units, so `p^k | r` with probability "
         "`1/phi(p^k)` rather than `1/p^k` -- at `p = 2`, `r` is always even. The "
         "fair baseline therefore draws integers within 10% of `|r|` from exactly "
         "this local law at every prime below 60, coprime to every prime the "
         "vector used.")
report.p()

report.p("## Measurement")
report.p()
rows = []
yields = {}
parity_totals = []
for bits, n, lattices in [(32, 20, 300), (32, 30, 300), (36, 30, 300),
                          (40, 30, 300), (44, 30, 300)]:
    while True:
        p = rng.randrange(1 << (bits // 2 - 1), 1 << (bits // 2)) | 1
        q = rng.randrange(1 << (bits // 2 - 1), 1 << (bits // 2)) | 1
        if p != q and is_prime(p) and is_prime(q):
            break
    N = p * q
    primes = sieve(1000)[:n]
    seen, rels = {}, []
    exp_a = exp_b = exp_par = 0.0
    split, used = None, 0
    sizes = []
    t0 = time.time()
    for _ in range(lattices):
        used += 1
        for u, v, r in prime_lattice_relations(N, primes, 1.0, rng):
            if (u, v) in seen:
                continue
            seen[(u, v)] = r
            sizes.append(abs(r).bit_length())
            a, b = residue_baselines(r, u, v, primes, rng)
            exp_a += a
            exp_b += b
            exp_par += parity_baseline(r, u, v, primes)
            if smooth(r, primes):
                rels.append((u, r))
        if len(rels) > n + 3:
            split = split_from_relations(N, primes, rels)
            if split:
                break
    yields[bits, n] = len(rels) / used
    parity_totals.append(exp_par)
    rows.append([bits, n, used, f"{len(seen):,}", statistics.median(sizes), len(rels),
                 f"{exp_a:.1f}", f"{exp_b:.1f}",
                 f"`{split}`" if split else f"no ({len(rels)} of {n + 2})",
                 f"{time.time() - t0:.0f} s"])
report.table(["bits of N", "primes n", "lattices", "distinct (u, v)",
              "median bits of r", "smooth r", "expected: same size",
              "expected: same size and local law", "factor found", "time"], rows)

obs = sum(r[5] for r in rows)
fair = sum(float(r[7]) for r in rows)
naive = sum(float(r[6]) for r in rows)
z = (obs - fair) / math.sqrt(fair) if fair else float("nan")
verdict = ("as smooth as integers of their size with the same local law, and no "
           "smoother -- the lattice buys a smaller `r`, and nothing about its "
           "factorisation"
           if abs(z) < 2.5 else
           "measurably different from integers of their size with the same local "
           "law -- worth a closer look")
report.p(f"Across the table: **{obs}** smooth residues against **{fair:.1f}** "
         f"expected from the fair baseline (Poisson z = {z:+.1f}) and {naive:.1f} "
         f"from the naive one. The lattice-found residues are {verdict}.")
report.p()
par = sum(parity_totals)
zp = (obs - par) / math.sqrt(par) if par else float("nan")
report.p(f"The baseline took three attempts, and the first two were wrong in "
         f"instructive ways. Integers of the same size predict {naive / obs:.1f} "
         f"times too many smooth residues -- they ignore the handicap. Adding "
         f"only coprimality and parity predicts {par:.1f} (z = {zp:+.1f}): it "
         f"omits the unit effect at 3, 5, 7, ..., which biases it "
         f"{(1 - par / fair):.0%} low -- enough, in a larger sample, to pass for "
         f"a lattice advantage. With the full local law the difference is "
         f"z = {z:+.1f}.")
report.p()
y32, y40 = yields[32, 30], yields[40, 30]
bigger = sum(1 for r in rows if r[4] > r[0])
factored = [r[0] for r in rows if not str(r[8]).startswith("no")]
worked = (f"The method works -- it factored the {', '.join(map(str, factored))}-bit "
          f"number{'s' if len(factored) > 1 else ''} -- " if factored else
          "Within these budgets it factored nothing, ")
report.p(f"And the size it buys is not small: the median residue is larger than "
         f"`N` in {bigger} of {len(rows)} rows. The relation yield falls from "
         f"{y32:.3f} per lattice at 32 bits to {y40:.3f} at 40 bits and "
         f"{yields[44, 30]:.3f} at 44, while the number of relations needed grows "
         f"with the factor base. {worked}at a cost set by the smoothness probability "
         f"of the residues it finds; Dixon's random squares, by comparison, have "
         f"residues below `N`.")
report.p()
report.p("This is a test at toy sizes, with LLL and Kannan's embedding rather than "
         "the BKZ and pruned enumeration Schnorr proposed; stronger reduction finds "
         "closer vectors and smaller residues. What it cannot change is the "
         "structural handicap, which is a property of the residues and not of the "
         "solver" + (", nor the absence of a smoothness bias beyond it."
                     if abs(z) < 2.5 else "."))
report.write()
