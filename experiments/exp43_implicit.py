"""Round 37: how many bits two moduli may share.

Rounds 12 and 13 priced information handed to you about one `N`.  This round
prices information that is not handed to anyone: bits that several moduli
*share*, whose value nobody knows.  The lattice of May and Ritzenhofen turns
`t` shared low bits of the `p_i` into a factorisation of every `N_i` at once,
as soon as `t` exceeds `alpha k / (k - 1)`, where `alpha` is the size of the
unshared cofactor `q`.

The threshold is the whole point: it never goes below `alpha`, so the shared
part must be larger than the secret it is used to find, and balanced moduli --
where `alpha` is the length of `p` itself -- are out of reach at any `k`.
"""

import random
import statistics
import time

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.implicit import (
    implicit_factor,
    shared_low_primes,
    shared_lsb_lattice,
    threshold_bits,
)

report = Report(
    "exp43_implicit",
    "Bits two moduli share, and what they are worth",
    "Round 37. The implicit-factoring threshold, measured against its prediction.",
)

rng = random.Random(43)
PBITS, ALPHA, SHARED = 200, 40, 100
KS = (2, 3, 5, 8)
TRIALS = 20
SCAN = range(40, 101, 4)


def instance(k, pbits=PBITS, alpha=ALPHA, shared=SHARED):
    pairs = shared_low_primes(k, pbits, shared, alpha, rng)
    return [p * q for p, q in pairs], [p for p, _ in pairs]


# ---------------------------------------------------------------- section 1
report.p("## The threshold")
report.p()
report.p(f"Each cell is the share of {TRIALS} instances factored, for `k` moduli "
         f"of {PBITS + ALPHA} bits whose `p_i` (of {PBITS} bits) agree on their "
         f"`t` low bits, with cofactors `q_i` of {ALPHA} bits. The attack is given "
         f"the moduli and `t` -- never the shared bits themselves.")
report.p()

cases = {k: [instance(k) for _ in range(TRIALS)] for k in KS}
grid = {}
for k in KS:
    for t in SCAN:
        wins = 0
        for moduli, primes in cases[k]:
            got = implicit_factor(moduli, t)
            if got is not None:
                assert sorted(got) == sorted(primes), "returned a wrong factorisation"
                wins += 1
        grid[k, t] = wins / TRIALS
report.table(["t"] + [f"k = {k}" for k in KS],
             [[t] + [f"{grid[k, t]:.0%}" for k in KS] for t in SCAN])

rows = []
for k in KS:
    hits = [t for t in SCAN if grid[k, t] >= 0.5]
    measured = min(hits) if hits else None
    predicted = threshold_bits(k, ALPHA)
    rows.append([k, f"{ALPHA} * {k}/{k - 1} = {predicted}",
                 measured if measured else "not reached",
                 f"{measured - predicted:+d}" if measured else "--"])
report.table(["k", "predicted threshold `alpha k / (k - 1)`",
              f"smallest scanned `t` with >= 50% (step {SCAN.step})", "difference"], rows)
report.p("The prediction is the Gaussian heuristic and nothing more: the wanted "
         "vector has norm about `2^alpha sqrt k`, a generic shortest vector of "
         "this lattice has norm `det^(1/k) = 2^(t (k-1)/k)`, and the attack works "
         "exactly when the first is the smaller. More moduli lower the bar "
         "towards `alpha` and never past it.")
report.p()

# ---------------------------------------------------------------- section 2
report.p("## Two controls")
report.p()
report.p("The first asks whether the lattice is doing anything at all: the same "
         "code on moduli whose primes share nothing. The second is the case that "
         "matters cryptographically -- balanced `p` and `q`, where the cofactor is "
         "as long as `p`, so the threshold exceeds the length of `p` itself and "
         "no amount of sharing can reach it.")
report.p()

unrelated_wins = 0
for _ in range(TRIALS):
    moduli = []
    for _ in range(3):
        while True:
            p = rng.randrange(1 << (PBITS - 1), 1 << PBITS) | 1
            if is_prime(p):
                break
        while True:
            q = rng.randrange(1 << (ALPHA - 1), 1 << ALPHA) | 1
            if is_prime(q):
                break
        moduli.append(p * q)
    unrelated_wins += implicit_factor(moduli, SHARED) is not None

bal_bits, bal_shared = 60, 48
balanced_wins = 0
for _ in range(TRIALS):
    moduli, _ = instance(8, pbits=bal_bits, alpha=bal_bits, shared=bal_shared)
    balanced_wins += implicit_factor(moduli, bal_shared) is not None

report.table(["control", "what the heuristic says", f"factored (of {TRIALS})"],
             [[f"3 unrelated moduli, `t = {SHARED}` claimed",
               "no short vector: the wanted one is not in the lattice",
               unrelated_wins],
              [f"8 balanced moduli ({bal_bits}-bit `p` and `q`), sharing "
               f"{bal_shared} of {bal_bits} bits",
               f"threshold {threshold_bits(8, bal_bits)} > {bal_bits} bits of `p`",
               balanced_wins]])
report.p(f"The second row is the reason this is not an attack on RSA. "
         f"{bal_shared} of {bal_bits} bits of `p` shared between eight moduli is "
         f"already an absurd amount of "
         f"correlation, and it is still short of the threshold, because with `q` "
         f"as long as `p` that threshold ({threshold_bits(8, bal_bits)} bits) is "
         f"longer than `p` itself ({bal_bits} bits) -- at every `k`, since "
         f"`alpha k/(k-1)` is strictly greater than `alpha`.")
report.p()

# ---------------------------------------------------------------- section 3
report.p("## What it costs")
report.p()
rows = []
for k in KS:
    moduli, _ = cases[k][0]
    t0 = time.perf_counter()
    for _ in range(5):
        implicit_factor(moduli, SHARED)
    dt = (time.perf_counter() - t0) / 5
    rows.append([k, f"{k} x {k}", f"{SHARED}", f"{dt * 1000:.1f} ms"])
report.table(["k", "lattice", "entry size (bits)", "one attack"], rows)
report.p("Polynomial in everything, because the lattice is small and its entries "
         "are `t` bits, not `log N` bits. When the condition holds the answer is "
         "immediate; when it fails, no budget helps. That is the same shape as "
         "Coppersmith in round 12 -- a lattice that either has the vector or "
         "does not.")
report.p()

# ---------------------------------------------------------------- section 4
report.p("## What a shared bit is worth")
report.p()
report.table(["what you are given", "how much of it is needed", "who finishes"],
             [["nothing", "--", f"`N^(1/4)` work (rounds 33-36)"],
              ["known bits of `p`, value and all", "a quarter of `log N` (round 12)",
               "Coppersmith"],
              ["bits of `p` shared with one other modulus, value unknown",
               f"more than `2 alpha` -- {threshold_bits(2, ALPHA)} here",
               "this lattice"],
              ["the same, shared with `k` moduli",
               f"more than `alpha k/(k-1)`, tending to `alpha`",
               "this lattice"]])
report.p("A shared bit is worth less than a known bit, and the exchange rate is "
         "the interesting number. Coppersmith finishes from a quarter of `log N`; "
         "this lattice needs the shared part to be longer than the whole cofactor "
         "it is solving for -- for a balanced modulus, more than half of `log N`, "
         "which is more than the secret itself. Correlation between instances is "
         "real information and it is measurably cheaper to exploit than to "
         "acquire, but it is priced above the same wall: what has to be supplied "
         "is always at least the size of what is being hidden.")
report.write()
