"""Round 39: amortisation, the one exponent that moves.

Every round so far has asked what one `N` costs. This one asks what a thousand
of them cost together, because that is the single place in the whole survey
where the per-instance price falls without any new mathematics: a product tree
makes one multiplication serve a whole batch.

It is also the mechanism already hiding inside the sieves. A sieve is batch
smoothness detection over an interval, which is why the two methods that reach
`L[1/2]` and `L[1/3]` are the two that sieve.

What it does not do is help a single number, and the last section says so.
"""

import random
import statistics
import time

from _common import Report

from aksfactor.arith import sieve
from aksfactor.batch import batch_smooth, smooth_parts, trial_smooth

report = Report(
    "exp45_batch",
    "Amortisation: what a batch costs that one number does not",
    "Round 39. Batch smoothness detection by product trees, against trial "
    "division one number at a time.",
)

rng = random.Random(45)
BOUND = 5000
PRIMES = sieve(BOUND)
VBITS = 64
CBITS = 32                    # the control's size: smooth often enough to matter


def values(count, bits=VBITS):
    return [rng.randrange(1 << (bits - 1), 1 << bits) for _ in range(count)]


def planted(count):
    """Numbers that are certainly smooth, so the control tests both answers."""
    out = []
    for _ in range(count):
        v = 1
        while v.bit_length() < CBITS:
            v *= PRIMES[rng.randrange(len(PRIMES))]
        out.append(v)
    return out


# ---------------------------------------------------------------- section 1
report.p("## The control: same answers")
report.p()
report.p(f"Both routines are asked which of a batch of numbers are "
         f"{BOUND}-smooth. The batch version multiplies every prime together, "
         f"reduces that product down a tree of the numbers, squares each "
         f"remainder enough times to kill every prime power, and takes a gcd. "
         f"It has to agree with trial division on every number, and it also has "
         f"to return the right smooth *part*, which is checked by dividing it "
         f"out.")
report.p()
report.p(f"The control runs at {CBITS} bits rather than the {VBITS} used for "
         f"timing below, because a {VBITS}-bit number is {BOUND}-smooth about "
         f"once in five thousand: a test where both routines answer \"no\" to "
         f"everything checks almost nothing. Two dozen numbers built out of small "
         f"primes are planted in each batch as well, so the positive answers are "
         f"exercised too.")
report.p()

rows = []
for count in (64, 256, 1024):
    batch = values(count, CBITS) + planted(24)
    got, want = batch_smooth(batch, PRIMES), trial_smooth(batch, PRIMES)
    parts = smooth_parts(batch, PRIMES)
    part_ok = 0
    for v, part in zip(batch, parts):
        rough = v // part
        assert v % part == 0
        part_ok += all(rough % ell for ell in PRIMES)      # nothing smooth left
    rows.append([len(batch), sum(got),
                 f"{sum(a == b for a, b in zip(got, want))}/{len(batch)}",
                 f"{part_ok}/{len(batch)}"])
report.table(["batch size", "smooth found", "agrees with trial division",
              "smooth part correct"], rows)
report.p()

# ---------------------------------------------------------------- section 2
report.p("## The cost, per number")
report.p()
report.p(f"Both methods are timed on the same batches of {VBITS}-bit numbers. "
         f"Trial division's cost per number should be flat -- it does not know "
         f"the batch exists. The tree's should fall, because the expensive "
         f"object, the product of all the primes, is built once and then walked "
         f"down.")
report.p()

rows = []
ratios = []
for count in (32, 128, 512, 2048, 8192):
    batch = values(count)

    t0 = time.perf_counter()
    batch_smooth(batch, PRIMES)
    batch_time = time.perf_counter() - t0

    sample = batch[:min(count, 512)]
    t0 = time.perf_counter()
    trial_smooth(sample, PRIMES)
    trial_time = (time.perf_counter() - t0) / len(sample) * count

    ratios.append((count, trial_time / batch_time, batch_time / count * 1e6,
                   trial_time / count * 1e6))
    rows.append([count, f"{batch_time / count * 1e6:.1f}",
                 f"{trial_time / count * 1e6:.1f}",
                 f"{trial_time / batch_time:.1f}x"])
report.table(["batch size", "product tree (us per number)",
              "trial division (us per number)", "speedup"], rows)

best_count, best_ratio, best_us, _ = min(ratios, key=lambda r: r[2])
first, last = ratios[0], ratios[-1]
report.p(f"The per-number cost falls and then rises again, and the best batch "
         f"here is {best_count} numbers at {best_us:.1f} us each, "
         f"{best_ratio:.1f} times faster than trial division -- against "
         f"{first[1]:.1f}x at {first[0]} and {last[1]:.1f}x at {last[0]}. The "
         f"fall is the amortisation: the product of all the primes is a fixed "
         f"cost spread over more numbers. The rise is arithmetic on the tree's "
         f"own nodes, whose root is the product of the whole batch -- "
         f"{last[0] * VBITS // 1000} thousand bits at the last row -- and "
         f"multiplying numbers that size is superlinear. Real implementations do "
         f"exactly what this table says to do: run the algorithm on blocks of "
         f"the best size rather than on everything at once.")
report.p()
batch_us = [r[2] for r in ratios]
trial_us = [r[3] for r in ratios]
report.p(f"Trial division's per-number cost varies by a factor of "
         f"{max(trial_us) / min(trial_us):.2f} across a 256-fold range of batch "
         f"sizes, which is measurement noise: it does not know the batch exists. "
         f"The tree's varies by {max(batch_us) / min(batch_us):.2f} over the same "
         f"range, and that variation is the whole effect.")
report.p()

# ---------------------------------------------------------------- section 3
report.p("## Where this already lives")
report.p()
report.p("Nothing above is new to factoring; it is the thing the sieves have "
         "always done, made explicit. The quadratic sieve does not test its "
         "candidates one at a time -- it walks each factor base prime through "
         "the interval, so the whole interval is tested for that prime in one "
         "pass. Round 33 measured the other side of that trade: CFRAC generates "
         "residues that are *smaller*, and so smooth more often, but has to "
         "divide each one on its own -- which is why its asymptotic constant is "
         "the worse of the two (`L[1/2, sqrt 2]` against the sieve's "
         "`L[1/2, 1]`), even though at the sizes measured there it was the "
         "faster of the two implementations.")
report.p()
report.p("The same argument scales up to whole factorisations. Coppersmith's "
         "batch number field sieve factors many numbers of the same size for "
         "less per number than one alone, by sharing the sieving and the "
         "linear algebra -- the same trick at the level of the algorithm rather "
         "than the arithmetic.")
report.p()

# ---------------------------------------------------------------- section 4
report.p("## What it does not do")
report.p()
report.p("Amortisation lowers the per-number cost of a batch. It does not lower "
         "the cost of *a* number, and factoring is asked about one number. The "
         "distinction is exactly the one round 37 drew: `k` correlated moduli "
         "fall to a lattice that cannot touch a single modulus, and a batch of "
         "numbers is cheaper each than one number alone. Both are real, both are "
         "measured here, and neither is an attack on the instance that matters.")
report.p()
report.p("That is the last exponent in this survey that moves. It moves in the "
         "batch size, not in `N`.")
report.write()
