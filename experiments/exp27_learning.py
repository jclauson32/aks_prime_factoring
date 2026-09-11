"""Round 19: can a learner read the factors off the bits of N?

A question that comes up whenever factoring meets machine learning.  The
learner here is the low-degree Fourier one: correlate every parity of up to two
bits of N (three among the top bits) with a target bit of p, keep the
significant ones, predict.  Low-degree structure is what statistical learners
find first; if the bits of p had any, this would see it.

Three controls keep the null result honest: a planted parity it must find, a
bit of p that really is predictable (from the size of N), and a quantity N
really does determine ((p+q) mod 4).
"""

import random

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.leakage import walsh_scan

report = Report(
    "exp27_learning",
    "Can a learner read the factors off the bits of N?",
    "Round 19. A low-degree Fourier learner, three positive controls, and the null.",
)

rng = random.Random(19)
BITS = 40


def randprime(lo, hi):
    while True:
        x = rng.randrange(lo, hi) | 1
        if is_prime(x):
            return x


data = []
while len(data) < 60000:
    p = randprime(1 << 19, 1 << 20)
    q = randprime(p + 2, min(2 * p, 1 << 21))
    if (p * q).bit_length() <= BITS:
        data.append((p * q, p, q))
train, test = data[:45000], data[45000:]

targets = [
    ("control: planted parity bit 7 xor bit 23 of N", "planted",
     lambda n, p, q: ((n >> 7) ^ (n >> 23)) & 1),
    ("control: (p+q)/2 odd, i.e. (p+q) mod 4", "N determines it",
     lambda n, p, q: ((p + q) // 2) & 1),
    ("control: bit 18 of p (second highest)", "size of N",
     lambda n, p, q: (p >> 18) & 1),
    ("bit 17 of p", "size of N", lambda n, p, q: (p >> 17) & 1),
    ("bit 15 of p", "size of N", lambda n, p, q: (p >> 15) & 1),
]
for j in (10, 5, 3, 2, 1):
    targets.append((f"bit {j} of p", "", lambda n, p, q, j=j: (p >> j) & 1))
targets.append(("p = 1 (mod 3)", "", lambda n, p, q: p % 3 == 1))
targets.append(("p = 1 (mod 4)", "", lambda n, p, q: p % 4 == 1))

rows = []
results = {}
for name, why, fn in targets:
    got = walsh_scan(train, test, fn, BITS, triple_top=12)
    results[name] = got
    rows.append([name, f"{got['accuracy']:.4f}", f"{got['majority']:.4f}",
                 f"{got['accuracy'] - got['majority']:+.4f}",
                 f"{got['max_corr']:.4f}", f"{got['null_max']:.4f}", why or "--"])
chars = results[targets[0][0]]["characters"]
report.p(f"{len(data):,} balanced semiprimes of {BITS} bits (`2^19 < p < q < 2p`), "
         f"{len(train):,} to train and {len(test):,} to test; {chars:,} Walsh "
         f"characters (bit 0 dropped: `N` is odd, so it is constant).")
report.p()
report.table(["target", "test accuracy", "majority", "gain", "max |corr|",
              "null max", "source of signal"], rows)

planted = results[targets[0][0]]["accuracy"]
mod4 = results[targets[1][0]]["accuracy"]
top = results[targets[2][0]]
low_names = [t[0] for t in targets[5:]]
worst_gain = max(results[n]["accuracy"] - results[n]["majority"] for n in low_names)
worst_ratio = max(results[n]["max_corr"] / results[n]["null_max"] for n in low_names)
report.p(f"The learner is not blind. It recovers the planted parity "
         f"({planted:.0%}), the residue `N` genuinely determines ({mod4:.0%}: "
         f"`N = 1 (mod 4)` exactly when `p = q (mod 4)`, the `ell = 4` case of "
         f"round 13), and a real but partial signal in the top bits of `p` "
         f"({top['accuracy']:.3f} against {top['majority']:.3f}), which is the "
         f"size of `N` seen through `p ~ sqrt(N/r)`.")
report.p()
report.p(f"Below the top bits there is nothing. For bits 10 down to 1 of `p` and "
         f"its residues mod 3 and 4, the best test-set gain over the majority "
         f"guess is {worst_gain:+.4f}, and the strongest single correlation is "
         f"{worst_ratio:.2f} times the maximum expected from noise alone over "
         f"{chars:,} characters.")
report.p()
report.p("## Degree four, all bits")
report.p()
rows4 = []
deg4 = {}
for name, why, fn in targets[5:]:
    got = walsh_scan(train, test, fn, BITS, max_degree=4)
    deg4[name] = got
    rows4.append([name, f"{got['accuracy']:.4f}", f"{got['majority']:.4f}",
                  f"{got['max_corr']:.4f}", f"{got['null_max']:.4f}"])
chars4 = deg4[targets[5][0]]["characters"]
report.p(f"The same targets against every character of degree at most four over all "
         f"39 non-constant bits of `N` -- {chars4:,} characters:")
report.p()
report.table(["target", "test accuracy", "majority", "max |corr|", "null max"], rows4)
ratio4 = max(g["max_corr"] / g["null_max"] for g in deg4.values())
gain4 = max(g["accuracy"] - g["majority"] for g in deg4.values())
se4 = (0.25 / len(test)) ** 0.5
report.p(f"The strongest correlation is {ratio4:.2f} times the typical noise maximum "
         f"-- maxima over {chars4:,} characters land on either side of it -- and the "
         f"best gain over the majority guess is {gain4:+.4f}, "
         f"{'within' if gain4 < 2 * se4 else 'beyond'} two standard errors "
         f"(+-{2 * se4:.4f}) of the test accuracy. Characters that cross the "
         f"selection threshold by chance do not generalise: where they are used, "
         f"test accuracy does not rise above the majority guess.")
report.p()
report.p("What this rules out is narrow and worth stating exactly: structure of "
         "degree at most four in the bits of `N` that predicts the low or middle "
         "bits of `p`. A learner with more capacity could in principle find "
         "higher-degree structure that this one would miss; the controls "
         "establish only that the pipeline -- sampling, splitting, scoring -- "
         "registers signal when there is some. The picture matches every other "
         "round: `N` carries its size and its symmetric residues, and nothing "
         "that tells `p` from `q`.")
report.write()
