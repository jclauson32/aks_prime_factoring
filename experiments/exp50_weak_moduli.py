"""Round 44: the ways a modulus can be weak, each priced.

Every other round in this project attacks a modulus whose primes were drawn at
random, and every one of them ends at the same wall. This round attacks the
other kind: moduli whose *generator* went wrong. Four failures, each a real
published one, each implemented here and measured -- and each with a control
showing the same attack finding nothing on a properly generated modulus.

The point is not that these work. It is what they say about where the hardness
lives: not in multiplication, but in the randomness of what was multiplied.
"""

import random
import time
from math import gcd, isqrt, log2

from _common import Report

from aksfactor.arith import is_prime, sieve
from aksfactor.batch import batch_gcd
from aksfactor.ecm import is_b1b2_smooth, pminus1
from aksfactor.lattice import factor_with_congruence

report = Report(
    "exp50_weak_moduli",
    "Four ways a modulus can be weak, each priced",
    "Round 44. Close primes, smooth `p - 1`, structured primes, a shared prime "
    "in a corpus -- all generator failures, all cheap, none of them factoring.",
)

rng = random.Random(50)


def randprime(bits, condition=None):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v) and (condition is None or condition(v)):
            return v


def next_prime(v):
    v |= 1
    while not is_prime(v):
        v += 2
    return v


# ---------------------------------------------------------------- section 1
report.p("## 1. The primes are too close: Fermat")
report.p()
report.p("`N = ((p+q)/2)^2 - ((q-p)/2)^2`, so a search upward from `sqrt N` for "
         "a perfect square finds the factorisation after about `(q-p)^2 / "
         "(8 sqrt N)` steps. The gap is the whole cost.")
report.p()


def fermat(n, max_steps=2_000_000):
    a = isqrt(n)
    if a * a < n:
        a += 1
    for step in range(1, max_steps + 1):
        b2 = a * a - n
        b = isqrt(b2)
        if b * b == b2:
            return a - b, step
        a += 1
    return None, max_steps


rows = []
for exponent in (8, 12, 16, 20, 24):
    p = randprime(32)
    q = next_prime(p + (1 << exponent))
    n = p * q
    got, steps = fermat(n)
    predicted = (q - p) ** 2 / (8 * isqrt(n))
    rows.append([f"2^{exponent}", f"{q - p:,}", f"{steps:,}", f"{predicted:,.0f}",
                 "yes" if got in (p, q) else "no"])
report.table(["gap asked for", "gap obtained", "Fermat steps", "predicted steps",
              "factored"], rows)

p, q = randprime(32), randprime(32)
_, steps_random = fermat(p * q, max_steps=200_000)
report.p(f"The control: two independent 32-bit primes, gap about `2^31`, and "
         f"Fermat is still searching after {steps_random:,} steps. Predicted for "
         f"that gap is about {(abs(q - p) ** 2) / (8 * isqrt(p * q)):,.0f}. "
         f"Nothing about `N` has changed -- only how its primes were chosen.")
report.p()

# ---------------------------------------------------------------- section 2
report.p("## 2. `p - 1` is smooth: Pollard")
report.p()


def powersmooth(m, bound):
    """What stage 1 actually kills: every prime power at most `bound`."""
    return is_b1b2_smooth(m, bound, bound)


rows = []
for bound in (100, 1000):
    p = randprime(32, lambda v: powersmooth(v - 1, bound))
    q = randprime(32)
    n = p * q
    t0 = time.perf_counter()
    got = pminus1(n, bound)
    dt = time.perf_counter() - t0
    rows.append([bound, "yes" if got in (p, q) else "no", f"{dt * 1000:.1f} ms"])
p, q = randprime(32), randprime(32)
t0 = time.perf_counter()
control = pminus1(p * q, 2000)
rows.append(["control: `p` drawn at random",
             "yes" if control in (p, q) else "no",
             f"{(time.perf_counter() - t0) * 1000:.1f} ms"])
report.table(["`p - 1` powersmooth to", "factored", "time"], rows)
report.p("Round 35 priced the draw itself: about 40% of 15-bit primes have a "
         "200-powersmooth `p - 1`, and the share falls with the size of `p`, "
         "which is why this is a generator's job and not an attack.")
report.p()

# ---------------------------------------------------------------- section 3
report.p("## 3. The primes carry structure: ROCA in miniature")
report.p()
report.p("Infineon's RSALib generated primes as `p = k M + 65537^a mod M` with "
         "`M` a primorial. That makes `p mod M` known up to the `a`, and `a` "
         "ranges only over the order of 65537 mod `M` -- so if `M >= N^(1/4)`, "
         "round 12's lattice finishes from each guess. CVE-2017-15361; the "
         "version here is the same shape with a small `M`.")
report.p()

M = 1
for ell in (2, 3, 5, 7, 11, 13, 17):
    M *= ell
order = 1
value = 65537 % M
while value != 1:
    value = value * 65537 % M
    order += 1

structured_a = rng.randrange(1, order)
base = pow(65537, structured_a, M)
while True:
    k = rng.randrange(1 << 12, 1 << 13)
    p_struct = k * M + base
    if is_prime(p_struct):
        break
q_struct = randprime(p_struct.bit_length())
n_struct = p_struct * q_struct

LATTICE_M = 5                                # round 43 says how much lattice this needs
window = 2 * isqrt(n_struct) // M + 1        # `p` may be either side of sqrt(N)
t0 = time.perf_counter()
tries, found = 0, None
for a in range(order):
    tries += 1
    got = factor_with_congruence(n_struct, pow(65537, a, M), M, m=LATTICE_M,
                                 bound=window)
    if got:
        found = got
        break
dt = time.perf_counter() - t0

p_rand, q_rand = randprime(p_struct.bit_length()), randprime(p_struct.bit_length())
t0 = time.perf_counter()
control_hit = None
for a in range(order):
    got = factor_with_congruence(p_rand * q_rand, pow(65537, a, M), M,
                                 m=LATTICE_M,
                                 bound=2 * isqrt(p_rand * q_rand) // M + 1)
    if got:
        control_hit = got
        break
control_dt = time.perf_counter() - t0

report.table(["modulus", "bits", "`M`", "lattice `m`", "guesses for `a`",
              "lattice calls", "factored", "time"],
             [["`p = kM + 65537^a mod M`", n_struct.bit_length(), f"{M:,}",
               LATTICE_M, order, tries, "yes" if found else "no", f"{dt:.2f} s"],
              ["control: both primes random", (p_rand * q_rand).bit_length(),
               f"{M:,}", LATTICE_M, order, order,
               "yes" if control_hit else "no", f"{control_dt:.2f} s"]])
report.p(f"The window handed to the lattice is `2^{log2(window):.1f}` -- twice "
         f"`sqrt(N)/M`, because the structured prime may be either side of "
         f"`sqrt N` and the attacker knows the key size but not which. Round 43 "
         f"measured what that costs in lattice: `m = 3` reaches only `N^0.203` "
         f"and fails on most instances here, `m = {LATTICE_M}` reaches "
         f"`N^0.234` and does not. That is a prediction from one round deciding "
         f"a parameter in another, and it is why the table above says "
         f"`m = {LATTICE_M}`.")
report.p()
report.p(f"`M = {M:,}` is `2^{log2(M):.1f}`, above `N^(1/4)` for a "
         f"{n_struct.bit_length()}-bit modulus, and 65537 has order {order} "
         f"there -- so the entire keyspace of this generator is {order} lattice "
         f"calls wide. The real flaw used a much larger primorial and needed "
         f"Coppersmith's trick of substituting a divisor of `M`, but the "
         f"accounting is this one: a generator that leaks `p mod M` for "
         f"`M >= N^(1/4)` has given away exactly what round 12 priced at "
         f"`N^(1/4)` of work.")
report.p()

# ---------------------------------------------------------------- section 4
report.p("## 4. Two moduli share a prime: batch gcd")
report.p()

corpus_size, mbits = 2000, 64
corpus = []
for _ in range(corpus_size - 2):
    corpus.append(randprime(mbits // 2) * randprime(mbits // 2))
shared = randprime(mbits // 2)
victims = [shared * randprime(mbits // 2), shared * randprime(mbits // 2)]
corpus += victims
rng.shuffle(corpus)

t0 = time.perf_counter()
gcds = batch_gcd(corpus)
batch_time = time.perf_counter() - t0
hits = [(v, g) for v, g in zip(corpus, gcds) if g > 1]

sample = corpus[:200]
t0 = time.perf_counter()
for i, v in enumerate(sample):
    for w in sample[i + 1:]:
        gcd(v, w)
pair_time = (time.perf_counter() - t0)
pairs_done = len(sample) * (len(sample) - 1) / 2
projected = pair_time / pairs_done * (corpus_size * (corpus_size - 1) / 2)

report.table(["method", "gcds", "time"],
             [["batch gcd (one product tree)", f"{corpus_size:,}",
               f"{batch_time:.2f} s"],
              ["every pair", f"{corpus_size * (corpus_size - 1) // 2:,}",
               f"{projected:.2f} s (projected from 200)"]])
found_primes = sorted({g for _, g in hits})
report.p(f"{len(hits)} of the {corpus_size:,} moduli came back with a nontrivial "
         f"gcd, sharing {len(found_primes)} distinct "
         f"{'prime' if len(found_primes) == 1 else 'primes'}, and the planted "
         f"prime is {'among them' if shared in found_primes else 'NOT among them'}"
         f" -- both victims are factored by their own gcd, without anyone "
         f"knowing in advance which two moduli to compare. This is Heninger, "
         f"Durumeric, Wustrow and Halderman's survey of the "
         f"live web, which found tens of thousands of hosts whose keys fell to "
         f"exactly this, for exactly this reason: entropy failures at boot.")
report.p()

# ---------------------------------------------------------------- section 5
report.p("## What the four have in common")
report.p()
report.table(["failure", "what leaked", "cost of the attack", "cost without it"],
             [["primes too close", "`(q-p)` is small", "`(q-p)^2 / sqrt N` steps",
               "`N^(1/4)`"],
              ["`p - 1` smooth", "the order of `(Z/p)*` factors small",
               "one exponentiation", "`L[1/2]` with redraws (round 35)"],
              ["structured primes (ROCA)", "`p mod M` for a known `M >= N^(1/4)`",
               "`ord_M(65537)` lattice calls", "`N^(1/4)` (round 12)"],
              ["a shared prime", "one prime, twice, in a public corpus",
               "one gcd, or `O~(k)` for the whole corpus", "`L[1/3]` per modulus"],
              ["shared low bits (round 37)", "`t > alpha k/(k-1)` bits of `p`",
               "one LLL", "`L[1/3]` per modulus"]])
report.p("Every row is a property of the generator, and every attack reads "
         "something the generator was supposed to keep random. None of them is "
         "an attack on factoring: give any of these methods a modulus whose "
         "primes were drawn independently and uniformly, and the controls above "
         "are what happens -- nothing.")
report.p()
report.p("Which is the same sentence this project has arrived at from every "
         "other direction. The hardness is not in the multiplication. It is in "
         "having nothing to say about `p` except that `p` divides `N`.")
report.write()
