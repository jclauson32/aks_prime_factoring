"""Round 16: hunting for anything that sees p rather than N.

Round 13 reduced polynomial-time factoring to one question: can a polynomial
computation produce (p+q) mod ell -- equivalently phi(N) mod ell -- for the
small primes ell?  N itself pins that residue only to about (ell-1)/2 values.
This round looks for the missing half-bit everywhere it might hide: inside
Pascal's row, in the low bits, in cheap functions of N, in modular forms, and
in lattice-point counts.
"""

import random
from math import comb, isqrt

from _common import Report

from aksfactor.arith import is_prime, sieve
from aksfactor.central import jacobi
from aksfactor.lattice import residue_candidates
from aksfactor.leakage import (
    bayes_sum_accuracy,
    feature_accuracy,
    r4,
    tau_prime_mod_691,
    tau_sum_candidates,
    tau_table,
    trace_divisor_term,
)
from aksfactor.pascal import row_entry

report = Report(
    "exp24_leakage",
    "Hunting for anything that sees p rather than N",
    "Round 16. Five places the missing residue might hide, and what each costs.",
)

rng = random.Random(16)


def randprime(lo, hi):
    while True:
        x = rng.randrange(lo, hi) | 1
        if is_prime(x):
            return x


def next_prime(x):
    x |= 1
    while not is_prime(x):
        x += 2
    return x


# ---------------------------------------------------------------------------
report.p("## 1. Inside the row: Theorem 21")
report.p()
report.p("Row `pq` contains row `p` at stride `q`: `C(pq, mq) = C(p, m) (mod pq)`. "
         "Multiples of the smaller prime carry `C(q, m) mod p` instead, which is "
         "zero unless `p` fails to divide `C(q, m)`.")
report.p()
primes = [x for x in sieve(200) if x > 2]
tot = bad = 0
for i, p in enumerate(primes):
    for q in primes[i + 1:]:
        n = p * q
        for m in range(p + 1):
            tot += 1
            bad += row_entry(n, m * q) != comb(p, m) % n
rows = []
for p, q in [(101, 103), (1009, 1013), (1009, 1511), (4001, 4003)]:
    n = p * q
    at_q = sum(1 for m in range(1, p) if row_entry(n, m * q))
    at_p = sum(1 for m in range(1, q) if row_entry(n, m * p))
    rows.append([p, q, f"{at_q} of {p - 1}", f"{at_p} of {q - 1}"])
report.p(f"Stride-`q` identity: **{tot - bad} of {tot}** entries over all odd "
         f"prime pairs below 200.")
report.p()
report.table(["p", "q", "non-zero at multiples of q", "non-zero at multiples of p"], rows)
report.p("The second column follows Fine's theorem: `prod(d_i + 1) - 2` over the "
         "base-`p` digits of `q`. When `q` is just above `p` it is a handful and the "
         "support is almost exactly a copy of row `p`; when `q mod p` is large, as "
         "for `1009 x 1511`, it approaches `q`. (An earlier write-up quoted the "
         "density as `~1/p + 1/q`; that is the upper end of this range.) Either "
         "way the row *contains* `p` -- at stride `q`, which is exactly as hard to "
         "find as `q`.")
report.p()

# ---------------------------------------------------------------------------
report.p("## 2. The low bits never prune")
report.p()
n = 1009 * 1013
rows = []
for k in (4, 8, 12, 16):
    mod = 1 << k
    survivors = sum(1 for a in range(1, mod, 2)
                    if (a * (n * pow(a, -1, mod) % mod)) % mod == n % mod)
    rows.append([k, survivors, 1 << (k - 1)])
report.table(["bits k", "candidates for p mod 2^k", "2^(k-1)"], rows)
report.p("Every odd residue of `p` has a partner `q = N/p mod 2^k`, so Hensel "
         "lifting from the bottom keeps all `2^(k-1)` branches -- after `b/2` bits "
         "that is trial division's count. The same holds for any modulus coprime "
         "to `N`: **residues of `N` never constrain `p` beyond the `p <-> q` swap.**")
report.p()

# ---------------------------------------------------------------------------
report.p("## 3. Cheap functions of N")
report.p()
report.p("A held-out lookup table predicts `(p+q) mod ell` from `(N mod ell, f(N))`. "
         "The right yardstick is the **Bayes baseline** `2/(ell-1)`: the most "
         "likely allowed sum has two preimages `a, N/a`. (A first version of this "
         "scan compared against a uniform guess over the allowed set, about "
         "`1/((ell-1)/2 + 1)`; every feature 'beat' it by the same margin, which "
         "is the tell that the baseline was wrong, not the features right.)")
report.p()


def cf_second(n):
    a0 = isqrt(n)
    m, d, a = 0, 1, a0
    out = []
    for _ in range(2):
        m = d * a - m
        d = (n - m * m) // d
        if d == 0:
            return 0
        a = (a0 + m) // d
        out.append(a)
    return out[1]


FEATURES = {
    "floor(sqrt N)": lambda n, ell: isqrt(n) % ell,
    "Fermat gap ceil(sqrt N)^2 - N": lambda n, ell: ((isqrt(n) + 1) ** 2 - n) % ell,
    "second CF digit of sqrt N": lambda n, ell: cf_second(n) % ell,
    "next base-ell digit of N": lambda n, ell: (n // ell) % ell,
    "jacobi(3, N)": lambda n, ell: jacobi(3, n),
    "jacobi(ell, N)": lambda n, ell: jacobi(ell, n),
    "popcount(N)": lambda n, ell: bin(n).count("1") % ell,
}
ELLS = (5, 7, 11, 13)


def sample(count, close):
    out = []
    while len(out) < count:
        p = randprime(1 << 19, 1 << 20)
        q = next_prime(p + 2) if close else randprime(p + 2, p + int(0.6 * p))
        out.append((p * q, p, q))
    return out


results = {}
for label, close in (("close primes (positive control)", True),
                     ("random balanced semiprimes", False)):
    data = sample(6000, close)
    rows = []
    for name, f in FEATURES.items():
        cells = [name]
        for ell in ELLS:
            acc = feature_accuracy(data, ell, f)
            results[(label, name, ell)] = acc
            cells.append(f"{acc:.3f}")
        rows.append(cells)
    rows.append(["**Bayes baseline 2/(ell-1)**"] +
                [f"{bayes_sum_accuracy(ell):.3f}" for ell in ELLS])
    bits = sorted({n.bit_length() for n, _, _ in data})
    report.p(f"**{label}**, {len(data):,} semiprimes of {bits[0]}-{bits[-1]} bits:")
    report.p()
    report.table(["feature"] + [f"ell = {ell}" for ell in ELLS], rows)

control = min(results[("close primes (positive control)", "floor(sqrt N)", ell)]
              for ell in ELLS)
worst = max(results[("random balanced semiprimes", name, ell)] - bayes_sum_accuracy(ell)
            for name in FEATURES for ell in ELLS)
report.p(f"The detector works: on close primes `floor(sqrt N)` predicts the sum "
         f"with accuracy at least **{control:.3f}** at every `ell`, because there "
         f"`p + q` is within a few units of `2 sqrt N`. On random balanced "
         f"semiprimes the largest excess of any feature over the Bayes baseline is "
         f"**{worst:+.3f}** -- {'sampling noise' if worst < 0.03 else 'worth a second look'}. "
         f"A feature that sees the sum only when the primes are close is a size "
         f"method (Fermat), not a new source.")
report.p()

# ---------------------------------------------------------------------------
report.p("## 4. Ramanujan's tau")
report.p()
tau = tau_table(400)
ok_cong = all((tau[m] - sum(d ** 11 for d in range(1, m + 1) if m % d == 0)) % 691 == 0
              for m in range(1, 400))
report.p(f"`tau(n) = sigma_11(n) (mod 691)` for all `n < 400`: **{ok_cong}**. With "
         f"multiplicativity, `tau(N) = (1 + p^11)(1 + q^11) (mod 691)`, and because "
         f"`gcd(11, 690) = 1` the eleventh-power map is a bijection on `F_691` -- so "
         f"`tau(N) mod 691` pins `{{p, q}} mod 691`.")
report.p()
before = []
after = []
for _ in range(300):
    p = randprime(10 ** 4, 10 ** 5)
    q = randprime(p + 2, 2 * p)
    if p % 691 == 0 or q % 691 == 0:
        continue
    t = tau_prime_mod_691(p) * tau_prime_mod_691(q) % 691
    got = tau_sum_candidates(p * q, t)
    assert (p + q) % 691 in got
    before.append(residue_candidates(p * q, 691)["candidates"])
    after.append(len(got))
report.table(["knowing", "candidates for (p+q) mod 691 (mean of 300)"],
             [["N", f"{sum(before) / len(before):.1f}"],
              ["N and tau(N) mod 691", f"{sum(after) / len(after):.2f}"]])
report.p("It is worth a full residue per prime -- exactly the currency round 13 "
         "asked for. Every route to it passes through the factorisation:")
report.p()
rows = []
for p, q in [(5, 7), (11, 13), (29, 31)]:
    rows.append([p * q, trace_divisor_term(p * q), 2 + 2 * p ** 11])
report.table(["N", "Eichler-Selberg divisor term", "2 + 2 p^11"], rows)
parity = all((tau[m] % 2 == 1) == (isqrt(m) ** 2 == m and isqrt(m) % 2 == 1)
             for m in range(1, 400))
report.p(f"- the **trace formula** for `tau(N)` contains `sum_(d | N) min(d, N/d)^11`, "
         f"which for `N = pq` is `2 + 2 p^11` -- it names the smaller factor;")
report.p("- the **q-expansion** reaches `tau(N)` only after `N` coefficients;")
report.p(f"- **mod 2** the expansion collapses to `sum q^((2k+1)^2)` (checked below "
         f"400: {parity}) -- cheap, and it says only whether `N` is an odd square.")
report.p()
report.p("This is a known wall: Bach and Charles (*The hardness of computing an "
         "eigenform*) showed that evaluating eigenform coefficients such as "
         "`tau(N)` at RSA moduli is as hard as factoring them. Edixhoven and "
         "Couveignes compute `tau(p)` in polynomial time for **prime** `p`; their "
         "method uses the primality of `p` in an essential way.")
report.p()

# ---------------------------------------------------------------------------
report.p("## 5. Lattice points on a sphere")
report.p()
rows = []
for p, q in [(3, 5), (3, 7), (5, 7), (7, 11), (3, 13)]:
    n = p * q
    count = r4(n)
    rows.append([n, count, 8 * (n + 1 + p + q), count // 8 - n - 1, p + q])
report.table(["N", "r_4(N)", "8 sigma(N)", "recovered p+q", "true p+q"], rows)
report.p("Jacobi's four-square theorem: counting the integer points on the "
         "3-sphere of radius `sqrt N` gives `sigma(N) = N + 1 + p + q`. The volume "
         "term `pi^2 N / 2` averages over `N`; `p + q` lives entirely in the "
         "fluctuation around the average -- the same place the divisor function "
         "keeps it. The direct count (`r_4 = r_2 * r_2`, a convolution over all "
         "`k <= N`) costs about `N` operations, and no method is known that computes "
         "it without, in effect, factoring `N`.")
report.p()

# ---------------------------------------------------------------------------
report.p("## Verdict")
report.p()
report.p("| source | worth per small prime | cost |")
report.p("|---|---|---|")
report.p(f"| row `n` mod `n` at stride `q` | the whole of `p` | finding the stride is finding `q` |")
report.p(f"| low bits / residues of `N` | nothing beyond `p <-> q` | free |")
report.p(f"| cheap functions of `N` | nothing measurable (max excess {worst:+.3f}) | free |")
report.p(f"| `tau(N) mod 691` | an exact residue ({sum(after) / len(after):.2f} candidates) | as hard as factoring (Bach-Charles) |")
report.p(f"| `r_4(N)` | `p + q` exactly | about `N` directly; otherwise factoring |")
report.p()
report.p("The pattern from round 13 holds without exception: every quantity that "
         "carries the residue is a function of the divisors, and every quantity "
         "computable from `N` in polynomial time is symmetric under `p <-> q` and "
         "carries nothing more than `N` itself.")
report.write()
