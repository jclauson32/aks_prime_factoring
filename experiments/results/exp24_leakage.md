# Hunting for anything that sees p rather than N

Round 16. Five places the missing residue might hide, and what each costs.

## 1. Inside the row: Theorem 21

Row `pq` contains row `p` at stride `q`: `C(pq, mq) = C(p, m) (mod pq)`. Multiples of the smaller prime carry `C(q, m) mod p` instead, which is zero unless `p` fails to divide `C(q, m)`.

Stride-`q` identity: **58300 of 58300** entries over all odd prime pairs below 200.

| p | q | non-zero at multiples of q | non-zero at multiples of p |
|---|---|---|---|
| 101 | 103 | 100 of 100 | 4 of 102 |
| 1009 | 1013 | 1008 of 1008 | 8 of 1012 |
| 1009 | 1511 | 1008 of 1008 | 1004 of 1510 |
| 4001 | 4003 | 4000 of 4000 | 4 of 4002 |

The second column follows Fine's theorem: `prod(d_i + 1) - 2` over the base-`p` digits of `q`. When `q` is just above `p` it is a handful and the support is almost exactly a copy of row `p`; when `q mod p` is large, as for `1009 x 1511`, it approaches `q`. (An earlier write-up quoted the density as `~1/p + 1/q`; that is the upper end of this range.) Either way the row *contains* `p` -- at stride `q`, which is exactly as hard to find as `q`.

## 2. The low bits never prune

| bits k | candidates for p mod 2^k | 2^(k-1) |
|---|---|---|
| 4 | 8 | 8 |
| 8 | 128 | 128 |
| 12 | 2048 | 2048 |
| 16 | 32768 | 32768 |

Every odd residue of `p` has a partner `q = N/p mod 2^k`, so Hensel lifting from the bottom keeps all `2^(k-1)` branches -- after `b/2` bits that is trial division's count. The same holds for any modulus coprime to `N`: **residues of `N` never constrain `p` beyond the `p <-> q` swap.**

## 3. Cheap functions of N

A held-out lookup table predicts `(p+q) mod ell` from `(N mod ell, f(N))`. The right yardstick is the **Bayes baseline** `2/(ell-1)`: the most likely allowed sum has two preimages `a, N/a`. (A first version of this scan compared against a uniform guess over the allowed set, about `1/((ell-1)/2 + 1)`; every feature 'beat' it by the same margin, which is the tell that the baseline was wrong, not the features right.)

**close primes (positive control)**, 6,000 semiprimes of 39-40 bits:

| feature | ell = 5 | ell = 7 | ell = 11 | ell = 13 |
|---|---|---|---|---|
| floor(sqrt N) | 1.000 | 1.000 | 1.000 | 1.000 |
| Fermat gap ceil(sqrt N)^2 - N | 0.631 | 0.595 | 0.554 | 0.546 |
| second CF digit of sqrt N | 0.641 | 0.506 | 0.461 | 0.457 |
| next base-ell digit of N | 0.547 | 0.449 | 0.297 | 0.250 |
| jacobi(3, N) | 0.554 | 0.439 | 0.357 | 0.333 |
| jacobi(ell, N) | 0.558 | 0.386 | 0.257 | 0.210 |
| popcount(N) | 0.554 | 0.382 | 0.232 | 0.191 |
| **Bayes baseline 2/(ell-1)** | 0.500 | 0.333 | 0.200 | 0.167 |

**random balanced semiprimes**, 6,000 semiprimes of 39-41 bits:

| feature | ell = 5 | ell = 7 | ell = 11 | ell = 13 |
|---|---|---|---|---|
| floor(sqrt N) | 0.505 | 0.330 | 0.200 | 0.152 |
| Fermat gap ceil(sqrt N)^2 - N | 0.495 | 0.328 | 0.196 | 0.160 |
| second CF digit of sqrt N | 0.496 | 0.342 | 0.207 | 0.158 |
| next base-ell digit of N | 0.491 | 0.331 | 0.201 | 0.148 |
| jacobi(3, N) | 0.497 | 0.325 | 0.192 | 0.160 |
| jacobi(ell, N) | 0.483 | 0.331 | 0.210 | 0.155 |
| popcount(N) | 0.497 | 0.345 | 0.188 | 0.175 |
| **Bayes baseline 2/(ell-1)** | 0.500 | 0.333 | 0.200 | 0.167 |

The detector works: on close primes `floor(sqrt N)` predicts the sum with accuracy at least **1.000** at every `ell`, because there `p + q` is within a few units of `2 sqrt N`. On random balanced semiprimes the largest excess of any feature over the Bayes baseline is **+0.011** -- sampling noise. A feature that sees the sum only when the primes are close is a size method (Fermat), not a new source.

## 4. Ramanujan's tau

`tau(n) = sigma_11(n) (mod 691)` for all `n < 400`: **True**. With multiplicativity, `tau(N) = (1 + p^11)(1 + q^11) (mod 691)`, and because `gcd(11, 690) = 1` the eleventh-power map is a bijection on `F_691` -- so `tau(N) mod 691` pins `{p, q} mod 691`.

| knowing | candidates for (p+q) mod 691 (mean of 300) |
|---|---|
| N | 345.5 |
| N and tau(N) mod 691 | 1.00 |

It is worth a full residue per prime -- exactly the currency round 13 asked for. Every route to it passes through the factorisation:

| N | Eichler-Selberg divisor term | 2 + 2 p^11 |
|---|---|---|
| 35 | 97656252 | 97656252 |
| 143 | 570623341224 | 570623341224 |
| 899 | 24401019531411660 | 24401019531411660 |

- the **trace formula** for `tau(N)` contains `sum_(d | N) min(d, N/d)^11`, which for `N = pq` is `2 + 2 p^11` -- it names the smaller factor;
- the **q-expansion** reaches `tau(N)` only after `N` coefficients;
- **mod 2** the expansion collapses to `sum q^((2k+1)^2)` (checked below 400: True) -- cheap, and it says only whether `N` is an odd square.

This is a known wall: Bach and Charles (*The hardness of computing an eigenform*) showed that evaluating eigenform coefficients such as `tau(N)` at RSA moduli is as hard as factoring them. Edixhoven and Couveignes compute `tau(p)` in polynomial time for **prime** `p`; their method uses the primality of `p` in an essential way.

## 5. Lattice points on a sphere

| N | r_4(N) | 8 sigma(N) | recovered p+q | true p+q |
|---|---|---|---|---|
| 15 | 192 | 192 | 8 | 8 |
| 21 | 256 | 256 | 10 | 10 |
| 35 | 384 | 384 | 12 | 12 |
| 77 | 768 | 768 | 18 | 18 |
| 39 | 448 | 448 | 16 | 16 |

Jacobi's four-square theorem: counting the integer points on the 3-sphere of radius `sqrt N` gives `sigma(N) = N + 1 + p + q`. The volume term `pi^2 N / 2` averages over `N`; `p + q` lives entirely in the fluctuation around the average -- the same place the divisor function keeps it. The direct count (`r_4 = r_2 * r_2`, a convolution over all `k <= N`) costs about `N` operations, and no method is known that computes it without, in effect, factoring `N`.

## Verdict

| source | worth per small prime | cost |
|---|---|---|
| row `n` mod `n` at stride `q` | the whole of `p` | finding the stride is finding `q` |
| low bits / residues of `N` | nothing beyond `p <-> q` | free |
| cheap functions of `N` | nothing measurable (max excess +0.011) | free |
| `tau(N) mod 691` | an exact residue (1.00 candidates) | as hard as factoring (Bach-Charles) |
| `r_4(N)` | `p + q` exactly | about `N` directly; otherwise factoring |

The pattern from round 13 holds without exception: every quantity that carries the residue is a function of the divisors, and every quantity computable from `N` in polynomial time is symmetric under `p <-> q` and carries nothing more than `N` itself.
_Generated by `experiments/exp24_leakage.py` in 3.4s._
