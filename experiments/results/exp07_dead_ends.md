# Dead ends, and the mechanism they share

Four attempts to beat `Theta(spf(n))`, each tested and each killed.

## A. Symmetric functions of the fold (norms, resultants)

A single fold coefficient leaks a factor with probability `~1/p`. The *norm* `prod_t F(w^t)` aggregates all of them at once, so it looked like a way to buy every ticket with one purchase.

| n | r | a | Norm((x+a)^n) | (a^r + 1)^n mod n | equal |
|---|---|---|---|---|---|
| 15 | 5 | 2 | 12 | 12 | True |
| 35 | 5 | 2 | 17 | 17 | True |
| 77 | 5 | 2 | 66 | 66 | True |
| 143 | 5 | 2 | 132 | 132 | True |
| 105 | 5 | 2 | 48 | 48 | True |
| 1155 | 5 | 2 | 132 | 132 | True |
| 27 | 5 | 2 | 0 | 0 | True |
| 49 | 5 | 2 | 19 | 19 | True |
| 221 | 5 | 2 | 50 | 50 | True |

Verified for every `(n, r, a)` tested (odd `r`): **True**.

**Why it dies.** The norm equals `(a^r + 1)^n mod n` -- a function of `n`, `r`, `a` alone. It takes the same value modulo *every* prime factor of `n`, so `gcd(Norm - expected, n) = n` identically. Aggregating the tickets symmetrically destroys exactly the asymmetry that a factor consists of.

**The general lesson.** Any quantity invariant under the Galois action that permutes the `r`-th roots of unity is Frobenius-invariant mod every `p | n`, hence determined by `n`. Factoring needs a statistic that distinguishes `p` from `q`; symmetric ones cannot.

## B. Twisted moduli `x^r - c`

If the plain fold gives each coefficient a `1/p` chance, maybe weighting the wrapped terms by powers of `c` biases some coefficient toward vanishing.

| p | q | plain `x^r-1` rate | twisted `x^r-c` rate | 1/p |
|---|---|---|---|---|
| 251 | 443 | 0.01440 | 0.00849 | 0.00398 |
| 811 | 2909 | 0.00206 | 0.00206 | 0.00123 |
| 2011 | 7369 | 0.00000 | 0.00026 | 0.00050 |
| 7027 | 14957 | 0.00000 | 0.00000 | 0.00014 |

Both track `1/p`. The twist changes which tickets you hold, not the price. Dead.

## C. The AKS-ring generalisation of Pollard `p-1`  *(partly retracted)*

In `F_p[x]/(x^r - 1) = prod_i F_{p^{k_i}}`, an element's order divides `lcm_i (p^{k_i} - 1)`. Raising `x + a` to `lcm(1..B)` and taking a gcd generalises Pollard `p-1`, which is the `k = 1` case. More `k` values, more chances -- so it seemed strictly better.

| p | smoothness of p-1 | of p^2-1 | of p^3-1 | larger k helps? |
|---|---|---|---|---|
| 101 | 5 | 17 | 10303 | no |
| 1009 | 7 | 101 | 9181 | no |
| 10007 | 5003 | 5003 | 461521 | no |
| 65537 | 2 | 331 | 116085511 | no |
| 99991 | 101 | 431 | 372751 | no |

**Why it dies:** `p - 1` divides `p^k - 1` for every `k`, so `p^k - 1` is `B`-smooth only if `p - 1` already was. Working with `x + a` in the full unit group is dominated by the classical `p-1` method it generalises.

> **Correction (round 3).** The sentence that used to stand here said *strictly dominated*, full stop. That was too strong, and [exp10](exp10_norm_one.md) refutes it. The argument above is about the **full unit group**. The norm-one subgroup of `F_{p^d}*` has order `(p^d - 1)/(p - 1)`, which `p - 1` does **not** divide -- and you can land in it without knowing `p`, by choosing a monic polynomial whose roots multiply to `1`. For `d = 2` that is `x^2 - a x + 1`, i.e. Lucas sequences, i.e. Williams `p+1`. Primes with `p-1` rough and `p+1` smooth are plentiful, and the norm-one method factors them while `p-1` cannot. This entry is a dead end only for the unconstrained element `x + a`.

## D. Reading position out of the fold

If some residue class mod `r` were *occupied differently* by the support of the row, the class index would leak `p mod r`, and CRT over several `r` would reconstruct `p` in polylog time. This is the one that would actually have been polynomial.

| n | r | support occupancy per class mod r | spread | even |
|---|---|---|---|---|
| 35 | 9 | [0, 2, 1, 2, 0, 2, 1, 2, 0] | 2 | True |
| 77 | 5 | [3, 4, 3, 3, 3] | 1 | True |
| 77 | 9 | [1, 2, 2, 2, 2, 1, 2, 2, 2] | 1 | True |
| 143 | 5 | [4, 5, 5, 4, 4] | 1 | True |
| 143 | 7 | [2, 3, 3, 2, 4, 4, 4] | 2 | True |
| 143 | 9 | [2, 2, 3, 2, 4, 2, 3, 2, 2] | 2 | True |
| 221 | 5 | [5, 5, 6, 6, 6] | 1 | True |
| 221 | 7 | [3, 4, 4, 4, 3, 5, 5] | 2 | True |
| 221 | 9 | [2, 3, 3, 3, 3, 2, 4, 4, 4] | 2 | True |
| 1001 | 5 | [56, 56, 56, 56, 56] | 0 | True |
| 1001 | 9 | [31, 30, 31, 32, 32, 30, 30, 32, 32] | 2 | True |

**Why it dies.** If `gcd(r, p) = 1` then `i -> i*p mod r` is a bijection of `Z/r`, so the multiples of `p` fall into every class equally. The occupancy is flat by *theorem*, not by luck. Positional information about `p` requires `gcd(r, p) > 1`, i.e. `r >= p` -- an aliasing barrier exactly like Nyquist: you cannot resolve a period-`p` signal by sampling it into fewer than `p` bins.

## The shared mechanism

All four die the same way. Modulo `p` the AKS object is `(x^(p^v) + a)^(n/p^v)`, whose only distinguishing feature is a period of `p^v`. Every cheap thing you can compute from it is either

- **symmetric** in the prime factors (A, and C for the unconstrained element) -- same value mod every `p`, so the gcd is `n`; or
- **aliased** below the period (B, D) -- the period-`p` structure spreads evenly over all `r < p` buckets, leaving only the accidental vanishing of a bucket sum, a `1/p` event.

Extracting `p` is *period finding*. Classically that needs `Omega(p)` samples; it is the same problem Shor's algorithm solves in polylog time quantumly, by taking a Fourier transform of size `n` rather than size `r << p`. That is a sharp statement of what this framework is missing, and it is not something a cleverer choice of `r` or `a` can supply.
_Generated by `experiments/exp07_dead_ends.py` in 1.1s._
