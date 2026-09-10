# Theory

Everything below is elementary and self-contained. Each statement has a machine
checker in [`aksfactor/theorems.py`](../aksfactor/theorems.py) and is verified
exhaustively by `tests/test_theorems.py`.

## Setup

For an integer `n >= 2` write the reduced Pascal row

```
R(n) = [ C(n,0) mod n, C(n,1) mod n, ..., C(n,n) mod n ]
```

and call the positions `0 < k < n` the **interior**. Both ends are `1`.

`R(n)` is the AKS identity in coefficient form. AKS rests on

```
(x + a)^n  ==  x^n + a   (mod n)      for gcd(a,n) = 1, and n prime
```

whose `a = 1` case says precisely that the interior of `R(n)` vanishes. The
question this repository asks is what the interior says when it *doesn't*
vanish.

Throughout, for `0 < k < n` put

```
d = gcd(n, k),      x = n/d,      k' = k/d.
```

## Lemma 0 (absorption)

For `1 <= k <= n`,  `k * C(n,k) = n * C(n-1,k-1)`.

*Proof.* Both sides count pairs `(S, s)` with `S` a `k`-subset of `[n]` and
`s in S`: pick the subset then the element, or pick the element then the rest. ∎

## Theorem 1 — every residue is `x * y` with `x` a proper divisor of `n`

For `0 < k < n`, `x = n/gcd(n,k)` divides `C(n,k)`. Hence

```
C(n,k) mod n  =  x * y     with   0 <= y < d,
```

where `x` is a divisor of `n` with `2 <= x <= n/2`.

*Proof.* Lemma 0 divided by `d` gives `k' * C(n,k) = x * C(n-1,k-1)`. So
`x | k' * C(n,k)`, and `gcd(x,k') = 1`, so `x | C(n,k)`. Write `C(n,k) = x*m`.
Since `n = x*d`, reducing mod `n` keeps the factor `x`:
`x*m mod x*d = x*(m mod d)`. Set `y = m mod d`. Finally `d = gcd(n,k) <= k < n`,
so `d` is a proper divisor and `x = n/d >= 2`. ∎

**This is the conjecture, and it is a theorem: the split always exists, and `x`
is always a factor of `n`.**

### Corollary 1.1 — every non-zero residue hands you a factor

If `C(n,k) mod n != 0` then `gcd(C(n,k) mod n, n)` is a non-trivial proper
divisor of `n`.

*Proof.* `x` divides both the residue and `n`, so `x | gcd`, and `x >= 2`. The
residue is a non-zero value below `n`, so the gcd is at most it, hence `< n`. ∎

## Theorem 2 — the support sits only on positions sharing a factor with `n`

If `gcd(k,n) = 1` and `0 < k < n`, then `n | C(n,k)`.

*Proof.* `d = 1` in Theorem 1, so `x = n` divides `C(n,k)`. ∎

Equivalently `supp(interior of R(n)) ⊆ { k : gcd(k,n) > 1 }`. For prime `n` that
set is empty, recovering the classical fact that the interior of row `p` is
divisible by `p`.

## Theorem 3 — the exact value at a prime factor

If `p` is prime and `p | n`, then

```
C(n, p) mod n  =  n/p       exactly  (i.e. y = 1).
```

*Proof.* Take `k = p`, so `d = p` and `k' = 1`. Lemma 0 gives
`C(n,p) = (n/p) * C(n-1,p-1)`, so in the notation of Theorem 1,
`m = C(n-1,p-1)` and the residue is `(n/p) * (C(n-1,p-1) mod p)`.

Now evaluate `C(n-1, p-1) mod p` by Lucas' theorem. Because `p | n` we have
`n - 1 ≡ p - 1 (mod p)`, so the least significant base-`p` digit of `n-1` is
`p-1`. The base-`p` digits of `p-1` are `(p-1, 0, 0, ...)`. Lucas multiplies
digitwise:

```
C(n-1, p-1) ≡ C(p-1, p-1) * prod_{i>=1} C(digit_i(n-1), 0) = 1 * 1 * ... = 1  (mod p).
```

Hence the residue is `(n/p) * 1 = n/p`. ∎

## Theorem 4 — the first failure names the smallest prime factor

Let `n >= 2`.

1. `n` is prime **iff** the interior of `R(n)` is identically zero.
2. If `n` is composite with smallest prime factor `p`, then

```
min { k > 0 : C(n,k) mod n != 0 }  =  p        and the residue there is  n/p.
```

*Proof.* (1) For prime `n` every interior `k` is coprime to `n`, so Theorem 2
applies. Conversely a composite `n` has the non-zero entry supplied by (2).

(2) By Theorem 2 a non-zero interior entry forces `gcd(k,n) > 1`, so `k` shares
a prime factor with `n` and therefore `k >= p`. By Theorem 3 the entry at `k = p`
equals `n/p`, which is non-zero because `p < n`. So the minimum is attained
exactly at `p`. ∎

### Corollary 4.1 — one lookup, both factors

For composite `n`, if `k0` is the first non-zero interior position and `r0` its
value, then `k0 = spf(n)`, `r0 = n/spf(n)`, and

```
n = k0 * r0.
```

### Corollary 4.2 — semiprimes

If `n = p*q` with `p <= q` prime, the first non-zero interior entry of `R(n)` is
**exactly `q`**. So for the semiprime case the conjecture holds in its strongest
form: the remainder *is* a prime factor, with `y = 1`.

## Theorem 5 — closed form for an arbitrary entry

For `0 < k < n`, with `d = gcd(n,k)` and `k' = k/d`:

```
C(n,k) mod n  =  (n/d) * ( ( C(n-1,k-1) / k' )  mod d ).
```

*Proof.* From `k' C(n,k) = (n/d) C(n-1,k-1)` and `gcd(n/d, k') = 1` we get
`k' | C(n-1,k-1)`. Put `m = C(n-1,k-1)/k'`; then `C(n,k) = (n/d) m` and
Theorem 1's reduction gives the claim. ∎

### Corollary 5.1 — random access

Because `k = k' * d`, the identity `C(n-1,k-1) = k' * m` gives

```
C(n-1,k-1) mod k  =  k' * (m mod d),
```

so a single entry of `R(n)` needs `C(n-1,k-1)` **only modulo the small number
`k`**. Lucas (prime modulus) and Granville (prime powers) evaluate that in time
polynomial in `log n`, with a table whose size is the largest prime power
dividing `k`.

This is what `aksfactor.pascal.row_entry` does: `C(n,k) mod n` for `n` with
thousands of digits, in microseconds, with no row ever materialised.

## Theorem 6 — the folded object

Let `p` be prime with `p^v || n`, and `s = n / p^v`. In
`(Z/n)[x] / (x^r - 1)` write `(x + a)^n = sum_j S_j x^j`. Then for every `j`:

```
S_j  ≡  sum over { i : i*p^v ≡ j (mod r) }  of  C(s, i) * a^(s-i)     (mod p).
```

*Proof.* In `F_p[x]`, Frobenius gives `(x+a)^(p^v) = x^(p^v) + a^(p^v)`, and
`a^(p^v) = a` by Fermat. So `(x+a)^n = (x^(p^v) + a)^s`. Expanding by the
binomial theorem and folding exponents modulo `r` gives the stated sum. ∎

## Theorem 7 — the aliasing barrier

Let `p` be a prime factor of `n` with `p^v || n`, and let `r >= 1` with
`gcd(r, p) = 1`. Then the multiples of `p^v` are equidistributed among the `r`
residue classes modulo `r`.

*Proof.* The multiples of `p^v` below `N` are `i * p^v` for `i = 0, 1, 2, ...`.
Since `gcd(p^v, r) = 1`, the map `i -> i * p^v (mod r)` is a bijection of `Z/r`,
so consecutive blocks of `r` multiples hit each class exactly once. Occupancies
therefore differ by at most one partial period. ∎

**Consequence.** Folding the row at level `r` sorts positions into classes by
`k mod r`. If `r < p` then `gcd(r, p) = 1`, so every class contains the same
number of support positions: the classes are *indistinguishable by occupancy*.
No positional information about `p` survives the fold. The only signal left is
the numerical value of the class sums — the `1/p` event of Theorem 6.

To get positional information one needs `gcd(r, p) > 1`, hence `r >= p`. This is
an aliasing bound of exactly the Nyquist kind: a period-`p` structure cannot be
resolved by sampling it into fewer than `p` bins.

## Theorem 8 — the symmetry obstruction

For odd `r` and any `a`, in `(Z/n)[x] / (x^r - 1)`,

```
Norm( (x + a)^n )  =  prod over r-th roots of unity w of  (w + a)^n
                   =  (a^r + 1)^n     (mod n).
```

*Proof.* `prod_{w^r = 1} (X - w) = X^r - 1`, so
`prod_w (a + w) = (-1)^r prod_w (-a - w) = (-1)^r ((-a)^r - 1) = a^r + 1` for
odd `r`. The norm is multiplicative, so the norm of the `n`-th power is the
`n`-th power of the norm. ∎

**Consequence.** The right-hand side depends only on `n`, `r`, `a`. It takes the
same value modulo *every* prime factor of `n`, so

```
gcd( Norm((x+a)^n) - (a^r + 1)^n ,  n )  =  n     identically.
```

The norm aggregates all `r` fold coefficients into one number and, in doing so,
destroys precisely the asymmetry a factor consists of.

**The general principle.** Any quantity invariant under the Galois action
permuting the `r`-th roots of unity is Frobenius-invariant modulo every `p | n`,
hence expressible in terms of `n` alone. Symmetric aggregation can never
factor. A factoring statistic must distinguish `p` from `q`, and by Theorem 6
the only such statistic available at level `r` is the accidental vanishing of a
class sum.

## Theorem 9 — the second n-adic digit

Let `0 < k < n` with `gcd(k,n) = 1`. Theorem 2 gives `n | C(n,k)`, so write

```
C(n,k) = n * m.
```

Then `m ≡ C(n-1,k-1) * k^(-1) (mod n)`, and for every prime `p | n`:

```
m ≡ 0 (mod p)
   <=>  C(n-1, k-1) ≡ 0 (mod p)
   <=>  some base-p digit of k-1 strictly exceeds the corresponding digit of n-1.
```

*Proof.* Lemma 0 gives `k * C(n,k) = n * C(n-1,k-1)`. Substituting `C(n,k) = n m`
and cancelling `n` gives `k * m = C(n-1,k-1)`. Since `gcd(k,n) = 1` we have
`p` does not divide `k`, so `m ≡ 0 (mod p)` iff `C(n-1,k-1) ≡ 0 (mod p)`. Lucas'
theorem expresses `C(n-1,k-1) mod p` as a product of digit binomials
`C(a_i, b_i)`, which vanishes iff `b_i > a_i` for some `i`. ∎

**Consequence.** `gcd(m mod n, n)` is a proper factor of `n` exactly when the
digit condition holds for some prime factors of `n` and not others.

## Proposition 10 — the second digit is factoring-hard, at constant density

Let `n = pq`. For `k` uniform among the residues coprime to `n`, let `A_p` be the
event of Theorem 9 for `p`. Then `gcd(m mod n, n)` splits `n` exactly on the
symmetric difference `A_p Δ A_q`, and

```
Pr[A_p] = 1 - prod_i (a_i + 1)/p,      a_i = base-p digits of n-1.
```

Each factor `(a_i + 1)/p` is uniform-ish in `(0,1]`, and `n` has only about
`log_p n` digits base `p` — three, for a balanced semiprime. So `Pr[A_p]` is
bounded away from `0` and `1` **independently of the size of `p`**, and so is
`Pr[A_p Δ A_q]`.

Measured over `p` from `37` to `8.1e5`, the splitting rate stays in roughly
`0.26 – 0.97` while `1/p` falls from `2.7e-2` to `1.2e-6`
([exp09](../experiments/results/exp09_second_digit.md)).

**Consequence.** An algorithm computing `C(n,k) mod n^2` for uniform random `k`
in time `poly(log n)` factors `n = pq` in `O(1)` expected iterations. Computing
binomial coefficients modulo `n^2` is therefore factoring-hard.

### Corollary 10.1 — random access cannot be deepened

Theorem 5 gives genuine random access to `C(n,k) mod n`, for `n` of any size, in
time polynomial in `log n` — because the closed form only ever needs
`C(n-1,k-1)` modulo the *small* number `k`, where Lucas and Granville apply.
Proposition 10 says extending that access by a single n-adic digit, to modulus
`n^2`, would break factoring.

So the requirement in Granville-style algorithms that the modulus be *factored*
before binomials can be reduced modulo it is intrinsic, not an artifact of how
those algorithms happen to be written.

### The barrier, in its sharpest form

The factors are not hidden in the Pascal row. They sit in the second n-adic
digit at **constant density** — one random probe in a few would reveal one. The
entire obstruction is that reading that digit is itself the problem being
solved.

The rest of this document — the `1/p` fold lottery, the aliasing of Theorem 7,
the symmetry of Theorem 8 — describes the price of working at the
first digit, where the information has been quotiented away. Proposition 10 says
the second digit is where the information lives, and that it is priced exactly
at factoring.

## Proposition 11 — why the first-digit family stalls at `Theta(spf(n))`

Two corollaries of Theorem 6 bound what any "read the AKS object" strategy can
do.

**(a) Truncation is blind below degree `p`.** Modulo `p`, the object is
`(x^(p^v) + a)^s`, supported only on multiples of `p^v`. Every coefficient of
degree `0 < k < p^v` vanishes mod `p`, and Theorem 2 upgrades that to vanishing
mod `n` at every `k` coprime to `n`. A prefix of length `L < spf(n)` is
identically zero: it certifies "no factor below `L`" and nothing more.

**(b) Folding buys compression but not information.** Theorem 7 says the fold's
classes are equidistributed, so nothing is learned from *which* class fires, and
Theorem 8 says symmetric aggregates of the classes are constant. Reducing mod `x^r - 1`
packs the whole row into `r` numbers for `O(log n)` ring multiplications. But by
Theorem 6 the coefficient `S_j` leaks a factor exactly when a binomial
progression sum vanishes mod `p`. Nothing forces it to; empirically it behaves
like a uniform residue, so each coefficient is a `~1/p` event, and
`Theta(p)` coefficients are needed. This is measured in
[`experiments/results/exp04_fold_scaling.md`](../experiments/results/exp04_fold_scaling.md):
the ratio *coefficients inspected / p* stays near 1 across three orders of
magnitude, and extra AKS bases `a` change the number of tickets but never the
price.

So both routes cost `Theta(spf(n))` — trial-division complexity.

**What a breakthrough would need.** Not a proof of impossibility; a statement of
where the wall is. One would need either

- a fold coefficient that vanishes mod `p` with probability `>> 1/p`, i.e. an
  algebraic reason for the progression sums of Theorem 6 to collapse, or
- an operator that isolates the degree-`p^v` structure of `(x^(p^v)+a)^s`
  without enumerating `Theta(p)` coefficients. Differentiation is the obvious
  candidate and it fails: `d/dx` does not commute with reduction mod `x^r - 1`,
  and under truncation it needs degree `> p` before it sees anything, by (a).

## Proposition 12 — the first digit is a period-finding problem

Modulo `p`, the AKS object is `(x^(p^v) + a)^(n/p^v)`: a series whose support is
the arithmetic progression `p^v * Z`. Recovering `p` from it is exactly the
problem of recovering the period of that progression.

- By Theorem 7, sampling into `r < p` buckets aliases the period away.
- By Theorem 8, symmetric aggregates of the buckets are period-blind.
- By Theorem 6, an individual bucket reveals the period only through a `1/p`
  coincidence.

Classical period finding over `Z/n` needs `Omega(p)` samples. This is the same
problem Shor's algorithm solves in polylog time, by taking a Fourier transform
of size `n` rather than of size `r << p`. What this framework lacks is not a
better `r`, a better base `a`, or a better aggregation rule — it is the ability
to transform at full resolution.

## Relationship to AKS

AKS verifies `(x+a)^n == x^n + a (mod n, x^r - 1)` for `r` of size `polylog(n)`
and many `a`, and needs only a **yes/no**: does the identity hold? Theorem 6
says the identity's failure mod `p` is governed by binomial progression sums of
the cofactor. Primality only asks whether those sums are all consistent with a
prime; factoring asks *which* `p` produced them, and by Proposition 11 that
answer is spread over `Theta(p)` coefficients rather than concentrated in
`polylog(n)` of them.

That gap — detecting failure is cheap, localising it is not — is the precise
sense in which AKS does not extend to factoring here.

And Proposition 10 says where that localisation cost actually lives: not in the
Pascal row, which carries the factors at constant density one n-adic digit down,
but in the price of reading that digit — which is factoring itself.
