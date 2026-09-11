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

## Theorem 13 — the norm-one order, and the gap in the domination argument

Let `p` be prime and `d >= 2`. The norm map `N : F_{p^d}* -> F_p*`,
`N(z) = z^(1 + p + ... + p^(d-1))`, is surjective, so its kernel has order

```
|ker N| = (p^d - 1)/(p - 1).
```

`p - 1` does **not** divide this in general. (For `d = 2` it is `p + 1`, which
shares only a factor of `2` with `p - 1`.)

*Consequence, and a correction.* An earlier draft of this project argued that
generalising Pollard `p-1` into the AKS ring is strictly dominated, because
`p - 1 | p^k - 1` forces `p^k - 1` to be `B`-smooth only when `p - 1` already is.
That argument is valid for the **full unit group**, where the element `x + a`
lives. It does not apply to `ker N`, whose order is the quotient above.

*Reaching the kernel without knowing `p`.* If `f` is monic of degree `d` with
constant term `(-1)^d`, its roots multiply to `1`, so a root `z` satisfies
`N(z) = 1` **by construction**. For `d = 2`, `f = x^2 - a x + 1` and the trace of
`z^k` is the Lucas sequence `V_k(a, 1)` — Williams' `p+1` method.

For `d = 2` the root lands in `F_{p^2} \ F_p` exactly when `a^2 - 4` is a
quadratic non-residue mod `p`; otherwise `z ∈ F_p` and the order divides `p - 1`
again. So the construction must be retried over several `a`, each a coin flip.

## Proposition 14 — rigid versus varying group orders

Each factoring method reachable from this framework succeeds on a prime `p`
exactly when some integer attached to `p` is smooth:

| method | group | order |
|---|---|---|
| Pollard `p-1` | `F_p*` | `p - 1` |
| norm-one, degree `d` (Theorem 13) | `ker N <= F_{p^d}*` | `(p^d - 1)/(p - 1)` |
| ECM | `E(F_p)` | `p + 1 - t`, with `t` bounded by `2 sqrt(p)` |

In the first two rows the order is a **fixed function of `p` and `d`**: once `p`
is given, so is the number whose smoothness decides the outcome. Small `d` gives
finitely many such numbers, so a prime for which all of them are rough is out of
reach, and no amount of extra work at that `p` helps.

The third row is different in kind. Each elliptic curve over `F_p` supplies a
*new* order in an interval of width `4 sqrt(p)` around `p + 1`, so one may keep
drawing fresh orders at fixed `p` until a smooth one appears. That single
structural difference — a family of group orders rather than one — is what turns
smoothness luck into the `L[1/2]` running time of ECM.

*Measured.* Over 400 random monic `f` of degree `d` at each of
`p = 211, 503, 1009, 2003`, the number of distinct values of `|(F_p[x]/f)*|` is
exactly `partitions(d)` — `2, 3, 5` for `d = 2, 3, 4` — and does not grow with
`p`. Over 400 random elliptic curves at the same primes, the number of distinct
`#E(F_p)` is `55, 84, 109, 140`, tracking the Hasse width `4 sqrt(p)`
([exp11](../experiments/results/exp11_rigid_vs_varying.md)).

The cap is structural. For squarefree `f`, `F_p[x]/f` is a product of finite
fields `F_{p^{d_1}} x ... x F_{p^{d_m}}`, so the unit group order is
`prod_i (p^{d_i} - 1)`: a function of the degree *partition* alone. Two
polynomials with the same pattern give the same order.

**Open.** Is there a family of AKS-ring-like objects whose group order varies
with a parameter at fixed `p`? Theorem 13 supplies one order per `(p, d)`;
elliptic curves supply unboundedly many per `p`. Nothing in this repository
bridges that gap, and closing it is a strictly stronger requirement than the
"asymmetric statistic" that Theorem 13 already provides.

The requirement is satisfiable, just not here: elliptic curves supply a fresh
order per curve (ECM), and class groups of imaginary quadratic orders supply a
fresh `h(-D)` per discriminant `D = -kn` (Schnorr-Lenstra). Both leave polynomial
quotients of `Z/n` behind entirely, and that departure is precisely what buys the
varying order.

## Theorem 15 — no row leaks below `spf(n)`

Let `k < spf(n)`. Then `C(N, k) mod n` depends only on `N mod n` and `k`.

*Proof.* Every `i <= k` is coprime to `n` (else `n` would have a prime factor
`<= k < spf(n)`), so `k!` is invertible mod `n` and

```
C(N,k) = (N)(N-1)...(N-k+1) / k!   ==   prod_i ((N mod n) - i) * (k!)^(-1)  (mod n).
```

Both sides depend on `N` only through `N mod n`. ∎

**Consequence.** Rows `N` and `N'` with `N ≡ N' (mod n)` are *identical* mod `n`
at every position below `spf(n)`. So "try a different row" — row `2n`, row
`n+1`, row `n^2`, any row at all — buys exactly nothing in the region where the
classical row is empty. Theorem 2 is the special case `N = n`, where the shared
value happens to be `0`.

This closes the row-shifting family completely, and it explains why: below
`spf(n)` the binomial coefficient is a *polynomial identity* in `N mod n`,
carrying no arithmetic information about how `n` factors.

## Theorem 16 — the q-analogue of Kummer's theorem

Let `q` be coprime to `n`, `p | n` prime, and `d = ord_p(q)`. For the Gaussian
binomial coefficient `[n,k]_q`:

```
p | [n,k]_q   <=>   (k mod d) > (n mod d)                  [clause 1]
               or    p | C(floor(n/d), floor(k/d))          [clause 2]
```

*Proof sketch.* `p | Phi_d(q)`, and q-Lucas gives
`[n,k]_q ≡ C(floor(n/d), floor(k/d)) * [n mod d, k mod d]_q (mod Phi_d(q))`.
The second factor vanishes identically iff `k mod d > n mod d`. ∎

Verified exhaustively against the explicit q-Pascal row in
`tests/test_qpascal.py`.

**What changes.** The classical row's period is `p`, rigidly (Theorem 7). The
q-row's period is `d = ord_p(q)`, which **moves with the base `q`** — the first
construction in this project whose governing parameter genuinely varies at fixed
`p`. Clause 1 fires at

```
k = (n mod d) + 1        (when n mod d <= d - 2),
```

which can be far below `spf(n)`, inside the region Theorem 2 seals off.

**What does not change.** `ord_p(q)` divides `p - 1` for every `q`, so the
achievable periods form the divisor lattice of `p - 1`. Clause 1 fires early
exactly when `ord_p(q) | n - j` for a small `j`, i.e. when
`gcd(q^(n-j) - 1, n) > 1`. Measured work stays `Theta(p)`
([exp12](../experiments/results/exp12_qdeformation.md)).

The reason is subtler than "it reduces to Pollard `p-1`", which was this
project's first guess and is **false**: with `p - 1` smooth a random base still
has order close to `p - 1`, so clause 1 needs a *small* order, not a smooth one,
and the raw q-deformation is strictly weaker than Pollard.

Tuning the base is what closes the gap, and it closes it onto a dichotomy. Put
`Q = q^M` for `M` a smooth prime-power ladder to bound `B`; then
`ord_p(Q) = ord_p(q)/gcd(ord_p(q), M)` is the `B`-rough part of the order, hence

```
ord_p(Q) = 1        (when ord_p(q) is B-smooth)      or      ord_p(Q) > B.
```

Measured over 120 random pairs, the middle regime `1 < ord_p(Q) <= B` was **empty**.
And `ord_p(Q) = 1` means `Q ≡ 1 (mod p)`, so `gcd(Q - 1, n)` has already split
`n` — that is Pollard `p-1` verbatim, with the q-row adding nothing. So tuning
the period either lands exactly on Pollard or leaves clause 1 as a `1/B`
lottery. There is no regime in which the tunable period pays for itself.

### The refined hierarchy

| construction | governing quantity | varies at fixed `p`? | over what set |
|---|---|---|---|
| classical Pascal row mod `n` | period `p` | no | -- |
| **q-Pascal row mod `n`** | period `ord_p(q)` | **yes** | divisors of `p - 1` |
| norm-one subgroup, degree `d` | order `Phi_d(p)` | no | -- |
| elliptic curve | order `p + 1 - t` | yes | an interval of width `4 sqrt(p)` |

Proposition 14 asked for a parameter that varies. Theorem 16 supplies one, and
in doing so sharpens the question: **it is not variation that matters, but the
*density* of the set varied over.**

The divisors of `p - 1` are a sparse, structured set, and the elements realising
the small ones are vanishingly rare — which is why tuning collapses to the
dichotomy above. An elliptic curve's order ranges over an *interval* of width
`4 sqrt(p)`: dense, so every fresh curve is a genuinely fresh number, and
sampling until one is smooth is a well-behaved process. That density, not the
mere existence of a parameter, is what separates `L[1/2]` from `Theta(p)`.

**Open, restated for round 5.** Find a construction over `Z/n` whose governing
quantity at fixed `p` ranges over a *dense* set — an interval, not a divisor
lattice — while remaining computable without knowing `p`.

## Proposition 17 — the smoothness ceiling

Every construction reached in this project succeeds on a prime `p` exactly when
some integer attached to `p` is `B`-smooth. Write `G(p, theta)` for the group
supplied by the construction at parameter `theta`, and `ord G(p, theta)` for the
integer whose smoothness is tested.

| construction | `ord G(p, theta)` | reach at fixed `p` |
|---|---|---|
| classical Pascal row (T4, T7) | `p` | a single value |
| q-Pascal row (T16) | `ord_p(q)` | divisors of `p - 1` |
| norm-one subgroup (T13) | `(p^d - 1)/(p - 1)` | one value per `(p, d)` |
| ring unit group (P14) | `prod_i (p^{d_i} - 1)` | `partitions(d)` values |
| elliptic curve | `p + 1 - t`, with `t` bounded by `2 sqrt(p)` | an interval of width `4 sqrt(p)` |
| class group of `-kn` | `h(-kn) ~ sqrt(kn)` | dense in `k` |

The first four rows are sparse and cap out at `Theta(p)` or at a smoothness
condition on one fixed number. The last two are dense, and dense sampling is
what makes `L[1/2]` achievable.

**But `L[1/2]` is also the ceiling of the whole table.** The probability that a
uniform integer of size `x` is `x^(1/u)`-smooth is `u^(-u + o(u))` (Dickman), so
the expected number of samples before a smooth order appears, balanced against
the cost of testing each, is minimised at `L[1/2]` regardless of how the samples
are drawn. Improving the *density* of the sampled set can move a construction up
to that ceiling; it cannot move it past.

This is the precise sense in which the number field sieve is a different kind of
algorithm rather than a better-sampled one: `L[1/3]` comes from testing
smoothness of *smaller* numbers (algebraic norms), not from sampling a denser
family of the same-sized ones.

**Consequence for this project.** Rounds 3, 4 and 5 supplied, in turn, an
asymmetric statistic, a varying parameter, and variation over a dense set. Each
was necessary and none sufficient, because all three refine the same mechanism.
A polynomial-time method must abandon smoothness sampling altogether, and
nothing in the Pascal/AKS setting suggests what would replace it.

## Theorem 18 — the geometry: Lucas drawn

Plotted, Pascal's triangle mod a prime `p` is a Sierpinski gasket of scaling
ratio `p`. This is Lucas' theorem restated: the entry at `(i, j)` survives mod
`p` exactly when every base-`p` digit of `j` is dominated by the corresponding
digit of `i`, which is the gasket's construction rule.

Two exact consequences:

**(a) Row count.** The number of survivors in row `i` is `prod (d + 1)` over the
base-`p` digits of `i`. Hence row `i` contains a zero iff
`prod (d + 1) < i + 1`.

**(b) Box count.** The number of survivors in rows `0 .. N-1` is computable by a
digit dynamic program in `O(log_p N)` time — without drawing the triangle:

```
support_count(N, p) = sum over digit positions t of
       [ prod_{s>t} (n_s + 1) ] * [ n_t (n_t + 1)/2 ] * (p(p+1)/2)^t
```

for `N` with base-`p` digits `n_s`. Specialising to `N = p^k` gives the
self-similarity relation

```
support_count(p^k, p) = (p(p+1)/2)^k,
```

so the gasket has box dimension `log(p(p+1)/2) / log p`, which increases with
`p` toward `2`. Both statements are verified exhaustively against brute force in
`tests/test_fractal.py`.

For composite `n = pq` the picture is two gaskets superimposed at scales `p` and
`q`. An entry vanishes mod `n` only when it vanishes mod both, so the zeros of
the composite are exactly where the two hole systems coincide.

## Theorem 19 — the dual of Theorem 4

Theorem 4 reads *along* row `n` and returns the **smallest** prime factor.
Reading *down* the rows returns the **largest**.

Let `n = pq` with `p < q` primes. Then the least row of Pascal's triangle
containing an entry `== 0 (mod n)` is the least `i` with

```
prod (digits_p(i) + 1) < i + 1     and     prod (digits_q(i) + 1) < i + 1,
```

i.e. the least row carrying a hole in *both* gaskets. In particular:

1. no row below `q` contains a zero mod `n`;
2. the first zero row equals `q` exactly when row `q` carries a hole mod `p`.

*Proof.* An entry vanishes mod `n` iff it vanishes mod `p` and mod `q`, and by
(a) above each vanishing is the digit condition for that prime. For `i < q` the
index has a single base-`q` digit `i`, so `prod(digits_q + 1) = i + 1` and row
`i` has no hole mod `q` — giving (1). At `i = q` the digits are `(0, 1)`, so
`prod = 2 < q + 1` and the whole interior vanishes mod `q`; a zero mod `n` then
appears iff row `q` also has a hole mod `p`, which is (2). ∎

Verified exhaustively: the criterion predicts the first zero row on **127/127**
semiprimes below 2000, and the row equals `q` in 91% of them — precisely the
cases satisfying (2)
([exp14](../experiments/results/exp14_fractal.md)).

### Cost, and why this is not a shortcut

Theorem 4 costs `Theta(p)` — a walk along one row. Theorem 19 costs
`Theta(q^2)` — the rows must be built to reach row `q`, and each is `O(q)` long.
The two-dimensional picture carries strictly more information (the largest prime
factor as well as the smallest) and it is priced by **area** rather than length.

The same holds for the cheapest sideways probe. Reading down column `c`, Lucas
makes `p | C(i,c)` a condition on `i mod p^k`, so `gcd(C(i,c) mod n, n)` splits
`n` whenever exactly one prime divides. Widening `c` raises the hit rate and the
cost of one evaluation in exactly the same proportion: measured, the ratio
`hit-rate / cost` sits at `~1/p` across three orders of magnitude of `p` and two
of `c`. **The barrier is isotropic** — every direction through the triangle
costs `Theta(p)`.

### What the picture explains

The fractal reading is not a new attack. It is the reason the earlier attacks
failed, in a form one can see:

* **Theorem 2** — the row is empty below `spf(n)`: the top of the gasket is
  solid, holes begin only at scale `p`.
* **Theorem 7** — aliasing: one cannot see a scale-`p` fractal by sampling below
  scale `p`; folding at level `r < p` averages over whole self-similar cells.
* **Proposition 12** — period finding: the gasket *is* a periodic structure of
  period `p`.
* **Proposition 17** — the smoothness ceiling: every group order in the ladder is
  an arithmetic shadow of that same scale.

A barrier one can see is easier to attack than one that can only be computed,
which is the case for keeping the geometry in view even though it broke nothing.

## Theorem 20 — the factorial threshold, and Strassen from the gasket

For every `n >= 2` and `i >= 1`,

```
gcd(i! mod n, n)  =  prod over p^e || n  of  p^min(e, v_p(i!)),
```

with `v_p(i!) = (i - s_p(i))/(p - 1)` by Legendre. In particular

```
gcd(i! mod n, n) > 1     <=>     i >= spf(n),
```

a **monotone** threshold. Verified on 29,155 `(n, i)` pairs against Legendre's
formula, exactly.

*Proof.* `gcd(a mod n, n) = gcd(a, n)`, and `gcd(i!, n)` is the product over
primes of `p` to the least of the two valuations. The threshold statement is the
case `e >= 1`: some prime of `n` divides `i!` iff some prime of `n` is `<= i`. ∎

**The fractal reading.** By Theorem 18 the mod-`p` gasket first develops holes at
row `p`. So

> *which gaskets have started making holes by row `i`* = *which primes divide
> `i!`* = `gcd(i! mod n, n)`,

and the least row where anything at all has happened is `spf(n)`. The jump rows
of the gcd are the rows where some `v_p(i!)` crosses `v_p(n)`; for squarefree `n`
they are exactly the prime factors, each entering at its own row.

### Consequence: the geometry derives the best known deterministic bound

1. The predicate `gcd(i! mod n, n) > 1` is monotone with threshold `spf(n)`.
2. So `O(log n)` binary-search steps locate `spf(n)` exactly.
3. `i! mod n` costs `O~(sqrt(i))` (Bostan-Gaudry-Schost): build
   `f(X) = (X+1)...(X+c)` with `c = isqrt(i)` and multipoint-evaluate it at
   `0, c, 2c, ...`.

Total `O~(n^(1/4))` — Strassen's deterministic bound, obtained by asking the
gasket the cheapest possible question rather than by construction. Implemented as
`aksfactor.fast.threshold_spf`.

### And it does not beat it

Measured, this binary-search form runs about **20x slower** than the block-scan
form of the same bound in `exp08`: binary search pays for `O(log n)` separate
factorial computations where one product tree covers the whole range. Same
asymptotics, a `log n` factor apart, and the constant is real.

Beating `O~(n^(1/4))` along this route would require `i! mod n` in less than
`O~(sqrt(i))`, which is itself a known open problem equivalent to improving
deterministic factoring. The one place the literature does better,
`O~(n^(1/5))` (Hittmeir; Harvey), reaches it by combining the factorial with
extra sieving structure — not by asking the gasket a better question.

So the fractal reading's contribution is explanatory, and precisely so: it shows
`O~(n^(1/4))` is *the geometry's own answer* to the cheapest question one can
pose, not an artifact of one clever algorithm.

## Theorem 21 — row `pq` contains row `p` at stride `q`

For primes `p < q` and `n = pq`:

```
C(pq, m q)  =  C(p, m)                       (mod pq),   0 <= m <= p;
C(pq, m p)  =  0 (mod q)  and  C(q, m) (mod p),           0 <  m <  q.
```

*Proof.* Lucas in each prime, then CRT. Modulo `q`, `pq` has base-`q` digits
`(0, p)` and `mq` has `(0, m)`, giving `C(p, m)`. Modulo `p`, the low digit of
`pq` is `0` while that of `mq` is `mq mod p != 0` for `0 < m < p`, so both sides
vanish (and `C(p,m) = 0 (mod p)` there too); `m = 0, p` agree trivially. For the
second line, `mp mod q != 0` forces `q | C(pq, mp)`, and modulo `p` the pair
`(pq, mp)` is `(q, m)` shifted by one base-`p` digit. ∎

Verified on 58,300 stride-`q` entries and 7,322 stride-`p` entries.

**Consequence (exact support count).** Every multiple of the *larger* prime is
non-zero -- those entries reproduce row `p` -- while by the second line the
multiple `mp` is non-zero exactly when `p` does not divide `C(q, m)`. Fine's
theorem counts those: if `q = sum d_i p^i` in base `p`, then row `q` has
`prod (d_i + 1)` entries not divisible by `p`. So for `0 < k < n`,

```
#{k : C(pq, k) != 0 (mod pq)}  =  (p - 1)  +  prod_i (d_i + 1)  -  2,
```

checked on all 2,264 pairs of odd primes `p < q < 6p` below 400. The second
term swings with `q mod p`: for `n = 1009 x 1013` (`q = 1·p + 4`) it is `8`, so
the support is almost exactly row `p` at stride `q`, density `~1/q`; for
`n = 1009 x 1511` (`q = 1·p + 502`) it is `1004`, density near `1/p + 1/q`.
An earlier write-up gave `~1/p + 1/q` as the density; that is the upper end of
this range, reached when `q = -1 (mod p)`, not its typical value. In every case
the count is at least `p - 1`.

## Proposition 22 — Pascal rho

Iterating a column of the triangle, `x -> C(x, k) mod N`, is a collision
method, and for `k = 2` it is exactly Pollard's rho: with `x = 2y + 1/2`,

```
C(x, 2)  =  2 (y^2 - 5/16) + 1/2 .
```

For general `k`, the heuristic rho length of a polynomial map `f` on `F_p` is
`sqrt(pi p / (2 kappa))` with `kappa = (1/p) sum_v m_v (m_v - 1)` over fibre sizes;
`kappa` counts the orbits of the fibre's Galois group on ordered pairs of
distinct roots. The row reflection `C(k-1-x, k) = (-1)^k C(x, k)` decides it:

* `k = 2` and odd `k`: the reflection sends the fibre over `v` to the fibre over
  `-v`, imposes nothing inside a fibre, and `kappa = 1` (generic `S_k`);
* even `k >= 4`: `C(x, k)` is a polynomial in `(x - (k-1)/2)^2`, the fibre
  group sits inside the hyperoctahedral group, which has two orbits on ordered
  pairs (antipodal or not), and `kappa = 2`. For `k = 4` explicitly,
  `C(x, 4) = ((u^2 - 5/4)^2 - 1)/24` with `u = x - 3/2`.

Measured `kappa` is `1.000`/`2.000` to three decimals for `k = 2..9`, and 300-walk
rho constants match `sqrt(pi/(2 kappa))` within one standard error (`exp25`). The
even columns' `sqrt 2` fewer steps costs a second squaring per step, so per
collision they are no cheaper. This is the Brent-Pollard `x^(2^j) + c`
phenomenon, here forced by the symmetry of the triangle for every `p`.

## Proposition 23 — the divisor predicate under the hyperbola

For `X <= sqrt(N - 1)`,

```
#{d | N : d <= X}  =  sum_{y <= X} floor(N/y)  -  sum_{y <= X} floor((N-1)/y),
```

since `floor(N/y) - floor((N-1)/y) = [y | N]`. The left side is monotone in `X`
and jumps at `X = spf(N)`, so it is a binary-searchable predicate. Better, no
search is needed: `{xy > N-1}` is `{xy > N}` plus the divisor points
`(N/d, d)`, which lie on the strictly convex curve `xy = N` and are therefore
extreme points of `{xy >= N}` -- *vertices* of the hull of the lattice points of
`{xy > N-1}`. One hull walk, testing `xy = N` at each point, finds them.

Each sum is computed exactly from the convex hull of `{(x, y) in Z^2 : xy > m}`.
That region is convex, so no lattice point lies strictly between its hull and
the curve, and the number of points on or under the curve in row `y` is
`ceil(hull_x(y)) - 1`. An edge from `(x, y)` along the primitive vector
`(dx, -dy)` covers rows `y-dy .. y-1` and contributes
`dy x + (dy+1)(dx-1)/2` points. Walking from row `sqrt N` to row `N^(1/3)` takes
`O(N^(1/3) log N)` steps with a Stern-Brocot stack (Vinogradov's bound; the
walk is Sladkey's), and rows below `N^(1/3)` are summed directly. The factoring
walk stops at `p`, so it costs about `3 N^(1/3)` steps per factor of two in
`sqrt(N)/p` -- `O~(N^(1/3))` in the worst case (`exp26`).

At the divisor vertex `(q, p)` the two hull edges `(dx, dy)` are Farey
neighbours bracketing `p/q`, and each gives the identity
`(dy q + dx p)^2 - 4 (dx dy) N = (dy q - dx p)^2` that Lehman's method searches
for (60 of 60 vertices, 120 of 120 edges). The two methods enumerate these
certificates differently -- two thirds of the hull's lie outside Lehman's
search box -- but they are the same objects at the same exponent.

**Curved pieces (heuristic).** Each hull edge is an exact linear piece of
`floor(N/y)`. On `[x, 2x]`, a degree-`d` Taylor piece of length `h` misses
`N/y` by about `N h^(d+1) / x^(d+2)`; allowing `O(1)` stray lattice points per
piece gives `h ~ x N^(-1/(d+2))`, hence `N^(1/(d+2))` pieces per dyadic block.
`d = 2, 3, 4` would give `N^(1/4)`, `N^(1/5)`, `N^(1/6)`, provided (i) exact
degree-`d` floor sums over a piece cost polylog time and (ii) the stray points
near each arc can be located exactly. Already for `d = 2`, (i) meets class
numbers: for primes `p = 3 (mod 4)`, `p > 3`,

```
sum_{k<p} floor(k^2/p)  =  (p-1)(2p-1)/6  -  (p-1-2 h(-p))/2 ,
```

from `sum_k (k^2 mod p) = 2 sum_(QR r) r` and Dirichlet's
`h(-p) = -(1/p) sum a (a/p)` (checked for all 154 such primes below 2000). The
hyperbola needs short arcs rather than full periods, so this does not prove (i)
false; it places it next to a problem with no known polynomial-time algorithm.

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
