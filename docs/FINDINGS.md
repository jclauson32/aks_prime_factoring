# Findings

What was conjectured, what turned out to be true, and what it costs.

## The starting conjecture

> Take the rows of Pascal's triangle and reduce row `n` mod `n`. The remainders
> can very often — maybe always — be broken into `x * y` where `x` is a prime
> factor of `n`.

## Verdict

**The structural half is a theorem, and it is stronger than conjectured. The
"prime" half needs one correction. The computational payoff is real but bounded
by trial division.**

### 1. The split always exists (Theorem 1)

For every interior position, `C(n,k) mod n = x * y` with

```
x = n / gcd(n, k)
```

and `x` is *always* a divisor of `n`, never trivially `1` or `n`. Verified at
**204,464** interior positions across all rows `n < 700`, with zero
counterexamples ([exp01](../experiments/results/exp01_row_structure.md)).

Stronger still: whenever a residue is non-zero, `gcd(residue, n)` is a
non-trivial proper divisor of `n`. A non-zero remainder is *always* a factor,
handed over directly.

### 2. "`x` is prime" — the correction

`x` is a *divisor*, prime only sometimes: **36.0%** of non-zero residues below
`n = 600` ([exp03](../experiments/results/exp03_gist_stats.md)). For example
`n = 1024, k = 2` gives `x = 512`.

But at the position that matters — the **first** non-zero — `x = n/spf(n)`, and
that is prime exactly when `n` is a product of two primes. So for semiprimes,
the RSA-shaped case, the conjecture is exactly right:

> for `n = p*q`, the first non-zero remainder in row `n` is **exactly `q`**.

### 3. The sharpening nobody asked for (Theorems 3 and 4)

The residues are not merely "divisible by some factor". At a prime `p | n` the
value is pinned:

```
C(n, p) mod n  =  n/p     exactly.
```

And therefore:

> **The first non-zero entry of row `n` mod `n` sits at position `spf(n)` and
> equals `n / spf(n)`. Position times value is `n`.**

Prime rows have an identically zero interior, so this doubles as a primality
test. Verified for every `n < 3000`, composite and prime
([exp02](../experiments/results/exp02_first_nonzero.md)).

### 4. Random access works, and it is fast

Theorem 5's closed form needs `C(n-1,k-1)` only modulo the *small* number `k`,
which Lucas and Granville deliver in time polynomial in `log n`. So single
residues are reachable for `n` of any size:

| digits of `n` | time to evaluate `C(n,p) mod n` | result |
|---|---|---|
| 26 | 33 µs | the cofactor, exactly |
| 161 | 42 µs | the cofactor, exactly |
| 315 | 52 µs | the cofactor, exactly |

A 300-digit cofactor falls out of one binomial residue in microseconds —
**provided you already know which `k` to ask about.**

## The barrier

That caveat is the whole story. Finding the right `k` is the hard part, and
every route measured here costs `Theta(spf(n))`.

### Why scanning cannot win

Theorem 2 says a residue can only be non-zero where `gcd(k,n) > 1`. So
"scan positions until a residue fires" is, provably, the same search as "scan
positions until a gcd fires" — trial division. The theorem that makes the
method elegant is the same theorem that caps its speed.

### Why folding cannot win either

The one route that avoids reading positions one at a time is to fold the whole
row into `r` buckets: compute `(x+a)^n mod (n, x^r - 1)` in `O(log n)` ring
multiplications — the actual AKS object. This **does** leak factors. On a first
pass it factored 24 of 25 random semiprimes with `r < 40`, which looks like a
break.

It is not. Theorem 6 explains it: mod `p`, Frobenius collapses the object to
`(x^(p^v) + a)^s`, so the fold coefficient `S_j` is a binomial progression sum
of the cofactor mod `p`. It leaks a factor exactly when that sum vanishes mod
`p` — a `~1/p` event. The early success was an artifact of testing on small
`p`.

Measured ([exp04](../experiments/results/exp04_fold_scaling.md),
[exp05](../experiments/results/exp05_multibase.md)):

- coefficients inspected before a hit, divided by `p`: median **0.87**, no
  trend in `p` across three orders of magnitude;
- observed per-coefficient hit rate tracks `1/p` across twelve AKS bases `a`.

Extra bases buy more lottery tickets, never a better price per ticket.

### The unified statement

Modulo a prime factor `p`, the AKS object **is** `(x^(p^v) + a)^(n/p^v)`. Its
first non-trivial structure lives at degree `p^v`:

- **truncate** below that degree and every coefficient is zero — you learn only
  "no factor this small";
- **fold** below that degree and the structure smears across all `r` classes,
  leaving each a `1/p` chance of collapsing.

Either way, `Theta(p)`. Detecting that the AKS identity fails is cheap;
localising *where* it fails is what costs, and that is exactly the gap between
primality testing and factoring.

## Round 2: four more theories tried, and one that worked

### Dead: symmetric aggregates (norms, resultants)

A fold coefficient is a `1/p` lottery ticket. The **norm**
`prod_t F(w^t)` aggregates all `r` of them into one number — apparently buying
every ticket at once. It is worthless, and provably so:

```
Norm( (x+a)^n )  =  (a^r + 1)^n   (mod n)        [Theorem 8, odd r]
```

The right-hand side is a function of `n`, `r`, `a` alone, so it is identical
modulo every prime factor and the gcd is always `n`. Verified for every
`(n, r, a)` tested. The general form: **anything Galois-symmetric in the `r`-th
roots of unity is Frobenius-invariant mod every `p | n`, hence carries no
factoring information.** Symmetric aggregation destroys the very asymmetry a
factor consists of.

### Dead: twisted moduli `x^r - c`

Weighting the wrapped terms by powers of `c` gives fresh coefficients, so
perhaps some are biased toward vanishing. Measured hit rates for `c != 1` track
`1/p` exactly as `c = 1` does. Different tickets, same price.

### Dead — but only halfway: the AKS-ring generalisation of Pollard `p-1`

In `F_p[x]/(x^r - 1) = prod_i F_{p^{k_i}}`, the order of `x + a` divides
`lcm_i (p^{k_i} - 1)`. Classic `p-1` is the `k = 1` case, so more `k` values
looked like strictly more chances. For that element it is strictly *fewer*:

> `p - 1` divides `p^k - 1` for every `k`, so `p^k - 1` is `B`-smooth only if
> `p - 1` already was.

**This was written up as "strictly dominated". That was too strong — see round 3
below, which refutes it.** The argument is sound about the *full unit group*, and
`x + a` sits in the full unit group. It says nothing about subgroups.

### Dead, and this was the one that would have been polynomial: reading position out of the fold

If the residue classes mod `r` were *occupied differently* by the row's
support, the class index would leak `p mod r`, and CRT over a few small `r`
would reconstruct `p` in polylog time. Occupancy is flat — by theorem, not by
luck:

> **Theorem 7.** If `gcd(r, p) = 1` then `i -> i*p^v mod r` is a bijection of
> `Z/r`, so the multiples of `p^v` fall into every class equally.

Positional information needs `gcd(r, p) > 1`, i.e. `r >= p`. This is a Nyquist
aliasing bound: a period-`p` signal cannot be resolved into fewer than `p` bins.

### Worked: `O~(n^(1/4))` instead of `O~(n^(1/2))`

Theorem 2 turns "find the first non-zero residue" into "find the first `k` with
`gcd(k,n) > 1`". Read one position at a time, that is trial division. But a
block of `c` positions can be tested with a *single* gcd, using

```
f(X) = (X+1)(X+2)...(X+c),      gcd( f(i*c) mod n, n ) > 1
```

exactly when some position in `(i*c, i*c + c]` shares a factor with `n`. One
degree-`c` polynomial evaluated at `c` points is `O~(c)` ring operations via a
product tree and a remainder tree — so `c^2` positions cost `O~(c)`. Taking
`c = n^(1/4)` searches the whole range in `O~(n^(1/4))` operations.

Implemented in [`aksfactor/fast.py`](../aksfactor/fast.py), verified against
trial division exhaustively. Measured crossover against a tuned 2-3-5 wheel
lands near `p ~ 7e8` (`n ~ 1e18`), after which it pulls away
([exp08](../experiments/results/exp08_quartic.md)).

Two honest caveats. The wall-clock exponent is `~0.85`, not `0.5`, because
CPython multiplies big integers with Karatsuba rather than FFT — the `O~(n^(1/4))`
claim is about *ring operations*, which are counted directly. And `O~(n^(1/4))`
is Strassen's classical bound, not the deterministic record: later baby-step /
giant-step work (Hittmeir; Harvey) reaches `O~(n^(1/5))`.

### The one that changed the picture: the second n-adic digit

Every route above works at the **first** n-adic digit of the row, where
Theorem 2 has already quotiented the information away: `gcd(k,n) = 1` forces
`C(n,k) ≡ 0 (mod n)`. So look one digit deeper. Write

```
C(n,k) = n * m        (Theorem 2 guarantees n divides it)
```

From `k*C(n,k) = n*C(n-1,k-1)` we get `m = C(n-1,k-1)/k`, and since `k` is
invertible mod `n`, for each prime `p | n`

```
m ≡ 0 (mod p)   <=>   some base-p digit of k-1 exceeds that of n-1   [Lucas]
```

So `gcd(m mod n, n)` splits `n` whenever that digit condition holds for one
prime factor and not the other. For random `k` that is a **constant-probability**
event — `n` has only about three digits base `p`, so each prime independently
fails the digit test a constant fraction of the time. Measured:

| `p` | split rate | `1/p` | ratio |
|---|---|---|---|
| 37 | 0.815 | 2.7e-2 | 30x |
| 193 | 0.256 | 5.2e-3 | 49x |
| 1,571 | 0.895 | 6.4e-4 | 1,406x |
| 44,501 | 0.854 | 2.3e-5 | 37,982x |
| 131,687 | 0.312 | 7.6e-6 | 41,152x |
| 810,853 | 0.973 | 1.2e-6 | 789,230x |

The rate swings with the digit profile of `n-1` but does not decay with `p`. Every other route in this project produced a
`1/p` rate; this one is flat.

**Which makes it a hardness result, not an algorithm.**

> **Proposition 10.** If `C(n,k) mod n^2` can be computed in `poly(log n)` for
> uniform random `k`, then `n = pq` factors in `O(1)` expected iterations.
> Computing binomial coefficients modulo `n^2` is factoring-hard.

And that finally makes the barrier *sharp*. Theorem 5 gives genuine random
access to `C(n,k) mod n` — any size `n`, microseconds — precisely because the
closed form only ever needs `C(n-1,k-1)` modulo the **small** number `k`, where
Lucas and Granville apply. Proposition 10 says extending that by a single
n-adic digit would break factoring.

Two things follow:

1. The requirement in Granville-style algorithms that the modulus be *factored*
   before binomials can be reduced modulo it is **intrinsic**, not an artifact of
   how those algorithms are written.
2. The factors are not hidden in the Pascal row. They sit in the second n-adic
   digit at constant density — a handful of random probes would find one. The
   entire obstruction is that reading that digit *is* the problem being solved.

## The mechanism all the dead ends share

Modulo a prime factor `p`, the AKS object **is** `(x^(p^v) + a)^(n/p^v)`, a
series whose only distinguishing feature is a period of `p^v`. Everything cheap
you can compute from it is either

- **symmetric** in the prime factors — the same mod every `p`, so the gcd is `n`
  (norms, resultants, the `p^k - 1` family); or
- **aliased** below the period — the period-`p` structure spreads evenly over
  all `r < p` buckets (folds, twists, class occupancy).

So at the first digit, extracting `p` is **period finding**, and classically
that costs `Omega(p)` samples. It is precisely the problem Shor's algorithm solves in polylog time
quantumly, by transforming at size `n` instead of size `r << p`. What is missing
is not a cleverer `r`, `a`, or aggregation rule — it is resolution.

## What would actually constitute a breakthrough

Not proved impossible — just measured, and localised. Beating this needs one of:

1. **An algebraic reason for the progression sums of Theorem 6 to collapse** —
   some `(r, j, a)` chosen from `n` alone that forces `S_j ≡ 0 (mod p)` far more
   often than `1/p`.
2. **An operator that isolates the degree-`p^v` structure** of
   `(x^(p^v)+a)^s` without enumerating `Theta(p)` coefficients. Differentiation
   is the obvious candidate and it fails twice over: `d/dx` does not commute
   with reduction mod `x^r - 1`, and under truncation it needs degree `> p`
   before it sees anything at all.
3. **Sparse interpolation of the row.** The row mod `n` is sparse — its support
   is contained in the multiples of the prime factors, density `~1/p + 1/q` —
   and evaluations `(1+z)^n mod n` are available in `O(log n)` each. But
   Prony/Ben-Or–Tiwari needs about as many evaluations as there are non-zero
   terms, and for `n = pq` that count is `p + q - 1 ≈ sqrt(n)` at best. Same
   wall, different shape.

Each of these has a concrete experiment harness in `experiments/`. The
negative results are the contribution: they say precisely which door to try
next, and why the obvious ones are shut.

A fourth candidate, added after round 2: **break the symmetry deliberately.**
Theorem 8 says symmetric aggregates are useless and Theorem 7 says positions are
aliased, but neither forbids an *asymmetric* statistic computable from `n` alone
— something playing the role smoothness of `p-1` plays for Pollard.

Round 3 found one: the norm-one construction of Theorem 13. It is a genuine
asymmetric statistic and it factors numbers Pollard `p-1` cannot. It is also not
enough, and the counting in round 3 says exactly why: the statistic is *rigid*,
one order per `(p, d)`. So this door is open but narrow, and the requirement
behind it has been restated below in the form that matters — a group order that
**varies** at fixed `p`.

## Round 3: the correction, and a real asymmetric statistic

Round 2 closed by naming the open problem: *find a statistic of the first digit
that is asymmetric in the prime factors*. Round 3 found one — by noticing that
round 2's own domination argument had a gap.

### The gap

| group | order | divisible by `p-1`? |
|---|---|---|
| full unit group of `F_{p^d}` | `p^d - 1` | yes — hence dominated |
| **norm-one subgroup** `ker(N: F_{p^d}* -> F_p*)` | `(p^d - 1)/(p - 1)` | **no** |

The domination argument covers the full unit group, and `x + a` lives there. It
says nothing about the norm-one subgroup, whose order `p - 1` does not divide.

And you can land inside it **without knowing `p`**: choose a monic polynomial
whose roots multiply to `1`. For `d = 2` that is `x^2 - a x + 1`, whose orbit is
tracked by the Lucas sequence `V_k(a,1)` — which is Williams' `p+1` method,
recovered as a special case rather than imported.

### Measured

Primes with `p-1` rough and `p+1` smooth are not exotic; sampling smooth `p+1`
candidates turns them up readily. Against a cofactor rough on both sides, at an
identical budget of `bound = 500`:

| p | largest prime factor of `p-1` | of `p+1` | Pollard `p-1` | norm-one |
|---|---|---|---|---|
| 4,512,218,267 | 3,716,819 | 47 | fail | **found** |
| 7,205,934,499 | 171,569,869 | 41 | fail | **found** |
| 2,850,909,643 | 41,669 | 53 | fail | **found** |
| 1,350,354,433 | 125,591 | 41 | fail | **found** |
| 34,886,131,169 | 88,339 | 53 | fail | **found** |
| 2,783,306,129 | 17,623 | 53 | fail | **found** |

**6 of 6.** Round 2's claim is refuted, and `experiments/exp07_dead_ends.py`
carries the retraction inline rather than quietly dropping the entry.

One methodological note, because it nearly produced a second wrong conclusion:
base retries are mandatory. For `d = 2` the root only lands in `F_{p^2}` when
`base^2 - 4` is a non-residue mod `p`; otherwise it falls into `F_p` and the
method silently degenerates to `p-1`. That is a coin flip per base. The first
run of this experiment tried seven bases that all happened to be residues and
recorded a failure that was not real.

### The statistic is still a binomial sum

```
V_k(a,1) = sum_j (-1)^j * (k/(k-j)) * C(k-j, j) * a^(k-2j)
```

Verified on thousands of cases. So what finally separates `p` from `q` *is* a
binomial sum — just not one along row `n`. The change is not abandoning
binomials; it is constraining the roots to multiply to `1`, which pins the order
to `p+1` instead of `p-1`.

### Why it still is not polynomial — the next obstruction, stated

| method | group | order | varies at fixed `p`? |
|---|---|---|---|
| Pollard `p-1` | `F_p*` | `p - 1` | no |
| norm-one, degree `d` | `ker N` in `F_{p^d}*` | `Phi_d(p)`-ish | no |
| ECM | `E(F_p)` | `p + 1 - t` | **yes**, one per curve |

Every method reachable in this framework has a **rigid** group order: fix `p` and
the number whose smoothness decides success is fixed too. One ticket per `p`, and
if `Phi_d(p)` is rough for all small `d`, you are stuck. ECM's advantage is not a
better group but a *family* of groups — a fresh order per curve, at the same `p`.
That is what buys `L[1/2]` instead of dependence on smoothness luck.

### Counted, not asserted

That difference is countable, so round 3 counted it. Distinct group orders
reachable at one fixed `p`, over 400 random samples each (found / theoretical
cap):

| p | ring `d=2` | ring `d=3` | ring `d=4` | distinct `#E(F_p)` | Hasse width `4*sqrt(p)` |
|---|---|---|---|---|---|
| 211 | 2 / 2 | 3 / 3 | 5 / 5 | 55 | 56 |
| 503 | 2 / 2 | 3 / 3 | 5 / 5 | 84 | 88 |
| 1,009 | 2 / 2 | 3 / 3 | 5 / 5 | 109 | 124 |
| 2,003 | 2 / 2 | 3 / 3 | 5 / 5 | 140 | 176 |

The ring columns **saturate at `partitions(d)` and stop** — 2, 3, 5 — no matter
how large `p` grows or how many `f` are tried. The reason is structural, not a
sampling artifact: `F_p[x]/f` for squarefree `f` is a product of finite fields
`F_{p^{d_1}} x ... x F_{p^{d_m}}`, so its unit group order is
`prod_i (p^{d_i} - 1)` — a function of the degree *partition* alone. Two
polynomials with the same degree pattern give literally the same order. There is
nothing left to vary. Restricting to Theorem 13's norm-one subgroup changes
*which* rigid number you get, never that it is rigid.

The elliptic column grows with `sqrt(p)`, tracking the Hasse interval.

So if `Phi_e(p)` is rough for every `e <= d`, every construction in this
repository fails at that `p`, and no additional work at that `p` helps — the
finitely many available orders are exhausted. Raising `d` costs more while
supplying only `partitions(d)` new orders.

**The gap, as a requirement.** A polynomial-time method in this family would need
a construction over `Z/n` whose group order at fixed `p` varies with a parameter
we control. Quotients of `(Z/n)[x]` cannot supply it — their orders are pinned to
`prod (p^{d_i} - 1)` by the degree partition. Escaping needs a genuinely
different algebraic group, which is what an elliptic curve is, and what no
rearrangement of Pascal's triangle will produce.

Two known constructions do satisfy the requirement, which is worth saying so the
target does not read as a fantasy: **elliptic curves** (a fresh order `p + 1 - t`
per curve, giving ECM) and **class groups of imaginary quadratic orders** (a
fresh order `h(-D)` per discriminant `D = -kn`, giving Schnorr-Lenstra). Both
leave the world of polynomial quotients of `Z/n` entirely. That is the price of
a varying order, and nothing in the Pascal/AKS setting pays it.

## Round 4: deforming the period

Round 3 ended asking for a construction whose parameter *varies* at fixed `p`.
Round 4 found one, inside Pascal's triangle itself — and then found out why
varying is not enough.

### The q-deformation

Replace `C(n,k)` by the Gaussian binomial `[n,k]_q`. Divisibility by `p` is then
governed not by `p` but by `d = ord_p(q)`, which **moves as the base `q` moves**:

```
p | [n,k]_q   <=>   (k mod d) > (n mod d)                 [clause 1]
               or    p | C(floor(n/d), floor(k/d))         [clause 2]
```

a q-analogue of Kummer's theorem, verified exhaustively against the explicit
q-Pascal row. Clause 1 fires at `k = (n mod d) + 1`, which can sit far below
`spf(n)` — **inside the region Theorem 2 seals off**:

| n | factors | spf(n) | classical row below spf | q-row first hit | base |
|---|---|---|---|---|---|
| 143 | 11 x 13 | 11 | 0 everywhere | k = 4 | 2 |
| 221 | 13 x 17 | 13 | 0 everywhere | k = 2 | 4 |
| 437 | 19 x 23 | 19 | 0 everywhere | k = 6 | 2 |
| 899 | 29 x 31 | 29 | 0 everywhere | k = 4 | 2 |

(Bases with `gcd(q-1,n) > 1` are excluded — those have `ord_p(q) = 1` and a
single gcd would already have split `n`.)

That is a real deformation of the barrier. The classical period `p` is rigid;
this one is tunable.

### And a first guess that was wrong

The obvious conclusion — "`ord_p(q)` divides `p-1`, so this is just Pollard
`p-1`" — is **false**, and the measurement says so. With `p-1` smooth, a random
base *still* has order close to `p-1`:

| p (smooth p-1) | largest prime factor of p-1 | `ord_p(q)/(p-1)` for q = 2..8 |
|---|---|---|
| 1,193,011 | 23 | 1.00, 0.20, 0.50, 0.50, 0.17, 0.33, 0.33 |

Smooth `p-1` supplies many small *divisors* but almost no elements of small
*order*. Clause 1 needs a small order, not a smooth one — so the raw
q-deformation is strictly **weaker** than Pollard, not equal to it.

### The dichotomy that actually closes it

Tune the base: put `Q = q^M` with `M` a smooth prime-power ladder to bound `B`.
Then `ord_p(Q)` is the `B`-rough part of `ord_p(q)`, so it is `1` when `ord_p(q)`
is `B`-smooth and `> B` otherwise. Measured over 120 random pairs:

| ladder bound B | tuned period = 1 | 1 < period <= B | period > B |
|---|---|---|---|
| 50 | 0 | **0** | 60 |
| 200 | 2 | **0** | 58 |

**The middle column is empty.** And `period = 1` means `Q ≡ 1 (mod p)`, so
`gcd(Q-1, n)` has already split `n` — Pollard `p-1` verbatim, q-row contributing
nothing. So tuning either lands exactly on Pollard, or leaves clause 1 as a
`1/B` lottery. Measured work stays `Theta(p)`.

### Also closed: every other row

While here, one more family died in a single line.

> **Theorem 15.** For `k < spf(n)`, `C(N,k) mod n` depends only on `N mod n`.

Because `k!` is invertible there, `C(N,k) ≡ prod_i((N mod n) - i) * (k!)^(-1)`.
So rows `N` and `N'` congruent mod `n` are *identical* below `spf(n)`: row `2n`,
row `n+1`, row `n^2`, any row at all, all carry the same nothing. Theorem 2 is
just the case `N = n`, where the shared value happens to be zero. Verified on
thousands of (row, position) pairs.

### The lesson

| construction | governing quantity | varies at fixed `p`? | over what set |
|---|---|---|---|
| classical Pascal row | period `p` | no | — |
| **q-Pascal row** | period `ord_p(q)` | **yes** | divisors of `p-1` |
| norm-one subgroup, degree `d` | order `Phi_d(p)` | no | — |
| elliptic curve | order `p+1-t` | yes | interval of width `4 sqrt(p)` |

Round 3 asked for a parameter that varies. Round 4 supplied one and learned that
this was the wrong thing to ask for:

> **It is not variation that matters, but the density of the set varied over.**

The divisors of `p-1` are sparse and structured, and the elements realising the
small ones are vanishingly rare — hence the dichotomy. An elliptic curve's order
ranges over an *interval*: dense, so every fresh curve is a genuinely fresh
number and sampling until smooth is well behaved. That density is what separates
`L[1/2]` from `Theta(p)`.

**Open, for round 5.** A construction over `Z/n` whose governing quantity at
fixed `p` ranges over a *dense* set — an interval, not a divisor lattice — while
staying computable without knowing `p`.

## Round 5: the density requirement is met — and it caps at L[1/2]

Round 4's target was precise: a construction over `Z/n` whose governing quantity
at fixed `p` ranges over a **dense** set, while staying computable without
knowing `p`. Class groups of imaginary quadratic orders meet it.

For discriminant `D = −kn`, the class number `h(D)` sits near `√|D|` and moves
essentially arbitrarily with `k`. Ambiguous forms — those of order dividing 2 —
*are* factorizations of `D`, so a smooth `h(D)` hands back a factor of `n`. That
is Schnorr–Lenstra, implemented here in `aksfactor/classgroup.py` so the contrast
is measurable rather than asserted.

**Verified.** Identity, inverse, `f^h = 1` and associativity over 120
discriminants and 512 forms: zero failures. It splits 24/24 test semiprimes
including 1009·1013.

One bug is worth recording because casual testing misses it entirely: an
extended-gcd returning a **negative** gcd makes the change-of-basis matrix have
determinant `−1`. That is an *improper* equivalence — it lands in a different
class, so composition silently computes in the wrong group. Every structural
check failed until the sign was normalised.

**Dense, measured.**

| n | multipliers tried | distinct `h(−kn)` | distinct share | range of h |
|---|---|---|---|---|
| 143 | 40 | 17 | 42% | 4..124 |
| 1,147 | 40 | 22 | 55% | 6..392 |
| 5,183 | 40 | 26 | 65% | 24..504 |

against the rigid families at fixed `p = 1009`: ring unit orders give **2, 3, 5**
distinct values for degree 2, 3, 4 (the `partitions(d)` cap), while elliptic
curves give 100 and class numbers give a fresh number for most `k`.

### The ceiling

Meeting the density requirement does not give polynomial time. It gives
`L[1/2]` — the same class as ECM. The reason is structural and is the real
conclusion of this project. Every method in the ladder is a **smoothness
lottery**: build a group whose order is some integer attached to `p`, then hope
that integer is `B`-smooth.

| the set the parameter ranges over | example | what you get |
|---|---|---|
| a single value | classical Pascal row: period `p` | `Θ(p)` — trial division |
| divisor lattice of `p−1` | q-Pascal row; Pollard `p−1` | fast only when `p−1` is smooth |
| one value per `(p, d)` | norm-one subgroup: `Φ_d(p)` | Williams `p+1` and relatives |
| `partitions(d)` values | ring unit groups | a constant, independent of `p` |
| **dense interval near `p`** | ECM: `p + 1 − t` | `L[1/2]` |
| **dense, near `√(kn)`** | class groups: `h(−kn)` | `L[1/2]` |

Density separates the last two rows from the rest. It is still not enough,
because the probability that a number of size `p` is `B`-smooth is itself
governed by the Dickman function, and optimising `B` against sampling cost gives
`L[1/2]` however good the sampling is. **Denser sampling cannot beat the
smoothness density itself.**

That is also why the number field sieve is faster: `L[1/3]` comes not from
sampling better but from making the numbers *smaller* — testing smoothness of
algebraic norms rather than of numbers of size `p`. It changes *what* is tested,
not how it is sampled.

### Where this leaves it

Five rounds walked a ladder of requirements, each answered, each answer exposing
the next constraint:

1. *An asymmetric statistic* — supplied by the norm-one subgroup (round 3).
2. *A parameter that varies* — supplied by the q-deformation (round 4).
3. *Variation over a dense set* — supplied by class groups (round 5).

Each was necessary; none was sufficient. All three refine the **same** mechanism
— sampling for smoothness — whose ceiling is `L[1/2]` regardless of refinement.
The next requirement is not a better group but a different mechanism, and that
is where this line of attack, which began with a conjecture about remainders in
Pascal's triangle, honestly runs out.

What the Pascal framing did produce is exact and worth keeping: the first
non-zero entry of row `n` mod `n` sits at `spf(n)` and equals `n/spf(n)`; the
factors live in the second n-adic digit at constant density; and reading that
digit is provably as hard as factoring.

## Round 6: the picture

Pascal's triangle mod `n`, plotted, is a **Sierpiński gasket** — and that is not
an analogy. It *is* Lucas's theorem drawn: the entry at `(i,j)` survives mod `p`
exactly when every base-`p` digit of `j` is dominated by that of `i`, which is
the gasket's construction rule.

For composite `n = pq` the picture is **two gaskets superimposed**, at scales `p`
and `q`. An entry vanishes mod `n` only when it vanishes mod both — so the zeros
of the composite are exactly where the two hole systems *coincide*:

```
                       #
                      ##
                     ###
                    #pp#
                   ##p##
                  #qqqq#
                 #p0q0p#
                ##pqqp##
               ####q####
              #pppppppp#
             #q000p000q#
            ##q00pp00q##
           #ppq0ppp0qpp#
          ##p#qppppq#p##
         ######ppp######
        #00q0pq00qp0q00#
       ##0qqp#q0q#pqq0##
      ###qq###qq###qq###
     #ppp0ppppqpppp0ppp#
    ##ppppppp##ppppppp##
   #qq00p000q#q000p00qq#
  #p0q0pp00qppq00pp0q0p#
 ##pqqppp0q#p#q0pppqqp##
####q#pppq####qppp#q####
```

`p` = killed by 3 alone, `q` = killed by 5 alone, `0` = zero mod 15 (both),
`#` = survivor.

### What the geometry gives

**Exact box counting (T18).** The survivors in rows `0..N−1` are counted by a
digit dynamic program in `O(log N)` — without drawing the triangle. Specialising
to `N = p^k` gives the self-similarity relation `(p(p+1)/2)^k` exactly, so the
box dimension is `log(p(p+1)/2)/log p`, rising toward 2 with `p`. Verified
against brute force on every case tested.

**A dual to Theorem 4 (T19).** Theorem 4 reads *along* row `n` and returns the
**smallest** prime factor. Reading *down* the rows returns the **largest**: the
first row containing a zero mod `n` is the first row carrying a hole in *both*
gaskets. No row below `q` can contain one, and the first zero row equals `q`
exactly when row `q` has a hole mod `p`. The criterion predicts it on **127/127**
semiprimes below 2000; the row is `q` itself in 91% of them.

### Does it break the barrier? No — and now you can see why

The cheapest sideways probe reads *down a column*. Mod `p`, Lucas makes
`p | C(i,c)` a condition on `i mod p^k`, so `gcd(C(i,c) mod n, n)` splits `n`
whenever exactly one prime divides. Widening the column raises the hit rate — and
the cost of one evaluation by exactly as much:

| p | column c | hit rate | cost/probe | rate ÷ cost × p |
|---|---|---|---|---|
| 431 | 1 | 0.0025 | 1 | 1.08 |
| 431 | 256 | 0.5450 | 256 | 0.92 |
| 2,239 | 256 | 0.1943 | 256 | 1.70 |
| 39,239 | 256 | 0.0120 | 256 | 1.84 |

Flat at ~1 across three orders of magnitude. **The barrier is isotropic** — every
direction through the triangle costs `Θ(p)`.

And the picture explains every earlier barrier at a glance:

- **T2** (row empty below `spf(n)`): the top of the gasket is *solid*; holes
  begin only at scale `p`.
- **T7** (aliasing): you cannot see a scale-`p` fractal by sampling below scale
  `p` — folding at `r < p` averages over whole self-similar cells.
- **P12** (period finding): the gasket *is* a periodic structure of period `p`.
- **P17** (smoothness ceiling): every group order in the ladder is an arithmetic
  shadow of that same scale.

The 2-D picture holds strictly *more* information than row `n` — it contains the
largest prime factor as well as the smallest — but extracting it costs **area**
where the row cost length. The information is real and priced accordingly.

The fractal didn't break the wall. It made the wall visible, which is a better
place to stand than where round 5 left off.

## Round 7: the gasket derives Strassen — and doesn't beat him

Round 6 said the mod-`p` gasket first makes holes at row `p`. That has a sharp
arithmetic shadow:

> *which gaskets have started making holes by row `i`* = *which primes divide
> `i!`* = `gcd(i! mod n, n)`

**Theorem 20.** `gcd(i! mod n, n) = ∏ p^min(e, v_p(i!))` over `p^e ‖ n`. Hence
`gcd(i! mod n, n) > 1` **iff** `i ≥ spf(n)` — a *monotone* threshold. Verified
exactly on 29,155 `(n, i)` pairs.

The jump rows of that gcd are where some `v_p(i!)` crosses `v_p(n)`. For
squarefree `n` they are **exactly the prime factors**, each entering at its own
row:

| n | jump rows | distinct primes |
|---|---|---|
| 210 | 2, 3, 5, 7 | 2, 3, 5, 7 |
| 1,155 | 3, 5, 7, 11 | 3, 5, 7, 11 |
| 1,022,117 | 1009, 1013 | 1009, 1013 |
| 12 = 2²·3 | 2, 3, **4** | 2, 3 |

(the extra jump at 4 is where `v_2(i!)` reaches 2 — same statement, multiplicity
filling.)

### The derivation

1. Row `p` is where the mod-`p` gasket first has holes (T18).
2. So `gcd(i! mod n, n) > 1` iff `i ≥ spf(n)`.
3. Monotone ⟹ `O(log n)` binary-search steps locate `spf(n)`.
4. `i! mod n` costs `Õ(√i)` (Bostan–Gaudry–Schost).

Total `Õ(n^(1/4))` — **Strassen's deterministic bound, falling out of the
geometry** rather than being imposed on it. Implemented as `threshold_spf`.

### And the honest part

Measured, this pretty derivation is **~20× slower** than the block-scan form of
the *same* bound already in `exp08`: binary search pays for `O(log n)` separate
factorial computations where one product tree covers the whole range.

| p | binary search + BGS | block scan | ratio |
|---|---|---|---|
| 1,290,539 | 3.43 s | 0.17 s | 20× |
| 24,251,417 | 28.24 s | 1.14 s | 25× |
| 212,771,729 | 109.83 s | 7.75 s | 14× |

**No new bound.** Beating `Õ(n^(1/4))` this way needs `i! mod n` in less than
`Õ(√i)` — a known open problem, equivalent to improving deterministic factoring.
The one place the literature does better, `Õ(n^(1/5))` (Hittmeir; Harvey), gets
there by combining the factorial with extra sieving, not by asking the gasket a
better question.

So the contribution is explanatory, and precisely so: `Õ(n^(1/4))` is **the
geometry's own answer** to the cheapest question you can ask the picture — *has
any hole appeared by row `i`?* — not an artifact of one clever construction.

## Round 8: Harvey's N^(1/5), implemented

Round 7 named the live thread: Hittmeir and Harvey beat the `Õ(n^(1/4))` bound
this project kept rediscovering. So I read the paper
([arXiv:2010.05450](https://arxiv.org/abs/2010.05450)) and implemented it.

Three ingredients, **none of which come from Pascal's triangle**:

1. **Lehman's strategy** — for `N = pq` there are small `a, b` with `aq + bp` in
   a short, explicitly known interval; knowing `u = aq + bp` recovers `p, q` from
   the roots of `y² − uy + abN`.
2. **Hittmeir's congruence** — `α^(aq+bp) ≡ α^(aN+b) (mod p)` by Fermat, so a
   candidate `u` is testable *modulo p without knowing p*, via a gcd.
3. **One global BSGS sweep** — writing the offset as `i + jm` turns the search
   into matching `α^(−jm)·t_{a,b}` against a table of powers. Hittmeir sweeps
   chunks; Harvey's exponential gain is sweeping the whole space at once.

**Correctness.** Algorithm 4.2 on 256 adjacent-prime semiprimes: 256/256.
Algorithm 4.3 end to end on 25 semiprimes, balanced and unbalanced: 25/25.
Primes correctly reported prime: 5/5.

**Where the sweep takes over.** Algorithm 4.3 clears factors below `(N/r)^(1/2)`
with Strassen and hands the rest to the sweep; `r ~ N^0.2/lg^0.8 N` grows, so the
sweep's share grows with `N`. At 24 bits it handles under 1% of the range; by 512
bits, almost all of it. **It only becomes an `N^(1/5)` algorithm at scale.**

**And it loses at every size I can test:**

| p = spf(N) | digits | Harvey 4.3 | Strassen | ratio |
|---|---|---|---|---|
| 5,807 | 8 | 0.000 s | 0.000 s | 3.4× |
| 31,513 | 10 | 0.014 s | 0.015 s | 0.9× |
| 328,357 | 12 | 0.678 s | 0.072 s | 9.4× |

The stated bounds cross near **180 bits (`N ~ 10^54`)** — far past anything
reachable from Python, though well *below* cryptographic sizes, so at RSA scale
the `N^(1/5)` bound genuinely is the better one. (Both hopeless in absolute terms
there; it's a comparison of two astronomical numbers.)

### What it needed, measured against this project's ladder

The gasket supplies a **scale**, `p`, and every question I asked of it across
seven rounds was some form of *at what scale does the structure change?* Harvey
asks a different question:

- Lehman's relation `aq + bp` is **additive**; every construction in this
  repository was multiplicative. Pascal's triangle knows `p` as a scale — it
  knows nothing about `aq + bp` landing in a short interval.
- Hittmeir's test is Fermat — the one classical tool that survives from the AKS
  side of this project.
- Harvey's contribution is that the resulting search space, unlike a scale, is
  **flat enough to square-root with a single BSGS**.

That is the answer to round 7's question. **A scale cannot be square-rooted** —
Theorem 7's aliasing says sampling below it returns nothing. A *list of
candidates* can be, and Lehman's theorem is what turns factoring into a list. The
`N^(1/5)` speedup is baby-step/giant-step over that list, and no rearrangement of
Pascal's triangle produces the list in the first place.

## Round 9: an anatomy of Harvey's open question

Harvey's paper ends Remark 3.4 with an explicit open problem:

> *"whether it is possible to obtain a fully square-root speedup for Lehman's
> original choice `r ≍ N^(1/3)`. This would presumably lead to a factoring
> algorithm with complexity `N^(1/6+o(1))`."*

This round does not resolve it. It **locates** it, with two measurements.

### 1. Which term binds — and it binds by exactly 2

The candidate count is `s ~ √N·lg r/(2√r·m) + r·lg r`: a `j` sweep per pair, plus
one candidate per pair with `ab ≤ r`. Measured at the optimum `r = m = N^(1/5)`:

| bits of N | j-candidates | pair-candidates | ratio | pair share |
|---|---|---|---|---|
| 60 | 1.703e+04 | 3.407e+04 | **2.0000** | 66.7% |
| 140 | 2.605e+09 | 5.210e+09 | **2.0000** | 66.7% |
| 260 | 8.116e+16 | 1.623e+17 | **2.0000** | 66.7% |

**Exactly 2, at every size** — not an artifact of one `N`. At `r = m = N^(1/5)`
the terms are `N^(1/5)ln r/2` and `N^(1/5)ln r`. So two thirds of the work is the
pair enumeration, and *that* is what pins the exponent: `s ≥ r` forces `cost ≥ r`,
while `m ~ N^(1/4)/r^(1/4)` gives `r = N^(1/5)`. Reaching `N^(1/6)` means sweeping
the `Θ(r)` pairs themselves in `O(√r)`.

### 2. The structure a sweep would need is measurably absent

BSGS needs the candidates to be a geometric progression — the exponents
`e(a,b) = aN + b − ⌈2√(abN)⌉` to be arithmetic. `√(ab)` isn't additive, but it
could be *locally* additive. On the `a = 1` slice, the first difference holds
constant only while `f'(b) = √(N/b)` moves less than 1, and
`f''(b) = −√N/(2b^1.5)`, giving

```
run length  L(b)  ≤  2·b^1.5 / √N
```

| bits of N | b up to r | longest run | mean run | bound | run length BSGS needs |
|---|---|---|---|---|---|
| 30 | 1,024 | 3 | 1.20 | 2.00 | 32 |
| 40 | 10,321 | 3 | 1.20 | 2.00 | 102 |
| 50 | 20,000 | 1 | 1.00 | 0.17 | 141 |
| 60 | 20,000 | 1 | 1.00 | 0.01 | 141 |

Measurement tracks the bound; both are `O(1)`. A run of length `L` needs
`b ~ (L√N)^(2/3)`, so `L = √r = N^(1/6)` needs `b ~ N^(4/9)` — but `b ≤ r = N^(1/3)`,
and `1/3 < 4/9`. **The runs are never long enough, by a fixed margin in the
exponent.** That's not a shortage of cleverness; the structure isn't there.

### 3. Where that points

A square-root speedup over the pairs, if it exists, cannot come from group
structure in the candidates. It would have to come from the *arithmetic of the
pairs*: `{(a,b) : ab ≤ r}` with `a/b` a convergent to `p/q` is a Farey /
Stern–Brocot fan, not a geometric progression. Whether a fan admits a square-root
sweep is a question about continued fractions — and nothing in Pascal's triangle
speaks to it.

## Round 10: the Farey fan — Lehman's covering is not redundant

Round 9 reduced Harvey's Remark 3.4 to one question: the `Θ(r log r)` pairs are
two thirds of the cost, so `N^(1/6)` needs them handled in `O(√r)`. There are
exactly two ways — **use fewer pairs**, or **process them faster**. Round 9 closed
the second. This round closes the first.

Each pair certifies `p` only inside an interval. Lehman's condition
`0 ≤ aN/p + bp − 2√(abN) < W`, multiplied by `p`, becomes
`b·p² − (2√(abN)+W)·p + aN < 0` — so certified `p` lie between the roots, centred
at `p* = √(aN/b)` where AM–GM is tight, with half-width `√N/(2b√r)` that depends
**only on b**.

So it's a covering problem.

### Lehman's theorem, verified

The intervals must cover `[√(N/r), √N)`. They do, at every size — and the total
width is only about **1.8×** the range. The system is barely thicker than it has
to be.

### And it's essentially non-redundant

Greedy interval covering is optimal, so this is the exact minimum:

| r | pairs | minimum subcover | subcover/pairs | subcover/√r |
|---|---|---|---|---|
| 100 | 246 | 178 | 0.724 | 17.8 |
| 800 | 2,755 | 1,941 | 0.705 | 68.6 |
| 3,200 | 13,199 | 9,114 | 0.691 | 161.1 |
| 6,400 | 28,589 | 19,566 | **0.684** | **244.6** |

**Flat at ~0.70 across a 64× range of `r`.** At most 30% of pairs are droppable;
the minimum subcover is `Θ(r log r)` — same order as the full set. And
`subcover/√r` climbs 18 → 245, so it is emphatically *not* `O(√r)`.

The geometry says why: half-width `√N/(2b√r)` depends only on `b`, while centres
`√(aN/b)` are spaced `≈ √N/(2√(ab))`. Width beats spacing only when `4a ≥ br`,
which at `ab = r` needs `a ≥ r/2` — a vanishing corner of the fan. Everywhere
else the intervals sit edge to edge, and dropping one opens a gap.

### Both routes to N^(1/6) are now closed

| route | what it needs | status |
|---|---|---|
| **fewer pairs** | a subfamily of size `O(√r)` covering the range | closed here — min subcover is `0.70 × pairs` |
| **faster sweep** | local arithmetic progressions in `e(a,b)` | closed in round 9 — runs are `O(1)`, need `√r` |

Neither is a shortage of ingenuity. Both are *measurable absences of structure*.

### Where the fan actually leads

The natural remaining idea is the **Stern–Brocot tree**: convergents to a fixed
`ξ = p/q` form a path of length `O(log)`, not a set of size `Θ(r)`. Descending it
needs one comparison per step — *is `a/b < p/q`?* And that comparison is

```
a/b < p/q  ⟺  aq < bp  ⟺  aN < b·p²  ⟺  p > √(aN/b)
```

— exactly a threshold query *"is p bigger than this?"*, which is **Strassen's
problem**, costing `Õ(√threshold) = Õ(N^(1/4))` in Lehman's range. The tree walk
costs more than the answer it's looking for.

That's a clean place for the thread to end. The Farey structure genuinely *does*
compress the search — `Θ(r)` candidates down to `O(log)` tree steps — but each
step is priced at exactly the barrier the compression was meant to dodge. **The
fan and the threshold oracle are the same problem in different coordinates.**

## Where the search space stands now

After two rounds the picture is no longer a list of failed attempts; it is a
dichotomy, and it covers the natural approaches:

| | first n-adic digit (`C(n,k) mod n`) | second n-adic digit (`C(n,k) mod n^2`) |
|---|---|---|
| **is it computable?** | yes, `poly(log n)` for any `n` (Theorem 5) | not known, and Proposition 10 says not without factoring |
| **does it leak factors?** | only via a `1/p` accident (Theorem 6) | at **constant** density (Theorem 9) |
| **why** | Theorem 2 has already quotiented the information away; what survives is aliased (T7) or symmetric (T8) | Lucas' digit condition separates the primes |

The first column is cheap and nearly empty. The second is rich and priced at
factoring.

Round 3 added a third column that is neither — it sidesteps the digit question
entirely by leaving row `n` and going to a group whose order is a *different*
function of `p`:

| | norm-one subgroup, degree `d` |
|---|---|
| **is it computable?** | yes, `O(log M)` ring operations |
| **does it leak factors?** | when `Phi_d(p)` is smooth — a real, common, asymmetric event |
| **why it stops** | the order is rigid: one ticket per `(p, d)`, no way to redraw |

The open problem round 2 left behind was: **find a statistic of the first digit
that is asymmetric in the prime factors.** Round 3 answered it — the norm-one
construction is exactly such a statistic — and replaced it with a harder one:
**find a family of such objects whose group order varies at fixed `p`.** That is
the property separating everything in this repository from ECM.
