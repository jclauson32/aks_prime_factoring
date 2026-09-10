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

### Dead: the AKS-ring generalisation of Pollard `p-1`

In `F_p[x]/(x^r - 1) = prod_i F_{p^{k_i}}`, the order of `x + a` divides
`lcm_i (p^{k_i} - 1)`. Classic `p-1` is the `k = 1` case, so more `k` values
looked like strictly more chances. It is strictly *fewer*:

> `p - 1` divides `p^k - 1` for every `k`, so `p^k - 1` is `B`-smooth only if
> `p - 1` already was.

The generalisation is dominated by the method it generalises. Isolating the
genuinely new cyclotomic part `Phi_k(p)` requires the norm-one subgroup —
Lucas sequences, Williams `p+1` — which is not reachable inside this ring.

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
— something playing the role smoothness of `p-1` plays for Pollard. Every such
statistic tried here reduced to smoothness of `p^k - 1` and collapsed.

## Where the search space stands now

After two rounds the picture is no longer a list of failed attempts; it is a
dichotomy, and it covers the natural approaches:

| | first n-adic digit (`C(n,k) mod n`) | second n-adic digit (`C(n,k) mod n^2`) |
|---|---|---|
| **is it computable?** | yes, `poly(log n)` for any `n` (Theorem 5) | not known, and Proposition 10 says not without factoring |
| **does it leak factors?** | only via a `1/p` accident (Theorem 6) | at **constant** density (Theorem 9) |
| **why** | Theorem 2 has already quotiented the information away; what survives is aliased (T7) or symmetric (T8) | Lucas' digit condition separates the primes |

The first column is cheap and nearly empty. The second is rich and priced at
factoring. There is no third column in this framework — and that is a more
useful thing to know than another failed heuristic would have been.

The open problem this leaves behind is correspondingly precise: **find a
statistic of the first digit that is asymmetric in the prime factors.** Theorems
7 and 8 rule out the two obvious families (positional, symmetric). Nothing rules
out a third — but nothing tried here found one.
