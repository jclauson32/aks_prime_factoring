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
next, and why the obvious three are shut.
