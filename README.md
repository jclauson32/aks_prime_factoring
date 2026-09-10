# aks_prime_factoring

**Reduce row `n` of Pascal's triangle mod `n`. The first non-zero entry sits at
position `spf(n)` and its value is `n / spf(n)`. Position times value is `n`.**

That single sentence is a theorem, it is proved in [`docs/THEORY.md`](docs/THEORY.md),
and it turns the AKS primality identity into a factor extractor. It is also,
provably, no faster than trial division — and this repository measures exactly
why, which is the more useful half of the result.

## The picture

```
  n │ interior of row n mod n   (· = 0,  █ = non-zero,  ▲ = first non-zero)
────┼────────────────────────────────────────────────────────
  4 │ ·▲·                                          spf=2  n/spf=2
  5 │ ····                                         prime
  6 │ ·▲██·                                        spf=2  n/spf=3
  7 │ ······                                       prime
  8 │ ·▲·█·█·                                      spf=2  n/spf=4
  9 │ ··▲··█··                                     spf=3  n/spf=3
 10 │ ·▲··█··█·                                    spf=2  n/spf=5
 11 │ ··········                                   prime
 12 │ ·▲██···███·                                  spf=2  n/spf=6
 13 │ ············                                 prime
 14 │ ·▲·█·███·█·█·                                spf=2  n/spf=7
 15 │ ··▲·██··██·█··                               spf=3  n/spf=5
 16 │ ·▲·█·█·█·█·█·█·                              spf=2  n/spf=8
 17 │ ················                             prime
 18 │ ·▲█··█··█··█··██·                            spf=2  n/spf=9
 19 │ ··················                           prime
 20 │ ·▲·██··█·█·█··██·█·                          spf=2  n/spf=10
 21 │ ··▲···█·█··█·█···█··                         spf=3  n/spf=7
 22 │ ·▲·█·█····█····█·█·█·                        spf=2  n/spf=11
 23 │ ······················                       prime
 24 │ ·▲██·█·██··█··██·█·███·                      spf=2  n/spf=12
 25 │ ····▲····█····█····█····                     spf=5  n/spf=5
 26 │ ·▲·····█·█··█··█·█·····█·                    spf=2  n/spf=13
```

Prime rows are empty. Composite rows fire for the first time at their smallest
prime factor, every time.

## Quick start

No dependencies — Python 3.10+, standard library only.

```bash
git clone https://github.com/jclauson32/aks_prime_factoring
cd aks_prime_factoring

python3 -m aksfactor factor 1234567
python3 -m aksfactor row 35 --kmax 12
python3 -m aksfactor verify --upto 400
python3 run_tests.py
```

```
$ python3 -m aksfactor factor 1234567
1234567 = 127 * 9721
  certificate: C(1234567, 127) mod 1234567 = 9721 = 1234567/127
  first non-zero residue of row 1234567 is at k = 127
  smallest prime factor 127, cofactor 9721
```

```python
from aksfactor import row_entry, factor, certificate

factor(1234567)                  # {127: 1, 9721: 1}
row_entry(35, 5)                 # 7   -- position 5, value 7, and 5 * 7 = 35
certificate(1234567, 127)        # the C(n,p) mod n = n/p witness

# a 67-digit modulus, one residue, no row materialised:
n = 1000003 * (2**200 - 75)
row_entry(n, 1000003)            # returns the 61-digit cofactor exactly, in microseconds
```

## The claim that started this

> Take the rows of Pascal's triangle and mod `n` them; the remainder can very
> often — if not always — be broken into `x * y` where `x` is a prime factor.

**The split always exists.** For every interior `k`,
`C(n,k) mod n = x * y` with `x = n/gcd(n,k)`, and `x` is always a proper
divisor of `n`. Checked at 204,464 interior positions across every row
`n < 700`: zero counterexamples.

**`x` is a divisor; it is prime 36% of the time in general** — but at the
*first* non-zero position, `x = n/spf(n)`, which is prime exactly when `n` is a
semiprime. So for `n = p*q` the conjecture is exactly right: the first non-zero
remainder **is** the prime factor `q`, with `y = 1`.

## Results

| | statement | status |
|---|---|---|
| **T1** | `n/gcd(n,k)` divides `C(n,k)`; every residue is `x * y` with `x \| n` | proved, exhaustively verified |
| **T2** | `gcd(k,n) = 1` ⟹ `C(n,k) ≡ 0 (mod n)` | proved |
| **T3** | `C(n,p) mod n = n/p` **exactly**, for every prime `p \| n` | proved (Lucas) |
| **T4** | first non-zero interior entry is at `spf(n)`, with value `n/spf(n)` | proved |
| **T5** | closed form for any entry; gives `O(1)` random access for huge `n` | proved, implemented |
| **T6** | folded coefficients mod `p` are binomial progression sums of the cofactor | proved, verified |
| **T7** | with `gcd(r,p)=1` the row's support is *equidistributed* mod `r` — an aliasing barrier | proved, verified |
| **T8** | `Norm((x+a)^n) = (a^r+1)^n mod n` — symmetric aggregates carry zero information | proved, verified |
| **T9** | the *second* n-adic digit of the row splits `n` at **constant** rate, not `1/p` | proved, verified |

Every row is machine-checked in [`aksfactor/theorems.py`](aksfactor/theorems.py)
and exercised by `run_tests.py` (41 tests, all passing).

## The honest verdict

The method is exact, deterministic, certificate-producing — and
`Theta(spf(n))`, i.e. trial-division complexity.

**Scanning cannot beat it,** because T2 proves that "scan until a residue
fires" *is* "scan until a gcd fires". The theorem that makes the method pretty
is the theorem that caps its speed.

**Folding cannot beat it either.** Computing `(x+a)^n mod (n, x^r - 1)` reads
the entire row in `O(log n)` ring multiplications, and it genuinely leaks
factors — on the first pass it cracked 24 of 25 random semiprimes with `r < 40`.
That looked like a break. It was an artifact of small `p`. T6 explains it:
each fold coefficient leaks a factor with probability `~1/p`, so `Theta(p)`
coefficients are needed. Measured over three orders of magnitude, the ratio
*coefficients-inspected / p* has median **0.87** and no trend
([exp04](experiments/results/exp04_fold_scaling.md)), and the per-coefficient
hit rate tracks `1/p` across twelve AKS bases
([exp05](experiments/results/exp05_multibase.md)).

**The unified reason:** modulo a prime factor `p`, the AKS object *is*
`(x^(p^v) + a)^(n/p^v)`. Its first non-trivial structure lives at degree
`p^v`. Truncate below that and you see zeros; fold below that and it smears
across all `r` buckets. Either way you pay `Theta(p)`. Detecting that the AKS
identity fails is cheap; localising *where* it fails is not — and that gap is
precisely the gap between primality testing and factoring.

[`docs/FINDINGS.md`](docs/FINDINGS.md) lists the three doors that are shut and
what a fourth would have to look like.

## Pushing on it: round 2

Four more routes to polynomial time, tested and killed — plus one that moved the
exponent. Full scorecard in [exp07](experiments/results/exp07_dead_ends.md).

| theory | why it died |
|---|---|
| **Symmetric aggregates** (norms, resultants) — buy every lottery ticket at once | `Norm((x+a)^n) = (a^r+1)^n mod n`, a function of `n` alone. Same value mod *every* prime factor, so the gcd is always `n`. Symmetric aggregation destroys the asymmetry a factor consists of. |
| **Twisted moduli `x^r - c`** — bias a coefficient toward vanishing | Hit rate still tracks `1/p`. Different tickets, same price. |
| **AKS-ring Pollard `p-1`** — many `k`, many chances | `p-1` divides `p^k-1`, so `p^k-1` is smooth only if `p-1` already was. *Strictly dominated* by the method it generalises. |
| **Position from the fold** — leak `p mod r`, CRT over several `r`, done in polylog | **Theorem 7:** if `gcd(r,p)=1` then `i ↦ i·p^v mod r` is a bijection, so the classes are equally occupied *by theorem*. Positional info needs `r ≥ p`. Nyquist aliasing. |

### What did work: `Õ(n^(1/4))`

Theorem 2 makes the search "first `k` with `gcd(k,n) > 1`". One position per gcd
is trial division — but a *block* of `c` positions costs one gcd, since
`gcd(f(i·c), n) > 1` for `f(X) = (X+1)···(X+c)` exactly when the block
`(i·c, i·c+c]` contains a position sharing a factor with `n`. Evaluating one
degree-`c` polynomial at `c` points is `Õ(c)` ring operations via a product tree
and a remainder tree, so `c²` positions cost `Õ(c)`. Take `c = n^(1/4)`.

Implemented in [`aksfactor/fast.py`](aksfactor/fast.py), verified against trial
division exhaustively, and it beats a tuned 2-3-5 wheel from `p ≈ 7×10⁸` onward
([exp08](experiments/results/exp08_quartic.md)).

Two caveats, stated plainly: the wall-clock exponent is `~0.85` rather than
`0.5` because CPython uses Karatsuba rather than FFT for big integers — the
`Õ(n^(1/4))` claim is about *ring operations*, which are counted directly — and
`Õ(n^(1/4))` is Strassen's classical bound, not the deterministic record, which
is `Õ(n^(1/5))` (Hittmeir; Harvey).

### The one that changed the picture: `C(n,k) mod n²`

Every route above works at the **first** n-adic digit, where Theorem 2 has
already quotiented the information away. So go one digit deeper: write
`C(n,k) = n·m`. Then `m = C(n-1,k-1)/k`, and by Lucas, `m ≡ 0 (mod p)` exactly
when some base-`p` digit of `k-1` exceeds that of `n-1`. For random `k` that is a
**constant-probability** event — `n` has only ~3 digits base `p`.

| `p` | split rate of `gcd(m, n)` | `1/p` | ratio |
|---|---|---|---|
| 37 | 0.815 | 2.7e-2 | 30× |
| 1,571 | 0.895 | 6.4e-4 | 1,406× |
| 44,501 | 0.854 | 2.3e-5 | 37,982× |
| 810,853 | 0.973 | 1.2e-6 | 789,230× |

Flat in `p` — it depends on the base-`p` digit profile of `n-1`, not on the
size of `p`. Which makes it a **hardness result**, not an algorithm:

> **Proposition 10.** If `C(n,k) mod n²` can be computed in `poly(log n)` for
> random `k`, then `n = pq` factors in `O(1)` expected tries. Binomial
> coefficients mod `n²` are factoring-hard.

That is where the barrier becomes sharp. Theorem 5 gives real random access to
`C(n,k) mod n` for any size `n`, because the closed form only needs
`C(n-1,k-1)` modulo the **small** number `k`. Going one digit further would
break factoring. So the factors are *not hidden* in the Pascal row — they sit in
the second digit at constant density. The whole obstruction is that reading that
digit is the problem itself.

### The mechanism every dead end shares

Modulo a prime factor `p`, the AKS object **is** `(x^(p^v) + a)^(n/p^v)` — a
series whose only distinguishing feature is a **period** of `p^v`. Everything
cheap you can compute from it is either *symmetric* in the prime factors (same
value mod every `p`, gcd is `n`) or *aliased* below the period (spread evenly
over all `r < p` buckets).

So at the first digit, extracting `p` is **period finding**, which classically
costs `Ω(p)` samples. It is exactly the problem Shor's algorithm solves in
polylog time quantumly, by transforming at size `n` instead of size `r ≪ p`.
What the framework lacks is not a cleverer `r`, base `a`, or aggregation rule —
it is resolution.

## Benchmarks

Balanced semiprimes, finding `spf(n)`:

| `p = spf(n)` | naive row scan | residue scan | certified | trial division | Pollard rho |
|---|---|---|---|---|---|
| 631 | 0.009 s | 0.00008 s | 0.00003 s | 0.00001 s | 0.00004 s |
| 3,593 | 0.169 s | 0.00043 s | 0.00008 s | 0.00006 s | 0.00004 s |
| 65,119 | (too slow) | 0.00855 s | 0.00110 s | 0.00109 s | 0.00020 s |
| 360,817 | (too slow) | 0.04970 s | 0.00612 s | 0.00664 s | 0.00016 s |

Representative single-run timings (CPython 3.12, no native extensions); regenerate
with `./run_experiments.sh`. The *naive row scan* refuses to use Theorem 2 and
reads real residues; the *residue scan* calls `row_entry`, which applies Theorem 2
and so degenerates into a gcd loop -- that collapse is the theorem doing the work,
and it is why columns 2-4 converge.

Pollard rho is `O(n^(1/4))` and pulls away as `p` grows. That is the ceiling.

## Layout

```
aksfactor/
  arith.py      primality, sieving, Pollard rho, reference factoring
  binom.py      C(N,M) mod m for huge N,M and small m -- Lucas, Granville, and an
                independent falling-factorial reference; all three cross-checked
  pascal.py     row_entry (random access), row_series (Kronecker-packed AKS
                truncation), row_exact; the x*y decomposition
  theorems.py   machine-checkable form of T1-T5, T7, T9 and the conjecture
  factor.py     pascal_spf / pascal_split / factor, with certificates
  ring.py       (x+a)^n mod (n, x^r - 1), the fold identity, the fold attack,
                and the norm (symmetry obstruction)
  fast.py       O~(n^(1/4)) search: product tree, Newton division, remainder
                tree, fast multipoint evaluation
  cli.py        python -m aksfactor {factor,row,entry,verify,fold}
docs/           THEORY.md (proofs), FINDINGS.md (what it buys)
experiments/    nine reproducible scripts; results/ holds their generated reports
tests/          41 tests; run_tests.py needs no pytest
```

Regenerate every measurement (and the figure above) with `./run_experiments.sh`.

## License

MIT.
