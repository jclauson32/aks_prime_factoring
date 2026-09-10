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
| **T13** | the norm-one subgroup has order `(p^d−1)/(p−1)`, which `p−1` does not divide | proved; refutes an earlier claim |
| **P14** | ring group orders are capped at `partitions(d)` per `p`; elliptic orders grow with `√p` | proved, measured |
| **T15** | below `spf(n)`, `C(N,k) mod n` depends only on `N mod n` — *no* row leaks there | proved, verified |
| **T16** | q-analogue of Kummer: the q-row's period is `ord_p(q)`, which **varies** with `q` | proved, verified |

Every row is machine-checked in [`aksfactor/theorems.py`](aksfactor/theorems.py)
and exercised by `run_tests.py` (64 tests, all passing).

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

*(One conclusion in this section was later found to be too strong; the retraction
is in [round 3](#round-3-the-correction) and inline in the experiment that made
it.)*

Four more routes to polynomial time, tested and killed — plus one that moved the
exponent. Full scorecard in [exp07](experiments/results/exp07_dead_ends.md).

| theory | why it died |
|---|---|
| **Symmetric aggregates** (norms, resultants) — buy every lottery ticket at once | `Norm((x+a)^n) = (a^r+1)^n mod n`, a function of `n` alone. Same value mod *every* prime factor, so the gcd is always `n`. Symmetric aggregation destroys the asymmetry a factor consists of. |
| **Twisted moduli `x^r - c`** — bias a coefficient toward vanishing | Hit rate still tracks `1/p`. Different tickets, same price. |
| **AKS-ring Pollard `p-1`** — many `k`, many chances | `p-1` divides `p^k-1`, so for the element `x+a` it is dominated by the method it generalises. **Partly retracted in round 3** — the argument covers only the full unit group. |
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

## Round 3: the correction

Round 2 ended by naming the open problem — *find a statistic of the first digit
that is asymmetric in the prime factors*. Round 3 found one, by noticing that
round 2's own domination argument had a gap.

| group | order | divisible by `p-1`? |
|---|---|---|
| full unit group of `F_{p^d}` | `p^d − 1` | yes — hence dominated |
| **norm-one subgroup** `ker(N)` | `(p^d − 1)/(p − 1)` | **no** |

The domination argument is about the *full* unit group, and `x + a` lives there.
It says nothing about subgroups. And you can land in the norm-one subgroup
**without knowing `p`** — pick a monic polynomial whose roots multiply to `1`.
For `d = 2` that is `x² − ax + 1`, tracked by the Lucas sequence `V_k(a,1)`:
Williams' `p+1` method, recovered rather than imported.

Head to head at an identical budget, against a cofactor rough on both sides:

| p | largest prime factor of `p−1` | of `p+1` | Pollard `p−1` | norm-one |
|---|---|---|---|---|
| 4,512,218,267 | 3,716,819 | 47 | fail | **found** |
| 7,205,934,499 | 171,569,869 | 41 | fail | **found** |
| 2,850,909,643 | 41,669 | 53 | fail | **found** |
| 34,886,131,169 | 88,339 | 53 | fail | **found** |

**6 of 6** in the full run ([exp10](experiments/results/exp10_norm_one.md)). The
retraction is written into `exp07_dead_ends.py` inline, not quietly dropped.

And the statistic is *still a binomial sum* —
`V_k(a,1) = Σ_j (−1)^j (k/(k−j)) C(k−j, j) a^(k−2j)`, verified on thousands of
cases. What changes is not leaving binomials behind; it is constraining the roots
to multiply to `1`, which pins the order to `p+1` instead of `p−1`.

### Why it is still not polynomial

| method | group | order | varies at fixed `p`? |
|---|---|---|---|
| Pollard `p−1` | `F_p*` | `p − 1` | no |
| norm-one, degree `d` | `ker N` in `F_{p^d}*` | `Φ_d(p)`-ish | no |
| ECM | `E(F_p)` | `p + 1 − t` | **yes**, one per curve |

Everything reachable here has a **rigid** group order: fix `p` and the number
whose smoothness decides success is fixed too. One ticket per `p`; if `Φ_d(p)` is
rough for every small `d`, you are stuck. ECM's advantage is not a better group
but a *family* of them — a fresh order per curve at the same `p`.

### Counted, not asserted

Distinct group orders reachable at one fixed `p`, 400 random samples each
(found / theoretical cap):

| p | ring `d=2` | ring `d=3` | ring `d=4` | distinct `#E(F_p)` | Hasse width `4*sqrt(p)` |
|---|---|---|---|---|---|
| 211 | 2 / 2 | 3 / 3 | 5 / 5 | 55 | 56 |
| 503 | 2 / 2 | 3 / 3 | 5 / 5 | 84 | 88 |
| 1,009 | 2 / 2 | 3 / 3 | 5 / 5 | 109 | 124 |
| 2,003 | 2 / 2 | 3 / 3 | 5 / 5 | 140 | 176 |

The ring columns **saturate at `partitions(d)` and stop** — 2, 3, 5 — regardless
of `p`. That is structural: `F_p[x]/f` for squarefree `f` is a product of finite
fields, so its unit group order is `∏ (p^{d_i} − 1)`, a function of the degree
*partition* alone. Same pattern, same order — nothing left to vary. The elliptic
column grows with `√p`.

So the open problem, as a requirement: **a construction over `Z/n` whose group
order at fixed `p` varies with a parameter we control.** Quotients of `(Z/n)[x]`
cannot supply one — their orders are pinned by the degree partition. Escaping
needs a different algebraic group, which is what an elliptic curve is and what no
rearrangement of Pascal's triangle will produce.

## Round 4: deforming the period

Round 3 asked for a parameter that *varies* at fixed `p`. Round 4 found one
inside Pascal's triangle — then found out why that was the wrong thing to ask
for.

Replace `C(n,k)` with the Gaussian binomial `[n,k]_q`. Divisibility by `p` is
then governed by `d = ord_p(q)`, which moves as `q` moves (**T16**, a q-analogue
of Kummer verified exhaustively). It fires at `k = (n mod d) + 1` — which can sit
**below `spf(n)`, inside the region Theorem 2 seals off**:

| n | spf(n) | classical row below spf | q-row first hit | base |
|---|---|---|---|---|
| 143 = 11·13 | 11 | 0 everywhere | k = 4 | 2 |
| 221 = 13·17 | 13 | 0 everywhere | k = 2 | 4 |
| 899 = 29·31 | 29 | 0 everywhere | k = 4 | 2 |

A real deformation of the barrier. My first guess about why it still fails —
"`ord_p(q)` divides `p−1`, so it's just Pollard" — was **wrong**: with `p−1`
smooth a random base still has order near `p−1`, so the raw q-deformation is
strictly *weaker* than Pollard. What actually closes it is a dichotomy. Tune the
base to `Q = q^M` with `M` a smooth ladder to bound `B`; then `ord_p(Q)` is the
`B`-rough part of the order, so:

| ladder bound B | tuned period = 1 | 1 < period ≤ B | period > B |
|---|---|---|---|
| 50 | 0 | **0** | 60 |
| 200 | 2 | **0** | 58 |

The middle is **empty**. Period 1 means `Q ≡ 1 (mod p)` — `gcd(Q−1,n)` has
already split `n`, i.e. Pollard verbatim. Otherwise clause 1 is a `1/B` lottery.
Work stays `Θ(p)`.

Also closed, in one line: **Theorem 15** — for `k < spf(n)`, `C(N,k) mod n`
depends only on `N mod n`, so *every* row is identical below `spf(n)`. Row `2n`,
row `n+1`, row `n²` — all carry the same nothing. "Try another row" is dead.

### The lesson

Round 3 asked for a parameter that varies. Round 4 supplied one and learned:

> **It is not variation that matters, but the density of the set varied over.**

Divisors of `p−1` are sparse and structured, and elements realising the small
ones are vanishingly rare — hence the dichotomy. An elliptic curve's order ranges
over an *interval*: dense, so every fresh curve is a genuinely fresh number.

**Open, for round 5:** a construction over `Z/n` whose governing quantity at
fixed `p` ranges over a *dense* set — an interval, not a divisor lattice — while
staying computable without knowing `p`.

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
  theorems.py   machine-checkable form of T1-T5, T7, T9, T15, T16 and the conjecture
  factor.py     pascal_spf / pascal_split / factor, with certificates
  ring.py       (x+a)^n mod (n, x^r - 1), the fold identity, the fold attack,
                and the norm (symmetry obstruction)
  cyclo.py      norm-one cyclotomic factoring: Lucas sequences, Williams p+1,
                general degree-d via resultants
  qpascal.py    q-deformed Pascal row: q-Kummer criterion, tunable period
  grouporder.py counts reachable group orders: ring unit orders vs #E(F_p)
  fast.py       O~(n^(1/4)) search: product tree, Newton division, remainder
                tree, fast multipoint evaluation
  cli.py        python -m aksfactor {factor,row,entry,verify,fold}
docs/           THEORY.md (proofs), FINDINGS.md (what it buys)
experiments/    twelve reproducible scripts; results/ holds their generated reports
tests/          64 tests; run_tests.py needs no pytest
```

Regenerate every measurement (and the figure above) with `./run_experiments.sh`
(about 20 minutes; `exp09` dominates, because measuring the second-digit rate
needs Lucas at random base-`p` digits, which costs `O(p)` per call — the
barrier charges you even to observe it).

## License

MIT.
