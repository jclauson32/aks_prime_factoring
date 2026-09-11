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
| **R5** | class numbers `h(−kn)` are dense in `k` — the density requirement is satisfiable | implemented, measured |
| **T18** | the mod-`p` triangle is a Sierpiński gasket; exact `O(log N)` box count | proved, verified |
| **T19** | first row with a zero mod `n` = the **largest** prime factor (dual of T4) | proved, 127/127 |
| **T20** | `gcd(i! mod n, n) > 1` iff `i ≥ spf(n)` — monotone; derives Strassen's `Õ(n^(1/4))` | proved, 29,155 pairs |
| **R8** | Harvey's `N^(1/5)` algorithm (arXiv:2010.05450) implemented from the paper | correct, 281/281 |
| **R9** | anatomy of Harvey's Remark 3.4: the pair term binds by exactly 2×; the runs BSGS needs are absent | measured |
| **R10** | Lehman's interval covering is non-redundant — min subcover is `0.70 × pairs`, never `O(√r)` | measured |
| **R11** | Coppersmith implemented: **poly-time** factoring given `N^(1/4)` of `p`; guessing it is `Θ(N^(1/4))` by counting | implemented, proved |
| **R12** | poly-time factoring ⟺ poly-time approximation of `p` to `N^(1/4)`; guess+Coppersmith never beats Strassen | proved, measured |
| **R13** | the budget is `φ(N) mod ℓ` for primes `ℓ ≤ (1/4)ln N`; `N` reveals exactly the `p↔q` symmetry and nothing more | proved, measured |
| **R14** | generic ring programs achieve exactly `1−φ(N)/N`; every method's content is one *mechanism*; four are known (taxonomy corrected in round 15) | measured |
| **R15** | the central column computes Legendre symbols: a Pascal-native *sign* method, Strassen's exponent | implemented, measured |
| **T21** | row `pq` contains row `p` at stride `q`; exact support count via Fine's theorem (corrects an earlier density figure) | proved, 58,300 + 2,264 pairs |
| **R16** | no cheap function of `N` beats the Bayes baseline for `(p+q) mod ℓ`; `τ(N) mod 691` would pin it exactly, but is as hard as factoring | measured |
| **P22** | Pascal rho: `x → C(x,2)` *is* Pollard's rho; the row palindrome doubles the fibre statistic for even `k`, √2 fewer steps, no cheaper | proved (k=2), measured |

Every row is machine-checked in [`aksfactor/theorems.py`](aksfactor/theorems.py)
and exercised by `run_tests.py` (152 tests, all passing).

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
| 431 | 64 | 0.2015 | 64 | 1.36 |
| 431 | 256 | 0.5450 | 256 | 0.92 |
| 2,239 | 256 | 0.1943 | 256 | 1.70 |
| 39,239 | 64 | 0.0025 | 64 | 1.53 |
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

## Round 11: there *is* a polynomial-time algorithm — and its price is exactly N^(1/4)

You asked for polynomial time. It exists, and it has been outside this project's
scope the whole way: **Coppersmith's method**. Given `N = pq` and `p` known to
within about `N^(1/4)`, it recovers `p` in time **polynomial in log N** — no
search, no smoothness, no group order. It builds a lattice whose short vectors
are polynomials sharing `f`'s small root *over the integers*, LLL-reduces, and
reads the root off.

Implemented from scratch in `aksfactor/lattice.py` (exact-rational LLL +
Howgrave-Graham). It works:

| bits of N | unknown low bits of p | factored | seconds |
|---|---|---|---|
| 26 | 4 | yes | 0.02 |
| 34 | 6 | yes | 0.05 |
| 38 | 6 | yes | 0.06 |

Two bugs worth recording: my first root-finder enumerated divisors of the
constant term — which in a Coppersmith lattice is `~N^m`, so it hung; and my
first LLL recomputed Gram–Schmidt from scratch each step. Both were fatal, both
fixed.

**How much information does it need?** The bound approaches `N^(1/4)` as the
lattice grows, exactly as theory predicts (computed live by the experiment, not
transcribed):

| m | lattice dim | max unknown bits recovered | (1/4)·log₂N | fraction of limit |
|---|---|---|---|---|
| 2 | 4 | 8 | 10 | 0.80 |
| 3 | 6 | 9 | 10 | 0.90 |
| 5 | 10 | 9 | 10 | **0.90** |

The limit is `N^(1/4)` — **a quarter of the bits of `N`, half the bits of `p`.**

### Why it can't be bootstrapped

Guess the approximation and run Coppersmith on each guess. The cost is fixed by
pure counting:

> covering a range of length `L` with windows of width `w` needs at least `L/w`
> of them.

`√N / N^(1/4) = N^(1/4)`, always. Guess-and-Coppersmith is **`Θ(N^(1/4))`** — the
exponent Strassen reached in 1977 — and cannot reach `N^(1/5)`, let alone
`N^(1/6)`. This is not a measured tendency like rounds 9 and 10; it is a counting
bound with no escape. And no structure on the guesses changes it, because it
depends on only two numbers: range length and window width. Lehman's fan
reorganises the guesses; it does not reduce them.

### What that clarifies

Factoring **is** polynomial time given a quarter of the bits of `N`. So the
entire difficulty is *the cost of those bits* — and every route in this project
pays at least `N^(1/4)`:

| route | cost |
|---|---|
| trial division / Pascal row scan | `Θ(p)` |
| Strassen / the gasket threshold | `Õ(N^(1/4))` |
| Harvey's BSGS sweep | `Õ(N^(1/5))` |
| guess + Coppersmith | `Θ(N^(1/4))`, by counting |
| **Coppersmith alone** | **`poly(log N)`** — given `N^(1/4)` of it |

Two ways out, both famous open problems: **improve Coppersmith's exponent
`β²/d`** (widening the window past `N^(1/4)` would immediately beat everything
above — a known barrier in lattice cryptanalysis), or **find a source of
high-order bits of `p` cheaper than `N^(1/4)`**, which is precisely what eleven
rounds here failed to do.

## Round 12: the search reduces to an *equivalent* problem, not an easier one

Round 11 left a target: get `N^(1/4)` bits of `p` cheaply. The hope was that a
quarter of the bits might cost less than all of them. It doesn't — and now I can
say that precisely.

### The window is `N^(β²)`, and `β ≤ 1/2` always

Coppersmith recovers roots below `N^(β²/d)`, with `d = 1` for factoring. Measured:

| bits of N | β | bits recovered | β²·log₂N | ratio |
|---|---|---|---|---|
| 27 | 0.502 | 7 | 6.8 | 1.03 |
| 29 | 0.393 | 4 | 4.5 | 0.89 |
| 32 | 0.310 | 4 | 3.1 | 1.30 |
| 38 | 0.210 | 2 | 1.7 | 1.19 |

`p` is the *smaller* factor, so `β ≤ 1/2`: **the widest possible window is
`N^(1/4)`, at exactly the balanced semiprime.** No regime lets the lattice see
further.

### Guessing never beats Strassen

Covering `N^β` with windows of `N^(β²)` needs `N^(β−β²)` guesses; Strassen finds
`p` in `Õ(N^(β/2))`:

| β | window `β²` | guess+Coppersmith `β(1−β)` | Strassen `β/2` | winner |
|---|---|---|---|---|
| 0.50 | 0.2500 | 0.2500 | 0.2500 | **tie** |
| 0.40 | 0.1600 | 0.2400 | 0.2000 | Strassen |
| 0.25 | 0.0625 | 0.1875 | 0.1250 | Strassen |
| 0.10 | 0.0100 | 0.0900 | 0.0500 | Strassen |

`β(1−β) ≤ β/2` iff `β ≥ 1/2`, and `β ≤ 1/2` always. **Guess+Coppersmith is never
strictly better than Strassen**, and strictly worse for every unbalanced
semiprime. They meet at `N^(1/4)` exactly where `β(1−β)` peaks — the hardest
input for the lattice route is the hardest input for everything else.

### The equivalence

> **Theorem.** For balanced semiprimes `N = pq` with `p ≤ q ≤ 2p`, these are
> equivalent:
>
> **(a)** a deterministic `poly(log N)` algorithm factoring `N`;
> **(b)** a deterministic `poly(log N)` algorithm outputting `p̃` with
> `|p − p̃| ≤ N^(1/4)`.
>
> *Proof.* (a)⇒(b): factor and output `p`. (b)⇒(a): feed `p̃` to Coppersmith,
> whose window at `β = 1/2` is exactly `N^(1/4)`. ∎

Twelve rounds hunted for *partial* information on the theory that a quarter of
the bits might be cheaper than all of them. **They are the same problem.** Any
algorithm producing a quarter of the bits in polynomial time *is* a
polynomial-time factoring algorithm, and every barrier measured here applies to
it unchanged.

### What this settles, and what it doesn't

It settles the **strategy**: cheap partial information is not a way around the
barrier, it's the barrier restated. Two doors remain, both long known and both
shut:

- **widen the window** — improve Coppersmith's `β²/d`. The bound is tight for
  this lattice family, so it needs a genuinely different construction.
- **shrink the search** — find structure in the *location* of `p` that no
  algebraic, geometric, group-theoretic or fractal property examined across
  twelve rounds provided.

It does **not** settle the question. Factoring isn't known to be hard — it's in
NP ∩ co-NP, so it's very unlikely to be NP-complete, and nothing proved rules out
a polynomial-time algorithm. What this repository can offer is a map: eleven
distinct routes, each closed by a measurement rather than an intuition, and a
precise statement of what a twelfth must do that none of them did.

## Round 13: taking the Coppersmith door

The remaining lever was Coppersmith's `β²/d` bound. I didn't widen it — that
bound is tight for this lattice family. I attacked the *shape* of what fills the
window instead, and it compressed the problem further than anything else here.

### The budget can be spent as a congruence

Coppersmith is usually stated with the high bits of `p`. It works equally well
with `p mod M`: writing `p = r + Mx`, the unknown is bounded by `√N/M`, so it
succeeds once `M ≥ N^(1/4)`. (`r + Mx` isn't monic, so multiply by `M⁻¹ mod N` —
doesn't move the roots mod `p`.) Implemented as `factor_with_congruence`, verified.

**Which matters, because `M` can be built from small primes.** `p mod M` follows
by CRT from `p mod ℓ` for primes up to `~(1/4)ln N`.

### What N actually tells you about `p mod ℓ`

Over `F_ℓ`, `p` and `q` are the roots of `z² − sz + N` with `s = (p+q) mod ℓ`.
`N mod ℓ` is free; **`s` is the entire unknown** — and since `p+q = N+1−φ(N)`,
knowing `s` is knowing `φ(N) mod ℓ`.

| ℓ | mean candidates for `p mod ℓ` | (ℓ−1)/2 | ratio |
|---|---|---|---|
| 13 | 6.52 | 6 | 1.085 |
| 23 | 11.61 | 11 | 1.045 |
| 37 | 18.49 | 18 | 1.027 |

The map `a ↦ a + N/a` is two-to-one, since `a` and `N/a` collide. **N gives away
exactly the `p ↔ q` symmetry — a factor of two — and nothing else.**

Is any of it free? Scanning every `ℓ ≤ 13`: exactly **one** case. `N ≡ 2 (mod 3)`
forces `p+q ≡ 0 (mod 3)`, since `{1,2}` is the only unordered pair with product 2.
Nothing else, anywhere.

| bits of N | primes needed | largest ℓ | search `∏(ℓ−1)/2` | N^(1/4) | symmetry saving |
|---|---|---|---|---|---|
| 512 | 27 | 103 | 2^105.1 | 2^128 | 2^23 |
| 2048 | 76 | 383 | 2^439.2 | 2^512 | **2^73** |

`∏(ℓ−1)/2 = N^(1/4)/2^π(y)`, and `π(y) = O(log N/log log N)`, so the saving is
`N^o(1)` — enormous in absolute terms (2^73!), invisible in the exponent.

### Where the door leads

> **Factoring a balanced semiprime in polynomial time is equivalent to computing
> `φ(N) mod ℓ` for every prime `ℓ` up to about `(1/4)·ln N`.**

Each is a *single element of a field with `O(log N)` elements*. `N mod ℓ` is
handed to you free; `φ(N) mod ℓ` is the whole problem. For a **2048-bit modulus
the entire difficulty is 76 residues, modulo primes no larger than 383.**

That's the smallest this project has made factoring. It is also — the third time
now — an equivalence rather than a reduction: enough `φ(N) mod ℓ` gives `φ(N)`,
which factors `N` outright.

The pattern is consistent enough to name:

> **Every reformulation of factoring this repository found *computable* turned
> out to be symmetric in `p` and `q`; every *asymmetric* one turned out to be
> equivalent to factoring.**

The `p ↔ q` symmetry that Theorem 8 first found in the AKS ring is the same
obstruction that survives here — in a setting with no binomial coefficients
anywhere in it.

## Round 14: generic algebra is worth *nothing*; one alignment is worth everything

You're right that shaving exponents solves nothing. So this round stops proposing
algorithms and measures the **space** of them.

Factoring is the manufacture of a **zero divisor**. Under CRT `Z/N = F_p × F_q`,
and a nontrivial `gcd(x,N)` is an `x` vanishing in exactly one coordinate. Ring
operations act *diagonally* — the same `+`, `−`, `×` on both coordinates — so no
arrangement of them can tell `p` from `q`.

### The generic baseline

Random straight-line programs over `Z/N`: random base, arbitrary `+ − ×`. Since
`×` can square a register, depth `d` reaches exponents up to `2^d` — **this model
contains Pollard's `p−1`.**

| p | q | depth | measured split rate | `1 − φ(N)/N` | ratio |
|---|---|---|---|---|---|
| 179 | 419 | 6 | 0.005617 | 0.007960 | 0.706 |
| 179 | 419 | 24 | 0.007792 | 0.007960 | **0.979** |
| 733 | 761 | 24 | 0.002717 | 0.002677 | **1.015** |
| 5,849 | 7,927 | 24 | 0.000283 | 0.000297 | **0.954** |

It converges to `1/p + 1/q` — **exactly the chance a uniformly random element of
`Z/N` happens to be a zero divisor.** Depth buys astronomical exponents and buys
nothing else. Generic algebra is worth precisely as much as reaching into a hat.

### What one alignment is worth

| lcm bound B | Pollard `p−1` rate | generic baseline | gain |
|---|---|---|---|
| 50 | 0.0787 | 0.000427 | 184× |
| 200 | 0.9947 | 0.000427 | **2,332×** |

Same ring, same operations, same depth budget. The **only** difference: the
exponent is chosen divisible by everything small — aligned with the *order* of
the group mod `p`. That one choice is the entire method.

### The taxonomy — *corrected in round 15*

> **Correction.** This section originally said "only two alignments are known,"
> order and size. That was wrong: **Pollard rho** (collisions) and the
> **quadratic/number field sieves** (square-root signs) are neither, and NFS's
> `L[1/3]` beats both caps the original table listed.

| mechanism | how the zero divisor arises | methods | best known |
|---|---|---|---|
| **order** | exponent divisible by `\|G\|` for a group attached to `p` | Pollard, Williams, ECM, class groups | `L[1/2]` — smooth-number density |
| **size** | a window that contains `p` | trial division, Fermat, Lehman, Strassen, Coppersmith, Harvey | `N^(1/5)` deterministic |
| **collision** | two iterates agreeing mod `p` only | Pollard rho | `N^(1/4)` — birthday bound |
| **sign** | a square root with mixed CRT signs | Dixon, QS, NFS; round 15's central column | `L[1/3]` (NFS) |

Every round up to 14 was *order* or *size*. None was collision or sign — part of
why none approached `L[1/3]`.

### So what a polynomial-time algorithm has to do

In every mechanism the cost sits in the same place: **finding** the structure
(smooth order, window, collision, mixed-sign root), never using it — once found,
the factor is one gcd away. A polynomial-time algorithm needs one of those
finding steps to be polynomial, or a mechanism not on the list. The barriers in
this repository are each statements about a finding step:

- symmetric constructions (T8, `N mod ℓ`) find nothing — they see `N`, not `p`;
- aliased constructions (T7) can't locate a window below the scale `p`;
- rigid orders (P14) give one order per `p`, no redraw;
- dense families (round 5) redraw freely and hit smoothness density.

Not a proof that no polynomial mechanism exists — a measurement that fourteen
rounds found none, and a test for the next attempt: **which mechanism does it
use, and why would its finding step be cheaper than the best known one?**

## Round 15 onward: the idea ledger

From round 15 the attempts are logged in [`docs/LAB.md`](docs/LAB.md) — one row
per idea, with its mechanism, its verdict, and the experiment behind it — so this
README stops growing a section per round.

Round 15 also **corrects round 14**: its taxonomy claimed only two mechanisms
(order, size) were known. Pollard rho (collision) and the quadratic/number field
sieves (square-root sign) are two more, and NFS's `L[1/3]` beats every cap round
14 listed. The round-14 section above now carries the corrected table.

The headline of round 15 itself: the **central column** of Pascal's triangle,
`Σ C(2k,k) xᵏ = (1−4x)^(−1/2)`, truncated anywhere in `[(p−1)/2, p−1]`, equals the
Legendre symbol `((1−4a)/p)` mod `p` (3,168/3,168 verified). That gives a
Pascal-native factoring method on the *sign* mechanism, which detects `p` at
`T ≈ p/2` — **below** `p`, where Strassen's factorial test is blind — at
Strassen's exponent and about 15× his constant.

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
  theorems.py   machine-checkable form of T1-T5, T7, T9, T15, T16, T20, T21 and the conjecture
  factor.py     pascal_spf / pascal_split / factor, with certificates
  ring.py       (x+a)^n mod (n, x^r - 1), the fold identity, the fold attack,
                and the norm (symmetry obstruction)
  cyclo.py      norm-one cyclotomic factoring: Lucas sequences, Williams p+1,
                general degree-d via resultants
  qpascal.py    q-deformed Pascal row: q-Kummer criterion, tunable period
  classgroup.py binary quadratic forms, Gauss composition, Schnorr-Lenstra
  fractal.py    the gasket: Lucas geometry, digit-DP box counts, row/column duals
  harvey.py     Harvey's N^(1/5) deterministic factoring: Lehman + Fermat + BSGS
  lattice.py    exact LLL and Coppersmith: polynomial-time factoring given a hint
  generic.py    the generic-ring baseline: what algebra buys without an alignment
  leakage.py    round 16: feature scans, Ramanujan tau mod 691, four-square counts
  collision.py  round 17: Pascal rho, x -> C(x, k) mod N, and its fibre statistic
  central.py    the central column: Legendre symbols from Pascal's triangle
  grouporder.py counts reachable group orders: ring unit orders vs #E(F_p)
  fast.py       O~(n^(1/4)) search: product tree, Newton division, remainder
                tree, multipoint evaluation, BGS factorial, threshold search
  cli.py        python -m aksfactor {factor,row,entry,verify,fold}
docs/           THEORY.md (proofs), FINDINGS.md (what it buys), LAB.md (idea ledger)
experiments/    twenty-five reproducible scripts; results/ holds their generated reports
tests/          152 tests; run_tests.py needs no pytest
```

Regenerate every measurement (and the figure above) with `./run_experiments.sh`
(about 20 minutes; `exp09` dominates, because measuring the second-digit rate
needs Lucas at random base-`p` digits, which costs `O(p)` per call — the
barrier charges you even to observe it).

## License

MIT.
