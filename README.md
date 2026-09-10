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

Every row is machine-checked in [`aksfactor/theorems.py`](aksfactor/theorems.py)
and exercised by `run_tests.py` (97 tests, all passing).

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
  classgroup.py binary quadratic forms, Gauss composition, Schnorr-Lenstra
  fractal.py    the gasket: Lucas geometry, digit-DP box counts, row/column duals
  harvey.py     Harvey's N^(1/5) deterministic factoring: Lehman + Fermat + BSGS
  grouporder.py counts reachable group orders: ring unit orders vs #E(F_p)
  fast.py       O~(n^(1/4)) search: product tree, Newton division, remainder
                tree, multipoint evaluation, BGS factorial, threshold search
  cli.py        python -m aksfactor {factor,row,entry,verify,fold}
docs/           THEORY.md (proofs), FINDINGS.md (what it buys)
experiments/    sixteen reproducible scripts; results/ holds their generated reports
tests/          97 tests; run_tests.py needs no pytest
```

Regenerate every measurement (and the figure above) with `./run_experiments.sh`
(about 20 minutes; `exp09` dominates, because measuring the second-digit rate
needs Lucas at random base-`p` digits, which costs `O(p)` per call — the
barrier charges you even to observe it).

## License

MIT.
