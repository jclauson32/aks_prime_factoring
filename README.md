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

Every row is machine-checked in [`aksfactor/theorems.py`](aksfactor/theorems.py)
and exercised by `run_tests.py` (30 tests, all passing).

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
  theorems.py   machine-checkable form of T1-T5 and the original conjecture
  factor.py     pascal_spf / pascal_split / factor, with certificates
  ring.py       (x+a)^n mod (n, x^r - 1), the fold identity and the fold attack
  cli.py        python -m aksfactor {factor,row,entry,verify,fold}
docs/           THEORY.md (proofs), FINDINGS.md (what it buys)
experiments/    six reproducible scripts; results/ holds their generated reports
tests/          30 tests; run_tests.py needs no pytest
```

Regenerate every measurement (and the figure above) with `./run_experiments.sh`.

## License

MIT.
