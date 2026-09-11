# The idea ledger

Every attempt at a faster factoring method made in this repository, with the
mechanism it uses and what killed it. The test applied to each entry is the one
from round 14 (as corrected in round 15):

> **Which mechanism does it use -- order, size, collision, sign, or something
> new -- and why would its finding step be cheaper than the best known one?**

| # | idea | mechanism | verdict | evidence |
|---|---|---|---|---|
| 1 | first non-zero entry of row `n` mod `n` | size | exact theorem (sits at `spf(n)`, value `n/spf(n)`); `Theta(p)` to find | `exp02` |
| 2 | AKS fold `(x+a)^n mod (n, x^r - 1)` | order | leaks at rate `1/p` per coefficient | `exp04`, `exp05` |
| 3 | norm of the fold | symmetric | identically `(a^r+1)^n` -- zero information | `exp07` |
| 4 | twisted moduli `x^r - c` | order | still `1/p` | `exp07` |
| 5 | reading position from the fold | size | aliased: support equidistributed mod `r` (T7) | `exp07` |
| 6 | Strassen batching of the row scan | size | `O~(n^(1/4))` -- works, and is the ceiling of this family | `exp08` |
| 7 | second n-adic digit `C(n,k) mod n^2` | sign-like | leaks at constant rate, but computing it is factoring-hard (P10) | `exp09` |
| 8 | norm-one subgroup, Williams `p+1` | order | real method; rigid order, one ticket per `p` | `exp10` |
| 9 | ring unit groups over varying `f` | order | capped at `partitions(d)` orders | `exp11` |
| 10 | q-deformed Pascal row | order | tunable period, but tuned period is `1` or `> B` (dichotomy) | `exp12` |
| 11 | other rows `N` | -- | identical below `spf(n)` (T15) | `exp12` |
| 12 | class groups (Schnorr-Lenstra) | order | dense orders, capped at `L[1/2]` | `exp13` |
| 13 | the Sierpinski gasket, row/column probes | size | isotropic `Theta(p)`; dual theorem T19 | `exp14` |
| 14 | factorial threshold + binary search | size | derives Strassen, loses to it by 20x | `exp15` |
| 15 | Harvey's `N^(1/5)` | size | implemented; both routes to `N^(1/6)` measurably closed | `exp16`-`exp18` |
| 16 | Coppersmith, any shape of hint | size | poly time given `N^(1/4)` of `p`; getting it is `Theta(N^(1/4))` | `exp19`-`exp21` |
| 17 | generic straight-line programs | none | exactly the zero-divisor density `1 - phi(N)/N` | `exp22` |
| 18 | central column `sum C(2k,k) x^k` | **sign** | Legendre symbol at `T ~ p/2`, below `p`; Strassen's exponent, ~15x the constant | `exp23` |
| 19 | `d`-th power characters | sign -> order | threshold `p/d`, but needs `d \| p-1` -- walks back into order | `exp23` |
| 20 | row `pq` at stride `q` (T21) | size | row `p` copied verbatim; locating the stride is locating `q` | `exp24` |
| 21 | Hensel-lifting `p` from the low bits | -- | all `2^(k-1)` branches survive; never prunes | `exp24` |
| 22 | cheap functions of `N` vs `(p+q) mod ell` | -- | max excess over Bayes baseline `+0.011` (noise); control 100% | `exp24` |
| 23 | Ramanujan `tau(N) mod 691` | order-free, divisor-bound | pins `(p+q) mod 691` exactly; as hard as factoring (Bach-Charles) | `exp24` |
| 24 | four-square count `r_4(N) = 8 sigma(N)` | divisor-bound | gives `p+q`; the count costs about `N` | `exp24` |
| 25 | Pascal rho `x -> C(x, k)` | **collision** | `k=2` is Pollard's rho (`c = -5/16`); even `k` collide `sqrt 2` sooner (P22), paid back in squarings | `exp25` |

## Round 15: the central column

The one mechanism no earlier round had used here was the **sign** of a square
root, and Pascal's triangle has it built in: the central column is the power
series of `(1-4x)^(-1/2)`. Truncated anywhere in `[(p-1)/2, p-1]`, it equals the
Legendre symbol `((1-4a)/p)` modulo `p` -- 3,168 of 3,168 triples verified -- so
`gcd(P_T(a) -+ 1, N)` splits `N`.

Two traps along the way, both caught by the harness rather than by luck:

- the first test "factored" 12 of 12 semiprimes, all at `T = 2` or `4`, where
  `P_2(a) - 1 = 2a(1+3a)` simply hit `p` -- trial division by a polynomial, with
  primes too small to rule out coincidence. With `p` in the thousands, splits
  are **100%** inside the window and **0%** below it;
- the splitter labeled every split "Legendre symbol", including one at
  `T = 512` for `p = 16249`, where the sum hit `-1 mod p` by chance. It now
  classifies after the fact, and the test requires the label to be honest.

The mechanism fires **below `p`**, where Strassen's factorial test is blind --
but it costs `O~(sqrt(p))` like Strassen, with a constant about 15 times worse
(2x2 polynomial matrices instead of scalars). Its `d`-th power generalisation
lowers the threshold to `(p-1)/d` (7 of 8 exactly; the eighth a below-window
coincidence), but only when `d | p-1`, which is the order mechanism's smoothness
condition wearing a different hat.

## Round 16: hunting for anything that sees `p` rather than `N`

Round 13 set the price of polynomial-time factoring at one residue
`(p+q) mod ell` per small prime. This round looked for that residue in five
places (`exp24`):

- **Inside the row.** Theorem 21: row `pq` contains row `p` verbatim at stride
  `q` (58,300 of 58,300 entries). The row *has* `p` in it -- at a stride that is
  `q` itself. Checking the theorem also exposed a wrong figure in the round-3
  write-up: the support density is not `~1/p + 1/q` in general. Fine's theorem
  gives the exact count, `(p-1) + prod(d_i+1) - 2` over the base-`p` digits of
  `q`, which runs from about `p + 2(q-p)` for close primes up to `p + q` when
  `q = -1 (mod p)`.
- **Low bits.** Every odd residue of `p` mod `2^k` has a partner; lifting never
  prunes.
- **Cheap functions of `N`.** Seven features, four moduli, 6,000 semiprimes. The
  positive control (close primes, `floor(sqrt N)`) scores 100%; on random
  semiprimes nothing exceeds the Bayes baseline `2/(ell-1)` by more than 0.011.
  The first scan used a too-weak baseline (uniform over the allowed set) and
  every feature "beat" it by the same margin -- a uniform excess across
  unrelated features is a baseline bug, and was.
- **Ramanujan's tau.** `tau(N) = (1+p^11)(1+q^11) (mod 691)`, and the eleventh
  power is a bijection on `F_691`, so `tau(N) mod 691` cuts `(p+q) mod 691`
  from 345.5 candidates to 1.00. It is the first quantity in the project worth
  exactly the currency round 13 asked for -- and it is known to be as hard as
  factoring (Bach-Charles). The trace formula names `p^11` in its divisor term;
  the q-expansion needs `N` coefficients; mod 2 it is cheap and says only
  whether `N` is an odd square.
- **Four squares.** `r_4(N) = 8 sigma(N)` gives `p + q`. The volume term is the
  average; `p + q` is the fluctuation.

Every source that carries the residue is a function of the divisors. Every
source computable from `N` in polynomial time is symmetric and carries nothing
beyond `N`.

## Round 17: the fourth mechanism

With the central column, Pascal's triangle had been placed on three of the four
mechanisms. The fourth is collision, and the triangle has it too: iterate the
second column. `x -> C(x, 2)` is conjugate to `y -> y^2 - 5/16` by
`x = 2y + 1/2`, so "Pascal rho" is Pollard's rho with a particular constant.

The higher columns give a small, real, and useless surprise. Each row of the
triangle is a palindrome, so `C(x, k)` is a polynomial in `(x - (k-1)/2)^2` for
even `k`. That doubles the fibre statistic `kappa` (measured `2.000` against
`1.000`), and the walks collide `sqrt 2` times sooner -- 300-walk constants agree
with `sqrt(pi/(2 kappa))` to within a standard error for every `k` from 2 to 9.
The saving is exactly consumed by the extra squaring each step needs.

Two bookkeeping errors were caught before anything was committed: the first
step counter skipped Brent's advance loop, which made the constants look better
than Pollard's; and a timing column compared a generic `k`-multiplication
implementation, which measured the implementation, not the method. Evaluations
are now counted in full, and cost is priced in multiplications.

With all four mechanisms placed, the map is complete, and it says one thing:
the triangle hosts every mechanism at that mechanism's baseline. ECM's speed
comes from redrawing groups, Harvey's from Lehman's geometry, and NFS's from
smoothness in number fields. None of them comes from anything the triangle
supplies.
