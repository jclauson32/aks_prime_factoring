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
