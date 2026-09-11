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
| 26 | divisor count under the hyperbola; the divisor is a hull vertex (P23) | size | one exact walk; `O~(N^(1/3))`; its edges are Lehman certificates | `exp26` |
| 27 | curved hull pieces (degree `d`) | size | would give `N^(1/(d+2))`; `d = 2` floor sums compute class numbers | `exp26` |
| 27b | jump to the divisor vertex by slope (2D integer program) | size | divisor's normal cone is `Theta(N^(-1/3))`, ~1.4x its neighbours'; aiming needs `p` to `N^(1/6)` | `exp26` |
| 28 | low-degree Fourier learner on the bits of `N` | -- | finds planted parity, `(p+q) mod 4`, top bits of `p` (size); nothing below | `exp27` |
| 29 | Schnorr's prime-number lattice (2021) | sign (relations) | residues avoid the primes the vector used; as smooth as same-size integers with the same local law, no more; yield collapses by 44 bits | `exp28` |
| 30 | quadratic sieve (reference) | **sign + smoothness** | pure Python, 110-bit `N` in seconds; overtakes rho at 80 bits; `L[1/2]` | `exp29` |
| 31 | Lenstra's ECM (reference) | **order, redrawn** | first curve 60% on `p` where `p - 1` fails 100%; 64-bit `p` of a 160-bit `N` in seconds; `L_p[1/2]` | `exp30` |
| 32 | Pascal's triangle over divisibility sequences (P24) | order | gasket cell = rank of apparition; elliptic cell = point order, redrawable; row `lcm(1..B)` is ECM | `exp31` |
| 33 | Shor's algorithm, simulated classically | order, **read directly** | factors to `N = 1517` on a 22-qubit register; simulation cost `N^2`, the amplitudes a quantum computer holds for free | `exp32` |
| 34 | dequantising Shor by heavy Fourier coefficients | order | every computable function of `a^x mod N` tested is flat (weight `2-11/r`); Jacobi is heavy but reveals only `r mod 2` | `exp33` |
| 35 | a toy number field sieve (reference) | **sign, two rings** | factors 40-80-bit `N`; at degree 3 its numbers are 9-15 bits larger than the QS's, so the QS wins here | `exp34` |
| 36 | every method on the same numbers | all four | size stops first, then collision; ECM and the sieves go on | `exp35` |
| 37 | evolving straight-line programs from random | order | selection rediscovers Pollard `p - 1` (smooth exponents), matching an optimised `p - 1`; nothing else | `exp36` |
| 38 | counting points on `E(Z/N)` | order | one count factors `N` (40 of 40); equivalent to factoring (Kunihiro-Koyama) | `exp31` |

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

## Round 18: the binary search, done with lattice points

The question from before the restart -- *is there any way to binary search for
the factor?* -- has a factorial answer (Theorem 20, `O~(N^(1/4))`) and, it turns
out, a geometric one. `#{d | N : d <= X}` is the difference of two lattice-point
counts under the hyperbolas `xy = N` and `xy = N - 1`, and each count is exact
from the convex hull of the points above the curve: no lattice point hides
between a convex region's hull and its boundary. A Stern-Brocot stack walks the
hull in `O(N^(1/3) log N)` steps. Then the search disappears altogether: the
divisor points lie on the strictly convex curve `xy = N`, so they are vertices
of the hull of `{xy > N-1}`, and one walk that checks `xy = N` at each point
finds `p` (60 of 60). The counts are correct on 400 of 400 random cases against
an independent reference, and the method factors every test semiprime -- at
`N^(1/3)` in the worst case, which is Lehman's exponent, and slower than rho.

It is not a coincidence of exponents. The two hull edges at the divisor vertex
bracket `p/q`, and each one *is* a Fermat-Lehman certificate:
`(dy q + dx p)^2 - 4 (dx dy) N` is the square `(dy q - dx p)^2`. The hull and
Lehman enumerate these certificates in different orders (two thirds of the
hull's fall outside Lehman's search box), but they are the same Farey
fractions round 10 found, at the same cost.

The interesting part is why it stops there. Each hull edge is an exact linear
piece of `floor(N/y)`; exact curved pieces of degree `d` would need only
`N^(1/(d+2))` of them. At `d = 4` that would beat every known deterministic
factoring method. But the first step up, exact quadratic floor sums, contains
the class number: `sum_{k<p} floor(k^2/p)` is a linear function of `h(-p)`,
checked on 154 of 154 primes. Short arcs are not full periods, so that is a
signpost rather than a proof -- and it points into known hard territory.

## Round 19: a learner on the bits of `N`

The machine-learning version of the question: is there structure in the bits of
`N` that predicts bits of `p`? The learner correlates all 1,000 Walsh characters
of degree at most two (three among the top bits) with a target bit, on 45,000
balanced 40-bit semiprimes, and predicts on 15,000 more.

It is calibrated both ways. It finds a planted parity at 100%; it finds
`(p+q) mod 4` at 100%, because `N mod 4` genuinely determines it (the `ell = 4`
case of round 13: symmetric information, not a leak); and it finds the top bits
of `p` partially (0.82 against a 0.75 majority), because they are the size of
`N` seen through a square root. For bits 10 down to 1, and for `p mod 3` and
`p mod 4`, the gain over the majority is zero and the strongest correlation sits
at the maximum expected from noise.

One bookkeeping trap caught on the way: bit 0 of `N` is always 1, so its
"character" is a constant and its "correlation" is the target's bias. It
reported `0.73` for an unbalanced target before it was excluded.

## Round 20: Schnorr's lattice, tested

In 2021 Schnorr claimed that lattice reduction finds the smooth relations of a
sieve fast enough to break RSA. The claim did not survive scrutiny at the time;
this round rebuilds its core with an exact integer LLL and Kannan's embedding
and asks the one question that decides it: are the residues `r = u - vN` from
close vectors in the prime-number lattice smooth more often than chance?

They are not, and the reason is structural. `u` and `v` use disjoint primes, and
`r` is divisible by none of them, so a vector that uses half the factor base
leaves `r` the other half. At the unused primes `u` and `v` are units, which
makes `p | r` slightly *more* likely than for a random integer (`1/(p-1)`
against `1/p`; at `p = 2`, always). Against integers of the same size drawn
from exactly that local law, the lattice's smooth count is within a standard
deviation. The method does factor a 32-bit number; its yield per lattice falls
more than tenfold by 40 bits and to nothing at 44, because its residues are
larger than `N`.

The baseline took three tries: same-size integers predict three times too many
smooth residues; adding only coprimality and parity undershoots by about a
fifth, which in one run looked like a z = +3 lattice advantage. The full
local law removes it. The experiment now computes all three live.

## Round 21: the frontier, in one codebase

A single-polynomial quadratic sieve in pure Python factors 110-bit semiprimes
in seconds. On the same numbers it overtakes rho at 80 bits, and its running
time grows about 2.7x per ten bits against rho's 6.0x (theory: 5.7). The sieve
and the central column of round 15 use the same mechanism, the sign of a
square root. The difference is that the sieve never asks for a Legendre symbol
mod `p`: it lets smooth numbers and linear algebra over `F_2` produce the
square. That is where all the sub-exponential methods get their power, and it
is exactly what Pascal's triangle does not supply.

## Round 22: the order mechanism, redrawn

The Pascal row, the AKS fold, the norm-one torus and the q-deformation all live
in groups whose order mod `p` is `p^d - 1` or a divisor of it: fixed by `p`.
Proposition 14 called this "one ticket per prime". Lenstra's ECM is the same
mechanism with the ticket redrawn: each curve has its own order `p + 1 - t`.

On 200 semiprimes with a 28-bit `p`, Pollard's `p - 1` succeeds on 109 of the
110 primes whose `p - 1` is smooth to its bounds and on none of the other 90,
and four fresh bases rescue none of its 91 failures. With the same bounds,
ECM's first curve succeeds on 60% and 62% of the two classes: the curve does not
know `p - 1` exists. With `N` fixed at 160 bits, ECM found 64-bit factors in
about five seconds of pure Python, its cost growing 3.5x per 8 bits of `p`
where rho's grows 16x.

Two bugs were caught on the way. Stage 2 of `p - 1` stepped `a^(q - B1 + 1)`
instead of `a^q`, which undercounted its successes by forty percent. And the
first scaling table stopped at 52-bit `p`, where ECM mostly wins on the first
curve, so the times were flat and demonstrated nothing about the exponent. It
now runs to 64 bits, and the conclusion is computed from the fit.

With this round every mechanism of round 14 has its best practical method in
the repository, run on the same numbers: ECM (order), Harvey (size), rho
(collision) and the quadratic sieve (sign). Pascal's triangle reaches each
mechanism, but never these methods, because what makes them fast (redrawn
groups, Lehman's geometry, smooth numbers) is not in the triangle.

## Round 23: the triangle that is ECM

Pascal's triangle is the case `a_n = n` of a general construction: for any
strong divisibility sequence, `a_n .. a_(n-k+1) / (a_1 .. a_k)` is an integer.
With `a_n = (q^n - 1)/(q - 1)` it is the Gaussian triangle of round 12; with
Fibonacci numbers, the Fibonomials; with an elliptic divisibility sequence --
the division polynomials of a point on a curve -- an elliptic triangle.

Modulo a prime `p` every one of them is a Sierpinski gasket, and the Kummer-type
carry rule in the mixed radix `(r, p, p, ...)` predicts every zero we tested.
Only the width `r` of the first cell changes. In Pascal's triangle it is `p`,
which is why round 1 found the first non-zero entry at `spf(n)`: the cell is the
answer. In the Gaussian and Fibonomial triangles it divides `p - 1` or `p + 1`.
In the elliptic triangle it is the order of the point mod `p` -- 268 different
widths from 300 random curves at `p = 10007`, where every Gaussian cell is stuck
with the factor 5003 of `p - 1`.

Jumping to row `lcm(1..B)` of an elliptic triangle mod `N`, with the sequence's
own double-and-add, and taking a gcd is exactly ECM stage 1; it factors 32-bit
`p` out of 96-bit `N` in tens of curves. The opening observation of this
project -- look along a row of the triangle mod `N` for entries sharing a
factor with `N` -- was a real mechanism. It needed a triangle whose cell width
could be redrawn.

## Round 24: Shor, simulated

The only known polynomial-time factoring algorithm is quantum, and it uses the
order mechanism. The difference from every classical order method -- `p - 1`,
ECM, the Pascal and elliptic triangles -- is that it *reads* the order `r`
instead of exponentiating blindly by a multiple of it, so it never needs `r` to
be smooth. Simulated exactly, with the first register held as `Q = 2^t >= N^2`
amplitudes, a Fourier transform and continued fractions, it factors every test
number up to `1517 = 37 x 41` through the period finding alone.

The simulation pays for every amplitude: its time per run grows like `N^2`,
which is worse than trial division. That is the whole shape of the gap. A
quantum computer holds the `Q` values of `a^x` in `2 log2 N` qubits and extracts
`r` by interference; a classical one has to write them down. Twenty-three
classical rounds found no way to read `r` without them, and by Miller's
reduction any such way would factor `N`.

## Round 25: can the interference be done classically?

Shor's measurement reads `r` from the Fourier spectrum of a function of
`a^x mod N`. Classical algorithms can find heavy Fourier coefficients too --
Goldreich-Levin, Kushilevitz-Mansour -- in about `1/tau` queries for weight
`tau`, so the question is whether any computable function of `a^x mod N` puts
`1/poly` weight on a frequency that reveals `r`.

None does. Additive and quadratic characters, low, middle and high bits, and a
residue symbol mod 7 all have their heaviest revealing coefficient at 2 to 11
times the uniform weight `1/r`, the size of the maximum of a random spectrum,
with orders from a few hundred to thirty thousand. The Jacobi symbol is the
exception, with all its weight at `j = r/2`, and it reveals only that `r` is
even. Even Shor's own collapsed register -- the indicator of one coset of
`<a>`, which *is* classically computable -- is perfectly flat over its `r`
peaks. What the quantum measurement provides is a sample from a spread-out
distribution, not a heavy coefficient a classical search could find.

Two claims were corrected before commit. The collapsed register's indicator had
been called uncomputable; it is computable, and the experiment now measures its
flat spectrum instead. And a test function `e(c/v)` turned out to be the
additive character with time reversed, so it was replaced.

## Round 26: the number field sieve

The last best-known method missing from the repository was the number field
sieve, the sign mechanism at `L[1/3]`. The toy here is complete: base-`m`
polynomials ranked by a crude root score, a two-sided line sieve over `(a, b)`,
rational and algebraic factor bases, quadratic characters, GF(2) elimination,
and the square root in `Z[alpha]` by Newton lifting from an inert prime. It
factors every test semiprime from 40 to 80 bits.

It is slower than the quadratic sieve at every size here, and the experiment
shows why rather than asserting it: at degree 3 the two numbers it needs
smooth are together 9 to 15 bits larger than the sieve's `Q(x)`. Its advantage
needs the degree to grow with `N`, and appears in practice near a hundred
decimal digits.

The square root had a real bug on the way: the Newton lift estimated its
`p`-adic precision from bit lengths, undershot by one power of `p` once the
exponent reached 16, and then silently failed on every dependency of one
40-bit number. Planted squares with 60-bit coefficients exposed it; the lift
now tracks the exponent exactly and recovers 2000-bit square roots. Two
poorly chosen polynomials also starved the sieve of relations, which is why
polynomial selection now ranks candidates by their roots modulo small primes.

## Round 28: evolution rediscovers `p - 1`

Round 14 showed that *random* straight-line programs over `Z/N` are worth
exactly the chance of hitting a zero divisor. This round let selection loose on
them: populations of short programs (multiplications, additions, subtractions,
powers by small primes, under a fixed multiplication budget) scored on how often
their output shares a factor with `N`, trained on 500 semiprimes and tested on
500 others, starting from random programs with nothing seeded.

Every run converged on the same idea -- the base raised to a product of small
primes, compared with 1 -- sometimes directly as `a^E`, once as a difference of
two powers `a^e (a^D - 1)`. That is Pollard's `p - 1`, rediscovered by selection
alone. The best evolved program scores 0.190 on the test set against 0.192 for
`p - 1` with its exponent optimised on the same training data: a match, 70
times better than random programs, and nothing outside the order mechanism.

One comparison had to be corrected first. Against a hand-written `p - 1` that
spent its budget on a poor exponent, one evolved program seemed to win (0.190
against 0.160). Decoded symbolically, it *was* `p - 1`, with a better exponent;
the experiment now optimises the baseline's exponent before comparing.

## Round 27: every method, the same numbers

Thirteen methods, three balanced semiprimes per size from 40 to 100 bits, each
method stopped once it projects past fifteen seconds (`exp35`). The size methods
stop first -- the Pascal row scan and Harvey by 50 bits here, Strassen by 60,
the hyperbola walk by 70 -- then collision, with rho reaching 90. The rigid
order methods (`p - 1`, `p + 1`) turn into fractions as soon as `p` is too large
for `p +- 1` to be smooth by luck. What keeps going is what is powered by
smooth numbers: ECM, which redraws its groups, and the quadratic sieve; the toy
number field sieve reaches 70 bits. Among the Pascal-native constructions the
two that go furthest are Pascal rho, which is Pollard's rho, and the elliptic
triangle, which is ECM with its sequence changed.

The first run of the benchmark had no hard timeout, and one call of Harvey's
method at 60 bits ran for almost five minutes; every call now has a sixty-second
alarm.

## Round 29: the cells, counted

The cell widths of the elliptic triangles divide the group orders `#E(F_p)` and
`#E(F_q)`, and `#E(Z/N)` is their product. An oracle for that single number
factors `N`: `[#E(Z/N)] P = O` modulo both primes, and dividing out a small prime
that divides only one of the two orders leaves a multiple that kills `P` modulo
one prime but not the other. One count each factored 40 of 40 semiprimes. With
Schoof's algorithm for the other direction, point counting over `Z/N` is
equivalent to factoring (Kunihiro and Koyama, 1998) -- another entry on the list
that began with `phi(N)`: quantities that carry the factorisation in plain
sight and cost as much as the factorisation to compute.
