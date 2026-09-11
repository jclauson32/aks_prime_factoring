# Factoring by counting lattice points under a hyperbola

Round 18. A binary-searchable predicate with no factorial in it, and where it stops.

## 1. The count, exactly

The lattice points above `xy = N` form a convex set, so none lies between its hull and the curve, and the number of points under the curve in row `y` is read off the hull edge crossing that row. A Stern-Brocot stack finds the next hull vertex in amortised `O(1)` steps.

Against an independent `O(sqrt N)` reference: **400 of 400** random `(N, X)` pairs agree.

| N | hull vertices | steps | steps / N^(1/3) | steps / (N^(1/3) ln N) | rows covered | time |
|---|---|---|---|---|---|---|
| ~1e8 | 3268 | 6815 | 13.39 | 0.716 | 11,483 | 0.002 s |
| ~1e10 | 19072 | 38983 | 17.00 | 0.732 | 109,816 | 0.011 s |
| ~1e12 | 117635 | 232592 | 21.04 | 0.753 | 1,162,497 | 0.075 s |
| ~1e14 | 613689 | 1203451 | 24.78 | 0.766 | 10,701,944 | 0.403 s |
| ~1e16 | 3560903 | 6943468 | 28.86 | 0.776 | 118,025,975 | 2.484 s |

The step count grows like `N^(1/3) log N`: the fifth column moves by +8% across eight orders of magnitude. Approximating the hyperbola by segments of rational slope is also how Voronoi (1903) got his `x^(1/3)` error term in the divisor problem; here the segments are exact rather than approximate.

## 2. The predicate, and one pass instead of a binary search

| X | #{d | 2021027 : d <= X} |
|---|---|
| 1007 | 1 |
| 1008 | 1 |
| 1009 | 2 |
| 1010 | 2 |
| 1421 | 2 |

Binary search on `X` would cost `log N` pairs of walks. It is not needed. The hull of `{xy > N - 1}` contains every lattice point of `{xy > N}` plus the divisor points `(N/d, d)`, and those lie on the strictly convex curve `xy = N`, so each is an extreme point -- a *vertex* of the hull. One walk of that hull, testing `x * y = N` at each point, finds `p`.

| bits of N | p found | steps / N^(1/3) | hyperbola walk | Pollard rho |
|---|---|---|---|---|
| 32 | 4 of 4 | 0.8 | 0.000 s | 0.0001 s |
| 38 | 4 of 4 | 0.8 | 0.001 s | 0.0002 s |
| 44 | 4 of 4 | 1.2 | 0.007 s | 0.0006 s |
| 50 | 4 of 4 | 1.4 | 0.034 s | 0.0021 s |

The walk stops at the divisor, so its cost is the number of hull vertices between rows `p` and `sqrt N` -- a constant times `N^(1/3)` per factor of two in `sqrt(N)/p`. Balanced semiprimes (`q < 2p` above) are its best case; the worst is `p` near `N^(1/3)`:

| bits of N | bits of p (target) | steps / N^(1/3) | log2(sqrt(N)/p) |
|---|---|---|---|
| 44 | 22 | 0.46 | 0.1 |
| 44 | 18 | 14.45 | 4.4 |
| 44 | 16 | 19.36 | 6.2 |
| 50 | 25 | 1.13 | 0.3 |
| 50 | 20 | 17.89 | 5.3 |
| 50 | 18 | 23.37 | 7.4 |

Away from `sqrt N` the cost is 3.2 `N^(1/3)` steps per factor of two in `sqrt(N)/p`, as the curvature count predicts.

Exactly and deterministically, then, at `O~(N^(1/3))` in the worst case -- Lehman's exponent, and worse than the `N^(1/4)` of Strassen and rho, which pull away in the last column of the table before.

### The divisor vertex and Lehman

| check | count |
|---|---|
| divisor point is a hull vertex | 60 of 60 |
| its two edges bracket the slope p/q | 60 of 60 |
| edge (dx, dy) satisfies (dy q + dx p)^2 - 4(dx dy)N = (dy q - dx p)^2 | 120 of 120 |
| ... with (k, c) = (dx dy, dy q + dx p) inside Lehman's search box | 43 of 120 |

The edges at the divisor are rational approximations of `p/q` from both sides, and each is literally a Fermat-Lehman certificate. But the hull and Lehman enumerate differently: Lehman walks `k <= N^(1/3)` and a short `c`-window for each; the hull walks the curve's own best approximations at the local curvature, and its certificates fall outside Lehman's box 77 times in 120. Same objects -- Farey fractions near `p/q`, the fan round 10 found in Lehman -- and the same exponent.

### Random access by slope does not help

The hull vertex supporting a given slope can be found without walking -- it minimises a linear form over the lattice points of a convex region, a two-dimensional integer program. So if the divisor vertex had a wide normal cone, random slopes would hit it quickly. It does not:

| bits of N | divisor cone / median nearby cone | percentile | divisor cone x N^(1/3) |
|---|---|---|---|
| 32 | 1.31 | 70% | 1.45 |
| 40 | 1.29 | 68% | 1.71 |
| 48 | 1.16 | 64% | 1.54 |

The divisor's cone is a little wider than its neighbours' -- it sits exactly on the curve -- but it scales as `N^(-1/3)` like theirs, so random slopes need `Theta(N^(1/3))` tries. Aiming instead requires `p/q` to within `N^(-1/3)`, i.e. `p` to within `N^(1/6)`, and Coppersmith finishes from far less.

## 3. Where it stops: curved pieces and class numbers

Each hull edge is an exact *linear* piece of `floor(N/y)`, summed in `O(1)`. If exact degree-`d` pieces could be summed as cheaply, the piece count per dyadic block would fall from `N^(1/3)` to `N^(1/(d+2))`: a degree-`d` Taylor piece of length `h` at `x` misses `N/y` by `N h^(d+1) / x^(d+2)`, and `O(1)` stray lattice points per piece allows `h ~ x N^(-1/(d+2))`.

| piece degree d | pieces | exponent | same exponent as |
|---|---|---|---|
| 1 | N^(1/3) | 0.333 | Lehman; this walk |
| 2 | N^(1/4) | 0.250 | Strassen, Pollard rho |
| 3 | N^(1/5) | 0.200 | Harvey (deterministic) |
| 4 | N^(1/6) | 0.167 | nothing known -- below every deterministic method |
| 5 | N^(1/7) | 0.143 | nothing known |

The obstruction arrives at `d = 2`. By Dirichlet's class number formula, the full-period quadratic floor sum is a class number in disguise:

```
sum_{k<p} floor(k^2/p)  =  (p-1)(2p-1)/6  -  (p-1-2 h(-p))/2      (p = 3 mod 4)
```

(checked for **154 of 154** primes below 2000). No polynomial-time algorithm for `h(-p)` is known, so exact quadratic floor sums in polylogarithmic time would be news in their own right. The pieces the hyperbola needs are short arcs, not full periods, so this is a warning rather than a proof -- but it puts the next step in a known hard neighbourhood.

A second requirement hides in "`O(1)` stray points per piece": the linear walk never has strays, because hull edges are chosen so that no lattice point lies between edge and curve. A curved analogue needs the lattice points *near* a short parabolic arc located exactly -- small fractional parts of a quadratic sequence, a Diophantine problem with no known Euclid-like algorithm.
_Generated by `experiments/exp26_hyperbola.py` in 13.5s._
