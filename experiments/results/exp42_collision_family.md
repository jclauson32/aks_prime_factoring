# The collision mechanism, and the price of a hidden group

Round 36. Fibre statistics, cycle finding, and why m machines buy sqrt(m).

## What the map contributes: the fibre statistic

A walk on `x -> C(x, k) mod p` collides after about `sqrt(pi p / (2 kappa))` steps, where `kappa = (1/p) sum_v m_v (m_v - 1)` counts how much the map folds. Round 17 predicted `kappa = 2` for even `k >= 4` (Pascal's reflection closes every fibre) and `kappa = 1` otherwise. Both are recomputed here, and the prediction is checked against measured rho lengths.

The average has to be taken across *maps*, not across starting points within one map. A fixed map has a fixed functional graph, and every start in its giant component reports much the same rho length, so a thousand starts on one map estimate that map's own value to high precision and say nothing about how far it sits from the law. The scatter that matters is between maps, and it is what the last column reports.

| k | kappa predicted | kappa measured (mean over maps) | measured / predicted rho length, over 12 maps | spread between maps |
|---|---|---|---|---|
| 2 | 1 | 1.000 | 0.91 +- 0.08 | 0.27 |
| 3 | 1 | 1.000 | 0.90 +- 0.07 | 0.25 |
| 4 | 2 | 2.000 | 1.06 +- 0.07 | 0.24 |
| 5 | 1 | 1.007 | 0.86 +- 0.07 | 0.23 |
| 6 | 2 | 1.990 | 0.78 +- 0.05 | 0.16 |

`kappa` is barely statistical: averaged over the twelve maps it lands within 0.7% of the predicted integer for every column, including the factor of two for even `k >= 4` that round 17 derived from Pascal's reflection. The rho length follows the law in shape and sits about 10% below it in constant, averaging 0.90 across the columns.

That residual is not all noise: `k = 4` and `k = 6` differ by 3.3 standard errors (1.06 against 0.78). So `kappa` captures the big effect -- which columns fold and which do not -- without capturing everything about these particular maps. The last column is the reason the averaging is done over maps: one map's own functional graph scatters by tens of percent whatever the law says.

Folding the map harder shortens the walk by `sqrt(kappa)` and lengthens each step by the extra multiplications, which is why no column of Pascal's triangle beats the `k = 2` column -- that column being Pollard's rho (round 17).

## What the cycle finder contributes

Same walk, two ways of noticing it has closed. The measure is evaluations of the map, not seconds, so the comparison does not depend on how the gcds are batched.

| cycle finder | median evaluations | numbers split |
|---|---|---|
| Floyd (tortoise and hare) | 1,641 | 40/40 |
| Brent | 1,598 | 40/40 |

Brent's evaluations come to 0.97 of Floyd's here. Floyd pays three evaluations per step to Brent's one, and Brent overshoots the collision by walking in powers of two; the saving is a constant, and a constant is all a cycle finder can ever be worth.

## What `m` machines contribute -- twice

Now the point of the round. Fix a prime `p` and walk `x -> x^2 + c mod p` (the collision is a fact about the walk mod `p`; how it is detected is the variable). Two ways to use `m` machines:

- **Honest**: `m` walks with independent `c`, each waiting to close up on itself. That is what Pollard's rho can detect from `Z/N` alone, because a self-collision is found by the cycle structure, not by comparing values.
- **Oracle**: all `m` walks on one map, with distinguished points stored in a shared table keyed by `x mod p`. Any two walks that meet are caught at the next distinguished point. This is van Oorschot-Wiener, and it needs `p`.

| machines `m` | honest: steps per machine | speedup | honest: total steps | oracle DP: steps per machine | speedup |
|---|---|---|---|---|---|
| 1 | 9,816 | 1.00x | 9,816 | 10,072 | 1.00x |
| 2 | 6,616 | 1.48x | 13,232 | 5,020 | 2.01x |
| 4 | 4,154 | 2.36x | 16,617 | 2,077 | 4.85x |
| 8 | 2,993 | 3.28x | 23,944 | 799 | 12.60x |
| 16 | 2,054 | 4.78x | 32,858 | 571 | 17.65x |

Walks that went 20 gaps without a distinguished point were restarted from a fresh point, as van Oorschot and Wiener prescribe -- without that rule a walk can fall into a short cycle carrying no distinguished point and never report. 0 of 150 runs still hit the step cap and were dropped.

Both columns count steps to the moment a collision exists or is caught; the honest column omits the constant a cycle finder adds and the oracle column includes the distance to the next distinguished point, so the constants are not comparable between columns. The exponent in `m` is, and it is the claim.

Fitted exponents in `m`: honest `m^(-0.57)`, oracle `m^(-1.09)`, both within the noise of a five-point fit of the theoretical `-1/2` and `-1`: the minimum of `m` independent rho lengths is `sqrt(pi p / 2m)`, while pooling the walks makes one birthday problem over all `m` trajectories, whose total is flat in `m`.

The honest column's *total* work grows: 9,816 steps at `m = 1` against 32,858 at `m = 16`. Parallel rho for factoring is not work-efficient; parallel rho for a discrete logarithm is.

## Why the hidden group costs exactly that

The distinguished-point table is a hash table on the collision value. Without `p` there is no value to hash: two iterates `x` and `y` can only be compared by `gcd(x - y, N)`, one pair at a time. So `K` stored points cost `K^2 / 2` gcds to search, while hashing costs `K`. A single walk escapes this because its self-collision is found by cycle structure -- Floyd and Brent compare `O(K)` pairs, not `K^2` -- and that escape is available once per walk, not once per pair of walks.

| how a collision is noticed | comparisons among K points | needs p? |
|---|---|---|
| cycle structure, one walk (Floyd/Brent) | O(K) = 4,096 | no |
| pairwise gcd, any two points | K^2/2 = 8,388,608 | no |
| distinguished points, hashed | O(K) = 4,096 | yes |

That table is the whole difference between `sqrt(m)` and `m`, and it is the same shape as every other result in this project: the mechanism is not the obstacle, the missing residue is.

## Where the four mechanisms land

Rounds 33, 34 and 35 did sign, size and order; with collision the survey is complete. Stripped of smooth numbers, every mechanism costs the same:

| mechanism | best exponent without smoothness | in this repository | with smooth numbers |
|---|---|---|---|
| size | `O~(N^(1/4))` | factorial threshold (T20), hull walk `N^(1/3)` (P23) | -- |
| collision | `N^(1/4)` | Pascal rho = Pollard rho (P22) | -- |
| sign | `N^(1/4)` | central column (round 15), SQUFOF (round 32) | `L[1/2]` quadratic sieve, `L[1/3]` NFS |
| order | `N^(1/4)` | baby-step giant-step over the unknown order, with the differences batched into one gcd -- the same trick as T20, and the one entry in this table not implemented separately here | `L[1/2]` ECM, Pollard `p - 1`, the elliptic triangle |

The two escapes from `N^(1/4)` are both visible in round 33-36's data. Smoothness buys sub-exponential time and is the only thing that ever has. Combining mechanisms buys a little: Harvey's `N^(1/5)` is the hull walk's geometry plus a Fermat congruence plus a collision sweep, and it is the best deterministic exponent known. Neither escape is polynomial, and nothing in this repository suggests a third.
_Generated by `experiments/exp42_collision_family.py` in 2.0s._
