# Generic algebra is worth exactly nothing; one alignment is worth everything

Round 14. A measurement of the algorithm space, not another algorithm.

## The generic baseline

A straight-line program over `Z/N`: start from a random base, apply `+`, `-`, `*` in any arrangement. Because `*` can square a register, depth `d` reaches exponents up to `2^d` -- **this model contains Pollard's `p-1`**. The question is what a *random* such program achieves.

| p | q | program depth | measured split rate | 1 - phi(N)/N | ratio |
|---|---|---|---|---|---|
| 179 | 419 | 6 | 0.005617 | 0.007960 | 0.706 |
| 179 | 419 | 12 | 0.006742 | 0.007960 | 0.847 |
| 179 | 419 | 24 | 0.007792 | 0.007960 | 0.979 |
| 733 | 761 | 6 | 0.002083 | 0.002677 | 0.778 |
| 733 | 761 | 12 | 0.002300 | 0.002677 | 0.859 |
| 733 | 761 | 24 | 0.002717 | 0.002677 | 1.015 |
| 5849 | 7927 | 6 | 0.000233 | 0.000297 | 0.785 |
| 5849 | 7927 | 12 | 0.000200 | 0.000297 | 0.673 |
| 5849 | 7927 | 24 | 0.000283 | 0.000297 | 0.954 |

The rate converges to `1 - phi(N)/N = 1/p + 1/q`: **exactly the chance that an element drawn uniformly at random from `Z/N` happens to be a zero divisor.** Depth buys enormous exponents and buys nothing else. Generic algebra is worth precisely as much as reaching into a hat.

## What one alignment is worth

| lcm bound B | Pollard p-1 split rate | generic baseline | gain |
|---|---|---|---|
| 50 | 0.0787 | 0.000427 | 184x |
| 100 | 0.0787 | 0.000427 | 184x |
| 200 | 0.9947 | 0.000427 | 2,332x |

Same ring. Same operations. Same depth budget -- `lcm(1..B)` is reached by repeated squaring like anything else. The **only** difference is that the exponent is chosen to be divisible by everything small, i.e. to align with the *order* of the group mod `p`.

That single choice is worth three to four orders of magnitude, and it is the entire content of the method. Generic algebra cannot make it: the order of the group mod `p` is not visible to `+`, `-`, `*` on `Z/N`.

## The taxonomy

Every working factoring method manufactures its zero divisor through one of a small number of mechanisms. **An earlier version of this experiment claimed there were only two, order and size. That was wrong** -- it omitted Pollard's rho and the whole quadratic/number field sieve family, and the latter beats both caps it listed. The corrected table:

| mechanism | how the zero divisor arises | methods | best known |
|---|---|---|---|
| **order** | exponent divisible by `\|G\|` for a group `G` attached to `p` | Pollard `p-1`, Williams `p+1`, ECM, class groups | `L[1/2]` (ECM) -- capped by smooth-number density |
| **size** | a window or interval that contains `p` | trial division, Fermat, Lehman, Strassen, Coppersmith, Harvey | `N^(1/5)` deterministic (Harvey) |
| **collision** | two iterates agreeing mod `p` but not mod `q` | Pollard rho | `N^(1/4)` -- the birthday bound on a set of size `p` |
| **sign** | a square root whose CRT signs differ between `p` and `q` | Dixon, CFRAC, QS, NFS; the central column of round 15 | `L[1/3]` (NFS, heuristic) |

The sign mechanism is the strongest known, and it is worth being precise about *where* its cost lives. Given a congruence `x^2 = y^2 (mod N)` with `x != +-y`, the split is a single gcd. Finding the congruence is the hard part, and the sieves do it through smooth relations -- so their running time is set by smoothness density, just in a smaller-number setting than the order methods. By Rabin's reduction, producing square roots mod `N` on demand is exactly as hard as factoring.

Mapped onto the fourteen rounds: the Pascal row, the AKS fold, the q-deformation, the norm-one subgroup and the class group are **order**; the gasket threshold, the factorial search, Lehman's fan and Coppersmith's window are **size**. None of them was collision or sign -- which is part of why none of them approached `L[1/3]`.

## What a polynomial-time algorithm has to do

Every mechanism above has a super-polynomial cost, and in each one the cost sits in the same place: **finding** the structure that makes the zero divisor, not using it. A smooth group order, a window containing `p`, a collision, a mixed-sign square root -- once any of these is in hand, the factor falls out in polynomial time. A polynomial-time algorithm needs one of them to be *findable* in polynomial time, or a mechanism not on the list.

Every barrier in this repository is a statement about one of those finding steps:

- symmetric constructions (Theorem 8, the norm; round 13, `N mod ell`) find nothing, because they see `N` and not `p`;
- aliased constructions (Theorem 7) cannot locate a window below the scale `p`;
- rigid group orders (Proposition 14) offer one order per `p` and no way to redraw;
- dense families (round 5) redraw freely and run into smoothness density instead.

This is not a proof that no polynomial-time mechanism exists. It is a measurement that fourteen rounds of algebraic, geometric, group-theoretic, fractal, lattice and analytic attempts produced none, and a question precise enough to judge the next attempt in a sentence: *which mechanism does it use, and why would its finding step be cheaper than the best known one?*
_Generated by `experiments/exp22_generic.py` in 12.8s._
