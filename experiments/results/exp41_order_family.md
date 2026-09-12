# The order mechanism, priced per draw

Round 35. How many groups each method can draw, how often a draw wins, what a draw costs -- and whether the product predicts the run.

## How many draws each family gets

| family | groups available per prime | if the draw loses |
|---|---|---|
| Pascal row, AKS fold, q-deformation | 1 (the order divides `p^d - 1`) | nothing to do |
| Pollard p - 1 | 1 | nothing to do |
| Williams p + 1 | 2 (the base's discriminant decides whether the order divides `p - 1` or `p + 1`) | nothing to do |
| class groups (Schnorr-Lenstra) | one per multiplier `k` | change `k` |
| ECM / the elliptic triangle | one per curve, filling the Hasse interval of about `4 sqrt p` orders | change the curve |

The split in the last column is the whole distinction between the rigid families and the elliptic ones. It is not about the group being elliptic -- class groups are not -- but about whether a fresh group is available at a price.

## How often one draw wins

All four samples are orders attached to primes of 15 bits, and the test is the one stage 1 actually applies: every prime power at most `B1 = 200`. The last row is the control -- a random integer from the Hasse interval, which is what a curve order would be if it carried no bias of its own. It does carry one, so the control is what measures it.

| what has to be smooth | chance it is 200-powersmooth | samples |
|---|---|---|
| `p - 1` (Pollard) | 0.403 +- 0.009 | 3000 |
| `p + 1` (Williams) | 0.409 +- 0.009 | 3000 |
| `#E` for a Suyama curve (ECM) | 0.596 +- 0.014 | 1187 |
| `#E` for a uniform curve | 0.378 +- 0.014 | 1200 |
| random integer in the Hasse interval (control) | 0.283 +- 0.013 | 1200 |

Ranked by the measured rate: `#E` for a Suyama curve (ECM) (0.596), `p + 1` (Williams) (0.409), `p - 1` (Pollard) (0.403), `#E` for a uniform curve (0.378), random integer in the Hasse interval (control) (0.283).

The comparisons that matter, each computed from the samples above:

- `p - 1` (Pollard) loses to `#E` for a Suyama curve (ECM) (0.403 against 0.596, 11.5 sigma). The two draws are not interchangeable, and the direction is worth pinning down: the curve ECM actually runs wins more often than `p - 1` does, because Suyama's parametrisation hands it a factor of 12 for free (next section).
- `#E` for a uniform curve beats random integer in the Hasse interval (control) (0.378 against 0.283, 5.0 sigma). A uniform curve order is *not* quite a random integer of its size -- it carries small factors more often, the bias Galbraith and McKee describe. Round 31 measures the same gap against a plain random integer at higher power (+5.4 sigma there).
- `#E` for a Suyama curve (ECM) beats `#E` for a uniform curve (0.596 against 0.378, 10.9 sigma). Suyama's parametrisation is not neutral -- see the next section.

## What Suyama's parametrisation is worth

| curves | orders resolved | divisible by 12 |
|---|---|---|
| Suyama (what ECM runs) | 1187 | 1187/1187 |
| uniform `y^2 = x^3 + ax + b` | 1200 | 225/1200 |

The order was counted over `F_p` and the twist resolved by the starting point (13 draws where the point's order divided both were discarded). Suyama's curves carry a rational 12-torsion structure, so a factor of 12 is free on every draw and the cofactor that still has to be smooth is smaller by that much, against 19% of uniform curves getting the same discount. On the smoothness rate that discount was worth +0.218 here (10.9 sigma). Either way it is bought once and never improved on: no parametrisation makes the remaining cofactor smooth more than a constant more often.

## What a draw costs

| one draw | time | what it is |
|---|---|---|
| Pollard p - 1 | 0.22 ms | one exponentiation |
| Williams p + 1 | 0.24 ms | one Lucas chain |
| ECM, one curve | 0.43 ms | a Montgomery ladder |
| elliptic triangle, row lcm(1..B1) | 2.11 ms | the sequence's own double-and-add |

The elliptic draw costs more than the rigid ones -- a curve is a bigger object than a residue -- but it is a draw, and the rigid families have only the one.

## Does the model predict the run?

The table above says a curve wins with probability 0.596 and `p - 1` with probability 0.403. Both are now run on fresh semiprimes, one draw each.

| method, one draw | predicted from the order sample | measured on 400 semiprimes | agree? |
|---|---|---|---|
| ECM, one Suyama curve | 0.596 +- 0.014 | 0.632 +- 0.024 | yes |
| Pollard p - 1 | 0.403 +- 0.009 | 0.400 +- 0.024 | yes |

Per instance, `p - 1` succeeded exactly when `p - 1` was 200-powersmooth on 394 of 400 numbers: the method is barely probabilistic, it is close to a lookup on a property of `p` that was decided before the method started. ECM's draw is the probabilistic one.

The 6 disagreements were chased down rather than waved at. 6 were wins without `p - 1` being powersmooth, and in 6 of those the order of the base `2` modulo `p` was powersmooth even though `p - 1` was not -- what stage 1 needs is the order of the base, and a divisor can be smoother than the number it divides.

## Multiplying it out

| family | chance a draw wins | cost per draw | expected cost of a split |
|---|---|---|---|
| Pollard p - 1 | 0.403 | 0.22 ms | succeeds on 40% of primes; on the rest, never |
| ECM (redrawn per curve) | 0.596 | 0.43 ms | 1.7 curves, 0.7 ms |

That quotient is the whole of ECM's analysis. Raising `B1` raises the cost of a draw roughly linearly and raises the chance it wins by the Dickman factor (round 31); minimising the quotient over `B1` is what produces `L_p[1/2]`, and it is a minimisation over a function whose smooth-number term falls like `u^(-u)`. Nothing in the order mechanism is sub-exponential except through that term.

Shor's algorithm belongs in this table too, in the row that does not exist classically: one group, one draw, no smoothness at all, because it reads the order instead of guessing a multiple of it (round 24). Every classical entry here is paying for not being able to do that.
_Generated by `experiments/exp41_order_family.py` in 8.7s._
