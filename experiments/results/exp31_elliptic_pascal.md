# Pascal's triangle over divisibility sequences, and the one that is ECM

Round 23. Four integer triangles, one Sierpinski rule, and a cell width you can redraw.

## 1. Four integer triangles

For a sequence with `gcd(a_m, a_n) = a_gcd(m,n)`, the generalised binomial `[n, k]_a = a_n a_(n-1) .. a_(n-k+1) / (a_1 .. a_k)` is an integer. The fourth sequence is an *elliptic divisibility sequence*: the division polynomials of a point `P` of infinite order, `1, 10, 456, 9440, ...`.

| triangle | integral entries, rows 0-39 |
|---|---|
| Pascal, a_n = n | 820 of 820 |
| Gaussian, a_n = (2^n - 1) | 820 of 820 |
| Fibonomial, a_n = F_n | 820 of 820 |
| elliptic, a_n = psi_n(P) on y^2 = x^3 + 17 | 820 of 820 |

## 2. One Sierpinski rule

Modulo a prime `p`, `[n, k]_a` vanishes exactly when adding `k` and `n - k` carries in the mixed radix `(r, p, p, ...)`, where `r = r(p)` is the rank of apparition -- the first `n` with `p | a_n` (Kummer's theorem, as generalised by Knuth and Wilf). Each triangle mod `p` is a gasket whose first-level cell is `r` wide.

| triangle | primes 5..47 used | entries where the carry rule is right |
|---|---|---|
| Pascal, a_n = n | 10 | 10350 of 10350 |
| Gaussian, a_n = (2^n - 1) | 13 | 13455 of 13455 |
| Fibonomial, a_n = F_n | 12 | 12420 of 12420 |
| elliptic, a_n = psi_n(P) on y^2 = x^3 + 17 | 8 | 8280 of 8280 |

The same prime, `p = 7`, four cells:

**Pascal, a_n = n** mod 7, cell `r = 7`:

```
#
##
###
####
#####
######
#######
#......#
##.....##
###....###
####...####
#####..#####
######.######
##############
#......#......#
##.....##.....##
###....###....###
####...####...####
#####..#####..#####
######.######.######
#####################
#......#......#......#
##.....##.....##.....##
###....###....###....###
####...####...####...####
#####..#####..#####..#####
######.######.######.######
############################
#......#......#......#......#
##.....##.....##.....##.....##
```

**Gaussian, a_n = (2^n - 1)** mod 7, cell `r = 3`:

```
#
##
###
#..#
##.##
######
#..#..#
##.##.##
#########
#..#..#..#
##.##.##.##
############
#..#..#..#..#
##.##.##.##.##
###############
#..#..#..#..#..#
##.##.##.##.##.##
##################
#..#..#..#..#..#..#
##.##.##.##.##.##.##
#####################
#....................#
##...................##
###..................###
#..#.................#..#
##.##................##.##
######...............######
#..#..#..............#..#..#
##.##.##.............##.##.##
#########............#########
```

**Fibonomial, a_n = F_n** mod 7, cell `r = 8`:

```
#
##
###
####
#####
######
#######
########
#.......#
##......##
###.....###
####....####
#####...#####
######..######
#######.#######
################
#.......#.......#
##......##......##
###.....###.....###
####....####....####
#####...#####...#####
######..######..######
#######.#######.#######
########################
#.......#.......#.......#
##......##......##......##
###.....###.....###.....###
####....####....####....####
#####...#####...#####...#####
######..######..######..######
```

**elliptic, a_n = psi_n(P) on y^2 = x^3 + 17** mod 7, cell `r = 13`:

```
#
##
###
####
#####
######
#######
########
#########
##########
###########
############
#############
#............#
##...........##
###..........###
####.........####
#####........#####
######.......######
#######......#######
########.....########
#########....#########
##########...##########
###########..###########
############.############
##########################
#............#............#
##...........##...........##
###..........###..........###
####.........####.........####
```

## 3. The cell is an order

| triangle | cell r(p) | checked on primes 7..397 |
|---|---|---|
| Pascal | `p` itself | 74 of 74 |
| Gaussian (q = 2) | the order of 2 mod p, a divisor of `p - 1` | 74 of 74 |
| Fibonomial | a divisor of `p - (5/p)` | 74 of 74 |
| elliptic | the order of `P` in `E(F_p)`, anywhere in `p + 1 +- 2 sqrt p` | 74 of 74 |

The Pascal triangle's cell *is* `p`: to see the cell is to have factored. The Gaussian and Fibonomial cells are orders, but orders fixed by `p` -- divisors of `p - 1` or `p + 1`, the one ticket of Proposition 14. The elliptic cell is the order of a point, and a different curve is a different triangle.

## 4. Redrawing the triangle

| triangle family | possible cells at p = 10007 | a 50-smooth cell available? |
|---|---|---|
| Gaussian, any q | divisors of p - 1 = 2 x 5003 | no |
| Fibonomial | divisors of p - (5/p) = 2^3 x 3^2 x 139 | only small divisors (139 blocks the rest) |
| elliptic, 300 random curves | 268 distinct orders in [30, 10199] | 74 of the 300 curves |

`p - 1` has the prime factor 5003, so every Gaussian triangle mod 10007 has a cell that is a multiple of 5003 or a divisor of 2; the Fibonomial cells are stuck with the divisors of 10008. The elliptic triangles take 268 different cell widths across the Hasse interval, and 74 of them are `50`-smooth.

How many different cells can each family offer at one prime? Sampling many members of each family at `p = 10007`:

| family | members sampled | distinct group orders (cells divide these) |
|---|---|---|
| Pascal | -- | 1: the cell is `p` |
| Gaussian, q = 2..399 | 398 | 2 cells, all dividing `p - 1`: one order |
| Fibonomial / Lucas | -- | 2: `p - 1` or `p + 1` |
| elliptic with CM, `y^2 = x^3 + b` at p = 10009 | 199 | 6 |
| elliptic, generic | 60 | 57 |

The family with complex multiplication offers only 6 orders -- the six twists of a curve with `j = 0` -- so it is barely more redrawable than `p +- 1`. Generic curves offer a new order almost every time: the Hasse interval holds about `4 sqrt p` of them. That count of tickets per prime is the whole difference between the rigid order methods and ECM.

## 5. Jumping to a row: ECM

Row `M` of an elliptic triangle mod `N` begins with `a_M = W_M`, and `p | W_M` exactly when the cell `r_p` divides `M`. Take `M = lcm(1..B)` -- reached in `O(log M)` steps by the sequence's own double-and-add -- and `gcd(W_M, N)` splits `N` whenever `r_p` is `B`-smooth. If not, draw another curve.

| bits of p | bits of N | curves needed (three semiprimes) |
|---|---|---|
| 20 | 84 | 9 / 4 / 4 |
| 24 | 88 | 3 / 1 / 1 |
| 28 | 92 | 2 / 11 / 5 |
| 32 | 96 | 16 / 2 / 9 |

That is stage 1 of Lenstra's elliptic curve method, read as the original observation of this project: *look along a row of the triangle mod `N` for entries that share a factor with `N`.* In Pascal's triangle the row that works is `N` itself and the entries that work sit at multiples of `p`, so finding them is finding `p`. In the elliptic triangle the row that works is any multiple of the cell width, and the cell width can be redrawn until it is smooth. The gasket the project started from was the right picture; it needed a different sequence under it.
_Generated by `experiments/exp31_elliptic_pascal.py` in 81.6s._
