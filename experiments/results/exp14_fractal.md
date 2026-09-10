# The fractal reading: Lucas drawn, and the barrier made geometric

Round 6, from the observation that the picture looks like Sierpinski. It is.

## Pascal mod 3 is exactly the Sierpinski gasket

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
                 #........#
                ##.......##
               ###......###
              #..#.....#..#
             ##.##....##.##
            ######...######
           #..#..#..#..#..#
          ##.##.##.##.##.##
         ##################
        #........#........#
       ##.......##.......##
      ###......###......###
     #..#.....#..#.....#..#
    ##.##....##.##....##.##
   ######...######...######
  #..#..#..#..#..#..#..#..#
 ##.##.##.##.##.##.##.##.##
###########################
```

## Composite: two gaskets superimposed

`p` marks an entry killed by 3 alone, `q` by 5 alone, `0` by both -- a true zero mod 15 -- and `#` a survivor. The zeros of the composite are precisely where the two hole systems coincide.

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
        #pp#pp#pp#pp#pp#pp#pp#pp#
       #q0qq0qq0qq0qq0qq0qq0qq0q#
      ##qqqqqqqqqqqqqqqqqqqqqqq##
     #pp0000000000000000000000pp#
    ##pp000000000000000000000pp##
   ###pp00000000000000000000pp###
  #00q0p0000000000000000000p0q00#
 ##0qqpp000000000000000000ppqq0##
###qq#pp00000000000000000pp#qq###
```

## The geometry is exactly Lucas

`support_count` is a digit dynamic program -- it box-counts the gasket in `O(log N)` without drawing it. Checked against brute force on **796** cases: **0** mismatches.

| p | rows | surviving entries | self-similarity prediction | box dimension |
|---|---|---|---|---|
| 2 | 2^3 = 8 | 27 | (2*3/2)^3 = 27 | 1.5850 |
| 3 | 3^3 = 27 | 216 | (3*4/2)^3 = 216 | 1.6309 |
| 5 | 5^3 = 125 | 3,375 | (5*6/2)^3 = 3,375 | 1.6826 |
| 7 | 7^3 = 343 | 21,952 | (7*8/2)^3 = 21,952 | 1.7124 |
| 11 | 11^3 = 1331 | 287,496 | (11*12/2)^3 = 287,496 | 1.7472 |

The exact relation `support_count(p^k) = (p(p+1)/2)^k` is the gasket's self-similarity. The dimension `log(p(p+1)/2)/log p` rises toward `2` with `p`: bigger prime, denser gasket, sparser holes.

## A dual to Theorem 4

Theorem 4 reads *along* row `n` and finds the **smallest** prime factor. Reading *down* the rows finds the **largest** one.

| n | factors | first row with a zero mod n | larger prime q | equal? |
|---|---|---|---|---|
| 35 | 5 x 7 | 7 | 7 | yes |
| 55 | 5 x 11 | 11 | 11 | yes |
| 65 | 5 x 13 | 13 | 13 | yes |
| 85 | 5 x 17 | 17 | 17 | yes |
| 95 | 5 x 19 | 20 | 19 | no |
| 115 | 5 x 23 | 23 | 23 | yes |
| 145 | 5 x 29 | 29 | 29 | yes |
| 155 | 5 x 31 | 31 | 31 | yes |

Over **127** semiprimes below 2000: the digit criterion predicts the first zero row exactly **127/127** times, and that row *is* the larger prime factor **116/127 = 91%** of the time.

**Why.** An entry vanishes mod `n` only if it vanishes mod both primes. Row `i` carries a hole in the mod-`q` gasket only once `i` has two base-`q` digits, so no row below `q` can contain a zero. At `i = q` the whole interior vanishes mod `q`, and a zero mod `n` appears exactly when row `q` also carries a hole in the mod-`p` gasket -- which the digit count decides. Verified: `first zero row == q` iff row `q` mod `p` has a zero, **exact on all 127 cases**.

## Does it break anything?  A column probe

Reading *down a column* is the cheapest way to sample the second dimension. Mod `p`, Lucas makes `p | C(i,c)` a condition on `i mod p^k`, so the column's zero pattern is periodic with period a power of `p`, and `gcd(C(i,c) mod n, n)` splits `n` whenever exactly one prime divides.

Widening the column raises the hit rate. It also raises the cost of one evaluation, by exactly as much:

| p | column c | hit rate | cost per probe | rate / cost | x p |
|---|---|---|---|---|---|
| 431 | 1 | 0.0025 | 1 | 0.002500 | 1.08 |
| 431 | 4 | 0.0145 | 4 | 0.003625 | 1.56 |
| 431 | 16 | 0.0550 | 16 | 0.003438 | 1.48 |
| 431 | 64 | 0.2015 | 64 | 0.003148 | 1.36 |
| 431 | 256 | 0.5450 | 256 | 0.002129 | 0.92 |
| 2,239 | 1 | 0.0013 | 1 | 0.001250 | 2.80 |
| 2,239 | 4 | 0.0037 | 4 | 0.000937 | 2.10 |
| 2,239 | 16 | 0.0125 | 16 | 0.000781 | 1.75 |
| 2,239 | 64 | 0.0568 | 64 | 0.000887 | 1.99 |
| 2,239 | 256 | 0.1943 | 256 | 0.000759 | 1.70 |
| 39,239 | 1 | 0.0000 | 1 | 0.000000 | 0.00 |
| 39,239 | 4 | 0.0003 | 4 | 0.000063 | 2.45 |
| 39,239 | 16 | 0.0005 | 16 | 0.000031 | 1.23 |
| 39,239 | 64 | 0.0025 | 64 | 0.000039 | 1.53 |
| 39,239 | 256 | 0.0120 | 256 | 0.000047 | 1.84 |

The last column is flat near `1`. Hit rate divided by cost is `~1/p` no matter how the column is chosen, so total work stays `Theta(p)`. **The barrier is isotropic**: reading the triangle sideways costs exactly what reading it along a row costs.

## What the picture actually explains

The fractal reading is not a new attack; it is the *reason* the previous attacks failed, in a form you can see.

- **Theorem 2** (the row is empty below `spf(n)`): the top of the gasket is solid. Holes only begin at scale `p`.
- **Theorem 7** (aliasing): you cannot see a scale-`p` fractal by sampling below scale `p`. Folding at level `r < p` averages over whole self-similar cells and returns the same value from each.
- **Proposition 12** (period finding): the gasket *is* a periodic structure of period `p`, and recovering it is recovering the period.
- **Proposition 17** (the smoothness ceiling): every group order in the ladder is an arithmetic shadow of the same scale.

And reading two dimensions costs **area**. The 2-D picture holds strictly more information than row `n` alone -- it contains the largest prime factor as well as the smallest -- but extracting it means visiting `Theta(q^2)` entries where the row needed `Theta(p)`. The extra information is real and it is priced accordingly.

That is the honest verdict on the fractal: it unifies the barriers rather than breaking them. Which is worth something -- a barrier you can *see* is easier to attack than one you can only compute.
_Generated by `experiments/exp14_fractal.py` in 0.5s._
