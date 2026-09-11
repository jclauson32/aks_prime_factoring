# Evolving factoring programs

Round 28. Selection over straight-line programs, starting from nothing.

Programs of 16 instructions over `Z/N` (`mul`, `add`, `sub`, and `pow` by a prime below 50), held to 40 multiplications. A program scores on `N` when `gcd(out, N)` or `gcd(out - 1, N)` is a proper factor. Training and test sets: 500 semiprimes each, with 12-bit primes, where smooth group orders are common enough for selection to have a gradient.

| program | multiplications | test success | what it computes |
|---|---|---|---|
| random programs (mean of 339) | -- | 0.0026 | -- |
| Pollard p - 1, written by hand | 38.5 | 0.160 | -- |
| Pollard p - 1, exponent optimised on the training set | 38.8 | 0.192 | `a^E - 1`, `E = 2^3 x 3^2 x 5 x 7 x 11 x 13 x 17 x 37` |
| evolved, run 1 (36 s) | 39.4 | 0.190 | `a^19275025 (a^D - 1)`, `D = 2^2 x 3^2 x 5^2 x 7 x 11 x 17 x 19 x 31` |
| evolved, run 2 (35 s) | 38.6 | 0.110 | `a^E`, `E = 2^2 x 3^2 x 5 x 7^2 x 11 x 17` |
| evolved, run 3 (33 s) | 39.8 | 0.124 | `a^E`, `E = 2^3 x 3 x 5 x 7 x 11 x 19 x 37` |

All three runs converged on the same idea: the base raised to an exponent that is a product of small primes (largest prime factor 37), compared with 1 -- either as `a^E` directly or as a difference of two powers, `a^e (a^D - 1)`, whose gcd with `N` is that of `a^D - 1`. That is Pollard's `p - 1`, rediscovered from random programs by selection alone.

Against random programs (0.0026) selection wins by a wide margin. The fair comparison is `p - 1` with its exponent optimised on the same training data (0.192 on the test set): the best evolved program (0.190) matches it (a difference within two standard errors, +-0.035, counts as a match). An early version compared evolution with a hand-written `p - 1` whose exponent spent the budget badly, and one evolved program appeared to beat it; decoded, that program was `p - 1` with a better exponent. Nothing any run found steps outside the order mechanism.

## With Pascal's collision primitive

Round 17 showed that iterating `x -> C(x, 2)` is Pollard's rho. Giving the programs that primitive, and a larger budget, asks whether selection assembles the collision mechanism too: iterate, then compare two iterates.

| program | multiplications | test success | what reaches the output |
|---|---|---|---|
| Floyd's rho on `C(x, 2)`, written by hand | 60 | 0.066 | iterates `C(x, 2)`, multiplies `x_2i - x_i` |
| Pollard p - 1, optimised exponent | 60 | 0.380 | `a^E - 1`, `E = 2^5 x 3^3 x 5 x 7^2 x 11 x 13 x 17 x 19 x 23 x 37` |
| evolved, arithmetic only | 60 | 0.052 | instructions reaching the output: add, mul, pow |
| evolved, with C(x, 2) | 57 | 0.126 | `a^E`, `E = 2^2 x 3^2 x 5^2 x 7 x 17 x 37` |

At a budget of 60 multiplications the collision mechanism is nearly worthless -- a hand-written Floyd walk succeeds on 6.6% of the test numbers, because a collision mod a 12-bit prime needs about `sqrt p` iterations -- while `p - 1` succeeds on 38.0%. The evolved program that had `C(x, 2)` available does not use it: the instructions its output depends on are powers of the base, and it is `p - 1` again. Differences between the two evolved scores are run-to-run variance of the search, not the primitive. The mechanism selection finds first is the one that pays first.

## When `p - 1` cannot pay

Safe primes `p = 2p' + 1` with `p'` above every exponent prime on offer make every order-mechanism program useless: `a^E = 1 (mod p)` needs `p' | E`. Training and test semiprimes come from disjoint pools of 39 and 39 such primes, so no program can memorise them; budget 100 multiplications.

| program | test success |
|---|---|
| Floyd's rho on `C(x, 2)`, written by hand | 0.090 |
| Pollard p - 1, optimised exponent | 0.000 |
| evolved, with C(x, 2) | 0.000 |
| evolved, arithmetic only | 0.003 |

The collision mechanism is now the only one that works -- the hand-written Floyd walk scores 0.090 -- and selection does not find it: the best evolved program scores 0.003. The difference is the shape of the search. Every small prime added to a `p - 1` exponent raises its success a little, so selection can climb to it one mutation at a time. A collision pays nothing until a long, coordinated chain of iterations and differences is in place, and random mutation does not assemble one.

## Given a loop

Straight-line programs cannot say *repeat this*. Loop programs can: a short body updates three state registers in place and is repeated as often as the budget allows, an `acc` instruction multiplies an accumulator by a difference of two registers, and the output is tested against `N` after every repetition, as a real rho implementation does. Rho is four instructions in this language -- `x <- C(x,2)`, `y <- C(C(y,2),2)`, `acc *= x - y`.

| program | test success, safe primes | body |
|---|---|---|
| rho body, written by hand | 0.287 | `[('cb2', 0, 0, 0), ('cb2', 1, 1, 0), ('cb2', 1, 1, 0), ('acc', 0, 1, 0)]` |
| evolved loop, run 1 | 0.420 | `[('sub', 1, 2, 1), ('cb2', 2, 1, 2), ('add', 2, 0, 2), ('sub', 0, 0, 1), ('add', 2, 2, 1)]` |
| evolved loop, run 2 | 0.470 | `[('sub', 1, 1, 2), ('sub', 2, 1, 0), ('cb2', 0, 2, 1)]` |

Now selection finds something that works where `p - 1` cannot. Does it scale? Give each body `3 sqrt p` repetitions on semiprimes of safe primes of growing size:

| primes in | rho body | best evolved body |
|---|---|---|
| 1000-2000 | 1.000 | 0.650 |
| 4000-8000 | 1.000 | 0.467 |
| 16000-32000 | 1.000 | 0.242 |
| 64000-128000 | 1.000 | 0.150 |

With repetitions growing as `sqrt p`, the hand-written rho succeeds on 100% to 100% of the numbers at every size -- the birthday bound at work. The best evolved body goes from 65% at the smallest primes to 15% at the largest. Selection found a trick that pays on the primes it trained on and fades on larger ones; the collision mechanism, which does not fade, is the one it did not assemble. The evolved body has no accumulator: it iterates a map and asks whether the output *hits* 0 or 1 mod `p`. A hit is found with probability growing like `R/p` in `R` steps; a collision, like `R^2/p` -- the birthday paradox. On small primes the two look alike, and selection took the one it could reach.
_Generated by `experiments/exp36_evolve.py` in 259.4s._
