# The 2-adic search tree, and what it costs to finish it

Round 38. Bitwise lifting has no pruning; with Coppersmith it is an `N^(1/4)` algorithm.

## The tree

A node at depth `k` is a pair `(x, y)` of odd residues mod `2^k` with `x y == N`. Its children are the four ways to add a bit to each; the congruence keeps some of them. Round 24 counted the survivors from the residue side -- every odd `r` is a legal `p mod 2^k` -- and the counts below enumerate the tree itself.

| depth k | nodes | children per node | distinct `x` values | odd residues mod 2^k |
|---|---|---|---|---|
| 2 | 2 | 2.00 | 2 | 2 |
| 3 | 4 | 2.00 | 4 | 4 |
| 4 | 8 | 2.00 | 8 | 8 |
| 5 | 16 | 2.00 | 16 | 16 |
| 6 | 32 | 2.00 | 32 | 32 |
| 7 | 64 | 2.00 | 64 | 64 |
| 8 | 128 | 2.00 | 128 | 128 |
| 9 | 256 | 2.00 | 256 | 256 |
| 10 | 512 | 2.00 | 512 | 512 |
| 11 | 1024 | 2.00 | 1024 | 1024 |
| 12 | 2048 | 2.00 | 2048 | 2048 |
| 13 | 4096 | 2.00 | 4096 | 4096 |
| 14 | 8192 | 2.00 | 8192 | 8192 |
| 15 | 16384 | 2.00 | 16384 | 16384 |
| 16 | 32768 | 2.00 | 32768 | 32768 |

Two of the four children survive at every node, so the tree doubles with depth, and the distinct `x` values are *all* the odd residues mod `2^k`. That is the whole difficulty in one line: for every odd `r` there is a `y` with `r y == N (mod 2^k)`, so the congruence says nothing whatever about `p` on its own. There is no branch to cut.

## Every prune available from `N` alone

The prune a bitwise search would actually use is the product bound: a node whose literal `x y` already exceeds `N` cannot be completed downward. Below, the share of nodes it removes, and the share removed by testing `x | N` (which is trial division wearing a hat).

| depth k | nodes | removed by the product bound | nodes with `x | N` |
|---|---|---|---|
| 2 | 2 | 0% | 0 |
| 3 | 4 | 0% | 0 |
| 4 | 8 | 0% | 0 |
| 5 | 16 | 0% | 0 |
| 6 | 32 | 0% | 0 |
| 7 | 64 | 0% | 0 |
| 8 | 128 | 0% | 0 |
| 9 | 256 | 0% | 0 |
| 10 | 512 | 0% | 0 |
| 11 | 1024 | 0% | 0 |
| 12 | 2048 | 0% | 0 |
| 13 | 4096 | 0% | 0 |
| 14 | 8192 | 22% | 2 |
| 15 | 16384 | 66% | 2 |
| 16 | 32768 | 88% | 2 |
| 17 | 65536 | 96% | 2 |
| 18 | 131072 | 99% | 2 |
| 19 | 262144 | 100% | 2 |

The product bound removes nothing until the two halves together exceed `N`, which happens at depth `log2(N)/2` = 13 here -- and the next section needs the tree only to depth `log2(N)/4` = 6. In the range that matters the prune is inactive. The `x | N` column is the other way to end the search: it finds `p` once the tree is deep enough to contain it, having enumerated about `p` nodes to get there, which is trial division with extra steps.

## Finishing the tree with Coppersmith

A node at depth `t` is exactly a claim about `p mod 2^t`, and round 12 showed Coppersmith finishes from such a claim once `2^t` is about `N^(1/4)`. So the tree does not have to be walked to the bottom: it has to be walked to depth `log2(N)/4`, and every leaf handed to a lattice. First, how much depth the lattice really needs, measured on instances where `p` is known.

| bits of N | `log2(N^(1/4))` | extra bits this instance needed | depth used |
|---|---|---|---|
| 28 | 7 | 5 | 12 |
| 34 | 9 | 5 | 14 |
| 40 | 10 | 1 | 11 |

Coppersmith needs a little more than the bare quarter, and how much more is a property of the instance rather than of its size -- 5, 5, 1 bits for these three, because the bound it has to beat is asymptotic and 28 bits is not asymptotic. Each row below walks the same instance whose depth was just measured. Measuring it needed `p`, which the search does not have: a real run would deepen the tree one bit at a time until a leaf worked, which costs about twice the final depth's leaves and changes nothing below.

| bits of N | depth t | leaves in the tree | `log(leaves)/log(N)` | leaves tried | time | how |
|---|---|---|---|---|---|---|
| 28 | 12 | 2,048 | 0.409 | 239 | 0.85 s | run to a split |
| 34 | 14 | 8,192 | 0.396 | 2,145 | 8.11 s | run to a split |
| 40 | 11 | 1,024 | 0.252 | 338 | 0.99 s | run to a split |

The rows at 28 bits, 34 bits and 40 bits were walked to the end: every odd residue mod `2^t` handed to the lattice until one of them factored the number. A full walk stops at whichever of `p` and `q` has the earlier residue, so it runs about four times short of the whole tree on average.

The algorithm is real -- every row factored its number -- and its exponent is the one the size mechanism already had. The `log(leaves)/log(N)` column runs 0.409 to 0.252: the tree is `N^(1/4)` leaves plus the lattice's fixed margin, and the margin is what keeps the measured column above 0.25 at sizes this small. The lattice itself is polynomial, so the whole thing is `N^(1/4)` times a polynomial, and round 34's table lists four other ways to pay exactly that.

## What would have to change

The tree is prunable the moment anything else is known. One fact -- the true `p mod 2^t` -- collapses it from `N^(1/4)` leaves to one, which is the measurement in the first table of section 3. That is the shape of Heninger and Shacham's attack on RSA keys with corrupted bits: their tree is the same tree, and what prunes it is the redundancy between `p`, `q`, `d`, `d_p` and `d_q`, all of which have to agree. Given only `N` there is no second equation, so nothing disagrees with anything, and every branch survives.

So the honest answer to 'can the bits of `p` be searched one at a time': yes, and it terminates, and it costs `N^(1/4)`. A binary search needs a comparison that eliminates half the space; the 2-adic congruence eliminates nothing, because it is satisfiable for every odd residue. What makes a search fast is not that the space is searchable in order -- it is that the order tells you which half to drop.
_Generated by `experiments/exp44_padic_tree.py` in 10.4s._
