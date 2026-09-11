"""Round 28: evolving factoring programs.

Round 14 measured *random* straight-line programs over ``Z/N`` and found them
worth exactly the zero-divisor density ``1 - phi(N)/N``.  Here the programs are
*evolved*: a population of short programs is scored on how often their output
shares a factor with ``N`` over a training set of semiprimes, and the better
ones breed.  If structure beyond the known mechanisms were within reach of
short programs, selection is the cheapest way to find it.

A program is a list of instructions over registers ``r0 = a`` (a random base
per run), ``r1 = 1`` and the results of earlier instructions:

* ``("mul", i, j)`` -- one multiplication;
* ``("add", i, j)``, ``("sub", i, j)`` -- free;
* ``("pow", i, c)`` -- ``r_i^c`` for a small prime ``c``, costing ``log2 c``
  multiplications.

The score of a program on ``N`` is 1 if ``gcd(out, N)`` or ``gcd(out - 1, N)``
is a proper factor, where ``out`` is its last register.  Programs are held to a
multiplication budget, so the comparison with Pollard's ``p - 1`` at the same
budget is fair.
"""

from __future__ import annotations

import math
import random
from math import gcd

SMALL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)


def cost(program) -> float:
    total = 0.0
    for op, _, arg in program:
        if op in ("mul", "cb2"):
            total += 1
        elif op == "pow":
            total += math.log2(arg) + 1
    return total


def run(program, n: int, a: int) -> int:
    regs = [a % n, 1]
    for op, i, arg in program:
        x = regs[i]
        if op == "mul":
            regs.append(x * regs[arg] % n)
        elif op == "add":
            regs.append((x + regs[arg]) % n)
        elif op == "sub":
            regs.append((x - regs[arg]) % n)
        elif op == "cb2":
            regs.append(x * (x - 1) * ((n + 1) // 2) % n)      # C(x, 2), n odd
        else:
            regs.append(pow(x, arg, n))
    return regs[-1]


def splits(program, n: int, a: int) -> bool:
    out = run(program, n, a)
    for g in (gcd(out, n), gcd(out - 1, n)):
        if 1 < g < n:
            return True
    return False


def score(program, cases) -> float:
    """Fraction of ``(n, a)`` cases the program splits."""
    return sum(splits(program, n, a) for n, a in cases) / len(cases)


OPS = ("mul", "mul", "add", "sub", "pow", "pow")
PASCAL_OPS = OPS + ("cb2", "cb2")


def random_instruction(k: int, rng, ops=OPS):
    op = rng.choice(ops)
    i = rng.randrange(k)
    if op == "pow":
        return (op, i, rng.choice(SMALL_PRIMES))
    if op == "cb2":
        return (op, i, 0)
    return (op, i, rng.randrange(k))


def random_program(length: int, rng, ops=OPS):
    prog = []
    for t in range(length):
        prog.append(random_instruction(t + 2, rng, ops))
    return prog


def repair(program):
    """Clamp register references so every instruction reads earlier registers."""
    out = []
    for t, (op, i, arg) in enumerate(program):
        k = t + 2
        i %= k
        if op not in ("pow", "cb2"):
            arg %= k
        out.append((op, i, arg))
    return out


def mutate(program, rng, ops=OPS):
    prog = list(program)
    t = rng.randrange(len(prog))
    choice = rng.random()
    if choice < 0.5:
        prog[t] = random_instruction(t + 2, rng, ops)
    elif choice < 0.8 and prog[t][0] == "pow":
        op, i, _ = prog[t]
        prog[t] = (op, i, rng.choice(SMALL_PRIMES))
    else:
        # point the last instruction at a different register
        op, i, arg = prog[-1]
        prog[-1] = (op, rng.randrange(len(prog) + 1), arg)
    return repair(prog)


def crossover(a, b, rng):
    cut = rng.randrange(1, min(len(a), len(b)))
    return repair(a[:cut] + b[cut:])


def pminus1_program(budget: float):
    """Pollard's ``p - 1`` as a program: ``a^(prod of small prime powers) - 1``
    within the multiplication budget (the output is compared with 1)."""
    prog, reg, spent = [], 0, 0.0
    for c in SMALL_PRIMES:
        e = 1
        while c ** (e + 1) <= 64:
            e += 1
        for _ in range(e):
            step = math.log2(c) + 1
            if spent + step > budget:
                return prog
            prog.append(("pow", reg, c))
            reg = len(prog) + 1
            spent += step
    # spend what is left on further powers of 2 and 3
    c = 2
    while spent + math.log2(c) + 1 <= budget:
        prog.append(("pow", reg, c))
        reg = len(prog) + 1
        spent += math.log2(c) + 1
        c = 3 if c == 2 else 2
    return prog


def evolve(cases, length: int, budget: float, generations: int, population: int, rng,
           history: list | None = None, ops=OPS):
    """Elitist genetic search; returns the best program found."""
    def fit(p):
        return score(p, cases) if cost(p) <= budget else -1.0

    pop = [random_program(length, rng, ops) for _ in range(population)]   # no seeding
    for g in range(generations):
        scored = sorted(((fit(p), p) for p in pop), key=lambda t: -t[0])
        if history is not None:
            history.append(scored[0][0])
        elite = [p for _, p in scored[: population // 5]]
        pop = list(elite)
        while len(pop) < population:
            if rng.random() < 0.5:
                child = crossover(rng.choice(elite), rng.choice(elite), rng)
            else:
                child = mutate(rng.choice(elite), rng, ops)
            pop.append(child)
    return max(pop, key=fit)


def exponent_of(program):
    """If the program only multiplies and powers ``r0``, the output is ``a^E``;
    return ``E`` (else ``None``)."""
    exps = [1, 0]
    for op, i, arg in program:
        if op == "mul":
            exps.append(exps[i] + exps[arg])
        elif op == "pow":
            exps.append(exps[i] * arg)
        else:
            return None
    return exps[-1]


def floyd_program(budget: float):
    """Pascal rho as a straight-line program: iterate ``x -> C(x, 2)`` from ``a``
    and multiply the differences ``x_(2i) - x_i`` (Floyd) within the budget."""
    prog = []
    xs = {0: 0}                  # iterate index -> register
    reg_of = lambda i: xs[i]
    spent, i = 0.0, 0
    acc = None
    while True:
        # extend the chain to 2(i+1)
        need = 2 * (i + 1)
        while max(xs) < need:
            if spent + 1 > budget:
                return prog
            prog.append(("cb2", xs[max(xs)], 0))
            xs[max(xs) + 1] = len(prog) + 1
            spent += 1
        i += 1
        prog.append(("sub", xs[2 * i], xs[i]))
        diff = len(prog) + 1
        if acc is None:
            acc = diff
        else:
            if spent + 1 > budget:
                return prog
            prog.append(("mul", acc, diff))
            acc = len(prog) + 1
            spent += 1


def as_polynomial(program):
    """The output as a sparse polynomial ``{exponent: coefficient}`` in ``a``."""
    regs = [{1: 1}, {0: 1}]

    def mul(p, q):
        out = {}
        for e1, c1 in p.items():
            for e2, c2 in q.items():
                out[e1 + e2] = out.get(e1 + e2, 0) + c1 * c2
        return {e: c for e, c in out.items() if c}

    for op, i, arg in program:
        x = regs[i]
        if op == "mul":
            regs.append(mul(x, regs[arg]))
        elif op in ("add", "sub"):
            y = regs[arg]
            s = 1 if op == "add" else -1
            out = dict(x)
            for e, c in y.items():
                out[e] = out.get(e, 0) + s * c
            regs.append({e: c for e, c in out.items() if c})
        else:
            r = {0: 1}
            for _ in range(arg):
                r = mul(r, x)
                if len(r) > 64:
                    return None
            regs.append(r)
        if regs[-1] is None or len(regs[-1]) > 64:
            return None
    return regs[-1]


def pminus1_success(exponent: int, cases) -> float:
    """Fraction of ``(n, a)`` with ``gcd(a^E - 1, n)`` a proper factor."""
    hits = 0
    for n, a in cases:
        g = gcd(pow(a, exponent, n) - 1, n)
        hits += 1 < g < n
    return hits / len(cases)


def optimised_pminus1(cases, budget: float):
    """Pollard ``p - 1`` with its exponent chosen greedily on ``cases``: add the
    prime whose extra power buys the most success per multiplication, while the
    cost ``sum k_l (log2 l + 1)`` stays within ``budget``."""
    exponent, spent = 1, 0.0
    current = pminus1_success(exponent, cases)
    while True:
        best = None
        for ell in SMALL_PRIMES:
            step = math.log2(ell) + 1
            if spent + step > budget:
                continue
            gain = (pminus1_success(exponent * ell, cases) - current) / step
            if best is None or gain > best[0]:
                best = (gain, ell, step)
        if best is None:
            return exponent, spent
        _, ell, step = best
        exponent *= ell
        spent += step
        current = pminus1_success(exponent, cases)


def output_slice(program):
    """The instructions the last register actually depends on."""
    need = {len(program) + 1}
    keep = []
    for t in range(len(program) - 1, -1, -1):
        reg = t + 2
        if reg in need:
            op, i, arg = program[t]
            keep.append((t, program[t]))
            need.add(i)
            if op not in ("pow", "cb2"):
                need.add(arg)
    return [ins for _, ins in reversed(keep)]


def output_exponent(program):
    """``E`` if the output is ``a^E`` (unused registers may be anything)."""
    exps = [1, 0]
    for op, i, arg in program:
        x = exps[i]
        if op == "mul":
            y = exps[arg]
            exps.append(None if x is None or y is None else x + y)
        elif op == "pow":
            exps.append(None if x is None else x * arg)
        else:
            exps.append(None)
    return exps[-1]


# ------------------------------------------------------------------ loop programs
LOOP_OPS = ("cb2", "mul", "add", "sub", "sq", "acc")


def loop_cost(body) -> float:
    return sum(1 for ins in body if ins[0] in ("cb2", "mul", "sq", "acc")) or 1


def run_loop(body, n: int, a: int, budget: float) -> int:
    """State ``s = [a, a, 1]`` updated in place by ``body`` repeated
    ``budget // cost(body)`` times; ``("acc", i, j)`` multiplies an accumulator
    by ``s_i - s_j``.  Output: the accumulator if the body has one, else ``s_0``."""
    s = [a % n, a % n, 1]
    acc, has_acc = 1, False
    half = (n + 1) // 2
    reps = int(budget // loop_cost(body))
    for _ in range(reps):
        for op, t, i, j in body:
            if op == "cb2":
                s[t] = s[i] * (s[i] - 1) * half % n
            elif op == "mul":
                s[t] = s[i] * s[j] % n
            elif op == "sq":
                s[t] = s[i] * s[i] % n
            elif op == "add":
                s[t] = (s[i] + s[j]) % n
            elif op == "sub":
                s[t] = (s[i] - s[j]) % n
            else:
                acc = acc * (s[i] - s[j]) % n
                has_acc = True
    return acc if has_acc else s[0]


def loop_splits(body, n: int, a: int, budget: float) -> bool:
    """Run ``body`` as in ``run_loop`` but test the output after every repetition,
    as a real rho implementation does -- a single final gcd fails as soon as
    both primes have collided, which would make success non-monotone in effort."""
    s = [a % n, a % n, 1]
    acc, has_acc = 1, False
    half = (n + 1) // 2
    reps = int(budget // loop_cost(body))
    for _ in range(reps):
        for op, t, i, j in body:
            if op == "cb2":
                s[t] = s[i] * (s[i] - 1) * half % n
            elif op == "mul":
                s[t] = s[i] * s[j] % n
            elif op == "sq":
                s[t] = s[i] * s[i] % n
            elif op == "add":
                s[t] = (s[i] + s[j]) % n
            elif op == "sub":
                s[t] = (s[i] - s[j]) % n
            else:
                acc = acc * (s[i] - s[j]) % n
                has_acc = True
        out = acc if has_acc else s[0]
        for g in (gcd(out, n), gcd(out - 1, n)):
            if 1 < g < n:
                return True
            if g == n and has_acc:
                return False                 # both primes collided in one step
    return False


def loop_score(body, cases, budget) -> float:
    return sum(loop_splits(body, n, a, budget) for n, a in cases) / len(cases)


def random_loop_instruction(rng):
    return (rng.choice(LOOP_OPS), rng.randrange(3), rng.randrange(3), rng.randrange(3))


def evolve_loops(cases, budget, generations, population, rng, max_body: int = 5):
    """Elitist search over loop bodies of up to ``max_body`` instructions."""
    def fit(b):
        return loop_score(b, cases, budget)

    pop = [[random_loop_instruction(rng) for _ in range(rng.randrange(1, max_body + 1))]
           for _ in range(population)]
    for _ in range(generations):
        scored = sorted(((fit(b), b) for b in pop), key=lambda t: -t[0])
        elite = [b for _, b in scored[: population // 5]]
        pop = list(elite)
        while len(pop) < population:
            parent = list(rng.choice(elite))
            r = rng.random()
            if r < 0.4:
                parent[rng.randrange(len(parent))] = random_loop_instruction(rng)
            elif r < 0.7 and len(parent) < max_body:
                parent.insert(rng.randrange(len(parent) + 1), random_loop_instruction(rng))
            elif len(parent) > 1:
                parent.pop(rng.randrange(len(parent)))
            pop.append(parent)
    return max(pop, key=fit)


RHO_BODY = [("cb2", 0, 0, 0), ("cb2", 1, 1, 0), ("cb2", 1, 1, 0), ("acc", 0, 1, 0)]
