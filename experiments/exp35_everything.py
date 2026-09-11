"""Round 27: every method in the repository, on the same numbers.

Thirty-four experiments built a zoo: the Pascal-native methods of the early
rounds, the reference methods of the late ones, and the elliptic triangle in
between.  This is the zoo in one table -- balanced semiprimes, three per size,
every method that finishes inside its budget -- with each method tagged by the
mechanism of round 14 it uses.
"""

import random
import signal
import statistics
import time

from _common import Report

from aksfactor.arith import is_prime
from aksfactor.central import central_split
from aksfactor.classgroup import schnorr_lenstra_split
from aksfactor.collision import pascal_rho
from aksfactor.cyclo import pollard_pminus1, williams_pplus1
from aksfactor.divseq import elliptic_triangle_factor
from aksfactor.ecm import ecm
from aksfactor.factor import pascal_split
from aksfactor.fast import fast_split
from aksfactor.harvey import harvey_factor
from aksfactor.hyperbola import hyperbola_factor
from aksfactor.nfs import number_field_sieve
from aksfactor.qs import quadratic_sieve

report = Report(
    "exp35_everything",
    "Every method, the same numbers",
    "Round 27. The whole repository in one table.",
)

rng = random.Random(35)
CAP = 15.0            # stop a method once its next size is projected past this
HARD = 60.0           # and never let a single call run longer than this
SIZES = (40, 50, 60, 70, 80, 90, 100)


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


signal.signal(signal.SIGALRM, _alarm)


def first(x):
    if x is None:
        return None
    if isinstance(x, tuple):
        return x[0]
    return x


def rho(n):
    for x0 in range(3, 300):
        g, _ = pascal_rho(n, 2, x0)
        if g:
            return g
    return None


NFS = {40: (400, 3000), 50: (700, 5000), 60: (1500, 8000), 70: (2500, 12000),
       80: (4000, 16000), 90: (6000, 20000), 100: (9000, 26000)}

METHODS = [
    ("Pascal row, residue scan", "size", lambda n, b: first(pascal_split(n, mode="scan"))),
    ("factorial threshold (Strassen)", "size", lambda n, b: first(fast_split(n))),
    ("Harvey N^(1/5)", "size", lambda n, b: first(harvey_factor(n))),
    ("hyperbola hull walk", "size", lambda n, b: hyperbola_factor(n)),
    ("central column", "sign", lambda n, b: first(central_split(n))),
    ("Pascal rho = Pollard rho", "collision", lambda n, b: rho(n)),
    ("Pollard p - 1, B = 10^5", "order, rigid", lambda n, b: pollard_pminus1(n, 100000)),
    ("Williams p + 1, B = 10^5", "order, rigid", lambda n, b: williams_pplus1(n, 100000)),
    ("class groups (Schnorr-Lenstra)", "order, redrawn", lambda n, b: first(schnorr_lenstra_split(n))),
    ("elliptic triangle, row lcm(1..B)", "order, redrawn",
     lambda n, b: elliptic_triangle_factor(n, 3000, 3000, random.Random(b))),
    ("ECM", "order, redrawn", lambda n, b: ecm(n, b1=3000, b2=150000, curves=3000)),
    ("quadratic sieve", "sign", lambda n, b: quadratic_sieve(n)),
    ("number field sieve", "sign", lambda n, b: number_field_sieve(n, 3, *(NFS[b][0], NFS[b][0]),
                                                                   24, NFS[b][1], 3000)),
]


def randprime(bits):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v):
            return v


numbers = {}
for bits in SIZES:
    numbers[bits] = []
    while len(numbers[bits]) < 3:
        p, q = randprime(bits // 2), randprime(bits // 2)
        if p != q:
            numbers[bits].append((p * q, p, q))

results = {}
skipped = set()
last = {}
for bits in SIZES:
    for name, mech, fn in METHODS:
        if name in skipped:
            results[name, bits] = "--"
            continue
        times, fails, timed_out = [], 0, False
        for n, p, q in numbers[bits]:
            t0 = time.perf_counter()
            signal.setitimer(signal.ITIMER_REAL, HARD)
            try:
                g = fn(n, bits)
            except Timeout:
                g, timed_out = None, True
            except Exception:
                g = None
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
            dt = time.perf_counter() - t0
            times.append(dt)
            if g not in (p, q):
                fails += 1
            if dt > CAP:
                break
        med = statistics.median(times)
        ok = len(times) - fails
        if timed_out:
            results[name, bits] = f"> {HARD:.0f} s"
        else:
            results[name, bits] = f"{med:.2f} s" + ("" if ok == 3 else f" ({ok}/{len(times)})")
        growth = med / last[name] if last.get(name) else 4.0
        last[name] = med
        if med * max(growth, 2.0) > CAP or max(times) > CAP:
            skipped.add(name)

rows = [[name, mech] + [results.get((name, b), "--") for b in SIZES] for name, mech, _ in METHODS]
report.p(f"Three balanced semiprimes per size (`p, q` of `bits/2` bits each). Median "
         f"time per number attempted; `(k/m)` where only `k` of the `m` attempted were "
         f"factored (a method stops attempting once one number takes over {CAP:.0f} s); "
         f"`> {HARD:.0f} s` where a single call hit the hard limit; `--` once a "
         f"method's next size is projected past {CAP:.0f} s.")
report.p()
report.table(["method", "mechanism"] + [f"{b} bits" for b in SIZES], rows)


def seconds(cell):
    try:
        return float(cell.split()[0]) if cell != "--" and "(" not in cell else None
    except ValueError:
        return None


winners = []
for b in SIZES:
    timed = [(seconds(results[name, b]), name) for name, _, _ in METHODS
             if seconds(results.get((name, b), "--")) is not None]
    if timed:
        winners.append((b, min(timed)[1]))
report.p("Fastest method that factored all three numbers, by size: " +
         "; ".join(f"{b} bits -- {w}" for b, w in winners) + ".")
report.p()


def reach(name):
    ok_sizes = [b for b in SIZES if seconds(results.get((name, b), "--")) is not None]
    return max(ok_sizes) if ok_sizes else None


lines = []
for name, mech, _ in METHODS:
    r = reach(name)
    lines.append(f"{name} ({mech}): " + (f"all three numbers up to {r} bits" if r else
                                          "never all three"))
report.p("Where each method stops -- the largest size at which it factored all "
         "three numbers within the budget:")
report.p()
for line in lines:
    report.p(f"- {line}")
report.p()
report.p("The size methods stop first, then collision; the methods that keep going "
         "are the ones powered by smooth numbers -- ECM, which redraws its groups, "
         "and the sieves. The rigid order methods succeed only on numbers whose "
         "`p +- 1` happens to be smooth -- common for small `p`, rare for large -- "
         "which is why their cells turn into fractions as `N` grows. Among the Pascal-native constructions, the one that goes "
         "furthest is Pascal rho, which is Pollard's rho, and the elliptic triangle, "
         "which is ECM once its sequence is changed.")
harvey = [results.get(("Harvey N^(1/5)", b), "--") for b in SIZES]
strassen = [results.get(("factorial threshold (Strassen)", b), "--") for b in SIZES]
report.p()
report.p(f"Harvey's `N^(1/5)` is the best *deterministic* exponent known; here it runs "
         f"{', '.join(f'{c} at {b} bits' for b, c in zip(SIZES, harvey) if c != '--')}, "
         f"against Strassen's "
         f"{', '.join(f'{c} at {b} bits' for b, c in zip(SIZES, strassen) if c != '--')}. "
         f"The jump is this implementation's parameter choice and pure-Python "
         f"constants, not the exponent.")
report.write()
