"""Round 42: the taxonomy against the literature, checked by running it.

The claim this project keeps leaning on is a claim of *coverage*: every
factoring method manufactures a zero divisor mod `N` through one of four
mechanisms -- order, size, collision, sign. A claim like that is only worth
what its counterexamples are worth, so this round lists every integer
factoring algorithm it can name, assigns each one, and -- for the ones this
repository implements -- runs it on an instance it should manage and checks
that it comes back with the right factor.

Two columns therefore mean different things. "Mechanism" is a judgement.
"Verified" is not: it is this script, calling that function, on that number.
"""

import random
import time

from _common import Report

from aksfactor.arith import is_prime, pollard_rho, spf_trial
from aksfactor.central import central_split
from aksfactor.cfrac import cfrac
from aksfactor.classgroup import schnorr_lenstra_split
from aksfactor.collision import pascal_rho
from aksfactor.cyclo import pollard_pminus1, williams_pplus1
from aksfactor.divseq import elliptic_triangle_factor
from aksfactor.ecm import ecm
from aksfactor.factor import pascal_split
from aksfactor.fast import fast_split, threshold_spf
from aksfactor.harvey import harvey_factor
from aksfactor.hyperbola import binary_search_factor, hyperbola_factor
from aksfactor.implicit import implicit_factor, shared_low_primes
from aksfactor.lattice import factor_with_hint
from aksfactor.nfs import number_field_sieve
from aksfactor.qs import quadratic_sieve
from aksfactor.shor import shor
from aksfactor.squfof import squfof

report = Report(
    "exp48_taxonomy",
    "Every algorithm named, assigned to a mechanism, and run",
    "Round 42. A coverage claim is worth its counterexamples: here is the list, "
    "and here is each implementation actually factoring something.",
)

rng = random.Random(48)


def randprime(bits, condition=None):
    while True:
        v = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if is_prime(v) and (condition is None or condition(v)):
            return v


def smooth(m, bound=200):
    for ell in range(2, bound + 1):
        while m % ell == 0:
            m //= ell
    return m == 1


BALANCED = None
while BALANCED is None:                      # 40 bits, two 20-bit primes
    a, b = randprime(20), randprime(20)
    if a != b:
        BALANCED = (a * b, min(a, b))
SMALL = None
while SMALL is None:                         # 24 bits, for the slow constants
    a, b = randprime(12), randprime(12)
    if a != b:
        SMALL = (a * b, min(a, b))
PM1 = (lambda p, q: (p * q, min(p, q)))(randprime(20, lambda v: smooth(v - 1)),
                                        randprime(20))
PP1 = (lambda p, q: (p * q, min(p, q)))(randprime(20, lambda v: smooth(v + 1)),
                                        randprime(20))
TINY = (3 * 7, 3)                            # Shor's simulator is a state vector


def _window(n):
    """`N^(1/4)`, rounded up to a power of two: what Coppersmith can close."""
    return 1 << (int(n ** 0.25).bit_length())


def _coppersmith(n, p):
    """Hand the lattice `p` truncated to a window, widest window first.

    The reachable window is `N^(1/4)` asymptotically and a little less at 40
    bits, which round 38 measured from the other side; trying four widths
    reports the one the lattice actually closed rather than assuming it.
    """
    for target in (n // p, p):          # the lattice's beta wants the larger one
        w = _window(n)
        while w > 1:
            got = factor_with_hint(n, target - target % w, w)
            if got:
                return got
            w //= 2
    return None


def attempt(fn, case):
    """Run one method; report whether it returned a true factor, and how long.

    The function is handed the true factor as well, because two rows need it to
    build their own input -- Coppersmith is given a hint by definition, and the
    hint has to come from somewhere.
    """
    n, want = case
    t0 = time.perf_counter()
    try:
        got = fn(n, want)
    except Exception as exc:                 # a method that throws has not verified
        return f"error: {type(exc).__name__}", time.perf_counter() - t0
    dt = time.perf_counter() - t0
    if isinstance(got, dict):
        got = tuple(got)
    if isinstance(got, tuple):
        # the return conventions differ: (p, q) here, (factor, multiplier)
        # there, so ask which of the returned values actually divides n
        real = [g for g in got if isinstance(g, int) and 1 < g < n and n % g == 0]
        got = min(real) if real else (got[0] if got else None)
    if got is None:
        return "no split", dt
    ok = isinstance(got, int) and 1 < got < n and n % got == 0
    return ("yes" if ok else f"wrong ({got})"), dt


ROWS = [
    ("trial division", "antiquity", "size", "`arith.spf_trial`",
     lambda n, _: spf_trial(n), SMALL),
    ("Pascal row residue scan (T1-T5)", "this project", "size",
     "`factor.pascal_split`", lambda n, _: pascal_split(n), SMALL),
    ("Pollard-Strassen product tree", "1974/76", "size", "`fast.fast_split`",
     lambda n, _: fast_split(n), BALANCED),
    ("factorial threshold search (T20)", "this project", "size",
     "`fast.threshold_spf`", lambda n, _: threshold_spf(n), SMALL),
    ("Fermat / Lehman", "1643/1974", "size", "in `harvey.lehman_recover`", None, None),
    ("Harvey's deterministic `N^(1/5)`", "2021", "size", "`harvey.harvey_factor`",
     lambda n, _: harvey_factor(n), BALANCED),
    ("hyperbola hull walk (P23)", "this project", "size",
     "`hyperbola.hyperbola_factor`", lambda n, _: hyperbola_factor(n), BALANCED),
    ("binary search by divisor count", "this project", "size",
     "`hyperbola.binary_search_factor`", lambda n, _: binary_search_factor(n), BALANCED),
    ("Coppersmith, given `p` to within about `N^(1/4)`", "1996", "size",
     "`lattice.factor_with_hint`", _coppersmith, BALANCED),
    ("Pollard rho", "1975", "collision", "`arith.pollard_rho`",
     lambda n, _: pollard_rho(n), BALANCED),
    ("Pascal rho `x -> C(x,k)` (P22)", "this project", "collision",
     "`collision.pascal_rho`", lambda n, _: pascal_rho(n)[0], BALANCED),
    ("Pollard `p - 1`", "1974", "order", "`cyclo.pollard_pminus1`",
     lambda n, _: pollard_pminus1(n, 1000), PM1),
    ("Williams `p + 1`", "1982", "order", "`cyclo.williams_pplus1`",
     lambda n, _: williams_pplus1(n, 1000), PP1),
    ("Schnorr-Lenstra class groups", "1984", "order",
     "`classgroup.schnorr_lenstra_split`, which returns `(factor, multiplier)`",
     lambda n, _: schnorr_lenstra_split(n), SMALL),
    ("Lenstra's ECM", "1987", "order", "`ecm.ecm`",
     lambda n, _: ecm(n, b1=2000, b2=100000, curves=200, rng=random.Random(1)), BALANCED),
    ("elliptic triangle, row `lcm(1..B)` (P24)", "this project", "order",
     "`divseq.elliptic_triangle_factor`",
     lambda n, _: elliptic_triangle_factor(n, 500, 40, random.Random(2)), BALANCED),
    ("Shor", "1994", "order (read, not guessed)", "`shor.shor`",
     lambda n, _: shor(n, random.Random(3)), TINY),
    ("central column Legendre symbols", "this project", "sign",
     "`central.central_split`", lambda n, _: central_split(n), SMALL),
    ("Shanks' SQUFOF", "1975", "sign", "`squfof.squfof`", lambda n, _: squfof(n), BALANCED),
    ("CFRAC", "1931/1975", "sign", "`cfrac.cfrac`", lambda n, _: cfrac(n), BALANCED),
    ("Dixon's random squares", "1981", "sign", "unsieved CFRAC/QS", None, None),
    ("quadratic sieve", "1981", "sign", "`qs.quadratic_sieve`",
     lambda n, _: quadratic_sieve(n), BALANCED),
    ("MPQS / SIQS", "1987", "sign", "engineering on top of the QS", None, None),
    ("number field sieve", "1993", "sign", "`nfs.number_field_sieve`",
     lambda n, _: number_field_sieve(n, 3, 1500, 1500, 24, 8000, 3000), BALANCED),
    ("Schnorr's prime-number lattice", "1991", "sign",
     "`schnorr.split_from_relations`", None, None),
    ("Regev's quantum variant", "2023", "order (read)", "not simulable here",
     None, None),
]

report.p("## The list")
report.p()
rows = []
for name, year, mechanism, where, fn, case in ROWS:
    if fn is None or case is None:
        rows.append([name, year, mechanism, where, "--", "--"])
        continue
    verdict, dt = attempt(fn, case)
    rows.append([name, year, mechanism, where, verdict, f"{dt:.2f} s"])
report.table(["algorithm", "year", "mechanism", "where it lives here",
              "factored its test number", "time"], rows)

ran = [r for r in rows if r[4] not in ("--",)]
good = [r for r in ran if r[4] == "yes"]
report.p(f"{len(good)} of {len(ran)} implementations returned a value that is a "
         f"nontrivial divisor of the number they were given. The instances differ by row on purpose: "
         f"`p - 1` and `p + 1` are handed a prime whose shifted value is smooth, "
         f"Shor is handed 21 because its simulator is a state vector, and the "
         f"methods with large constants are handed 24 bits rather than 40. This "
         f"table is a coverage check, not a race -- round 27 is the race.")
report.p()

report.p("## The rows with no implementation, and why")
report.p()
report.table(["algorithm", "why it is not here"], [
    ["Fermat / Lehman", "Lehman's certificates are produced by the hull walk "
     "(P23) and consumed by `harvey.lehman_recover`; Fermat is the `k = 1` case"],
    ["Dixon's random squares", "the quadratic sieve with the sieve removed -- "
     "same mechanism, strictly more work"],
    ["MPQS / SIQS", "polynomial switching and self-initialisation: better "
     "constants for the same `L[1/2, 1]`"],
    ["Schnorr's prime-number lattice", "the lattice and its relations are "
     "implemented (round 20); it never produced a usable relation, which is "
     "round 20's result"],
    ["Regev's variant", "quantum, and its circuit is not a state vector this "
     "machine can hold"],
])
report.p()

report.p("## What would falsify the taxonomy")
report.p()
report.p("A method that produces a nontrivial `gcd(x, N)` without `x` being a "
         "difference of two things congruent mod `p` (collision, sign), a "
         "candidate divisor or a batch of them (size), or a group element "
         "killed by an exponent (order). Nothing in the list above is such a "
         "method, and the two additions this project made -- the hyperbola hull "
         "walk and the elliptic triangle -- landed in size and order "
         "respectively rather than anywhere new.")
report.p()
report.p("The nearest thing to an escape is not a mechanism at all but a change "
         "of question. Implicit factoring (round 37) does not attack one `N`; "
         "it attacks `k` of them that share bits, and it is the one entry here "
         "whose input is not a single integer:")
report.p()

pairs = shared_low_primes(3, 200, 100, 40, rng)
moduli = [p * q for p, q in pairs]
t0 = time.perf_counter()
got = implicit_factor(moduli, 100)
dt = time.perf_counter() - t0
ok = got is not None and sorted(got) == sorted(p for p, _ in pairs)
report.table(["algorithm", "year", "input", "mechanism", "factored", "time"],
             [["implicit factoring (May-Ritzenhofen)", "2009",
               "3 moduli sharing 100 low bits of `p`", "lattice, between instances",
               "yes" if ok else "no", f"{dt:.2f} s"]])
report.p("Which is why round 37 is filed as a different currency rather than a "
         "fifth mechanism: it buys its zero divisor from a correlation between "
         "instances, and there is no correlation inside one instance to buy.")
report.write()
