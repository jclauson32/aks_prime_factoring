"""Pictures of triangles mod N, as PNG files, with the standard library only.

Each entry of a triangle mod ``N = p q`` is classified by which primes divide it:

* divisible by neither -- light;
* divisible by ``p`` only, or by ``q`` only -- coloured: these are exactly the
  entries whose gcd with ``N`` is a proper factor, the "x * y" entries of the
  observation this project started from;
* divisible by both -- dark.

For Pascal's triangle the two gaskets have cells ``p`` and ``q``.  For a
triangle over an elliptic divisibility sequence they have cells ``r_p`` and
``r_q``, the orders of the point mod ``p`` and mod ``q``.

Entries are computed exactly modulo each prime by tracking ``p``-adic
valuations: ``[n, k] = [n, k-1] a_(n-k+1) / a_k``, with each ``a_i`` written as
``p^v`` times a unit.
"""

from __future__ import annotations

import struct
import zlib


def write_png(path, pixels):
    """``pixels``: list of rows, each a list of ``(r, g, b)`` tuples."""
    height = len(pixels)
    width = len(pixels[0])
    raw = b"".join(b"\x00" + bytes(c for px in row for c in px) for row in pixels)

    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data +
                struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))

    png = (b"\x89PNG\r\n\x1a\n" +
           chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)) +
           chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))
    with open(path, "wb") as fh:
        fh.write(png)


def _units(seq_mod, p, e):
    """For ``a_i`` known mod ``p^e``: valuations and unit parts mod ``p``."""
    val, unit = [0], [1]
    pe = p ** e
    for i in range(1, len(seq_mod)):
        x = seq_mod[i] % pe
        if x == 0:
            raise ValueError("precision too low")
        v = 0
        while x % p == 0:
            x //= p
            v += 1
        val.append(v)
        unit.append(x % p)
    return val, unit


def zero_pattern(seq_mod, p, e, rows):
    """``zero[n][k]``: whether ``p`` divides the generalised binomial ``[n, k]``."""
    val, unit = _units(seq_mod, p, e)
    inv = [0] + [pow(u, -1, p) for u in unit[1:]]
    out = []
    for n in range(rows):
        v, u = 0, 1
        row = [False]
        for k in range(1, n + 1):
            v += val[n - k + 1] - val[k]
            u = u * unit[n - k + 1] * inv[k] % p
            row.append(v > 0)
        out.append(row)
    return out


LIGHT, DARK = (246, 244, 236), (40, 40, 48)
ONLY_P, ONLY_Q = (214, 72, 64), (52, 110, 200)


def classify(zero_p, zero_q, rows, scale=2):
    """Pixels for rows ``0..rows-1``, each entry a ``scale x scale`` block,
    the triangle drawn left-aligned (row ``n`` has ``n + 1`` entries)."""
    width = rows * scale
    pixels = []
    for n in range(rows):
        line = []
        for k in range(rows):
            if k > n:
                c = (255, 255, 255)
            else:
                zp, zq = zero_p[n][k], zero_q[n][k]
                c = DARK if zp and zq else ONLY_P if zp else ONLY_Q if zq else LIGHT
            line.extend([c] * scale)
        for _ in range(scale):
            pixels.append(list(line))
    return pixels
