#!/usr/bin/env python3
"""99 Bottles of Beer -- traditional counting song, public domain.

Not part of mytools; a standalone warm-up script.
"""


def bottles(n):
    """Count as a noun phrase: 0 -> "no more bottles", 1 -> "1 bottle", n -> "n bottles"."""
    if n == 0:
        return "no more bottles"
    return f"{n} bottle" if n == 1 else f"{n} bottles"


def capitalised(s):
    """Same, but capitalised for the start of a line."""
    return s[0].upper() + s[1:]


START = 99
out = []

for n in range(START, 0, -1):
    out += [
        f"{capitalised(bottles(n))} of beer on the wall, {bottles(n)} of beer.",
        f"Take one down and pass it around, {bottles(n - 1)} of beer on the wall.",
        "",
    ]

out += [
    f"{capitalised(bottles(0))} of beer on the wall, {bottles(0)} of beer.",
    f"Go to the store and buy some more, {bottles(START)} of beer on the wall.",
]

print("\n".join(out))
