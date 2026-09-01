"""Convert hexadecimal dyadic coefficients to exact integers at scale 2**96.

Only Python's standard library is used; no floating-point conversion occurs.
"""
from __future__ import annotations
from fractions import Fraction
from pathlib import Path
import re

BITS = 96
SCALE = 1 << BITS
HEX = re.compile(r'([+-]?)0x([0-9a-fA-F]+)(?:\.([0-9a-fA-F]*))?p([+-]?\d+)\Z')


def dyadic(token: str) -> Fraction:
    match = HEX.fullmatch(token)
    if match is None:
        raise ValueError('Invalid hexadecimal dyadic: ' + token)
    sign, whole, fractional, exponent = match.groups()
    fractional = fractional or ''
    mantissa = int(whole + fractional, 16)
    if sign == '-':
        mantissa = -mantissa
    exponent = int(exponent) - 4 * len(fractional)
    return (Fraction(mantissa * (1 << exponent)) if exponent >= 0
            else Fraction(mantissa, 1 << (-exponent)))


def generate(source: Path, destination: Path) -> None:
    """Create an integer witness, rejecting coefficients not exactly scaled."""
    lines = source.read_text(encoding='ascii').splitlines()
    if len(lines) < 2:
        raise ValueError('Missing witness header')
    header = [int(x) for x in lines[0].split()]
    if len(header) != 8:
        raise ValueError('Incorrect header length')
    _, components, length, _, _, _, _, maximum_large = header
    if len(lines) != length + 2:
        raise ValueError('Incorrect number of coefficient rows')
    if len(lines[1].split()) != maximum_large + 2:
        raise ValueError('Incorrect number of cutoffs')
    with destination.open('w', encoding='ascii') as stream:
        stream.write(str(BITS) + '\n' + lines[0] + '\n' + lines[1] + '\n')
        for line in lines[2:]:
            tokens = line.split()
            if len(tokens) != 2 * components:
                raise ValueError('Incorrect coefficient row length')
            values = []
            for token in tokens:
                value = dyadic(token) * SCALE
                if value.denominator != 1:
                    raise ArithmeticError('Coefficient is not exactly representable at scale 2**96')
                values.append(str(value.numerator))
            stream.write(' '.join(values) + '\n')
