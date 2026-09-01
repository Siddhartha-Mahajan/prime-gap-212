#!/usr/bin/env python3
"""Verify the fixed box witness with exact integer interval arithmetic."""
from __future__ import annotations
import argparse
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import runpy
import shutil
import subprocess
import tempfile
from fixed_inputs import BITS, generate

HERE = Path(__file__).resolve().parent
DIGEST = '4b82ce9794faeb80d1766dddc2884d6a19734490baf15a4200bb15bdfd5af718'


def decimal_bound(value: Fraction, rounding: str) -> str:
    with localcontext() as context:
        context.prec = 35
        context.rounding = rounding
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    os.chdir(HERE)
    digest = hashlib.sha256((HERE / 'witness.hex').read_bytes()).hexdigest()
    if digest != DIGEST:
        raise RuntimeError('Witness checksum mismatch')
    print('Witness SHA-256:', digest, flush=True)
    print('\nExact support and admissibility checks:', flush=True)
    runpy.run_path(str(HERE / 'verify_support.py'), run_name='__main__')
    compiler = os.environ.get('CXX', 'g++')
    if shutil.which(compiler) is None:
        raise SystemExit(f'{compiler} was not found. C++17, OpenMP, and Boost headers are required.')
    executable = HERE / 'certify_integer'
    flags = ['-O3', '-std=c++17', '-fopenmp']
    subprocess.run([compiler, *flags, str(HERE / 'certify_integer.cpp'), '-o', str(executable)], check=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC='FALSE')
    results = []
    totals = {key: Fraction(0) for key in ('Ilo', 'Ihi', 'Jlo', 'Jhi')}
    with tempfile.TemporaryDirectory(prefix='prime-gap-integer-') as tmp:
        integer_input = Path(tmp) / 'witness.fix'
        generate(HERE / 'witness.hex', integer_input)
        print('\nExact integer-interval calculation:', flush=True)
        for i in range(3):
            for j in range(i, 3):
                process = subprocess.run([str(executable), str(integer_input), str(i), str(j)],
                                         env=env, text=True, capture_output=True, check=True)
                item = json.loads(process.stdout)
                if item['bits'] != BITS or (item['i'], item['j']) != (i, j):
                    raise ArithmeticError('Incorrect component or integer scale')
                if int(item['Ilo']) > int(item['Ihi']) or int(item['Jlo']) > int(item['Jhi']):
                    raise ArithmeticError('Invalid component enclosure')
                results.append(item)
                multiplier = 1 if i == j else 2
                for key in totals:
                    totals[key] += multiplier * Fraction(int(item[key]), 1 << BITS)
                print(f'  pair ({i},{j}): integer enclosure computed', flush=True)
    if not (0 < totals['Ilo'] <= totals['Ihi'] and 0 < totals['Jlo'] <= totals['Jhi']):
        raise ArithmeticError('Aggregate enclosures are not valid positive intervals')
    lower = 45 * totals['Jlo'] / totals['Ihi']
    upper = 45 * totals['Jhi'] / totals['Ilo']
    # Verify all integral and quotient bounds printed in Proposition 7.1.
    printed_totals = {
        'Ilo': Fraction(2263523516527043, 2251799813685248),
        'Ihi': Fraction(1131761758374273, 1125899906842624),
        'Jlo': Fraction(402482637490103, 18014398509481984),
        'Jhi': Fraction(1609930550374301, 72057594037927936),
    }
    for integral in ('I', 'J'):
        lo, hi = integral + 'lo', integral + 'hi'
        if not (printed_totals[lo] <= totals[lo] <= totals[hi] <= printed_totals[hi]):
            raise ArithmeticError('A published integral enclosure was not certified')
    printed_lower = Fraction(6037239562351545, 6036062711329456)
    printed_upper = Fraction(72446874766843545, 72432752528865376)
    if not (1 < printed_lower <= lower <= upper <= printed_upper):
        raise ArithmeticError('The published ratio enclosure was not certified')
    output = dict(
        threads=args.threads, bits=BITS, witness_sha256=digest,
        environment=dict(
            compiler_version=subprocess.run([compiler, '--version'], check=True, capture_output=True,
                                            text=True).stdout.splitlines()[0],
            python_version=platform.python_version(), system=platform.system(), machine=platform.machine(),
            compiler_flags=flags),
        pairs=results, exact_totals={key: str(value) for key, value in totals.items()},
        exact_ratio_lower=str(lower), exact_ratio_upper=str(upper),
        ratio_lower=decimal_bound(lower, ROUND_FLOOR), ratio_upper=decimal_bound(upper, ROUND_CEILING),
        published_enclosure_verified=True, threshold_passed=True)
    (HERE / 'verified_integer_output.json').write_text(json.dumps(output, indent=2) + '\n')
    print('\nCERTIFIED: 45 J_boxes / I lies in')
    print('  [' + output['ratio_lower'] + ',')
    print('   ' + output['ratio_upper'] + ']')
    print('The published enclosure and the strict threshold are verified with integer arithmetic.')


if __name__ == '__main__':
    main()
