#!/usr/bin/env python3
"""Compile and run the independent certificate. Python standard library only."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
import runpy
import shutil
import subprocess
from fractions import Fraction
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING

HERE = Path(__file__).resolve().parent

def decimal_bound(x: Fraction, rounding: str) -> str:
    with localcontext() as ctx:
        ctx.prec = 30
        ctx.rounding = rounding
        return str(Decimal(x.numerator) / Decimal(x.denominator))

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')
    os.chdir(HERE)
    digest = hashlib.sha256((HERE / 'witness.hex').read_bytes()).hexdigest()
    expected = '4b82ce9794faeb80d1766dddc2884d6a19734490baf15a4200bb15bdfd5af718'
    if digest != expected:
        raise RuntimeError('Witness checksum mismatch.')
    print('Witness SHA-256:', digest, flush=True)
    print('\nExact support and admissibility checks:', flush=True)
    runpy.run_path(str(HERE / 'verify_support.py'), run_name='__main__')
    compiler = os.environ.get('CXX', 'g++')
    if shutil.which(compiler) is None:
        raise SystemExit(f'{compiler} was not found. A C++17 compiler with OpenMP is required.')
    executable = HERE / 'certify'
    command = [compiler, '-O3', '-std=c++17', '-fopenmp', '-frounding-math',
               '-ffp-contract=off', '-fno-fast-math', str(HERE / 'certify.cpp'),
               '-o', str(executable)]
    subprocess.run(command, check=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OMP_DYNAMIC='FALSE')
    totals = {name: Fraction(0) for name in ('Ilo', 'Ihi', 'Jlo', 'Jhi')}
    results = []
    print('\nDirected-rounding integral enclosures:', flush=True)
    for i in range(3):
        for j in range(i, 3):
            result = subprocess.run([str(executable), str(HERE / 'witness.hex'), str(i), str(j)],
                                    env=env, text=True, capture_output=True, check=True)
            item = json.loads(result.stdout)
            results.append(item)
            multiplier = 1 if i == j else 2
            for name in totals:
                # A hexadecimal binary64 representation is an exact dyadic rational.
                totals[name] += multiplier * Fraction.from_float(float.fromhex(item[name]))
            print(f'  pair ({i},{j}): enclosure computed', flush=True)
    if not (0 < totals['Ilo'] <= totals['Ihi'] and 0 < totals['Jlo'] <= totals['Jhi']):
        raise RuntimeError('The aggregate enclosures are not valid positive intervals.')
    low = 45 * totals['Jlo'] / totals['Ihi']
    high = 45 * totals['Jhi'] / totals['Ilo']
    if low <= 1:
        raise RuntimeError('The witness did NOT pass the sieve threshold.')
    environment = dict(compiler_version=subprocess.run(
        [compiler, '--version'], check=True, capture_output=True, text=True
    ).stdout.splitlines()[0], python_version=platform.python_version(),
        system=platform.system(), machine=platform.machine(),
        compiler_flags=command[1:-3])
    output = dict(threads=args.threads, witness_sha256=digest, environment=environment, pairs=results,
                  exact_totals={name: str(value) for name, value in totals.items()},
                  exact_ratio_lower=str(low), exact_ratio_upper=str(high),
                  ratio_lower=decimal_bound(low, ROUND_FLOOR),
                  ratio_upper=decimal_bound(high, ROUND_CEILING),
                  threshold_passed=True)
    (HERE / 'verified_output.json').write_text(json.dumps(output, indent=2) + '\n')
    print('\nCERTIFIED: 45 J_boxes / I is in')
    print('  [' + output['ratio_lower'] + ',')
    print('   ' + output['ratio_upper'] + ']')
    print('The lower endpoint is strictly greater than 1.')
    print('The exact support and admissibility checks also passed.')

if __name__ == '__main__':
    main()
