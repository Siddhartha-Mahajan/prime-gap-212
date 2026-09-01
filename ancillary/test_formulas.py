#!/usr/bin/env python3
"""Compare both interval implementations with exhaustive exact integration.

Run the two certificate scripts first to compile the C++ verifiers. These
six-variable tests use signed radial coefficients and only Python's standard
library. The default checks both implementations; --arithmetic selects one.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import json
import os
import subprocess
import tempfile
from fixed_inputs import BITS, generate

HERE = Path(__file__).resolve().parent


def test_case(seed: int, arithmetic: list[str]) -> dict:
    k, R, n, D, den, T, L, maxr = 6, 3, 19, 2, 64, 24, 18, 4
    B = [0, 5, 7, 9, 9, 9]
    g = [[Q(3 + ((a + 2*i + seed) % 5), 8) for a in range(n)] for i in range(R)]
    h = [[Q(((s*(i + 1) + 3*i + seed) % 11) - 5, 16) for s in range(n)] for i in range(R)]

    def value(indices: tuple[int, ...]) -> Q:
        s = sum(indices)
        r = sum(a >= D for a in indices)
        if s >= n or (r and (r > maxr or sum(a for a in indices if a >= D) + r > B[r])):
            return Q(0)
        out = Q(0)
        for i in range(R):
            v = h[i][s]
            for a in indices:
                v *= g[i][a]
            out += v
        return out

    exact_I = sum((value(a)**2 for a in product(range(B[1]), repeat=k)), Q(0)) / den**k
    exact_J = Q(0)
    for a in product(range(B[1]), repeat=k-1):
        if sum(a) > L - (k-1):
            continue
        fiber = sum((value(a + (u,)) for u in range(B[1])), Q(0))
        exact_J += fiber**2
    exact_J /= den**(k+1)
    intervals = {}
    with tempfile.TemporaryDirectory(prefix='gap-formula-test-') as tmp:
        path = Path(tmp) / 'toy.hex'
        with path.open('w') as f:
            f.write(' '.join(map(str, (k, R, n, D, den, T, L, maxr))) + '\n')
            f.write(' '.join(map(str, B)) + '\n')
            for a in range(n):
                # These toy coefficients are exactly representable in binary64.
                values = [*(row[a] for row in g), *(row[a] for row in h)]
                f.write(' '.join(float(v).hex() for v in values) + '\n')
        fixed = Path(tmp) / 'toy.fix'
        generate(path, fixed)
        env = dict(os.environ, OMP_NUM_THREADS='2', OMP_DYNAMIC='FALSE')
        for method in arithmetic:
            total = {name: Q(0) for name in ('Ilo', 'Ihi', 'Jlo', 'Jhi')}
            executable = HERE / ('certify_integer' if method == 'integer' else 'certify')
            source = fixed if method == 'integer' else path
            for i in range(R):
                for j in range(i, R):
                    run = subprocess.run([str(executable), str(source), str(i), str(j)],
                                         check=True, capture_output=True, text=True, env=env)
                    out = json.loads(run.stdout)
                    for name in total:
                        endpoint = (Q(int(out[name]), 1 << BITS) if method == 'integer'
                                    else Q.from_float(float.fromhex(out[name])))
                        total[name] += (1 if i == j else 2) * endpoint
            if not (total['Ilo'] <= exact_I <= total['Ihi'] and total['Jlo'] <= exact_J <= total['Jhi']):
                raise ArithmeticError(f'Exhaustive exact integral not enclosed by {method}')
            intervals[method] = {key: str(value) for key, value in total.items()}
    return dict(seed=seed, exact_I=str(exact_I), exact_J=str(exact_J), intervals=intervals, passed=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arithmetic', choices=('integer', 'binary64', 'both'), default='both')
    args = parser.parse_args()
    methods = ['integer', 'binary64'] if args.arithmetic == 'both' else [args.arithmetic]
    for method in methods:
        executable = HERE / ('certify_integer' if method == 'integer' else 'certify')
        if not executable.exists():
            raise SystemExit('Compile the verifiers first by running the certificate scripts.')
    results = [test_case(seed, methods) for seed in (0, 1, 7)]
    (HERE / 'test_output.json').write_text(json.dumps(results, indent=2) + '\n')
    print('Three exhaustive exact-integration tests passed for: ' + ', '.join(methods) + '.')


if __name__ == '__main__':
    main()
