#!/usr/bin/env python3
"""Check the article's support inequalities and admissible tuple exactly.

Every mathematical comparison uses fractions.Fraction.  The decimal bounds
reported with the certificate are parsed as terminating rational numbers, not binary floats.
The parameter record and the witness header are checked against the values
in the article before any inequalities are evaluated.
"""
from __future__ import annotations

import json
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent


def verify() -> dict:
    den = 65536
    A, eps, delta = (Q(v, den) for v in (16862, 557, 1019))
    b1, b2, b3, B = (Q(v, den) for v in (10171, 10378, 11321, 12082))
    xi1, xi2, xi3 = Q(19, 50), Q(2, 5), Q(2, 5)
    eta, tau, zeta = Q(1, 10**10), Q(1, 10**8), Q(1, 10**6)
    w = A - Q(1, 4)
    ds = delta + zeta
    G = Q(1, 3) + 8*w + Q(7, 3)*delta
    gm, gp = xi2 - eta, G + 3*eta
    alpha = gm - 2*delta - 8*w - 109*tau
    S = Q(1, 2) - 2*delta - 10*w - 116*tau
    C = max(Q(1, 2) - gp - 2*w - 7*tau, Q(4, 5)*(Q(1, 2) - gp))
    K = Q(2, 5) + Q(24, 5)*w + Q(7, 5)*delta
    Delta_Ih = (1 - 68*w)/14 - eta
    Delta_III = Q(1, 2) - Q(7, 2)*w - Q(9, 8)*(xi3 + eta) - eta
    U_III = Q(1, 3) + Q(4, 3)*Delta_III - Q(4, 3)*w
    L_III = U_III - Delta_III

    header = [45, 3, 17375, 1019, den, 17419, 16305, 11]
    caps = [0, 10171, 10378, 11321] + [12082]*9
    record = json.loads((HERE / 'parameters.json').read_text())
    if record['sha256'] != '4b82ce9794faeb80d1766dddc2884d6a19734490baf15a4200bb15bdfd5af718':
        raise ArithmeticError('The parameter record has an incorrect witness digest.')
    fields = ['k', 'R', 'n', 'D', 'den', 'T', 'L', 'maxr']
    if [record[key] for key in fields] != header or record['B'] != caps:
        raise ArithmeticError('The parameter record does not match the article.')
    if record['A_num'] != 16862 or record['eps_num'] != 557:
        raise ArithmeticError('The parameter record has an incorrect simplex cutoff.')
    with (HERE / 'witness.hex').open() as stream:
        if list(map(int, stream.readline().split())) != header:
            raise ArithmeticError('The witness header does not match the article.')
        if list(map(int, stream.readline().split())) != caps:
            raise ArithmeticError('The witness caps do not match the article.')

    # Each quantity must be positive, apart from the explicitly weak comparisons.
    checks = {
        'A > 1/4': w,
        'epsilon > 0': eps,
        'A > epsilon': A - eps,
        'A+epsilon < 1/2': Q(1, 2) - A - eps,
        'delta > 0': delta,
        'b1 > delta': b1 - delta,
        'b1 < 1': 1 - b1,
        'b2 > b1': b2 - b1,
        'b3 > b2': b3 - b2,
        'B > b3': B - b3,
        'b2 <= b1+delta': b1 + delta - b2,
        'b3 <= b2+delta': b2 + delta - b3,
        'B <= b3+delta': b3 + delta - B,
        '2xi1+3xi2 < 2': 2 - 2*xi1 - 3*xi2,
        'xi2 <= xi3': xi3 - xi2,
        'xi1+9xi2 < 4': 4 - xi1 - 9*xi2,
        '2xi1+xi2 > 1': 2*xi1 + xi2 - 1,
        '17xi2 < 7': 7 - 17*xi2,
        'roughness: b1 < 1-2xi2': 1 - 2*xi2 - b1,
        'partition: alpha - 2b1': alpha - 2*b1,
        'partition: S - (5/2)b2': S - Q(5, 2)*b2,
        'partition: 3C - b3': 3*C - b3,
        'partition: S - 2B': S - 2*B,
        'partition: gamma_min-2delta-4tau - 2B': gm - 2*delta - 4*tau - 2*B,
        'partition: S+alpha+2delta - 4B': S + alpha + 2*delta - 4*B,
        'IIc scalar 1': 1 - 8*w - 4*ds - 2*gp,
        'IIc scalar 2': gm - 32*w - 10*ds,
        'IIc scalar 3': 4*gm - 1 - 48*w - 16*ds,
        'factor-window -52 lower slack': zeta - 49*tau,
        'factor-window -52 upper slack': 49*tau,
        'factor-window -100 lower slack': zeta - tau,
        'factor-window -100 upper slack': tau,
        'gamma_max < 1/2': Q(1, 2) - gp,
        'gamma_max > gamma_min': gp - gm,
        'I low: Delta - delta': xi1 - 4*A + Q(2, 3) - 2*eta - delta,
        'I low: upper - 2B': xi1 - eta - 3*tau - 2*B,
        'I low: total - lower': Q(1, 6) - 4*w - eta,
        'I high: upper - 2B': Q(1, 2) - 2*w - 4*tau - 2*B,
        'I high: Delta - delta': Delta_Ih - delta,
        'IIa: minimum width surplus': Q(3, 7)*eta,
        'IIa: upper - 2B': K + 2*eta - 3*tau - 2*B,
        'IIa: total - lower': Q(1, 14) - Q(24, 7)*w - eta,
        'IIb: width surplus lower bound': Q(2, 7)*eta,
        'IIb: r upper - 2B': G + 3*eta - 3*tau - 2*B,
        'IIb: u upper positive': Q(1, 2) - K - 2*eta - 2*w - 6*tau,
        'IIb: smooth mass - u upper': G + 2*w - 2*B - 2*tau,
        'III: Delta - delta': Delta_III - delta,
        'III: upper - (2B-delta)': U_III - 2*B + delta,
        'III: available mass - lower': Q(1, 2) - 2*tau - b1 - L_III,
        'Type II cutoffs ordered': K + 2*eta - gp,
        'Type II upper cutoff < 1/2': Q(1, 2) - K - 2*eta,
    }
    weak = {'xi2 <= xi3', 'b2 <= b1+delta', 'b3 <= b2+delta', 'B <= b3+delta'}
    for name, value in checks.items():
        if value < 0 or (value == 0 and name not in weak):
            raise ArithmeticError(f'Failed inequality {name}: {value}')

    # Strict decimal lower bounds for the article and supporting certificate.
    printed_bounds = {
        'partition: alpha - 2b1': '0.0001576',
        'partition: S - (5/2)b2': '0.0000751',
        'partition: 3C - b3': '0.0001434',
        'partition: S - 2B': '0.02725',
        'partition: gamma_min-2delta-4tau - 2B': '0.0001891',
        'partition: S+alpha+2delta - 4B': '0.0001869',
        'IIc scalar 1': '0.02352',
        'IIc scalar 2': '0.01110',
        'IIc scalar 3': '0.00110',
        'I low: Delta - delta': '0.001943',
        'I low: upper - 2B': '0.01128',
        'I low: total - lower': '0.13749',
        'I high: upper - 2B': '0.11669',
        'I high: Delta - delta': '0.02045',
        'IIa: upper - 2B': '0.08806',
        'IIa: total - lower': '0.04642',
        'IIb: r upper - 2B': '0.05924',
        'IIb: u upper positive': '0.02863',
        'IIb: smooth mass - u upper': '0.07383',
        'III: Delta - delta': '0.00892',
        'III: upper - (2B-delta)': '0.00307',
        'III: available mass - lower': '0.01303',
    }
    for name, bound in printed_bounds.items():
        if checks[name] <= Q(bound):
            raise ArithmeticError(f'Incorrect printed lower bound for {name}.')

    data = json.loads((HERE / 'tuple.json').read_text())
    H = data['tuple']
    if (data['k'] != 45 or len(H) != 45 or len(set(H)) != 45
            or H != sorted(H) or H[-1] - H[0] != 212):
        raise ArithmeticError('Incorrect tuple dimensions or diameter.')
    primes = [p for p in range(2, 46) if all(p % d for d in range(2, isqrt(p) + 1))]
    missing = {p: sorted(set(range(p)) - {h % p for h in H}) for p in primes}
    if not all(missing.values()):
        raise ArithmeticError('The tuple is not admissible.')
    if data['primes'] != primes or len(data['excluded']) != len(primes):
        raise ArithmeticError('Incorrect list of small primes.')
    for p, residue in zip(primes, data['excluded']):
        if residue not in missing[p]:
            raise ArithmeticError(f'The displayed residue for {p} is not omitted.')

    values = dict(A=A, epsilon=eps, delta=delta, b1=b1, b2=b2, b3=b3, B=B,
                  xi1=xi1, xi2=xi2, xi3=xi3, eta=eta, tau_max=tau,
                  zeta=zeta, delta_star=ds, omega=w, G=G, gamma_min=gm,
                  gamma_max=gp, alpha=alpha, S=S, C=C, K=K,
                  Delta_Ih=Delta_Ih, Delta_III=Delta_III,
                  U_III=U_III, L_III=L_III)
    result = dict(parameters={n: str(v) for n, v in values.items()},
                  margins={n: str(v) for n, v in checks.items()},
                  printed_strict_lower_bounds=printed_bounds,
                  missing_residues=missing, cardinality=45, diameter=212,
                  all_checks_passed=True)
    (HERE / 'support_output.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    result = verify()
    for name, margin in result['margins'].items():
        print(f'{name}: {margin}')
    print('All exact parameter inequalities and printed bounds passed.')
    print('Admissible tuple: 45 entries, diameter 212.')
    print('Missing residues:', result['missing_residues'])
