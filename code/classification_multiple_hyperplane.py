#!/usr/bin/env python3
"""Exact small identities for the full resonant modular eigenspace.

The modular reduction and all-solution argument are in
MULTIPLE_HYPERPLANE_MODULAR_SUPPORT.md.  No sampled points are used here.
"""
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s


def matrix(entries):
    m = s.zeros(4)
    for (i, j), f in zip(itertools.combinations(range(4), 2), entries):
        m[i, j], m[j, i] = f, -f
    return m


def jacobi(m, xx):
    return [s.expand(sum(m[i, l] * s.diff(m[j, k], xx[l])
                        + m[j, l] * s.diff(m[k, i], xx[l])
                        + m[k, l] * s.diff(m[i, j], xx[l])
                        for l in range(4)))
            for i, j, k in itertools.combinations(range(4), 3)]


def pf(m):
    return s.expand(m[0, 1]*m[2, 3]-m[0, 2]*m[1, 3]+m[0, 3]*m[1, 2])


def run():
    started = time.monotonic()
    xx = s.symbols('x1:5')
    x1, x2, x3, x4 = xx
    a, b, y, c, d, p, q, t = s.symbols('a b y c d p q t')
    kk = (b, y, -a-b-y, -2*b-y, a+3*b+y, -2*(a+b))

    def tensor(k, cc=c, dd=d, pp=p, qq=q):
        return matrix([cc+k[0]*x1*x2, k[1]*x1*x3,
                       pp*x2+k[2]*x1*x4, qq*x1+k[3]*x2*x3,
                       k[4]*x2*x4, dd+k[5]*x3*x4])

    m = tensor(kk)
    mod = [s.expand(sum(s.diff(m[i, j], xx[j]) for j in range(4)))
           for i in range(4)]
    # The convention div(P)^i = -sum_j d_j P^{ij}.
    assert mod == [-a*x1, a*x2, -2*a*x3, 2*a*x4]
    assert pf(tensor(kk, 0, 0, 0, 0)) == 0
    jj = jacobi(m, xx)
    expected = [
        -2*b*c*x3 - q*(b+y)*x1**2,
        2*b*c*x4 + p*(a+2*b+y)*x2**2,
        (d*(a+b)+p*q)*x1 - 2*p*(a+2*b+y)*x2*x3,
        -(d*(a+b)+p*q)*x2 + 2*q*(b+y)*x1*x4,
    ]
    assert all(s.expand(f-g) == 0 for f, g in zip(jj, expected)), jj

    # Every p-only branch, and the p=q=0 constant branches, is included.
    families = {
        'p_only_C0': ((b, -a-2*b, b+t, a, b, -2*(a+b)), 0, 0, p, 0, a*t),
        'p_only_c_nonzero': ((t, -a, -t, a, t, -2*a), c, 0, p, 0, -2*a*t),
        'p_only_d_nonzero': ((-a, a+t, -a-t, a, -a, t), 0, d, p, 0, -a*t),
        'L0_c_nonzero': ((t, y, -a-y, -y, a+y, -2*a), c, 0, 0, 0, -2*a*t),
        'L0_d_nonzero': ((-a, y, -y, 2*a-y, y-2*a, t), 0, d, 0, 0, -a*t),
    }
    checks = {}
    for name, (k, cc, dd, pp, qq, kappa) in families.items():
        mt = tensor(k, cc, dd, pp, qq)
        assert jacobi(mt, xx) == [0]*4, name
        assert s.expand(pf(tensor(k, 0, 0, 0, 0))-kappa*x1*x2*x3*x4) == 0
        checks[name] = {'all_parameter_jacobi': True, 'highest_pfaffian': str(kappa)}
    return {
        'status': 'PASS',
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'ring': 'QQ[a,b,y,c,d,p,q,t,x1,x2,x3,x4]',
        'localization': 'a != 0; c or d nonzero only in their named coverage charts',
        'full_resonant_coefficients': list(map(str, kk)),
        'full_jacobi_equations': list(map(str, jj)),
        'source_modular_eigenvalues': ['a', '-a', '2*a', '-2*a'],
        'families': checks,
        'logical_scope': 'All invertible regular-semisimple modular paired resonant spectra; rank4 points only. Singular and nonsemisimple A remain separate.',
        'seconds': round(time.monotonic()-started, 4),
    }


if __name__ == '__main__':
    result = run()
    output = Path('reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RESONANT_20261107.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'report': str(output), 'seconds': result['seconds']}))
