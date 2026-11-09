#!/usr/bin/env python3
"""Independent exact identities for raw sources, before S61 normalization.

Reads coefficient formulae, never producer Jacobi routines or their PASS flags.
Checks universal projective Jacobi, not just a rational sample. This does not
prove source density, componenthood, coverage, or the serialized DAG identity.
"""
import argparse
import hashlib
import itertools
import json
import re
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[2]
X = s.symbols('x0:5')
PAIRS = list(itertools.combinations(range(5), 2))
TRIPLES = list(itertools.combinations(range(5), 3))


def reduce(v, root=False):
    numerator = s.cancel(v).as_numer_denom()[0]
    if root:
        r = s.Symbol('r')
        numerator = s.rem(numerator, r*r-3, r)
    return s.expand(numerator)


def matrix(values):
    P = s.zeros(5)
    for key, value in values.items():
        i, j = map(int, key.split(','))
        names = re.findall(r'\b[A-Za-z_]\w*\b', str(value))
        P[i, j] = s.sympify(value, locals={name: s.Symbol(name) for name in names})
        P[j, i] = -P[i, j]
    return P


def check(name, P, root=False):
    # Clearing parameter-only denominators preserves projective Jacobi up to
    # their square, and sparse polynomial arithmetic avoids expression swell.
    from sympy.polys.rings import ring
    denominator = s.S.One
    for v in P:
        denominator = s.lcm(denominator, s.cancel(v).as_numer_denom()[1])
    assert not (denominator.free_symbols & set(X))
    parameters = sorted(set().union(*(v.free_symbols for v in P))-set(X), key=str)
    R, *generators = ring([*parameters,*X], s.QQ)
    x = generators[-5:]
    entries = {(i,j):R.from_expr(s.cancel(P[i,j]*denominator)) for i in range(5) for j in range(5)}
    J = {}
    for i, j, k in TRIPLES:
        J[i, j, k] = sum((entries[i,l]*entries[j,k].diff(x[l])
                         + entries[j,l]*entries[k,i].diff(x[l])
                         + entries[k,l]*entries[i,j].diff(x[l])
                         for l in range(5)), R.zero)
    # E wedge [P,P]=0 is invariant under Euler gauge; strict cone Jacobi
    # need not vanish for the supplied homogeneous representative.
    for indices in itertools.combinations(range(5), 4):
        identity = sum(((-1)**a*x[i]*J[tuple(j for j in indices if j != i)]
                       for a, i in enumerate(indices)), R.zero)
        if root:
            r = generators[parameters.index(s.Symbol('r'))]
            identity = identity.rem(r*r-3)
        assert not identity, (name, indices)
    assert all(s.Poly(v, *X).is_homogeneous and s.Poly(v, *X).total_degree() == 2
               for v in P if v != 0), name
    print(name, 'universal projective Jacobi', flush=True)
    return {'id': name, 'universal_projective_Jacobi': True,
            'field': 'QQ(sqrt(3))' if root else 'QQ(parameters)'}


def last_five():
    h, x, y, z, w = X
    coefficients = s.symbols('r0:6')
    R = sum(a*m for a, m in zip(coefficients, [h*h, h*x, h*y, x*x, x*y, y*y]))
    tensors = [
        ('G', {(1,2):4*h*h-y*y, (1,3):-h*x+2*y*z+4*w*w,
               (1,4):4*y*(4*h-w), (2,3):h*y, (3,4):4*(y*y-h*w)}),
        ('X1', {(1,3):2*x*(x-h), (1,4):-2*x*y, (2,3):y*(3*x-h),
                (2,4):2*x*(h-x)-3*y*y, (3,4):-(x-h)*w+R}),
        ('X2', {(1,3):2*x*(x+h), (1,4):2*x*y, (2,3):3*x*y,
                (2,4):3*(y*y-(h+x)**2), (3,4):-x*w+R}),
        ('X3', {(1,3):-h*x, (1,4):3*x*(h-x), (2,3):3*(h-x)**2-3*h*y,
                (2,4):-9*x*y, (3,4):3*x*z+3*R})]
    # D uses the same index order, with variables (h,t,x,y,r).
    h, t, x, y, r = X
    tensors.append(('D', {(1,2):h*h+h*x-2*x*y, (1,3):y*(4*y-h),
                          (1,4):(h-6*y)*r+h*x, (2,4):-x*x, (3,4):h*h+2*x*y}))
    result = []
    for name, data in tensors:
        P = matrix({f'{i},{j}': value for (i,j), value in data.items()})
        row = check('last_five_'+name, P)
        pf = s.factor(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3])
        assert pf != 0
        assert not (pf.free_symbols & set(coefficients))
        row['actual_quintic'] = str(s.factor(X[0]*pf))
        row['all_six_magnetic_parameters_retained'] = name.startswith('X')
        result.append(row)
    return result


def spinor():
    t, z = X[1], X[4]
    r0,r1,r2,r3,r4,r5,r6,r7 = s.symbols('r0:8')
    a = r0*t-r4*z
    b = -r6*X[0]+r1*t-r4*X[2]+r2*z
    c = r3*t+r6*z
    d = r4*X[0]+r5*t+r6*X[3]+r7*z
    P = s.zeros(5)
    for (i,j), v in zip(PAIRS, [b*d,-b*c,-a*d,a*c,b*b,d*d,-a*b-c*d,a*b-c*d,c*c,a*a]):
        if 1 not in (i,j): P[i,j],P[j,i] = v,-v
    assert all(s.expand(sum(s.diff(P[i,j],X[i]) for i in range(5))) == 0 for j in range(5))
    Q1 = r0*X[0]*t-r3*t*X[3]-r4*X[0]*z-r5*t*z-r6*X[3]*z-r7*z*z/2
    Q2 = r0*t*X[2]+r1*t*z+r2*z*z/2-r3*X[0]*t-r4*X[2]*z-r6*X[0]*z
    gradients = s.Matrix([[s.diff(f,x) for x in X] for f in [t,Q1,Q2]])
    for i,j,k in TRIPLES:
        l,m = [v for v in range(5) if v not in (i,j,k)]
        perm = [l,m,i,j,k]
        sign = (-1)**sum(perm[u]>perm[v] for u in range(5) for v in range(u+1,5))
        assert s.expand(gradients[:,[i,j,k]].det()+sign*P[l,m]) == 0
    return {'universal_eight_parameter_spinor_Jacobian_identity': True}


def elliptic_graphs():
    """Decode the two small DAGs symbolically and reconstruct without S61 code."""
    import gzip
    cert = json.loads((ROOT/'legacy/v02/data/symplectic_slice_certificate.json').read_text())
    equations = json.loads((ROOT/'legacy/v02/data/symplectic_slice_equations.json').read_text())
    rows = []
    for name in ['T39', 'T40']:
        path = ROOT/f'reports/MPS_ELLIPTIC_GRAPH_CIRCUITS_20261103_v62/{name}.json.gz'
        data = json.loads(gzip.decompress(path.read_bytes()))
        variables = [s.Symbol(v) for v in data['input_names']]
        values = []
        for node in data['nodes']:
            op = node[0]
            if op == 'q': v = s.Rational(node[1])
            elif op == 'v': v = variables[node[1]]
            elif op == '+': v = sum(values[i] for i in node[1])
            elif op == '*': v = s.prod(values[i] for i in node[1])
            elif op == 'inv': v = 1/values[node[1]]
            else: raise ValueError(op)
            values.append(v)
        u = [s.expand(values[i]) for i in data['output_u1_through_u61']]
        monomials = [X[i]*X[j] for i,j in itertools.combinations_with_replacement(range(5),2)]
        coefficients = [s.S.Zero]*150
        coefficients[PAIRS.index((1,2))*15] = X[0]**0
        coefficients[PAIRS.index((3,4))*15] = X[0]**0
        for index, value in zip(cert['raw_free_indices'], u, strict=True): coefficients[index] = value
        for entry in equations['eliminated_raw_variables']:
            coefficients[entry['raw_index']] = s.expand(sum(s.Rational(a,b)*s.prod(u[j] for j in indices)
                for indices,a,b in entry['polynomial']))
        P = s.zeros(5)
        for a,(i,j) in enumerate(PAIRS):
            P[i,j] = s.expand(sum(c*m for c,m in zip(coefficients[15*a:15*(a+1)],monomials)))
            P[j,i] = -P[i,j]
        row = check(name, P)
        row.update(input=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                   symbolic_DAG_inputs=len(variables),source='serialized DAG, not source PASS flags')
        rows.append(row)
    return rows


def main(output):
    rows = []
    for n in range(2,39):
        raw = ROOT/f'reports/MPS_S61_RAW_CATALOGUE_20261103_v62/T{n:02d}.json'
        completed = ROOT/f'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T{n:02d}.json'
        path = completed if completed.exists() else raw
        data = json.loads(path.read_text())
        row = check(data['component_id'], matrix(data['homogeneous_bivector']))
        row.update(input=str(path.relative_to(ROOT)), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        rows.append(row)
    for path in sorted((ROOT/'reports/ADDITIONAL_RANK4_FAMILIES_20261103_v62').glob('R4*.json')):
        data = json.loads(path.read_text())
        row = check(data['component_id'], matrix(data['homogeneous_bivector']), data['component_id']=='R4_DPM')
        row.update(input=str(path.relative_to(ROOT)), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        rows.append(row)
    toric = s.zeros(5)
    for v,(i,j) in zip(s.symbols('b0:6'), itertools.combinations(range(1,5),2)):
        toric[i,j],toric[j,i] = v*X[i]*X[j],-v*X[i]*X[j]
    rows.append(check('T01', toric))
    rows.extend(elliptic_graphs())
    result = {'scope': __doc__, 'raw_families':rows, 'last_five':last_five(),
              'rank_two':spinor(), 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)
