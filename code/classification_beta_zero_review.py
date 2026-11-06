#!/usr/bin/env python3
"""Independent finite checks for E_zero; old support/component theorems remain inputs."""
from pathlib import Path
from itertools import combinations, combinations_with_replacement
import argparse
import hashlib
import json
import sys
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
h, x, y, z, w, eps = s.symbols('h x y z w epsilon')
A, B, D, E = s.symbols('A B D E')
t = s.symbols('t1:10')
V, C = [x, y, z, w], [h, x, y, z, w]
PAIRS = list(combinations(range(4), 2))
MONS = [V[i]*V[j] for i, j in combinations_with_replacement(range(4), 2)]
LOCALS = {str(v): v for v in C + [A, B, D, E] + list(t)}
J = s.zeros(4)
J[0, 1] = J[2, 3] = 1
J[1, 0] = J[3, 2] = -1

def parse(value):
    return s.sympify(value, locals=LOCALS)

def matrix(values):
    out = s.zeros(4)
    for (i, j), value in zip(PAIRS, values, strict=True):
        out[i, j], out[j, i] = value, -value
    return out

def normal(a, b, c, d, e, p=0, r=0):
    return matrix([a*x*y+b*z*w, e*x*z, (c-e)*x*w+p*y*y,
                   (-b-e)*y*z+r*x*x, (b-c+e)*y*w, c*x*y+d*z*w])

def jacobi(P):
    return [s.expand(sum(P[i,k]*s.diff(P[j,l],V[k])
                        +P[j,k]*s.diff(P[l,i],V[k])
                        +P[l,k]*s.diff(P[i,j],V[k]) for k in range(4)))
            for i,j,l in combinations(range(4),3)]

def determinant(mat, field):
    """Ordinary field Gaussian elimination, separate from producer DomainMatrix."""
    n = len(mat)
    a = [[field.from_sympy(s.expand(v)) for v in row] for row in mat]
    det = field.one
    for j in range(n):
        pivot = next((k for k in range(j,n) if a[k][j]), None)
        assert pivot is not None
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            det = -det
        d = a[j][j]
        det *= d
        for k in range(j+1,n):
            factor = a[k][j]/d
            for l in range(j+1,n):
                a[k][l] -= factor*a[j][l]
    return field.to_sympy(det)

def review(output):
    started = time.monotonic()
    job_path = 'work/classification_E_zero_review_job.json'
    job = json.loads((ROOT/job_path).read_text())
    for item in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest() == item['sha256']
    producer = 'reports/CLASSIFICATION_BETA_ZERO_FULL_COVERAGE_20261106.json'
    report = json.loads((ROOT/producer).read_text())
    # Reconstruct the ten normal forms from the old Q22 displayed formulas.
    Q = {'G1':normal(A,0,0,D,E), 'G2':normal(A,B,0,B,E),
         'T1_Q45':normal(1,0,0,3,1,1),
         'T2_Q45':normal(s.Rational(1,2),0,s.Rational(1,2),3,1,1),
         'T3_Q45':normal(1,s.Rational(3,2),0,s.Rational(3,2),-s.Rational(1,2),1),
         'T4_Q45':normal(-1,0,2,3,1,1,3),
         'E1':matrix([x*y,x*z,-2*x*y-x*w,-y*z,y*y+y*w,z*w]),
         'E2':matrix([x*y/2,x*z,-2*x*y-x*w/2,-y*z,y*y+y*w/2,x*y/2+z*w]),
         'E3':matrix([x*y,x*z/2,-x*w/2+y*z,-y*z/2,y*w/2,z*w]),
         'E4':matrix([x*y/2,x*z/2,y*z,-y*z/2,0,x*y/2+z*w])}
    kappas = [A*D,A*B,3,s.Rational(3,2),s.Rational(3,2),-3,1,s.Rational(1,2),1,s.Rational(1,2)]
    limits = []
    assert [row['Q45_component'] for row in report['nonnil_ten']] == list(Q)
    for row, kappa in zip(report['nonnil_ten'], kappas, strict=True):
        name = row['Q45_component']
        q, weights = Q[name], row['pullback_weights']
        assert all(v == 0 for v in jacobi(q) + jacobi(J+q))
        P, L = h*h*J+q, s.zeros(4)
        for i,j in PAIRS:
            transformed = s.expand(P[i,j].xreplace({c:eps**weight*c for c,weight in zip(C,weights)})
                                   /eps**(weights[i+1]+weights[j+1]))
            poly = s.Poly(transformed,eps)  # Rejects all negative powers.
            L[i,j] = poly.coeff_monomial(1)
            L[j,i] = -L[i,j]
            assert s.expand(L[i,j]-parse(row['limit'][f'{i+1},{j+1}'])) == 0
        f = s.expand(L[0,1]*L[2,3]-L[0,2]*L[1,3]+L[0,3]*L[1,2])
        assert s.expand(f-kappa*x*y*z*w) == 0
        assert s.expand(h*f-parse(row['projective_limit_degeneracy_divisor'])) == 0
        limits.append(dict(Q45_component=name, nonnegative_epsilon_powers=True,
                           Pfaffian=str(f), compatible_and_Jacobi_exact=True))
    print('PASS: ten weighted limits', flush=True)
    minimal = []
    for row in report['whole_minimal_two']:
        c = parse(row['c'])
        a = t[0]*x*x+(-2*c-t[5])*x*z-t[7]*z*z/2-2*t[8]*z*w+t[1]*x*w+t[2]*w*w
        b = t[3]*x*x+t[4]*x*z+t[5]*x*w+t[6]*z*z+t[7]*z*w+t[8]*w*w
        source = h*h*J+matrix([(s.Rational(1,2)-c)*x*x,0,0,a,b,c*x*x])
        f, H = parse(row['f']), parse(row['potential'])
        extra = 0 if c == 0 else -s.diff(f,x)/2
        target = matrix([f,0,0,extra*z-s.diff(H,w),extra*w+s.diff(H,z),h*h if c == 0 else f])
        assert all(s.expand(v) == 0 for v in source-target)
        minimal.append(row['Q45_component'])
    print('PASS: two full minimal families', flush=True)
    dense, tau = [], s.symbols('tau')
    for row in report['dense_orbits']:
        cert = row['orbit_certificate']
        q = matrix([parse(row['representative'][f'{i+1},{j+1}']) for i,j in PAIRS])
        assert all(v == 0 for v in jacobi(q)+jacobi(J+q))
        columns = []
        for k,l in combinations_with_replacement(range(4),2):
            S = s.zeros(4)
            S[k,l] = S[l,k] = 1
            T = J*S
            F, inverse = s.eye(4)+tau*T, s.eye(4)-tau*T
            sub = dict(zip(V,inverse*s.Matrix(V)))
            acted = F*q.applyfunc(lambda value:value.xreplace(sub))*F.T
            direction = acted.applyfunc(lambda value:s.expand(value).coeff(tau,1))
            columns.append([s.Poly(direction[i,j],*V).coeff_monomial(m) for i,j in PAIRS for m in MONS])
        columns.append([s.Poly(q[i,j],*V).coeff_monomial(m) for i,j in PAIRS for m in MONS])
        actual = s.Matrix.hstack(*(s.Matrix(col) for col in columns))
        saved = s.Matrix([[parse(v) for v in line] for line in cert['original60by11_matrix']])
        assert all(s.expand(v) == 0 for v in actual-saved)
        rows, cols = cert['rows'], cert['columns']
        assert len(rows) == len(cols) == 10
        field = s.QQ.algebraic_field(s.sqrt(3)) if row['Q45_component'] == 'Dpm' else s.QQ
        det = determinant([[actual[i,j] for j in cols] for i in rows],field)
        assert field.from_sympy(s.expand(det-parse(cert['determinant']))) == field.zero
        assert field.from_sympy(det) != field.zero
        dense.append(dict(Q45_component=row['Q45_component'], full_matrix_equal=True,
                          minor_size=10, determinant=str(det), coefficient_field=str(field)))
        print('PASS: dense orbit',row['Q45_component'],flush=True)
    # Match actual global catalogue families, without running their rank producers.
    sys.path.insert(0,str(ROOT/'code'))
    import additional_rank4_families_20261103_v62 as catalogue
    families = {record['component_id']:record for record in catalogue.families()}
    membership = []
    for row in report['dense_orbits']:
        record = families[row['known_target']]
        P = record['P'].subs(dict(zip(catalogue.X,C)),simultaneous=True)
        if row['Q45_component'] == 'Dpm':
            P = P.subs(s.Symbol('r'),s.sqrt(3))
        if row['Q45_component'] == 'Cplus':
            P = P.subs({s.Symbol('c'):-1,s.Symbol('beta'):0})
        q = matrix([parse(row['representative'][f'{i+1},{j+1}']) for i,j in PAIRS])
        assert all(s.expand(P[i+1,j+1]-h*h*J[i,j]-q[i,j]) == 0 for i,j in PAIRS)
        assert all(P[0,i] == 0 for i in range(5))
        membership.append(dict(Q45_component=row['Q45_component'],
                               global_family=row['known_target'], all_coefficients_equal=True,
                               fixed_sqrt3_embedding=row['Q45_component'] == 'Dpm'))
    NB = matrix([3*x*w/2,-3*x*x/2,3*x*z/2,3*x*y/2-z*w/2,
                 -3*y*z/2+w*w/2,-x*w/2+z*z/2])
    P = families['R4_TRIPLE_HYPERPLANE_QUADRIC']['P'].subs(dict(zip(catalogue.X,C)),simultaneous=True)
    assert all(s.expand(P[i+1,j+1]-h*h*J[i,j]-2*NB[i,j]) == 0 for i,j in PAIRS)
    a,b = s.symbols('a b')
    zero = matrix([-a*b*x*z,0,-b*b*x*x,-a*a*z*z,
                   b*(2*b-a)*x*y-a*(b-2*a)*z*w,a*b*x*z])
    P = families['R4_ZERO_LOG_5H']['P'].subs(dict(zip(catalogue.X,C)),simultaneous=True)
    assert all(s.expand(P[i+1,j+1]-h*h*J[i,j]-zero[i,j]) == 0 for i,j in PAIRS)
    print('PASS: five global representatives, NB scaling and entire L0 family',flush=True)
    body = dict(task_id='E_zero', all_checked=True, ten_weighted_limits=limits,
                full_nine_parameter_minimal_inclusions=minimal, five_dense_orbit_minor_checks=dense,
                five_global_representative_memberships=membership,
                NB_twice_representative_equals_global_quadratic_part=True,
                entire_L0_two_parameter_family_equals_global_family=True,
                scope='Finite new certificates only; old Q45_22 support and component-identification theorems remain explicit inputs.',
                input_sha256={item['path']:item['sha256'] for item in job['exact_inputs_and_hashes']},
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                elapsed_seconds=round(time.monotonic()-started,3),
                resource_limits=dict(timeout_seconds=60,virtual_memory_mib=1536))
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as stream:
        json.dump(body,stream,indent=2)
        stream.write('\n')
    print(json.dumps(dict(all_checked=True,elapsed_seconds=body['elapsed_seconds'])))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    review(parser.parse_args().output)
