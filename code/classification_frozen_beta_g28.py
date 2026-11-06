#!/usr/bin/env python3
"""One exact QQ G28 evaluation at the existing frozen T01 first jet only."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations, product
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT = Path(__file__).resolve().parents[1]
POINT = 'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'
BRIDGE = 'reports/CLASSIFICATION_SLICE_BRIDGE_20261106.json'
SPIN = 'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
INFO = 'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
EQUATION = 'reports/Q_NEW_FAMILY_DEGREE28_EXACT_EQUATION_20261104_v63/equation.json'

def read(path):
    raw = (ROOT/path).read_bytes()
    return json.loads(gzip.decompress(raw) if path.endswith('.gz') else raw)

def numeric(a):
    return [[F(c) for c in row] for row in a]

def multiply(a,b):
    columns = list(zip(*b))
    return [[sum((x*y for x,y in zip(row,col) if x and y),F(0)) for col in columns] for row in a]

def trace(a):
    return sum((row[i] for i,row in enumerate(a)),F(0))

def tracepair(a,b):
    return sum((a[i][j]*b[j][i] for i in range(len(a)) for j in range(len(a))),F(0))

def inverse(a):
    n = len(a)
    rows = [list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for j in range(n):
        k = next(i for i in range(j,n) if rows[i][j])
        rows[j],rows[k] = rows[k],rows[j]
        c = rows[j][j]
        rows[j] = [v/c for v in rows[j]]
        for i in range(n):
            if i!=j and rows[i][j]:
                c = rows[i][j]
                rows[i] = [v-c*w for v,w in zip(rows[i],rows[j])]
    return [row[n:] for row in rows]

def total_product(values):
    result = F(1)
    for value in values:
        result *= value
    return result

def evaluate(output,source_sha256):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    source_path = ROOT/'work/frozen-beta-g28-source.json'
    assert hashlib.sha256(source_path.read_bytes()).hexdigest() == source_sha256
    source = json.loads(source_path.read_text())
    for path,digest in source['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest,path
    point = read(POINT)
    original = next(row for row in point['points'] if row['component_id']=='T01')
    u61 = list(map(F,original['S61_coordinates_u1_through_u61']))
    bridge,spin,info,equation = map(read,[BRIDGE,SPIN,INFO,EQUATION])
    ids = bridge['original_weight1_parameter_indices0']
    assert ids == [3,4,11,12,13,24,25,26,33,34,35,36,47,48,49,50]
    C16 = numeric(bridge['primitive_basis_coordinate_bridge16'])
    b = [sum((c*u61[i] for c,i in zip(row,ids)),F(0)) for row in C16]
    # The invariant input uses the canonical primitive basis, not saved u16.
    from weight1_structure import constraints,sp_basis
    import sympy as s
    canonical = s.Matrix.hstack(*s.Matrix.vstack(*constraints()).nullspace())
    assert canonical == s.Matrix(info['primitive24_by16_basis'])
    assert canonical*s.Matrix(C16) == s.Matrix(bridge['primitive_embedding24x16'])
    mon = [tuple(v) for v in spin['quadratic_monomials']]
    assert mon == list(combinations_with_replacement(range(16),2))
    quad = [b[i]*b[j] for i,j in mon]
    cofs = [numeric(spin[key]) for key in
            ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']]
    CD = []
    for coefficient_rows in cofs:
        values = [sum((c*v for c,v in zip(row,quad)),F(0)) for row in coefficient_rows]
        CD.append([values[4*i:4*i+4] for i in range(4)])
    C,D = CD
    mu = [[-C[j][i]+2*D[j][i] for j in range(4)] for i in range(4)]
    Omega = numeric(info['symplectic_form16'])
    Oi = inverse(Omega)
    rho = [numeric(a) for a in info['all_ten_representation_matrices']]
    gs = [numeric(a.tolist()) for a in sp_basis()]
    Gram = [[tracepair(a,c) for c in gs] for a in gs]
    moments = []
    for action in rho:
        symmetric = multiply(Omega,action)
        moments.append(sum((b[i]*symmetric[i][j]*b[j] for i in range(16) for j in range(16)),F(0))/2)
    coefficients = [sum((c*v for c,v in zip(row,moments)),F(0)) for row in inverse(Gram)]
    actual_mu = [[sum((c*a[i][j] for c,a in zip(coefficients,gs)),F(0)) for j in range(4)] for i in range(4)]
    assert actual_mu == mu
    theta = [[sum((c*a[i][j] for c,a in zip(coefficients,rho)),F(0)) for j in range(16)] for i in range(16)]
    quartic = spin['quartic_trace_moment_square']
    F4 = tracepair(mu,mu)
    assert F4 == sum((F(c)*total_product(b[i] for i in word) for word,c in quartic),F(0))
    Hess = [[F(0) for _ in range(16)] for _ in range(16)]
    for word,c in quartic:
        counts = Counter(word)
        for i,ei in counts.items():
            for j,ej in counts.items():
                multiplicity = ei*(ej-int(i==j))
                if multiplicity:
                    rest = list(word)
                    rest.remove(i)
                    rest.remove(j)
                    Hess[i][j] += F(c)*multiplicity*total_product(b[k] for k in rest)/2
    # Independently recover Hess(F/2) from the global C,D quadratics.
    mu_cofs = [[-cofs[0][4*j+i][k]+2*cofs[1][4*j+i][k] for k in range(136)] for i in range(4) for j in range(4)]
    first = [[[F(0) for _ in range(4)] for _ in range(4)] for _ in range(16)]
    second = [[[[F(0) for _ in range(4)] for _ in range(4)] for _ in range(16)] for _ in range(16)]
    for pos,row in enumerate(mu_cofs):
        a,c = divmod(pos,4)
        for (i,j),value in zip(mon,row):
            first[i][a][c] += value*b[j]
            first[j][a][c] += value*b[i]
            second[i][j][a][c] += value
            second[j][i][a][c] += value
    assert Hess == [[tracepair(first[i],first[j])+tracepair(mu,second[i][j]) for j in range(16)] for i in range(16)]
    H16 = multiply(Oi,Hess)
    pencil = [F(0)]*5
    for perm in permutations(range(4)):
        sign = (-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        for letters in product(range(2),repeat=4):
            pencil[sum(letters)] += sign*total_product(CD[letters[i]][i][perm[i]] for i in range(4))
    C2,D2 = multiply(C,C),multiply(D,D)
    C3,D3 = multiply(C2,C),multiply(D2,D)
    CD2 = multiply(C,D2)
    Q0 = tracepair(CD2,CD2)-trace(D2)*tracepair(C2,D2)/2+pencil[4]*trace(C2)
    Q1 = trace(multiply(multiply(multiply(C3,D),C),D))
    Q2 = tracepair(C3,D3)
    level = H16
    Q3,T = None,None
    for degree in range(1,8):
        level = multiply(theta,level)
        if degree==5:
            Q3 = trace(level)
        if degree==7:
            T = trace(level)
    assert equation['name']=='G28' and len(equation['terms'])==78
    assert all(not any(row['exponents'][10:]) for row in equation['terms'])
    values = [F4,*pencil[1:],Q0,Q1,Q2,Q3,T,F(0),F(0),F(0)]
    assert len(values)==len(equation['variables'])==13
    value = sum((F(row['coefficient'])*total_product(v**e for v,e in zip(values,row['exponents'])) for row in equation['terms']),F(0))
    body = dict(task_id='G',coefficient_field='QQ',frozen_T01_beta_only=True,
                primitive_b16=list(map(str,b)),bridge_and_canonical_basis_checked=True,
                moment_C_D_and_Omega_rho_agree=True,Hessian_two_exact_recipes_agree=True,
                invariants11=dict(zip(['F','p0','p1','p2','p3','p4','Q0','Q1','Q2','Q3','T'],
                                      map(str,[F4,*pencil,Q0,Q1,Q2,Q3,T]))),
                A14_generators_not_evaluated_and_absent_from_all78_terms=True,
                G28_value_QQ=str(value),frozen_beta_on_G28=(value==0),
                generic_exhaustion_or_fiber_unit_claimed=False,
                input_sha256=source['input_sha256'],source_sha256=source_sha256,
                elapsed_seconds=round(time.monotonic()-started,3),
                limits=dict(seconds=60,virtual_memory_bytes=1024**3))
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as stream:
        json.dump(body,stream,indent=2)
        stream.write('\n')
    print(json.dumps({key:body[key] for key in ['G28_value_QQ','frozen_beta_on_G28','elapsed_seconds']}))

if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--source-sha256',required=True)
    args = parser.parse_args()
    evaluate(args.output,args.source_sha256)
