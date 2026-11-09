#!/usr/bin/env python3
"""Complete quadratic W sources for the eleven nilpotent cofactor normal forms."""
import argparse
import hashlib
import itertools
import json
import resource
import time
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU, (55, 56))
    resource.setrlimit(resource.RLIMIT_AS, (1536 * 1024**2,) * 2)
    start = time.monotonic()
    h, x, u, v, w = xx = s.symbols('h x u v w')
    mons = [xx[i] * xx[j] for i, j in itertools.combinations_with_replacement(range(5), 2)]
    unknowns = s.symbols('a0:45')
    forms = [('j3_zero', (v, w, 0)), ('j3_h', (v, w, h)),
             ('j3_hx', (v, w, h + x))]
    for name, L, M in [('zero', 0, 0), ('0_h', 0, h), ('0_hx', 0, h+x),
                       ('h_0', h, 0), ('hx_0', h+x, 0),
                       ('h_x', h, x), ('hx_h', h+x, h)]:
        forms.append(('j21_' + name, (v, L, M)))
    forms.append(('j111_rank2', (x, h, 0)))
    records = []
    for name, V in forms:
        P = s.zeros(5)
        for i, value in enumerate(V, 2):
            P[1, i] = x * value
        for k, (i, j) in enumerate([(2, 3), (2, 4), (3, 4)]):
            P[i, j] = sum(a*m for a, m in zip(unknowns[15*k:15*(k+1)], mons))
        for i, j in itertools.combinations(range(5), 2):
            P[j, i] = -P[i, j]

        def jac(i, j, k, matrix):
            return s.expand(sum(matrix[i, l]*s.diff(matrix[j, k], xx[l]) +
                                matrix[j, l]*s.diff(matrix[k, i], xx[l]) +
                                matrix[k, l]*s.diff(matrix[i, j], xx[l]) for l in range(5)))

        linear = []
        for i, j in itertools.combinations(range(2, 5), 2):
            linear += s.Poly(jac(1, i, j, P), *xx).coeffs()
        matrix, rhs = s.linear_eq_to_matrix(linear, unknowns)
        aug, pivots = matrix.row_join(rhs).rref()
        if 45 in pivots:
            records.append(dict(case=name, V=[str(f) for f in V], inconsistent=True))
            print(json.dumps(dict(case=name, inconsistent=True)), flush=True)
            continue
        free = [j for j in range(45) if j not in pivots]
        params = s.symbols('n0:' + str(len(free)))
        solution = {unknowns[j]: z for j, z in zip(free, params)}
        for i, pivot in enumerate(pivots):
            solution[unknowns[pivot]] = aug[i, 45] - sum(aug[i, j]*solution[unknowns[j]] for j in free)
        Q = P.subs(solution, simultaneous=True).applyfunc(s.expand)
        assert all(jac(1, i, j, Q) == 0 for i, j in itertools.combinations(range(2, 5), 2))
        residual = jac(2, 3, 4, Q)
        equations = list(dict.fromkeys(str(s.factor(f)) for f in s.Poly(residual, *xx).coeffs() if f))
        F = s.expand(V[0]*Q[3, 4]-V[1]*Q[2, 4]+V[2]*Q[2, 3])
        assert s.expand(sum(V[i]*s.diff(F, xx[i+2]) for i in range(3))) == 0
        records.append(dict(case=name, V=[str(f) for f in V], inconsistent=False,
            linear_matrix_shape=list(matrix.shape), linear_rank=len(pivots),
            free_parameters=list(map(str, params)), dimension=len(free),
            original_unknown_monomials=[str(m) for m in mons], free_columns=free,
            linear_rref=[[str(f) for f in aug.row(i)] for i in range(len(pivots))],
            brackets={','.join(map(str, ij)): str(Q[ij]) for ij in itertools.combinations(range(1,5),2)},
            actual_cubic=str(s.factor(F)), remaining_full_Jacobi=str(residual),
            remaining_coefficient_equations=equations))
        print(json.dumps(dict(case=name, dimension=len(free), equations=len(equations), F=str(s.factor(F)))), flush=True)
    source = Path(__file__).read_text()
    result = dict(task_id='E_reduced_nonnormal', cases=records,
        all_45_vertical_coefficients_retained=True,
        all_four_affine_Jacobi_equations_retained=True,
        normal_form_exhaustion_requires_the_written_centralizer_argument=True,
        irreducibility_and_MPS_membership_not_inferred_from_this_linear_calculation=True,
        seconds=time.monotonic()-start, producer_source=source,
        input_sha256={'code/classification_two_hyperplanes_cubic_nilpotent.py': hashlib.sha256(source.encode()).hexdigest()})
    with output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')


def components(input_path, output):
    """Use the existing full-ambient runner on analytically solved full sources."""
    from classification_modular_source_components import branch
    data = json.loads(input_path.read_text())
    source_path = ROOT / data['source_report']
    source = json.loads(source_path.read_text())
    adapted = {'cases': []}
    h, x, u, v, w = s.symbols('h x u v w')
    yy = s.symbols('x1:5')
    for case in source['cases']:
        if case.get('inconsistent'):
            adapted['cases'].append({})
            continue
        entries = [s.sympify(case['brackets'][key]).subs({h: 1, x: yy[0], u: yy[1], v: yy[2], w: yy[3]}, simultaneous=True)
                   for key in ['1,2', '1,3', '1,4', '2,3', '2,4', '3,4']]
        adapted['cases'].append(dict(full_free_parameter_count=case['dimension'],
            full_six_brackets=list(map(str, entries)),
            all_distinct_jacobi_coefficients=case['remaining_coefficient_equations']))
    results = [branch(adapted, spec) for spec in data['branches']]
    result = dict(status='PASS', task_id='E_reduced_nonnormal', branches=results,
                  source_report_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                  component_input=data, producer_source=Path(__file__).read_text(),
                  full_source_exhaustion_requires_written_proof=True)
    with output.open('x') as handle:
        json.dump(result, handle, indent=2, default=int)
        handle.write('\n')
    print(json.dumps([dict(name=r['name'], ranks=r['matrix_ranks_mod_prime'],
                         upper=r['projective_tangent_dimension_upper'],
                         lower=r['projective_family_dimension_lower']) for r in results]))


def normal_quadrics(report_path, case_index, output):
    """One finite genuine quadratic obstruction, with all actual source directions."""
    from sympy.polys.matrices import DomainMatrix
    from classification_rank3_mixed_component_review import COLUMNS, jacobi, sparse, vector
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU, (55, 56))
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2,)*2)
    record = json.loads(report_path.read_text())['branches'][case_index]
    x = s.symbols('h x1:5')
    def matrix(encoded):
        P = s.zeros(5)
        for key, value in encoded.items():
            i, j = map(int, key.split(',')); P[i,j] = s.sympify(value); P[j,i] = -P[i,j]
        return P
    P = matrix(record['witness_brackets']); ps = sparse(P, x)
    assert jacobi(ps) == s.zeros(350, 1)
    columns = []
    for ij, exponent in COLUMNS:
        variation = {key: dict(value) for key, value in ps.items()}
        variation.setdefault(ij, {})[exponent] = variation.get(ij, {}).get(exponent, 0) + 1
        columns.append(jacobi(variation))
    D = s.Matrix.hstack(*columns)
    columns = []; X = s.Matrix(x)
    for a,b in itertools.product(range(5), repeat=2):
        field=s.zeros(5,1); field[a]=x[b]
        columns.append(vector(sparse(X*field.T-field*X.T,x)))
    for a,b in itertools.product(range(5),repeat=2):
        E=s.zeros(5); E[a,b]=1
        columns.append(vector(sparse(E*P+P*E.T-x[b]*P.diff(x[a]),x)))
    columns.append(vector(ps))
    universal=matrix(record['universal_brackets'])
    point={s.Symbol(k):s.Rational(v) for k,v in record['full_source_specification']['witness'].items()}
    columns.extend(vector(sparse(universal.diff(s.Symbol(p)).subs(point),x))
                   for p in record['full_source_specification']['free_parameters'])
    F=s.Matrix.hstack(*columns); assert D*F == s.zeros(350,F.cols)
    rr,pivots=DomainMatrix.from_Matrix(D).convert_to(s.QQ).rref(); rr=rr.to_Matrix()
    free=[j for j in range(150) if j not in pivots]; kernel=[]
    for j in free:
        col=s.zeros(150,1); col[j]=1
        for i,p in enumerate(pivots): col[p]=-rr[i,j]
        kernel.append(col*s.ilcm(1,*[s.denom(c) for c in col]))
    K=s.Matrix.hstack(*kernel); assert D*K == s.zeros(350,len(free))
    basis=DomainMatrix.from_Matrix(F.row_join(K)).convert_to(s.QQ).rref()[1]
    normal_indices=[j-F.cols for j in basis if j>=F.cols]
    source_rank=len(basis)-len(normal_indices)
    normals=[]
    for j in normal_indices:
        R=s.zeros(5)
        for (ij,e),c in zip(COLUMNS,K[:,j]):
            R[ij]+=c*s.prod(z**power for z,power in zip(x,e))
        for i,k in itertools.combinations(range(5),2): R[k,i]=-R[i,k]
        normals.append(R)
    pairs=list(itertools.combinations_with_replacement(range(len(normals)),2))
    squares=[jacobi(sparse(R,x)) for R in normals]
    Q=s.Matrix.hstack(*[squares[i] if i==j else jacobi(sparse(normals[i]+normals[j],x))-squares[i]-squares[j] for i,j in pairs])
    dind=D[:,list(pivots)]
    rows=list(DomainMatrix.from_Matrix(dind.T).convert_to(s.QQ).rref()[1])
    block=DomainMatrix.from_Matrix(D.extract(rows,list(pivots))).convert_to(s.QQ)
    rhs=DomainMatrix.from_Matrix(Q.extract(rows,range(len(pairs)))).convert_to(s.QQ)
    correction=block.inv().matmul(rhs).to_Matrix()
    obstruction=Q-dind*correction
    qr,qp=DomainMatrix.from_Matrix(obstruction).convert_to(s.QQ).rref();qr=qr.to_Matrix()
    n=s.symbols('t0:'+str(len(normals)));monomials=[n[i]*n[j] for i,j in pairs]
    quadrics=[s.expand(sum(qr[i,j]*monomials[j] for j in range(len(pairs)))) for i in range(len(qp))]
    gb=s.groebner(quadrics,*n,order='grevlex') if quadrics else None
    result=dict(task_id='E_reduced_nonnormal',source_report=str(report_path),case=record['name'],
        full_D_rank=len(pivots),actual_source_rank=source_rank,normal_count=len(normals),
        all_D_times_directions_zero=True,normal_indices=normal_indices,
        normals=[{f'{i},{j}':str(R[i,j]) for i,j in itertools.combinations(range(5),2) if R[i,j]} for R in normals],
        normal_monomials=list(map(str,monomials)),quadratic_initial_generators=list(map(str,quadrics)),
        zero_dimensional_quadratic_support=bool(gb and gb.is_zero_dimensional),
        grevlex_basis=[str(g.as_expr()) for g in gb.polys] if gb else [],
        literal_Poisson_normal_axes=[i for i,q in enumerate(squares) if q==s.zeros(350,1)],
        full_initial_ideal_not_claimed=True,source_report_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),
        producer_source=Path(__file__).read_text(),seconds=time.monotonic()-start)
    with output.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:result[k] for k in ['case','normal_count','quadratic_initial_generators','zero_dimensional_quadratic_support','literal_Poisson_normal_axes','seconds']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--component-input', type=Path)
    parser.add_argument('--normal-report', type=Path)
    parser.add_argument('--case-index', type=int, default=0)
    args = parser.parse_args()
    if args.normal_report:
        normal_quadrics(args.normal_report,args.case_index,args.output)
    elif args.component_input:
        components(args.component_input, args.output)
    else:
        main(args.output)
