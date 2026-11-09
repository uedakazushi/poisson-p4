#!/usr/bin/env python3
"""Complete affine normal forms for nil211 with long-top translation.

This is a finite reduction of full source type10, not a component theorem.
The normal-form mode uses no producer or CAS imports. The optional ambient
mode reuses the independent sparse350 engine, exact QQ matrices and S61.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

from classification_minimal_quartic_review import LaurentRing
from classification_quadratic_cofactor_review import MONS, Ring

ROOT = Path(__file__).resolve().parents[1]
PAIRS = tuple(itertools.combinations(range(4),2))
PARAMETERS = ('alpha','A00','A10','A01','A20','A11','A02','d0','dz','dw','e0','ez',
              'delta','beta','gamma','eta','g10','g01','g20','g11','g02','g30','g21','g12','g03')


def matrix(ring, entries):
    P = [[{} for _ in range(4)] for _ in range(4)]
    for (i,j), value in zip(PAIRS,entries):
        P[i][j],P[j][i] = value,ring.scale(value,-1)
    return P


def jacobi(ring,P):
    return [ring.add(*(ring.mul(P[a][k],ring.derivative(P[b][d],k))
                        for a,b,d in [(i,j,l),(j,l,i),(l,i,j)] for k in range(4)))
            for i,j,l in itertools.combinations(range(4),3)]


def divergence(ring,P):
    return [ring.add(*(ring.derivative(P[i][j],i) for i in range(4))) for j in range(4)]


def pfaffian(ring,entries):
    return ring.add(ring.mul(entries[0],entries[5]),ring.scale(ring.mul(entries[1],entries[4]),-1),
                    ring.mul(entries[2],entries[3]))


def generic_data(ring,r,s,original_coordinates):
    A = ring.parse('A00+A10*z+A01*w+A20*z**2+A11*z*w+A02*w**2')
    d = ring.parse('d0+dz*z+dw*w')
    e = ring.parse('e0+ez*z+(alpha-1-dz)*w')
    L = ring.parse('eta+beta*w-gamma*z')
    G = ring.parse('g10*z+g01*w+g20*z**2+g11*z*w+g02*w**2+g30*z**3+g21*z**2*w+g12*z*w**2+g03*w**3')
    a = ring.add(ring.mul(ring.variable('alpha'),r),A)
    phi = ring.add(ring.scale(ring.mul(ring.variable('delta'),r,r),Q(1,2)),ring.mul(r,L),G)
    p = ring.add(ring.mul(ring.variable('beta'),r),ring.derivative(G,3))
    q = ring.add(ring.mul(ring.variable('gamma'),r),ring.scale(ring.derivative(G,2),-1))
    entries = [a,ring.add(p,ring.mul(s,d)) if original_coordinates else p,
               ring.add(q,ring.mul(s,e)) if original_coordinates else q,d,e,
               ring.add(ring.mul(ring.variable('delta'),r),L)]
    return dict(A=A,d=d,e=e,L=L,G=G,a=a,phi=phi,entries=entries)


def inverse(matrix_):
    n = len(matrix_)
    A = [[Q(value) for value in row]+[Q(i==j) for j in range(n)] for i,row in enumerate(matrix_)]
    determinant = Q(1)
    for j in range(n):
        pivot = next(i for i in range(j,n) if A[i][j])
        if pivot != j:
            A[pivot],A[j] = A[j],A[pivot]
            determinant *= -1
        value = A[j][j]
        determinant *= value
        A[j] = [t/value for t in A[j]]
        for i in range(n):
            if i != j and A[i][j]:
                value = A[i][j]
                A[i] = [a-value*b for a,b in zip(A[i],A[j])]
    assert [row[:n] for row in A] == [[Q(i==j) for j in range(n)] for i in range(n)]
    return [row[n:] for row in A],determinant


def affine_vectors(ring,entries):
    result = []
    for value in entries:
        for monomial in MONS:
            vector = [Q(0)]*(ring.n+1)
            for powers,coefficient in value.items():
                if powers[-4:] == monomial:
                    parameter_powers = powers[:-4]
                    assert min(parameter_powers)>=0 and sum(parameter_powers)<=1
                    index = 0 if not any(parameter_powers) else parameter_powers.index(1)+1
                    vector[index] += coefficient
            result.append(vector)
    return result


def substitute_constant(ring,polynomial,name,value):
    position = ring.names.index(name)
    result = {}
    for powers,coefficient in polynomial.items():
        coefficient *= value**powers[position]
        target = list(powers)
        target[position] = 0
        target = tuple(target)
        result[target] = result.get(target,Q(0))+coefficient
    return {e:c for e,c in result.items() if c}


def group_r_power(ring,polynomial,degree):
    result = {}
    for powers,coefficient in polynomial.items():
        if powers[ring.n] == degree:
            target = list(powers)
            target[ring.n] = 0
            result[tuple(target)] = coefficient
    return result


def verify_complete_source(source):
    ring = LaurentRing(PARAMETERS,('x','y','z','w'))
    r,s = ring.parse('x-y**2/2'),ring.variable('y')
    data = generic_data(ring,r,s,True)
    entries = data['entries']
    assert all(sum(e[-4:])<=2 for value in entries for e in value)
    P = matrix(ring,entries)
    assert divergence(ring,P) == [ring.variable('y'),ring.constant(1),{},{}]
    lie = []
    for i,j in PAIRS:
        value = ring.add(ring.mul(ring.variable('y'),ring.derivative(P[i][j],0)),ring.derivative(P[i][j],1))
        if i == 0:
            value = ring.add(value,ring.scale(P[1][j],-1))
        if j == 0:
            value = ring.add(value,ring.scale(P[i][1],-1))
        lie.append(value)
    assert lie == [{}]*6
    vectors = affine_vectors(ring,entries)
    free = source['free_coefficient_indices']
    assert len(free) == source['free_parameter_count'] == 25
    selected = [vectors[i][1:] for i in free]
    offsets = [vectors[i][0] for i in free]
    backwards,determinant = inverse(selected)
    old_ring = Ring(25)
    # Transport the entire ninety-coefficient affine solution via its free
    # coordinates. The inverse map is explicit, not a dimension-only claim.
    for actual,reported in zip(vectors,source['all_ninety_coefficients']):
        coefficients = [sum(actual[k+1]*backwards[k][j] for k in range(25)) for j in range(25)]
        constant = actual[0]-sum(a*b for a,b in zip(coefficients,offsets))
        polynomial = old_ring.constant(constant)
        for j,value in enumerate(coefficients):
            polynomial = old_ring.add(polynomial,old_ring.scale(old_ring.variable(f'n{j}'),value))
        assert polynomial == old_ring.parse(reported)
    Vphi = ring.add(ring.scale(ring.mul(data['a'],ring.derivative(data['phi'],0)),-1),
                    ring.mul(data['d'],ring.derivative(data['phi'],2)),
                    ring.mul(data['e'],ring.derivative(data['phi'],3)))
    assert pfaffian(ring,entries) == ring.scale(Vphi,-1)
    delta_zero = substitute_constant(ring,pfaffian(ring,entries),'delta',Q(0))
    assert not {e:c for e,c in delta_zero.items() if sum(e[-4:])==4}
    return dict(all_original_ninety_coefficients_matched=True,source_parameter_count=25,
                all_original_degree_two_constraints_preserved=True,
                full_divergence=['y','1','0','0'],full_cofactor_invariance_checked=True,
                free_indices=free,affine_parameter_matrix_determinant=str(determinant),
                inverse_parameter_matrix=[[str(a) for a in row] for row in backwards],
                original_free_coordinate_offsets=list(map(str,offsets)),
                raw_pfaffian_is_minus_V_phi=True,delta_zero_implies_highest_pfaffian_zero=True)


def verify_jacobi_and_square_completion():
    ring = LaurentRing(PARAMETERS,('r','s','z','w'),invert=('delta',))
    r,s = ring.variable('r'),ring.variable('s')
    data = generic_data(ring,r,s,False)
    P = matrix(ring,data['entries'])
    Vphi = ring.add(ring.scale(ring.mul(data['a'],ring.derivative(data['phi'],0)),-1),
                    ring.mul(data['d'],ring.derivative(data['phi'],2)),
                    ring.mul(data['e'],ring.derivative(data['phi'],3)))
    psi = ring.add(Vphi,data['phi'])
    assert jacobi(ring,P) == [ring.scale(ring.derivative(psi,3),-1),ring.derivative(psi,2),{},ring.derivative(psi,0)]
    assert group_r_power(ring,psi,2) == ring.parse('delta*(1/2-alpha)')
    half = lambda polynomial: substitute_constant(ring,polynomial,'alpha',Q(1,2))
    d,e = data['d'],half(data['e'])
    base_L = ring.add(ring.mul(d,ring.derivative(data['L'],2)),ring.mul(e,ring.derivative(data['L'],3)))
    assert group_r_power(ring,half(psi),1) == ring.add(ring.scale(data['L'],Q(1,2)),base_L,
                                                                      ring.scale(ring.mul(ring.variable('delta'),data['A']),-1))
    A = ring.mul(ring.parse('1/delta'),ring.add(ring.scale(data['L'],Q(1,2)),base_L))
    assert all(sum(exponents[-4:])<=1 for exponents in A)
    R = ring.add(r,ring.mul(data['L'],ring.parse('1/delta')))
    H = ring.add(data['G'],ring.scale(ring.mul(data['L'],data['L'],ring.parse('1/delta')),Q(-1,2)))
    Vr = ring.add(ring.scale(r,Q(-1,2)),ring.scale(A,-1))
    assert ring.add(Vr,ring.mul(base_L,ring.parse('1/delta'))) == ring.scale(R,Q(-1,2))
    assert ring.add(data['phi'],ring.scale(ring.mul(ring.variable('delta'),R,R),Q(-1,2))) == H
    normalized_psi = ring.add(ring.mul(Vr,ring.derivative(data['phi'],0)),
                               ring.mul(d,ring.derivative(data['phi'],2)),
                               ring.mul(e,ring.derivative(data['phi'],3)),data['phi'])
    assert normalized_psi == ring.add(ring.mul(d,ring.derivative(H,2)),ring.mul(e,ring.derivative(H,3)),H)
    return dict(full_Jacobi_equivalent_to_psi_constant=True,
                psi_definition='V(phi)+phi',r_squared_coefficient='delta*(1/2-alpha)',
                delta_nonzero_forces_alpha='1/2',
                r_coefficient='L/2+Vbase(L)-delta*A2',A2_is_affine=True,
                R='r+L/delta',H='G-L^2/(2delta)',V_R='-R/2',
                complete_square_identity_checked=True,
                remaining_equation='Vbase(H)+H=constant, trace(linear Vbase)=-1/2')


def weighted_table():
    monomials = [(i,j) for i in range(4) for j in range(4-i) if i+j]
    spectra = set()
    for i,j in monomials:
        if i != j:
            lam = Q(j-2,2*(i-j))
            spectra.add(tuple(sorted((lam,Q(-1,2)-lam))))
        else:
            assert -Q(i,2) != -1
    assert spectra == {(Q(-1),Q(1,2)),(Q(-1,2),Q(0)),(Q(-1,3),Q(-1,6))}
    table = []
    for lam,mu in sorted(spectra):
        allowed = [(i,j) for i,j in monomials if i*lam+j*mu == -1]
        table.append(dict(eigenvalues=[str(lam),str(mu)],all_degree_at_most_three_eigenmonomials=[list(t) for t in allowed]))
    assert not [(i,j) for i,j in monomials if Q(-1,4)*(i+j) == -1]
    zero_ring = LaurentRing(('a','b','kappa'),('z','w','p','q'))
    H = zero_ring.parse('z**2*(a+b*w)')
    assert zero_ring.mul(zero_ring.variable('kappa'),zero_ring.derivative(H,1)) == zero_ring.parse('kappa*b*z**2')
    return dict(all_trace_minus_half_spectra=table,
                repeated_eigenvalue_minus_quarter_has_no_eigenpolynomial=True,
                zero_eigen_translation_equation='kappa*b=0 for H=z^2(a+bw)',
                irreducible_F_forces_b_nonzero=True,
                resulting_normalized_H=['z','z^3','z^2w'])


def actual_models():
    ring = LaurentRing((),('x','y','z','w'))
    models = [dict(name='linear',H='z',lam='-1',mu='1/2',
                   brackets=['x/2-y**2/4','-y*z','-1+y*w/2','-z','w/2','x-y**2/2']),
              dict(name='cubic',H='z**3',lam='-1/3',mu='-1/6',
                   brackets=['x/2-y**2/4','-y*z/3','-3*z**2-y*w/6','-z/3','-w/6','x-y**2/2']),
              dict(name='mixed',H='z**2*w',lam='-1/2',mu='0',
                   brackets=['x/2-y**2/4','z**2-y*z/2','-2*z*w','-z/2','0','x-y**2/2'])]
    for model in models:
        entries = [ring.parse(t) for t in model['brackets']]
        P = matrix(ring,entries)
        assert jacobi(ring,P) == [{}]*4
        assert divergence(ring,P) == [ring.variable('y'),ring.constant(1),{},{}]
        F = ring.parse('(x-y**2/2)**2/2+'+model['H'])
        assert pfaffian(ring,entries) == F
        top = {e:c for e,c in F.items() if sum(e)==4}
        assert top == ring.parse('y**4/8')
        model.update(all_four_Jacobi_zero=True,full_divergence=['y','1','0','0'],
                     affine_pfaffian='(x-y^2/2)^2/2+'+model['H'],highest_pfaffian='y^4/8',
                     no_component_membership_or_componenthood_claim=True)
    return models


def ambient_audit(output):
    import sympy as S
    from sympy.polys.matrices import DomainMatrix
    import classification_rank3_mixed_component_review as independent
    from classification_degenerate_quadric_q4_nil31_I_graph import darboux
    from s61_model import S61, X, pushforward

    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun = time.monotonic()
    source_path = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_NORMAL_FORMS_20261107.json'
    source = json.loads(source_path.read_text())
    coordinates = S.symbols('h x y z w')
    h = coordinates[0]
    env = {str(v):v for v in coordinates}
    adapter = S61()

    def sparse(P):
        return {ij:dict(S.Poly(P[ij],*coordinates).terms()) for ij in independent.PAIRS if P[ij] != 0}

    def rational_matrix(M):
        return DomainMatrix.from_Matrix(M).convert_to(S.QQ)

    def certificate(M):
        A = rational_matrix(M)
        columns = list(A.rref()[1])
        rows = list(A.extract(range(A.shape[0]),columns).transpose().rref()[1])
        determinant = A.extract(rows,columns).det()
        assert determinant != 0
        return dict(shape=list(M.shape),rank_QQ=len(columns),rows=rows,columns=columns,
                    determinant=str(determinant))

    results = []
    for model in source['three_actual_Poisson_models']:
        started = time.monotonic()
        P = S.zeros(5)
        for ij,value in zip(itertools.combinations(range(1,5),2),model['brackets']):
            affine = S.sympify(value,locals=env)
            P[ij] = S.expand(h*h*affine.subs({x:x/h for x in coordinates[1:]},simultaneous=True))
            P[ij[::-1]] = -P[ij]
        scale = S.ilcm(*[coefficient.q for ij in independent.PAIRS for coefficient in S.Poly(P[ij],*coordinates).coeffs()])
        P *= scale
        raw = sparse(P)
        assert all(not value for value in independent.jacobi(raw,raw=True).values())
        assert independent.jacobi(raw) == S.zeros(350,1)
        derivative_columns = []
        for ij,exponent in independent.COLUMNS:
            variation = {pair:dict(values) for pair,values in raw.items()}
            variation.setdefault(ij,{})[exponent] = variation.get(ij,{}).get(exponent,0)+1
            derivative_columns.append(independent.jacobi(variation))
        D = S.Matrix.hstack(*derivative_columns)
        vertical,orbit = [],[]
        euler = S.Matrix(coordinates)
        for i,j in itertools.product(range(5),repeat=2):
            A = S.zeros(5)
            A[i,j] = 1
            field = A*euler
            vertical.append(independent.vector(sparse(euler*field.T-field*euler.T)))
            action = A*P+P*A.T-coordinates[j]*P.diff(coordinates[i])
            orbit.append(independent.vector(sparse(action)))
        VM = S.Matrix.hstack(*vertical)
        group = S.Matrix.hstack(*vertical,*orbit,independent.vector(raw))
        assert D*group == S.zeros(350,group.cols)
        certificates = {name:certificate(M) for name,M in [('D',D),('vertical',VM),('group',group)]}
        assert certificates['vertical']['rank_QQ'] == 24
        kernel = rational_matrix(D).nullspace().to_Matrix().T
        assert D*kernel == S.zeros(350,kernel.cols)
        assert kernel.cols == 150-certificates['D']['rank_QQ']
        combined = group.row_join(kernel)
        pivots = list(rational_matrix(combined).rref()[1])
        normals = [column-group.cols for column in pivots if column >= group.cols]
        assert len(normals) == kernel.cols-certificates['group']['rank_QQ']

        # Move a declared rank-four affine point to e0, then Darboux frame it.
        basepoint = (0,0,1,1 if model['name']=='mixed' else 0)
        translation = S.eye(5)
        for i,value in enumerate(basepoint,1):
            translation[i,0] = -value
        old = P.subs(dict(zip(coordinates,X)),simultaneous=True)
        translated = pushforward(old,translation)
        evaluation = translated.subs(dict(zip(X,(1,0,0,0,0))))[1:5,1:5]
        assert evaluation.det() != 0
        frame = S.eye(5)
        frame[1:5,1:5] = darboux(evaluation)
        normalized = adapter.normalize(pushforward(translated,frame))
        adapter.verify_point(normalized,check_equations=True)
        point = normalized['parameters']
        symbols = S.symbols('u1:62')
        graphs = {}
        for name,path in [('T39',ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json'),
                          ('T40',ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json')]:
            graph = json.loads(path.read_text())
            indices = graph['weight2_u_indices_zero_based']
            residuals = [S.cancel(point[index]-S.sympify(value).subs(dict(zip(symbols,point)),simultaneous=True))
                         for index,value in zip(indices,graph['quadratic_section_q'])]
            assert len(residuals) == 45
            graphs[name] = dict(all45_residuals=list(map(str,residuals)),
                                all45_graph_equations_zero=all(value==0 for value in residuals),
                                nonzero_residual_count=sum(value!=0 for value in residuals))
        rank_D = certificates['D']['rank_QQ']
        rank_G = certificates['group']['rank_QQ']
        result = dict(name=model['name'],integer_scale=str(scale),
                      full_homogeneous_witness={f'{i},{j}':str(P[i,j]) for i,j in independent.PAIRS if P[i,j] != 0},
                      all150_ambient_variables=True,all350_projective_equations=True,
                      all_strict_and_projective_Jacobi_zero=True,actual_orbit_and_vertical_columns_annihilated=True,
                      certificates=certificates,full_kernel_dimension=kernel.cols,
                      projective_tangent_dimension=125-rank_D,projective_actual_orbit_dimension=rank_G-25,
                      actual_orbit_equals_full_tangent=(150-rank_D==rank_G),
                      kernel_sparse_columns=[{str(i):str(kernel[i,j]) for i in range(150) if kernel[i,j]}
                                             for j in range(kernel.cols)],
                      orbit_quotient_normal_kernel_columns=normals,
                      orbit_quotient_normal_directions=[{str(i):str(kernel[i,j]) for i in range(150) if kernel[i,j]}
                                                       for j in normals],
                      matrix_sha256={name:hashlib.sha256(json.dumps(M.tolist(),default=str).encode()).hexdigest()
                                     for name,M in [('D',D),('group',group)]},
                      rank4_basepoint=list(basepoint),
                      translation=[[str(v) for v in row] for row in translation.tolist()],
                      nondegenerate_evaluation=[[str(v) for v in row] for row in evaluation.tolist()],
                      Darboux_frame=[[str(v) for v in row] for row in frame.tolist()],
                      perspective=list(map(str,normalized['perspective'])),all61_parameters=list(map(str,point)),
                      full150_quotient_reconstruction_checked=True,all134_equations_zero=True,
                      graph_membership=graphs,
                      T01_excluded_by_irreducible_quartic_factor=True,
                      component_recognition_requires_geometric_bridge_and_independent_review=True,
                      seconds=time.monotonic()-started)
        results.append(result)
        print(json.dumps(dict(name=model['name'],D=rank_D,group=rank_G,
                              tangent=125-rank_D,orbit=rank_G-25,
                              graphs={name:data['all45_graph_equations_zero'] for name,data in graphs.items()})),flush=True)
    paths = [Path(__file__),source_path,Path(independent.__file__),
             ROOT/'code/classification_degenerate_quadric_q4_nil31_I_graph.py',ROOT/'code/s61_model.py',
             ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json',ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json',
             *(ROOT/'legacy/v02/data'/name for name in ['symplectic_slice_certificate.json',
                 'symplectic_slice_equations.json','symplectic_slice_reduced_equations.json'])]
    report = dict(task_id='E_reduced_nonnormal',all_passed=True,full_model_count=3,results=results,
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  source_bytes=Path(__file__).read_text(),
                  scope='Exact full ambient ranks and declared S61 graph tests at each actual orbit representative. A graph hit recognizes the point using the accepted graph theorem; transporting to the whole type10 uses the complete normal-form proof. No independent replay of these new ambient matrices is claimed.',
                  seconds=time.monotonic()-begun)
    with output.open('x') as handle:
        json.dump(report,handle,indent=2)
        handle.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))


def quadratic_audit(output):
    import sympy as S
    from sympy.polys.matrices import DomainMatrix
    import classification_rank3_mixed_component_review as arithmetic

    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun = time.monotonic()
    source_path = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_AMBIENT_20261107.json'
    source = json.loads(source_path.read_text())
    coordinates = S.symbols('h x y z w')
    env = {str(x):x for x in coordinates}
    dm = lambda A:DomainMatrix.from_Matrix(A).convert_to(S.QQ)

    def sparse(P):
        return {ij:dict(S.Poly(P[ij],*coordinates).terms()) for ij in arithmetic.PAIRS if P[ij] != 0}

    results = []
    for saved in reversed(source['results']):
        started = time.monotonic()
        P = S.zeros(5)
        for key,value in saved['full_homogeneous_witness'].items():
            i,j = map(int,key.split(','))
            P[i,j] = S.sympify(value,locals=env)
            P[j,i] = -P[i,j]
        raw = sparse(P)
        assert arithmetic.jacobi(raw) == S.zeros(350,1)
        cols = []
        for ij,exponent in arithmetic.COLUMNS:
            variation = {pair:dict(values) for pair,values in raw.items()}
            variation.setdefault(ij,{})[exponent] = variation.get(ij,{}).get(exponent,0)+1
            cols.append(arithmetic.jacobi(variation))
        D = S.Matrix.hstack(*cols)
        assert hashlib.sha256(json.dumps(D.tolist(),default=str).encode()).hexdigest() == saved['matrix_sha256']['D']
        euler = S.Matrix(coordinates)
        group_columns = []
        for i,j in itertools.product(range(5),repeat=2):
            A = S.zeros(5)
            A[i,j] = 1
            field = A*euler
            group_columns.append(arithmetic.vector(sparse(euler*field.T-field*euler.T)))
        for i,j in itertools.product(range(5),repeat=2):
            A = S.zeros(5)
            A[i,j] = 1
            group_columns.append(arithmetic.vector(sparse(A*P+P*A.T-coordinates[j]*P.diff(coordinates[i]))))
        group_columns.append(arithmetic.vector(raw))
        group = S.Matrix.hstack(*group_columns)
        assert hashlib.sha256(json.dumps(group.tolist(),default=str).encode()).hexdigest() == saved['matrix_sha256']['group']
        normal_vectors = []
        normal_tensors = []
        scales = []
        for encoding in saved['orbit_quotient_normal_directions']:
            vector = S.zeros(150,1)
            for key,value in encoding.items():
                vector[int(key)] = S.Rational(value)
            scale = S.ilcm(1,*[S.denom(value) for value in vector])
            vector *= scale
            tensor = S.zeros(5)
            for (ij,exponents),value in zip(arithmetic.COLUMNS,vector):
                if value:
                    tensor[ij] += value*S.prod(x**e for x,e in zip(coordinates,exponents))
            for i,j in arithmetic.PAIRS:
                tensor[j,i] = -tensor[i,j]
            assert arithmetic.vector(sparse(tensor)) == vector
            normal_vectors.append(vector)
            normal_tensors.append(tensor)
            scales.append(str(scale))
        normals = S.Matrix.hstack(*normal_vectors)
        assert D*group == S.zeros(350,group.cols)
        assert D*normals == S.zeros(350,normals.cols)
        assert len(dm(group).rref()[1]) == saved['certificates']['group']['rank_QQ'] == 47
        assert len(dm(group.row_join(normals)).rref()[1]) == saved['full_kernel_dimension']
        pairs = list(itertools.combinations_with_replacement(range(normals.cols),2))
        squares = [arithmetic.jacobi(sparse(R)) for R in normal_tensors]
        Qmatrix = S.Matrix.hstack(*[squares[i] if i==j else
                    arithmetic.jacobi(sparse(normal_tensors[i]+normal_tensors[j]))-squares[i]-squares[j]
                    for i,j in pairs])
        rows = saved['certificates']['D']['rows']
        pivots = saved['certificates']['D']['columns']
        independent = D[:,pivots]
        block = dm(D.extract(rows,pivots))
        assert str(block.det()) == saved['certificates']['D']['determinant']
        inverse_block = block.inv()
        correction = inverse_block.matmul(dm(Qmatrix.extract(rows,range(len(pairs))))).to_Matrix()
        obstruction = Qmatrix-independent*correction
        assert obstruction.extract(rows,range(len(pairs))) == S.zeros(len(rows),len(pairs))
        reduced,pivot_quadrics = dm(obstruction).rref()
        reduced = reduced.to_Matrix()
        qrank = len(pivot_quadrics)
        variables = S.symbols(f't0:{normals.cols}')
        monomials = [variables[i]*variables[j] for i,j in pairs]
        quadrics = [S.expand(sum(reduced[i,j]*monomials[j] for j in range(len(pairs)))) for i in range(qrank)]
        gb = S.groebner(quadrics,*variables,order='grevlex') if quadrics else None
        isolated = bool(gb and gb.is_zero_dimensional)
        selected_functional_rows = list(dm(obstruction.T).rref()[1])
        functionals = []
        for row in selected_functional_rows:
            weights = dm(independent[row,:]).matmul(inverse_block).to_Matrix()
            functional = S.zeros(1,350)
            functional[row] = 1
            for j,selected in enumerate(rows):
                functional[selected] -= weights[j]
            assert functional*D == S.zeros(1,150)
            coefficients = functional*Qmatrix
            assert coefficients == obstruction[row,:]
            functionals.append(dict(row=row,left_coefficients={str(j):str(v) for j,v in enumerate(functional) if v},
                                    quadratic_coefficients=list(map(str,list(coefficients)))))
        result = dict(name=saved['name'],full_D_and_actual_group_hashes_reproduced=True,
                      full_D_rank=len(pivots),actual_group_vertical_scale_rank=47,
                      genuine_normal_count=normals.cols,normal_integer_scales=scales,
                      normal_bivectors=[{f'{i},{j}':str(R[i,j]) for i,j in arithmetic.PAIRS if R[i,j] != 0}
                                         for R in normal_tensors],
                      full_actual_and_normal_derivatives_zero=True,
                      group_and_normals_span_entire_kernel=True,
                      normal_monomial_order=list(map(str,monomials)),
                      projected_quadratic_rank=qrank,quadratic_initial_generators=list(map(str,quadrics)),
                      Groebner_basis=[str(p.as_expr()) for p in gb.polys] if gb else [],
                      quadratic_support_zero_dimensional=isolated,
                      left_cokernel_functionals=functionals,
                      implicit_minor_rows=rows,implicit_minor_columns=pivots,
                      implicit_minor_determinant=str(block.det()),
                      obstruction_matrix_sha256=hashlib.sha256(json.dumps(obstruction.tolist(),default=str).encode()).hexdigest(),
                      formal_slice_dimension_bound_if_isolated=0 if isolated else None,
                      orbit_component_dimension_if_isolated=22 if isolated else None,
                      componenthood_not_promoted_without_independent_review=True,
                      higher_jets_or_exact_nil_index_claimed=False,
                      seconds=time.monotonic()-started)
        results.append(result)
        print(json.dumps(dict(name=result['name'],normals=normals.cols,quadrics=list(map(str,quadrics)),
                              isolated=isolated,seconds=result['seconds'])),flush=True)
    paths = [Path(__file__),source_path,Path(arithmetic.__file__)]
    report = dict(task_id='E_reduced_nonnormal',all_passed=True,results=results,
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  source_bytes=Path(__file__).read_text(),
                  scope='Genuine quadratic initial support on kerD modulo actual orbit, scale and vertical directions. Positive-dimensional quadratic cones are retained; no higher-jet iteration or component count inferred.',
                  seconds=time.monotonic()-begun)
    with output.open('x') as handle:
        json.dump(report,handle,indent=2)
        handle.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))


def weighted_mixed_slice(output):
    """Exact finite equivariant IFT slice; the bound4 precedes the solve."""
    import sympy as S
    from sympy.polys.matrices import DomainMatrix
    import classification_rank3_mixed_component_review as arithmetic

    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun = time.monotonic()
    source_path = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_AMBIENT_20261107.json'
    source = next(v for v in json.loads(source_path.read_text())['results'] if v['name']=='mixed')
    coordinates = S.symbols('h x y z w')
    env = {str(v):v for v in coordinates}
    P = S.zeros(5)
    for key,value in source['full_homogeneous_witness'].items():
        i,j = map(int,key.split(','))
        P[i,j] = S.sympify(value,locals=env)
        P[j,i] = -P[i,j]
    dm = lambda A:DomainMatrix.from_Matrix(A).convert_to(S.QQ)

    def sparse(A):
        return {ij:dict(S.Poly(A[ij],*coordinates).terms()) for ij in arithmetic.PAIRS if A[ij] != 0}

    def from_vector(v):
        A = S.zeros(5)
        for (ij,powers),coefficient in zip(arithmetic.COLUMNS,v):
            if coefficient:
                A[ij] += coefficient*S.prod(x**e for x,e in zip(coordinates,powers))
        for i,j in arithmetic.PAIRS:
            A[j,i] = -A[i,j]
        return A.applyfunc(S.expand)

    raw = sparse(P)
    assert arithmetic.jacobi(raw) == S.zeros(350,1)
    columns = []
    for ij,exponents in arithmetic.COLUMNS:
        variation = {pair:dict(values) for pair,values in raw.items()}
        variation.setdefault(ij,{})[exponents] = variation.get(ij,{}).get(exponents,0)+1
        columns.append(arithmetic.jacobi(variation))
    D = S.Matrix.hstack(*columns)
    assert hashlib.sha256(json.dumps(D.tolist(),default=str).encode()).hexdigest() == source['matrix_sha256']['D']
    euler = S.Matrix(coordinates)
    actual = []
    for i,j in itertools.product(range(5),repeat=2):
        A = S.zeros(5)
        A[i,j] = 1
        field = A*euler
        actual.append(arithmetic.vector(sparse(euler*field.T-field*euler.T)))
    for i,j in itertools.product(range(5),repeat=2):
        A = S.zeros(5)
        A[i,j] = 1
        actual.append(arithmetic.vector(sparse(A*P+P*A.T-coordinates[j]*P.diff(coordinates[i]))))
    actual.append(arithmetic.vector(raw))
    group = S.Matrix.hstack(*actual)
    assert hashlib.sha256(json.dumps(group.tolist(),default=str).encode()).hexdigest() == source['matrix_sha256']['group']
    assert D*group == S.zeros(350,group.cols)
    weights = (0,2,1,1,2)
    character = 1
    relative = [weights[i]+weights[j]-sum(e*w for e,w in zip(powers,weights))-character
                for (i,j),powers in arithmetic.COLUMNS]
    row_weights = [sum(weights[i] for i in quad)-sum(e*w for e,w in zip(powers,weights))-2*character
                   for quad,powers in arithmetic.ROWS]
    pvector = arithmetic.vector(raw)
    assert {relative[i] for i,v in enumerate(pvector) if v} == {0}
    assert (min(relative),max(relative)) == (-4,3)
    assert all(len({relative[i] for i,value in enumerate(group[:,j]) if value})<=1 for j in range(group.cols))
    assert all(not D[i,j] or row_weights[i]==relative[j] for i in range(350) for j in range(150))

    group_columns = list(dm(group).rref()[1])
    assert len(group_columns)==47
    group_basis = group[:,group_columns]
    group_rows = list(dm(group_basis.T).rref()[1])
    assert len(group_rows)==47
    group_minor = dm(group_basis[group_rows,:]).det()
    assert group_minor != 0
    slice_indices = [i for i in range(150) if i not in group_rows]
    linear_slice = D[:,slice_indices]
    rr,pivots = dm(linear_slice).rref()
    rr = rr.to_Matrix()
    assert len(pivots)==101
    implicit = [slice_indices[i] for i in pivots]
    free_positions = [i for i in range(len(slice_indices)) if i not in pivots]
    assert len(free_positions)==2
    normal_data = []
    for free in free_positions:
        vector = S.zeros(150,1)
        vector[slice_indices[free]] = 1
        for row,pivot in enumerate(pivots):
            vector[slice_indices[pivot]] = -rr[row,free]
        assert D*vector == S.zeros(350,1)
        active_weights = {relative[i] for i,value in enumerate(vector) if value}
        assert len(active_weights)==1 and next(iter(active_weights))<0
        normal_data.append((-next(iter(active_weights)),slice_indices[free],vector))
    normal_data.sort(key=lambda value:value[0])
    assert [value[0] for value in normal_data]==[1,2]
    free_coordinates = [value[1] for value in normal_data]
    normal_vectors = S.Matrix.hstack(*(value[2] for value in normal_data))
    assert len(dm(group.row_join(normal_vectors)).rref()[1])==49
    rows = list(dm(D[:,implicit].T).rref()[1])
    block = dm(D.extract(rows,implicit))
    determinant = block.det()
    assert len(rows)==101 and determinant != 0
    inverse_block = block.inv()
    # Both the group cut and the selected equations are coordinate-weight
    # subspaces of the entire C150/350 representation, not just of kerD.
    assert sorted(group_rows+slice_indices)==list(range(150))
    assert all(not D[row,col] or row_weights[row]==relative[col] for row in rows for col in implicit)
    bound = max([0]+[-relative[i] for i in implicit])
    assert bound==4
    u,v = S.symbols('u v')
    parameters = (u,v)
    parameter_weights = (1,2)
    curve = P+from_vector(u*normal_vectors[:,0]+v*normal_vectors[:,1])

    def residual(A):
        return arithmetic.jacobi(sparse(A)).applyfunc(S.expand)

    def part(expression,degree):
        return S.expand(sum(coefficient*u**powers[0]*v**powers[1]
                            for powers,coefficient in S.Poly(expression,u,v).terms()
                            if sum(a*b for a,b in zip(powers,parameter_weights))==degree))

    def check_equivariance(values):
        for row,expression in enumerate(values):
            for powers,coefficient in S.Poly(expression,u,v).terms():
                assert not coefficient or sum(a*b for a,b in zip(powers,parameter_weights)) == -row_weights[row]

    corrections = []
    for degree in range(2,bound+1):
        values = residual(curve)
        check_equivariance(values)
        selected = S.Matrix([part(values[i],degree) for i in rows])
        correction = -inverse_block.matmul(dm(selected)).to_Matrix() if not selected.free_symbols else None
        if correction is None:
            correction = -(inverse_block.to_Matrix()*selected).applyfunc(S.expand)
        full_correction = S.zeros(150,1)
        for index,value in zip(implicit,correction):
            value = S.expand(value)
            if value:
                assert relative[index]==-degree
                assert all(sum(powers)>=2 for powers,coefficient in S.Poly(value,u,v).terms() if coefficient)
                full_correction[index] = value
        curve = (curve+from_vector(full_correction)).applyfunc(S.expand)
        after = residual(curve)
        assert all(part(after[i],degree)==0 for i in rows)
        corrections.append(dict(weighted_degree=degree,only_implicit_coordinates=True,
                                coefficients={str(i):str(value) for i,value in enumerate(full_correction) if value}))
    final = residual(curve)
    check_equivariance(final)
    assert final.extract(rows,[0]) == S.zeros(101,1)
    displacement = arithmetic.vector(sparse(curve-P))
    assert all(displacement[i]==0 for i in group_rows)
    assert [displacement[i] for i in free_coordinates]==[u,v]
    distinct = []
    for value in final:
        if value:
            normalized = S.Poly(value,u,v).monic().as_expr()
            if normalized not in distinct:
                distinct.append(normalized)
    gb = S.groebner(distinct,u,v,order='grevlex') if distinct else None
    isolated = bool(gb and gb.is_zero_dimensional)
    t = S.symbols('t')
    chart_polynomials = []
    for value in distinct:
        polynomial = S.Poly(value,u,v)
        degree = sum(a*b for a,b in zip(polynomial.monoms()[0],parameter_weights))
        substituted = S.cancel(value.subs(v,t*u*u)/u**degree)
        assert not substituted.has(u)
        chart_polynomials.append(S.Poly(substituted,t))
    chart_gcd = chart_polynomials[0] if chart_polynomials else S.Poly(0,t)
    for polynomial in chart_polynomials[1:]:
        chart_gcd = S.gcd(chart_gcd,polynomial)
    actual_axis = all(S.expand(value.subs(v,0))==0 for value in final)
    family = curve.subs(v,0) if actual_axis else None
    if actual_axis:
        assert residual(family)==S.zeros(350,1)
    paths = [Path(__file__),source_path,Path(arithmetic.__file__)]
    report = dict(task_id='E_reduced_nonnormal',all_passed=True,
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  source_bytes=Path(__file__).read_text(),
                  coordinate_weights=list(weights),P_character=character,
                  all150_coordinate_relative_weights=relative,all350_equation_relative_weights=row_weights,
                  entire_C150_weight_action_verified=True,actual_group_weight_stable=True,
                  actual_group_independent_columns=group_columns,actual_group_cut_rows=group_rows,
                  actual_group_cut_minor_determinant=str(group_minor),
                  weight_stable_slice_coordinate_indices=slice_indices,
                  implicit_coordinate_indices=implicit,selected_equation_rows=rows,
                  implicit_minor_determinant=str(determinant),free_normal_coordinate_indices=free_coordinates,
                  normal_relative_weights=[-1,-2],
                  normal_vectors=[{str(i):str(value) for i,value in enumerate(normal_vectors[:,j]) if value}
                                  for j in range(2)],
                  finite_weight_bound_proved_before_solving=bound,
                  finite_corrections=corrections,all_selected_equations_identically_zero=True,
                  full_polynomial_slice_brackets={f'{i},{j}':str(curve[i,j]) for i,j in arithmetic.PAIRS if curve[i,j] != 0},
                  all350_literal_residuals={str(i):str(value) for i,value in enumerate(final) if value},
                  complete_formal_slice_ideal=list(map(str,distinct)),
                  Groebner_basis=[str(p.as_expr()) for p in gb.polys] if gb else [],
                  isolated_support=isolated,weighted_chart_v_over_u_squared_gcd=str(chart_gcd.as_expr()),
                  v_zero_is_actual_Poisson_family=actual_axis,
                  actual_v_zero_family_brackets={f'{i},{j}':str(family[i,j]) for i,j in arithmetic.PAIRS if family[i,j] != 0} if actual_axis else None,
                  mathematical_dependency='The full C150 action, genuine actual-group cut and selected invertible derivative minor make the formal IFT map C*-equivariant. Normal weights -1,-2 force every implicit coordinate series to be a polynomial of weighted degree<=4; the displayed polynomial solves all selected equations identically. Therefore the full350 residuals are the actual complete formal slice equations, not a finite-jet guess.',
                  scope='Exact finite genuine formal slice at the mixed type10 orbit. Any componenthood/known-component recognition requires independent replay and the stated formal-slice argument. No unbounded jet iteration or exact nil index claim.',
                  seconds=time.monotonic()-begun)
    with output.open('x') as handle:
        json.dump(report,handle,indent=2)
        handle.write('\n')
    print(json.dumps(dict(all_passed=True,normal_weights=[-1,-2],bound=bound,
                          Groebner_basis=report['Groebner_basis'],isolated=isolated,
                          actual_axis=actual_axis,seconds=report['seconds'])))


def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(170,175))
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    start = time.monotonic()
    input_path = ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'
    source = json.loads(input_path.read_text())['full_nilpotent_sources'][9]
    assert source['variant'] == dict(jordan_blocks=[2,1,1],translation_top_length=2)
    bridge = verify_complete_source(source)
    reduction = verify_jacobi_and_square_completion()
    table = weighted_table()
    models = actual_models()
    inputs = [Path(__file__),Path(__file__).with_name('classification_minimal_quartic_review.py'),
              Path(__file__).with_name('classification_quadratic_cofactor_review.py'),input_path]
    result = dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_import=True,no_CAS_import=True,
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                  complete_original_source_bridge=bridge,full_Jacobi_reduction=reduction,
                  complete_weighted_table=table,three_actual_Poisson_models=models,
                  normalizing_coordinate_changes='r=x-y^2/2 is auxiliary only; r->r+L/delta is the actual affine x->x+L/delta. Base diagonalization/translations are affine. Final z->tz,w->uw gives delta_new=delta*t*u, b_new=b*t^(1-i)*u^(1-j).',
                  scale_normalization={'H=bz':['t=b/delta','u=1/b'],
                                       'H=bz^3':['choose t^3=b/delta','u=1/(delta*t)'],
                                       'H=bz^2w':['t=b','u=1/(delta*b)']},
                  inverted_elements=['delta!=0 forced by simpleH','nonzero H coefficientb after excluding reducibleF','nonzero determinant base affine coordinate transforms'],
                  complements=['delta=0 forces highestPfaffian0, hence H is not simple','repeated eigenvalues -1/4 have no degree<=3 eigenpolynomial','zero-eigen translation with kappa!=0 forces b=0 and a reducibleF','H coefficientb=0 makesF reducible'],
                  scope='All type10 simpleH/irreducible-quartic points reduce to three actual affine/PGL orbits. Ambient componenthood, known-component membership, pairwise/global distinctness remain open.',
                  seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  source_bytes=Path(__file__).read_text())
    with output.open('x') as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    print(json.dumps({key:result[key] for key in ('all_passed','seconds','maximum_RSS_KiB')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--ambient',action='store_true',help='Full150/350 ambient ranks and faithful S61 graph tests')
    parser.add_argument('--quadratic',action='store_true',help='Genuine quadratic initial support on all three actual normal quotients')
    parser.add_argument('--weighted-mixed',action='store_true',help='Exact weight-bounded genuine two-normal formal slice')
    args = parser.parse_args()
    assert sum([args.ambient,args.quadratic,args.weighted_mixed])<=1
    (weighted_mixed_slice if args.weighted_mixed else quadratic_audit if args.quadratic else ambient_audit if args.ambient else main)(args.output)
