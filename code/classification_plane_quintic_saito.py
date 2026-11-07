#!/usr/bin/env python3
"""Finite complete quadratic Saito kernels and all-parameter cotangent tensors."""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import resource
import time
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def candidate(record):
    x, y, z, u, v = variables = s.symbols('x y z u v')
    base = variables[:3]
    F = s.sympify(record['polynomial'], locals=dict(zip(map(str, base), base)))
    assert s.Poly(F, *base).total_degree() == 5
    assert all(sum(e) == 5 for e in s.Poly(F, *base).monoms())
    monomials = [x*x, x*y, x*z, y*y, y*z, z*z]
    target = sorted((e for e in itertools.product(range(7), repeat=3) if sum(e) == 6), reverse=True)
    gradient = [s.diff(F, a) for a in base]
    columns = [s.Poly(m*gradient[i], *base) for i in range(3) for m in monomials]
    matrix = s.Matrix([[p.coeff_monomial(e) for p in columns] for e in target])
    kernel = matrix.nullspace()
    assert len(kernel) == 2
    vectors = [s.Matrix([sum(k[6*i+j]*monomials[j] for j in range(6)) for i in range(3)]) for k in kernel]
    U0, V0 = vectors
    assert s.expand((s.Matrix(gradient).T*U0)[0]) == 0
    assert s.expand((s.Matrix(gradient).T*V0)[0]) == 0
    determinant = s.expand(s.Matrix.hstack(s.Matrix(base), U0, V0).det())
    coefficient = s.cancel(determinant/F)
    assert coefficient.is_Rational and coefficient != 0
    assert s.expand(determinant-coefficient*F) == 0

    # Complex absolute irreducibility, rather than factorization over QQ only.
    assert s.degree(F, z) == 1
    leading = s.expand(s.diff(F, z))
    constant = s.expand(F.subs(z, 0))
    assert leading == y**4 and constant.subs(y, 0) == x**5

    g00, g01, g10, g11 = frame = s.symbols('g00 g01 g10 g11')
    magnetic = s.symbols('r0:6')
    U = g00*U0+g01*V0
    V = g10*U0+g11*V0
    divU = sum(s.diff(U[i], base[i]) for i in range(3))
    divV = sum(s.diff(V[i], base[i]) for i in range(3))
    R = sum(c*m for c, m in zip(magnetic, monomials))
    B = s.expand(-divV*u+divU*v+R)
    P = s.zeros(5)
    for i in range(3):
        P[i, 3], P[3, i] = s.expand(U[i]), -s.expand(U[i])
        P[i, 4], P[4, i] = s.expand(V[i]), -s.expand(V[i])
    P[3, 4], P[4, 3] = B, -B
    jacobi = []
    for i, j, k in itertools.combinations(range(5), 3):
        expression = s.expand(sum(P[i, l]*s.diff(P[j, k], variables[l])
                                  + P[j, l]*s.diff(P[k, i], variables[l])
                                  + P[k, l]*s.diff(P[i, j], variables[l]) for l in range(5)))
        assert expression == 0
        jacobi.append(dict(triple=[i, j, k], polynomial=str(expression)))
    divergence = [s.expand(sum(s.diff(P[i, j], variables[j]) for j in range(5))) for i in range(5)]
    assert divergence == [0]*5
    casimir = [s.expand(sum(P[i, j]*s.diff(F, variables[j]) for j in range(5))) for i in range(5)]
    assert casimir == [0]*5
    pfaffians = []
    for missing in range(5):
        i, j, k, l = [a for a in range(5) if a != missing]
        pfaffians.append(s.expand(P[i, j]*P[k, l]-P[i, k]*P[j, l]+P[i, l]*P[j, k]))
    actualD = s.expand(sum((-1)**i*variables[i]*pfaffians[i] for i in range(5)))
    frame_det = g00*g11-g01*g10
    assert s.expand(actualD+coefficient*frame_det*F) == 0
    return dict(id=record['id'], polynomial=str(F), literature=record['literature'],
                coefficient_ring='QQ[g00,g01,g10,g11,r0,...,r5,x,y,z,u,v]',
                complete_quadratic_gradient_kernel=dict(
                    monomials=list(map(str, monomials)), target_exponents=target,
                    all_28_by_18_rows=[[str(a) for a in matrix.row(i)] for i in range(matrix.rows)],
                    rank=18-len(kernel), dimension=len(kernel),
                    full_18coefficient_basis=[[str(a) for a in k] for k in kernel],
                    vector_fields=[[str(a) for a in k] for k in vectors],
                    saito_determinant=str(determinant), determinant_over_F=str(coefficient)),
                complex_absolute_irreducibility=dict(
                    criterion='Primitive linear polynomial in z over C[x,y], irreducible over C(x,y)[z]; Gauss lemma.',
                    leading_coefficient=str(leading), constant_at_y0=str(constant.subs(y, 0)),
                    only_possible_common_prime='y', coprime_over_C=True),
                whole_cotangent_source=dict(
                    variables=list(map(str, variables)), full_frame=list(map(str, frame)),
                    frame_unit=str(frame_det), magnetic_parameters=list(map(str, magnetic)),
                    all_ten_brackets={f'{i},{j}':str(P[i, j]) for i, j in itertools.combinations(range(5), 2)},
                    all_ten_raw_strict_Jacobi=jacobi,
                    all_five_divergences=list(map(str, divergence)),
                    all_five_Casimir_contractions=list(map(str, casimir)),
                    all_five_raw_pfaffians=list(map(str, pfaffians)), actual_quintic=str(actualD),
                    exact_quintic_factor=str(-coefficient*frame_det),
                    complements='det(frame)=0 is excluded only from the rank-four fixed-F locus. All six magnetic coefficients, including zero, are retained.'),
                scope='One explicitly defined irreducible plane quintic has a complete Saito(2,2) kernel and full actual cotangent source. No exhaustive curve list, componenthood, ambient tangent/cut, or new component count.')


def ambient(input_path):
    import importlib.util
    from sympy.polys.matrices import DomainMatrix
    def domain(A):
        return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    def pivots(A):
        return list(domain(A).rref()[1])
    def sparse(v):
        return {str(i):str(a) for i,a in enumerate(v) if a}
    def minor(A):
        cc=pivots(A)
        rr=pivots(A[:,cc].T)
        return dict(rank=len(cc),rows=rr,columns=cc,
                    determinant=str(domain(A[rr,cc]).det()))
    helper=ROOT/'legacy/v03/code/poisson_p4_components_check.py'
    spec=importlib.util.spec_from_file_location('essential_three_ambient',helper)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    data=json.loads(input_path.read_text())
    assert data['all_passed'] and len(data['candidates'])==1
    row=data['candidates'][0]
    x,y,z,u,v=variables=s.symbols('x y z u v')
    env=dict(zip(map(str,variables),variables))
    env.update(dict(zip(('g00','g01','g10','g11'),(1,0,0,1))))
    env.update({f'r{i}':0 for i in range(6)})
    P=s.zeros(5)
    for pair,polynomial in row['whole_cotangent_source']['all_ten_brackets'].items():
        i,j=map(int,pair.split(','))
        P[i,j]=s.sympify(polynomial,locals=env)
        P[j,i]=-P[i,j]
    scale=s.ilcm(1,*[a.q for p in P for a in s.Poly(p,*variables).coeffs()])
    P*=scale
    bv={(i,j):{e:int(c) for e,c in s.Poly(P[i,j],*variables).terms() if c}
        for i,j in m.PAIRS if P[i,j]!=0}
    D=s.Matrix(m.linearized_jacobi_matrix(bv).tolist())
    assert D.shape==(350,150)
    def vector(A):
        return s.Matrix([s.Poly(A[i,j],*variables).coeff_monomial(e)
                         for (i,j),e in m.BIVECTOR_BASIS])
    directions=[]
    for i,j in itertools.product(range(5),repeat=2):
        A=s.zeros(5)
        A[i,j]=1
        X=A*s.Matrix(variables)
        directions.append(vector(sum((X[k]*P.diff(variables[k]) for k in range(5)),s.zeros(5))-A*P-P*A.T))
    vertical=s.Matrix(m.vertical_bivector_matrix().tolist())
    M=vertical.row_join(s.Matrix.hstack(*directions,vector(P)))
    assert vertical.shape==(150,25) and M.shape==(150,51)
    assert D*M==s.zeros(350,51)
    dcert,mcert=minor(D),minor(M)
    assert mcert['rank']==48
    sigma=(0,0,0,1,1)
    bw=[sigma[i]+sigma[j]-sum(a*b for a,b in zip(e,sigma))-1
        for (i,j),e in m.BIVECTOR_BASIS]
    ow=[sum(sigma[i] for i in quad)-sum(a*b for a,b in zip(e,sigma))-2
        for quad,e in m.FOUR_VECTOR_BASIS]
    assert all(not c or bw[j]==ow[i] for (i,j),c in D.todok().items())
    pointvector=vector(P)
    assert all(not c or bw[i]==0 for i,c in enumerate(pointvector))
    blocks=[]
    sourcebasis=[]
    normals=[]
    cutrows=[]
    for weight in sorted(set(bw)):
        cc=[i for i,w in enumerate(bw) if w==weight]
        rr=[i for i,w in enumerate(ow) if w==weight]
        blockD=D[rr,cc]
        kernel=blockD.nullspace()
        blockM=M[cc,:]
        selected=pivots(blockM)
        real=blockM[:,selected]
        assert blockD*real==s.zeros(len(rr),len(selected))
        assert len(kernel)>=len(selected)
        pivotrows=pivots(real.T)
        assert len(pivotrows)==len(selected)
        cutrows.extend(cc[i] for i in pivotrows)
        for k in range(len(selected)):
            col=s.zeros(150,1)
            for i,at in enumerate(cc):
                col[at]=real[i,k]
            sourcebasis.append(col)
        # The same coordinate cut gives a C*-stable genuine normal complement.
        available=[i for i in range(len(cc)) if i not in pivotrows]
        cutkernel=blockD[:,available].nullspace()
        assert len(cutkernel)==len(kernel)-len(selected)
        localnormals=[]
        for col0 in cutkernel:
            col=s.zeros(150,1)
            for i,at in enumerate(available):
                col[cc[at]]=col0[i]
            assert D*col==s.zeros(350,1)
            normals.append((weight,col))
            localnormals.append(sparse(col))
        blocks.append(dict(weight=weight,input_columns=cc,output_rows=rr,
                           full_D_block_rank=blockD.cols-len(kernel),kernel_dimension=len(kernel),
                           actual_action_rank=len(selected),normal_dimension=len(cutkernel),
                           actual_source_column_indices=selected,
                           source_cut_coordinate_rows=[cc[i] for i in pivotrows],
                           full150_true_normals=localnormals))
    F=s.Matrix.hstack(*sourcebasis)
    assert F.shape==(150,48) and len(cutrows)==48
    assert sum(b['actual_action_rank'] for b in blocks)==48
    assert len(normals)==150-dcert['rank']-48
    assert len(pivots(F.row_join(M)))==48
    assert domain(F[cutrows,:]).det()!=0
    remaining=[i for i in range(150) if i not in cutrows]
    Ucols=pivots(D[:,remaining])
    Uindices=[remaining[i] for i in Ucols]
    basis=s.eye(150)
    N=s.Matrix.hstack(*(c for _,c in normals)) if normals else s.zeros(150,0)
    splitting=F.row_join(N).row_join(basis[:,Uindices])
    assert splitting.shape==(150,150)
    splitting_det=domain(splitting).det()
    assert splitting_det!=0
    if dcert['rank']==102:
        implication='Full tangent rawdimension48 equals actual orbit/vertical/scalar source48. The point is smooth with unique projective23 component; the complete fixed-F orbit closure is that component. Prior-component distinctness not tested here.'
    elif dcert['rank']==101:
        implication='One true normal remains. Its relative character is saved for a separately gated weakest finite obstruction/model job; no componenthood conclusion yet.'
    else:
        implication='All true normals/weights saved. Design the weakest support statement before any new quadratic/cubic computation; no automatic jet continuation or component count.'
    return dict(candidate_id=row['id'],input_report=dict(path=str(input_path),sha256=hashlib.sha256(input_path.read_bytes()).hexdigest()),
                legacy_helper=dict(path=str(helper.relative_to(ROOT)),sha256=hashlib.sha256(helper.read_bytes()).hexdigest(),
                                   tracked_archive='sources/v03.zip',archive_member='code/poisson_p4_components_check.py'),
                canonical_point_integer_scale=str(scale),
                canonical_point={f'{i},{j}':str(P[i,j]) for i,j in m.PAIRS if P[i,j]!=0},
                full_350_by_150_D=dcert,full_D_sparse_rows=[sparse(D.row(i)) for i in range(350)],
                actual_Euler_GL_scalar_source=mcert,
                all51_actual_direction_columns=[sparse(M[:,i]) for i in range(51)],
                all_D_actual_action_products_zero=True,
                Cstar=dict(coordinate_weights=sigma,tensor_character=1,relative_input_weights=bw,
                           relative_output_weights=ow,full_D_weight_intertwining=True,
                           whole_cotangent_source_justifies_actual_stabilizer=True),
                full_weight_blocks=blocks,true_cut=dict(dimension=102,coordinate_rows=cutrows,
                    actual48_basis=[sparse(F[:,i]) for i in range(48)],
                    source_cut48_determinant=str(domain(F[cutrows,:]).det()),
                    independent_implicit_coordinate_indices=Uindices,
                    normal_weights=[w for w,_ in normals],all150_true_normals=[sparse(c) for _,c in normals],
                    full150_direct_sum_determinant=str(splitting_det)),
                projective_actual_orbit_dimension=23,projective_tangent_dimension=150-dcert['rank']-25,
                implication=implication,no_quadratic_cubic_or_higher_obstruction_run=True)


def component_support(input_path):
    import importlib.util
    from sympy.polys.matrices import DomainMatrix
    def domain(A):
        return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    def sparse(v):
        return {str(i):str(a) for i,a in enumerate(v) if a}
    def decode(data,n):
        v=s.zeros(n,1)
        for i,c in data.items():v[int(i)]=s.Rational(c)
        return v
    helper=ROOT/'legacy/v03/code/poisson_p4_components_check.py'
    spec=importlib.util.spec_from_file_location('essential_three_support',helper)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    old=json.loads(input_path.read_text())
    assert old['all_passed']
    row=old['full_ambient']
    assert row['full_350_by_150_D']['rank']==101
    assert row['actual_Euler_GL_scalar_source']['rank']==48
    assert row['true_cut']['normal_weights']==[-1]
    variables=s.symbols('x y z u v')
    D=s.Matrix.vstack(*[decode(r,150).T for r in row['full_D_sparse_rows']])
    R=decode(row['true_cut']['all150_true_normals'][0],150)
    def tensor(vec):
        result={}
        for c,(pair,e) in zip(vec,m.BIVECTOR_BASIS):
            if c:result.setdefault(pair,{})[e]=c
        return result
    def jac(vec):
        J=m.affine_jacobi_coefficients(tensor(vec))
        result=s.zeros(350,1)
        for triple,poly in J.items():
            for i in range(5):
                sign=m.permutation_sign((i,)+triple)
                if not sign:continue
                quad=tuple(sorted((i,)+triple))
                for e,c in poly.items():
                    ex=list(e);ex[i]+=1
                    result[m.FOUR_VECTOR_INDEX[(quad,tuple(ex))]]+=sign*c
        return result
    P=s.zeros(150,1)
    for pair,polynomial in row['canonical_point'].items():
        ij=tuple(map(int,pair.split(',')))
        for e,c in s.Poly(s.sympify(polynomial,locals=dict(zip(map(str,variables),variables))),*variables).terms():
            P[m.BIVECTOR_INDEX[(ij,e)]]=c
    assert jac(P)==jac(P+R)-jac(R)==s.zeros(350,1)
    dc=row['true_cut']['independent_implicit_coordinate_indices']
    assert len(dc)==101
    dr=list(domain(D[:,dc].T).rref()[1])
    assert len(dr)==101
    A=D[dr,dc]
    inv=domain(A).inv().to_Matrix()
    assert A*inv==s.eye(101)
    bw=row['Cstar']['relative_input_weights']
    ow=row['Cstar']['relative_output_weights']
    assert min(bw)==-3 and all(not c or bw[i]==-1 for i,c in enumerate(R))
    bound=3
    coefficients={0:P,1:R}
    cache={}
    def bilinear(a,b):
        key=tuple(sorted((a,b)))
        if key not in cache:
            cache[key]=jac(coefficients[a]) if a==b else jac(coefficients[a]+coefficients[b])-jac(coefficients[a])-jac(coefficients[b])
        return cache[key]
    def quadratic(degree):
        result=s.zeros(350,1)
        for a,b in itertools.combinations_with_replacement(sorted(k for k in coefficients if k),2):
            if a+b==degree:result+=bilinear(a,b)
        return result
    stages=[]
    certificate=None
    for degree in range(2,bound+1):
        Q=quadratic(degree)
        T=s.zeros(150,1)
        solution=-inv*Q[dr,:]
        for i,c in zip(dc,solution):T[i]=c
        assert all(not c or bw[i]==-degree for i,c in enumerate(T))
        residual=D*T+Q
        assert residual[dr,:]==s.zeros(101,1)
        assert all(not c or ow[i]==-degree for i,c in enumerate(residual))
        stages.append(dict(degree=degree,raw_quadratic_source=sparse(Q),
                           unique_selected_IFT_correction=sparse(T),full350_residual=sparse(residual)))
        if any(residual):
            at=next(i for i,c in enumerate(residual) if c)
            lam=s.zeros(1,350)
            lam[at]=1
            multipliers=D[at,dc]*inv
            for i,c in zip(dr,multipliers):lam[i]-=c
            assert lam*D==s.zeros(1,150)
            assert (lam*Q)[0]==residual[at]!=0
            certificate=dict(initial_degree=degree,residual_row=at,cokernel_functional=sparse(lam),
                             all150_D_columns_annihilated=True,nonzero_initial_coefficient=str(residual[at]))
            break
        coefficients[degree]=T
    full_residual={}
    if certificate is None:
        for degree in range(1,2*bound+1):
            Q=quadratic(degree)+(D*coefficients[degree] if degree in coefficients and degree else s.zeros(350,1))
            assert Q[dr,:]==s.zeros(101,1)
            assert all(not c or ow[i]==-degree for i,c in enumerate(Q))
            if any(Q):full_residual[str(degree)]=sparse(Q)
        if full_residual:
            first=min(map(int,full_residual))
            certificate=dict(initial_degree=first,full350_polynomial_residuals=full_residual,
                             implication='One-variable exact finite IFT residual has a nonzero initial coefficient; its reduced local support is the origin.')
    is_component=certificate is not None
    return dict(candidate_id=row['candidate_id'],input_report=dict(path=str(input_path),sha256=hashlib.sha256(input_path.read_bytes()).hexdigest()),
                reused_full_D_rank=101,actual_full_source_raw_rank=48,
                true150_direct_sum_determinant=row['true_cut']['full150_direct_sum_determinant'],
                same_true48_coordinate_cut=row['true_cut']['coordinate_rows'],
                implicit_columns=dc,implicit_rows=dr,implicit101minor=str(domain(A).det()),
                inverse_101_by_101=[[str(c) for c in inv.row(i)] for i in range(101)],
                coordinate_weights=row['Cstar']['coordinate_weights'],tensor_character=1,
                relative_normal_character=-1,prior_maximum_implicit_degree=bound,
                finite_bound_reason='The genuine one-normal IFT is C*-equivariant. Every implicit coordinate has relative character between -3 and 1; a power n^k has character -k. Hence all implicit series have degree at most3, before any coefficient calculation.',
                stages=stages,weakest_support_certificate=certificate,
                finite_complete_graph_computed=certificate is None or bool(full_residual),
                graph_coefficients={str(k):sparse(c) for k,c in coefficients.items() if k},
                full350_graph_residuals=full_residual,
                projective_orbit_closure_is_full_component=is_component,
                componenthood_reason='A nonzero one-variable true IFT residual makes its reduced support zero-dimensional. Krull gives raw local dimension<=48, while the full fixed-F orbit source has rawdimension48. Thus its projective23 orbitclosure is a full component. Prior distinctness is separate.' if is_component else 'All350 equations vanish on the complete polynomial IFT graph, giving a genuine extra normal curve and local projective dimension24; no new23 component.',
                initial_quadratic_was_sufficient=certificate is not None and certificate['initial_degree']==2,
                no_new_point_or_unbounded_jet_search=True)


def prior103_audit(input_path):
    snapshot=json.loads(input_path.read_text())
    rows=snapshot['all103_comparison_rows']
    assert len(rows)==103 and snapshot['immutable_catalogue']['count']==103
    xx=s.symbols('x0:5')
    def quintic(P):
        out=0
        for missing in range(5):
            i,j,k,l=[a for a in range(5) if a!=missing]
            out+=(-1)**missing*xx[missing]*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])
        return s.expand(out)
    def tensor(data,variables):
        P=s.zeros(5)
        env=dict(zip(map(str,variables),xx))
        for pair,text in data.items():
            i,j=map(int,pair.split(','));P[i,j]=s.sympify(text,locals=env);P[j,i]=-P[i,j]
        return P
    mps=[]
    paths=[]
    for row in rows:
        name=row['id']
        if not(name.startswith('T') and name[1:].isdigit() and row['projective_dimension']==23 and row['generic_rank']==4):continue
        path=ROOT/f'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/{name}.json'
        if not path.exists():path=ROOT/f'reports/MPS_S61_RAW_CATALOGUE_20261103_v62/{name}.json'
        data=json.loads(path.read_text())
        assert data.get('exact_projective_Jacobi_checked',data.get('projective_Jacobi_identically_zero',False))
        P=tensor(data['homogeneous_bivector'],xx)
        edges=set().union(*(p.free_symbols for p in P))-set(xx)
        edges={p for p in edges if str(p).startswith('t')}
        P=P.subs(dict.fromkeys(edges,0))
        D=quintic(P)
        kappa=s.factor(s.cancel(D/s.prod(xx)))
        assert kappa!=0 and not(set(xx)&kappa.free_symbols)
        assert s.expand(D-kappa*s.prod(xx))==0
        mps.append(dict(id=name,accepted_source_report=str(path.relative_to(ROOT)),
                        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        all_edges_zero=True,actual_SNC_quintic=str(D),kappa=str(kappa),
                        full_period_parameters=sorted(map(str,kappa.free_symbols)),
                        nonzero_gradient_five_minor=str(s.factor(kappa**5)),
                        minor_rows=[str(s.prod(xx[j] for j in range(5) if j!=i)) for i in range(5)],
                        essential_variable_rank=5,
                        membership='Edge-zero point is in the accepted full completed source; any nonzero period-chart point with kappa!=0 is an actual member, not merely a formal leading term.'))
        paths.append(path)
    assert len(mps)==29
    irred=[]
    p=ROOT/'reports/CLASSIFICATION_GA9_DISTINCT_COMPONENT_AUDIT_20261107.json'
    a=json.loads(p.read_text())
    p2=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_CURVE_RECOGNITION_20261107.json'
    b=json.loads(p2.read_text())['nil22_long_quartic_curve_recognition']
    r=next(row for row in b['all_five_actual_curve_recognitions'] if row['anchor']==5)
    specifications=[('R4_GA9_DISTINCT_IRREDUCIBLE_QUINTIC_23',p,a['canonical_eigen_coordinate_bivector'],s.symbols('u0:5'),a['canonical_eigen_coordinate_quintic']),
                    ('R4_NIL22_LONG_IRREDUCIBLE_QUINTIC_23',p2,r['primitive_eigen_tensor'],s.symbols('h x y z w'),r['primitive_eigen_quintic'])]
    for name,path,source,variables,knownF in specifications:
        P=tensor(source,variables)
        D=quintic(P)
        F=s.sympify(knownF,locals=dict(zip(map(str,variables),xx)))
        scalar=s.cancel(D/F)
        assert scalar.is_Rational and scalar!=0
        monomials=sorted({e for i in range(5) for e in s.Poly(s.diff(D,xx[i]),*xx).monoms()},reverse=True)
        G=s.Matrix([[s.Poly(s.diff(D,q),*xx).coeff_monomial(e) for q in xx] for e in monomials])
        piv=list(G.T.rref()[1])
        assert len(piv)==5
        determinant=G[piv,:].det()
        assert determinant!=0
        irred.append(dict(id=name,accepted_report=str(path.relative_to(ROOT)),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                          actual_quintic=str(D),quintic_scalar_to_accepted_normalform=str(scalar),
                          nonzero_gradient_five_minor=str(determinant),
                          minor_rows=[list(monomials[i]) for i in piv],essential_variable_rank=5,
                          membership='Actual tensor is the accepted canonical source of this full23 component; reconstructed full Euler/Pf quintic, not an arbitrary polynomial example.'))
        paths.append(path)
    lookup={r['id'] for r in mps+irred}
    reasons=[]
    for row in rows:
        name=row['id']
        if row['generic_rank']!=4:reason='Different generic rank of full component.'
        elif row['projective_dimension']!=23:reason='Different full component dimension.'
        elif name in lookup:reason='Whole new closure has closed essential-quintic-rank<=3, but this component has the saved actual essential-rank5 member.'
        else:
            divisor=row['generic_actual_degeneracy_divisor']
            assert divisor['irreducible_factor_degrees']!=[5] or not divisor.get('reduced',True)
            reason='Old accepted generic divisor is reducible/nonreduced; its whole component stays in that closed multiplication locus. New actual quintic is absolutely irreducible.'
        reasons.append(dict(id=name,reason=reason))
    return dict(immutable_prior103_catalogue=snapshot['immutable_catalogue'],
                input_projection=dict(path=str(input_path),sha256=hashlib.sha256(input_path.read_bytes()).hexdigest()),
                all103_comparison_rows=rows,all103_distinctness_reasons=reasons,
                all29_same23_MPS_actual_SNC_witnesses=mps,
                both_old_irreducible_quintic_actual_rank5_witnesses=irred,
                witness_input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                closed_invariant='The rank of the five partial derivatives of the actual quintic is <=3 on the entire new orbitclosure; this is a closed determinantal condition, including degenerations.',
                no_component_count_changed=True)


def eight_curves(input_path):
    """Complete finite eight-class source verification, without ambient point tests."""
    from sympy.polys.matrices import DomainMatrix
    def domain(A):
        return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    def certificate(A):
        cc=list(domain(A).rref()[1])
        rr=list(domain(A[:,cc].T).rref()[1])
        return dict(shape=list(A.shape),rank=len(cc),rows=rr,columns=cc,
                    determinant=str(domain(A[rr,cc]).det()),
                    all_rows=[[str(c) for c in A.row(i)] for i in range(A.rows)])
    def exponents(n,k):
        return sorted((a for a in itertools.product(range(n+1),repeat=k) if sum(a)==n),reverse=True)
    def polynomial(e,variables):
        return s.prod(v**a for a,v in zip(e,variables))
    raw=json.loads(input_path.read_text())
    assert len(raw['candidates'])==8
    ss,tt=s.symbols('s t')
    x,y,z,u,v=variables=s.symbols('x y z u v')
    base=variables[:3]
    quadratic=[x*x,x*y,x*z,y*y,y*z,z*z]
    rows=[]
    for record in raw['candidates']:
        phi=[s.expand(s.sympify(text,locals={'s':ss,'t':tt})) for text in record['parametrization']]
        assert all(all(sum(e)==5 for e in s.Poly(f,ss,tt).monoms()) for f in phi)
        binary_gcd=s.gcd(s.gcd(phi[0],phi[1]),phi[2])
        assert not binary_gcd.has(ss,tt)
        implicit={}
        for degree in [4,5]:
            source_exponents=exponents(degree,3)
            target_exponents=exponents(5*degree,2)
            polys=[s.Poly(s.expand(polynomial(e,phi)),ss,tt) for e in source_exponents]
            matrix=s.Matrix([[f.coeff_monomial(e) for f in polys] for e in target_exponents])
            ker=matrix.nullspace()
            assert len(ker)==(0 if degree==4 else 1)
            implicit[degree]=dict(coefficient_map=certificate(matrix),
                                  source_exponents=source_exponents,target_exponents=target_exponents,
                                  complete_kernel=[[str(c) for c in col] for col in ker])
        vector=s.Matrix([s.Rational(c) for c in implicit[5]['complete_kernel'][0]])
        den=s.ilcm(1,*[c.q for c in vector])
        integers=[int(c*den) for c in vector]
        gcd=s.igcd(*integers)
        integers=[c//gcd for c in integers]
        if next(c for c in integers if c)<0:integers=[-c for c in integers]
        F=s.expand(sum(c*polynomial(e,base) for c,e in zip(integers,implicit[5]['source_exponents'])))
        assert s.expand(F.subs(dict(zip(base,phi)),simultaneous=True))==0
        # The image of P1 is irreducible over C. A proper factor of a degree5
        # relation would give a degree<=4 relation, excluded by full rank above.
        gradient=[s.diff(F,a) for a in base]
        gmons=exponents(4,3)
        G=s.Matrix([[s.Poly(f,*base).coeff_monomial(e) for f in gradient] for e in gmons])
        gc=certificate(G)
        assert gc['rank']==3
        target=exponents(6,3)
        columns=[s.Poly(m*gradient[i],*base) for i in range(3) for m in quadratic]
        A=s.Matrix([[f.coeff_monomial(e) for f in columns] for e in target])
        kernel=A.nullspace()
        vectors=[s.Matrix([sum(k[6*i+j]*quadratic[j] for j in range(6)) for i in range(3)]) for k in kernel]
        fields=[]
        for W in vectors:
            assert s.expand(sum(W[i]*gradient[i] for i in range(3)))==0
            fields.append(list(map(str,W)))
        accepted=False
        det=None
        coefficient=None
        if len(kernel)==2:
            det=s.expand(s.Matrix.hstack(s.Matrix(base),*vectors).det())
            coefficient=s.cancel(det/F)
            accepted=bool(coefficient.is_Rational and coefficient!=0)
        derivative=s.Matrix(phi).diff(ss).cross(s.Matrix(phi).diff(tt))
        differential_gcd=s.gcd(s.gcd(derivative[0],derivative[1]),derivative[2])
        assert differential_gcd!=0
        sq=s.sqf_list(differential_gcd,ss,tt)[1]
        multiplicities=sorted(int(m) for factor,m in sq for _ in range(s.Poly(factor,ss,tt).total_degree()))
        matches=multiplicities==sorted(record['expected_complex_derivative_gcd_root_multiplicities'])
        current=dict(id=record['id'],expected_curve_type=record['expected_type'],
                     exact_degree5_parametrization=list(map(str,phi)),no_binary_basepoint=True,
                     complete_implicit_coefficient_maps={str(d):implicit[d] for d in [4,5]},
                     primitive_quintic=str(F),complex_absolute_irreducibility=dict(
                         all_degree4_relations_absent=True,unique_degree5_relation=True,
                         justification='Irreducible image of P1 over C; any proper factor would vanish on that image and contradict full degree4 coefficient rank. Polynomial degree<4 relations would give a degree4 relation by multiplication.'),
                     essential_gradient_certificate=gc,essential_gradient_monomials=gmons,
                     complete_quadratic_gradient_kernel=dict(coefficient_map=certificate(A),
                         all_degree6_target_exponents=target,quadratic_monomials=list(map(str,quadratic)),
                         dimension=len(kernel),full18_basis=[[str(c) for c in k] for k in kernel],
                         vector_fields=fields,saito_determinant=str(det),determinant_over_F=str(coefficient)),
                     derivative_gcd=str(differential_gcd),complex_differential_root_multiplicities=multiplicities,
                     differential_multiplicities_match_class=matches,
                     free_exponents_2_2=accepted)
        if accepted:
            U0,V0=vectors
            g00,g01,g10,g11=s.symbols('g00 g01 g10 g11')
            magnetic=s.symbols('r0:6')
            U=g00*U0+g01*V0
            V=g10*U0+g11*V0
            divU=sum(s.diff(U[i],base[i]) for i in range(3))
            divV=sum(s.diff(V[i],base[i]) for i in range(3))
            R=sum(c*m for c,m in zip(magnetic,quadratic))
            P=s.zeros(5)
            for i in range(3):
                P[i,3]=s.expand(U[i]);P[3,i]=-P[i,3]
                P[i,4]=s.expand(V[i]);P[4,i]=-P[i,4]
            P[3,4]=s.expand(-divV*u+divU*v+R);P[4,3]=-P[3,4]
            jacobi=[]
            for i,j,k in itertools.combinations(range(5),3):
                J=s.expand(sum(P[i,l]*s.diff(P[j,k],variables[l])+P[j,l]*s.diff(P[k,i],variables[l])+P[k,l]*s.diff(P[i,j],variables[l]) for l in range(5)))
                assert J==0
                jacobi.append(dict(triple=[i,j,k],raw_polynomial=str(J)))
            divergence=[s.expand(sum(s.diff(P[i,j],variables[j]) for j in range(5))) for i in range(5)]
            casimir=[s.expand(sum(P[i,j]*s.diff(F,variables[j]) for j in range(5))) for i in range(5)]
            assert divergence==casimir==[0]*5
            pfaffians=[]
            for missing in range(5):
                i,j,k,l=[a for a in range(5) if a!=missing]
                pfaffians.append(s.expand(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k]))
            actual=s.expand(sum((-1)**i*variables[i]*pfaffians[i] for i in range(5)))
            frame_det=g00*g11-g01*g10
            assert s.expand(actual+coefficient*frame_det*F)==0
            current['whole_GL2_frame_and_all6_magnetic_source']=dict(
                variables=list(map(str,variables)),frame_unit=str(frame_det),magnetic_parameters=list(map(str,magnetic)),
                all_ten_brackets={f'{i},{j}':str(P[i,j]) for i,j in itertools.combinations(range(5),2)},
                all_ten_raw_strict_Jacobi=jacobi,all_five_divergences=list(map(str,divergence)),
                all_five_Casimir_contractions=list(map(str,casimir)),
                all_five_raw_pfaffians=list(map(str,pfaffians)),actual_quintic=str(actual),
                exact_actual_quintic_factor=str(-coefficient*frame_det))
            current['total_Tjurina_number']=12
            current['tau_proof']='Complete Saito(2,2) gives saturated Hilbert-Burch resolution S(-6)^2 -> S(-4)^3 -> J_F, numerator1-3t^4+2t^6 and length12.'
            current['full_class_recognition']=bool(matches)
            current['full_class_reason']='Canonical parametrization is the complete pureproof class; all zero/nonzero magnetic/frame cases are covered by the accepted all-F actual inverse.' if record['canonical_class_from_pureproof'] else 'Free(2,2) forces rational cuspidal. Three simple differential roots give exactly three multiplicity2 cusps, so the complete pureproof list forces3A4; uniqueness identifies the algebraic canonical class, including both conjugate choices.'
        else:
            current['full_class_recognition']=bool(record['canonical_class_from_pureproof'])
            current['full_class_reason']='The canonical class has no complete Saito(2,2) basis, so the whole class is excluded by the canonical necessary theorem.' if record['canonical_class_from_pureproof'] else 'A failed rational alternative cannot exclude the algebraic3A4 canonical class; keep that class unresolved.'
        if not matches:current['full_class_recognition']=False
        rows.append(current)
        print(json.dumps(dict(curve=record['id'],implicit_degree=5,quadratic_kernel=len(kernel),free_2_2=accepted,identified=current['full_class_recognition'])),flush=True)
    return dict(all_eight_candidate_sources=rows,
                finite_curve_source_exhaustion_closed=all(r['full_class_recognition'] for r in rows),
                free_curve_source_count=sum(r['free_exponents_2_2'] for r in rows),
                no_new_component_count_changed=True,
                componenthood_scope='Each accepted fullsource is one actual projective23 orbit by the all-F magnetic theorem. No further ambient derivative or componenthood/known-parent recognition has been computed in this job.',
                no_new_ambient_points_or_unbounded_CAS=True)


def eight_ambient(input_path,prior104_path):
    """One finite seven-new-anchor ambient batch; FIRST104 is reused unchanged."""
    import importlib.util
    from sympy.polys.matrices import DomainMatrix
    def domain(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    def sparse(v):return {str(i):str(c) for i,c in enumerate(v) if c}
    def decode(row,n):
        out=s.zeros(n,1)
        for i,c in row.items():out[int(i)]=s.Rational(c)
        return out
    class InlineSource:
        def __init__(self,payload,pointer):
            self.payload=payload
            self.raw=json.dumps(payload,sort_keys=True).encode()
            self.pointer=pointer
        def read_text(self):return self.raw.decode()
        def read_bytes(self):return self.raw
        def __str__(self):return str(input_path)+'#'+self.pointer
    helper=ROOT/'legacy/v03/code/poisson_p4_components_check.py'
    spec=importlib.util.spec_from_file_location('essential_eight_ambient_support',helper)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    raw=json.loads(input_path.read_text())
    sources=raw['eight_curve_sources']
    assert raw['all_passed'] and sources['finite_curve_source_exhaustion_closed'] and sources['free_curve_source_count']==8
    prior=json.loads(prior104_path.read_text())
    assert prior['immutable_catalogue']['count']==104
    oldaudit_path=ROOT/'reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_FIRST_PRIOR103_20261107.json'
    oldaudit=json.loads(oldaudit_path.read_text())['prior103_audit']
    assert {r['id'] for r in prior['all104_comparison_rows']}=={r['id'] for r in oldaudit['all103_comparison_rows']}|{'R4_PLANE_UNICUSPIDAL_QUINTIC_45_23'}
    first_path=ROOT/'reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_FIRST_AMBIENT_20261107.json'
    first_support_path=ROOT/'reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_FIRST_COMPONENT_SUPPORT_20261107.json'
    assert json.loads(first_support_path.read_text())['component_support']['projective_orbit_closure_is_full_component']
    def jac(vec):
        bv={}
        for c,(pair,e) in zip(vec,m.BIVECTOR_BASIS):
            if c:bv.setdefault(pair,{})[e]=c
        J=m.affine_jacobi_coefficients(bv)
        out=s.zeros(350,1)
        for triple,polynomial in J.items():
            for i in range(5):
                sign=m.permutation_sign((i,)+triple)
                if not sign:continue
                quad=tuple(sorted((i,)+triple))
                for e,c in polynomial.items():
                    exponent=list(e);exponent[i]+=1
                    out[m.FOUR_VECTOR_INDEX[(quad,tuple(exponent))]]+=sign*c
        return out
    cases=[]
    for index,source in enumerate(sources['all_eight_candidate_sources']):
        assert source['free_exponents_2_2'] and source['full_class_recognition']
        if source['id']=='C45':
            cases.append(dict(id='C45',curve_type=source['expected_curve_type'],
                              status='known_first104_component23',
                              immutable_full_ambient=dict(path=str(first_path.relative_to(ROOT)),sha256=hashlib.sha256(first_path.read_bytes()).hexdigest()),
                              immutable_component_support=dict(path=str(first_support_path.relative_to(ROOT)),sha256=hashlib.sha256(first_support_path.read_bytes()).hexdigest()),
                              no_first104_derivative_or_obstruction_rerun=True))
            continue
        payload={'all_passed':True,'candidates':[{'id':source['id'],'whole_cotangent_source':source['whole_GL2_frame_and_all6_magnetic_source']}]}
        inline=InlineSource(payload,f'/eight_curve_sources/all_eight_candidate_sources/{index}')
        point=ambient(inline)
        point['input_report']=dict(path=str(input_path),sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
                                   json_pointer=inline.pointer,derived_payload_sha256=hashlib.sha256(inline.raw).hexdigest())
        rank=point['full_350_by_150_D']['rank']
        assert rank<=102
        result=dict(id=source['id'],curve_type=source['expected_curve_type'],full_ambient=point)
        if rank==102:
            result.update(status='smooth23_component_proved',
                support_reason='Full tangent raw48 equals complete orbit/vertical/scalar48; smooth unique projective23 component equals the whole orbitclosure.')
        elif rank==101:
            D=s.Matrix.vstack(*[decode(r,150).T for r in point['full_D_sparse_rows']])
            R=decode(point['true_cut']['all150_true_normals'][0],150)
            dc=point['true_cut']['independent_implicit_coordinate_indices']
            dr=list(domain(D[:,dc].T).rref()[1])
            A=D[dr,dc];inverse=domain(A).inv().to_Matrix()
            assert A*inverse==s.eye(101)
            Q=jac(R);T=s.zeros(150,1)
            for i,c in zip(dc,-inverse*Q[dr,:]):T[i]=c
            residual=D*T+Q
            assert residual[dr,:]==s.zeros(101,1)
            mu=point['true_cut']['normal_weights'][0]
            bw=point['Cstar']['relative_input_weights'];ow=point['Cstar']['relative_output_weights']
            assert all(not c or bw[i]==2*mu for i,c in enumerate(T))
            assert all(not c or ow[i]==2*mu for i,c in enumerate(residual))
            quadratic=dict(true_normal_character=mu,IFT_columns=dc,IFT_rows=dr,
                           nonzero101minor=str(domain(A).det()),
                           all101inverse_rows=[[str(c) for c in inverse.row(i)] for i in range(101)],
                           full350_Jacobi_quadratic_source=sparse(Q),unique_cut_quadratic_correction=sparse(T),
                           all350_true_quadratic_residual=sparse(residual))
            if any(residual):
                at=next(i for i,c in enumerate(residual) if c)
                lam=s.zeros(1,350);lam[at]=1
                for i,c in zip(dr,D[at,dc]*inverse):lam[i]-=c
                assert lam*D==s.zeros(1,150) and (lam*Q)[0]==residual[at]!=0
                quadratic.update(cokernel_functional=sparse(lam),all150_D_columns_annihilated=True,
                                 nonzero_pairing=str(residual[at]),residual_row=at)
                result.update(status='quadratic23_component_proved',
                    support_reason='One true normal has a nonzero quadratic cokernel coefficient, hence reduced truecut support is the origin. Krull upper48 and the complete source48 prove a full projective23 component.')
            else:
                result.update(status='open_one_normal_quadratic_zero',
                    support_reason='The single quadratic residual vanishes. No cubic or further coefficient was run; choose a finite weakest model from the saved character before continuing.')
            result['weakest_one_normal_quadratic']=quadratic
        else:
            result.update(status='open_multiple_true_normals',
                support_reason='At least two true normals remain. All frames/characters saved; no nonlinear coefficient was computed without a new mathematical support criterion.')
        result['prior104_distinctness_if_full_component23']=dict(
            all103_reasons_reused=oldaudit['all103_distinctness_reasons'],
            common_hypotheses='This full source has absolutely irreducible actualD and wholeclosed essentialrank<=3, exactly the two invariant hypotheses in the accepted all103 audit.',
            previous104_reason='The old104 curve has one multiplicity4 cusp; every remaining curve here has maximum cusp multiplicity<=3. Both fixed-F sources are actual orbit23, so different curve orbits cannot have the same23orbitclosure.',
            scope='These distinguish full components only once their local23 componenthood is proved; they do not exclude unresolved orbit points from larger known parents.')
        cases.append(result)
        print(json.dumps(dict(curve=source['id'],D_rank=rank,M_rank=point['actual_Euler_GL_scalar_source']['rank'],
                              normal_weights=point['true_cut']['normal_weights'],status=result['status'])),flush=True)
    pairs=[dict(first=a['id'],second=b['id'],reason='Different proved plane-cusp configurations or maximum multiplicity; source8 exact implicit/parametrization recognition and the pureproof uniqueness give different PGL3 curve orbits. Each fixed-F source is a dense actual orbit23, so its23closure cannot equal another such closure.')
           for a,b in itertools.combinations(cases,2)]
    return dict(all_eight_cases=cases,immutable_prior104_catalogue=prior['immutable_catalogue'],
                prior104_projection=dict(path=str(prior104_path),sha256=hashlib.sha256(prior104_path.read_bytes()).hexdigest()),
                accepted_prior103_audit=dict(path=str(oldaudit_path.relative_to(ROOT)),sha256=hashlib.sha256(oldaudit_path.read_bytes()).hexdigest()),
                all28_source_closure_distinctness_reasons=pairs,
                new_local_componenthood_count=sum(r['status'] in ['smooth23_component_proved','quadratic23_component_proved'] for r in cases),
                all_other_cases_remain_explicit=True,no_new_component_count_registered=True,
                no_cubic_or_higher_obstruction_computed=True,no_private_work_input=True)


def main():
    ap = argparse.ArgumentParser()
    selector = ap.add_mutually_exclusive_group(required=True)
    selector.add_argument('--input', type=Path)
    selector.add_argument('--ambient-input', type=Path)
    selector.add_argument('--support-input', type=Path)
    selector.add_argument('--prior103-input', type=Path)
    selector.add_argument('--eight-curves-input', type=Path)
    selector.add_argument('--eight-ambient-input', type=Path)
    ap.add_argument('--prior104', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    started = time.monotonic()
    selected=args.input or args.ambient_input or args.support_input or args.prior103_input or args.eight_curves_input or args.eight_ambient_input
    raw = selected.read_bytes()
    if args.input:
        rows = [candidate(r) for r in json.loads(raw)['candidates']]
        body=dict(candidates=rows)
    elif args.ambient_input:
        body=dict(full_ambient=ambient(selected))
    elif args.support_input:
        body=dict(component_support=component_support(selected))
    elif args.prior103_input:
        body=dict(prior103_audit=prior103_audit(selected))
    elif args.eight_curves_input:
        body=dict(eight_curve_sources=eight_curves(selected))
    else:
        assert args.prior104 is not None
        body=dict(eight_full_ambient=eight_ambient(selected,args.prior104))
    result = dict(task_id='E_reduced_nonnormal', all_passed=True, **body,
                  input=dict(path=str(selected), sha256=hashlib.sha256(raw).hexdigest()),
                  producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  producer_bytes=Path(__file__).read_text(), sympy_version=s.__version__,
                  seconds=time.monotonic()-started, max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(dict(all_passed=True, seconds=result['seconds'], max_rss_kib=result['max_rss_kib'])))


if __name__ == '__main__':
    main()
