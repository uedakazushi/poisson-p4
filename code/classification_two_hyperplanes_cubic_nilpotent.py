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


def normal_quadrics(report_path, case_index, output, record_override=None):
    """One finite genuine quadratic obstruction, with all actual source directions."""
    from sympy.polys.matrices import DomainMatrix
    from classification_rank3_mixed_component_review import COLUMNS, jacobi, sparse, vector
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU, (55, 56))
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2,)*2)
    record = record_override if record_override is not None else json.loads(report_path.read_text())['branches'][case_index]
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
    functional=None
    all_functionals=[]
    if qp:
        obstruction_rows=list(DomainMatrix.from_Matrix(obstruction.T).convert_to(s.QQ).rref()[1])
        inverse=block.inv().to_Matrix()
        for row in obstruction_rows:
            coefficients=D[row,list(pivots)]*inverse;ell=s.zeros(1,350);ell[row]=1
            for k,r in enumerate(rows):ell[r]-=coefficients[k]
            assert ell*D==s.zeros(1,150)
            all_functionals.append(dict(nonzero_coefficients={str(i):str(a)for i,a in enumerate(ell)if a},
                                        quadratic_pairing=list(map(str,list(ell*Q)))))
    if len(normals)==1 and any(obstruction):
        row=next(i for i,a in enumerate(obstruction[:,0])if a)
        inverse=block.inv().to_Matrix();coefficients=D[row,list(pivots)]*inverse
        ell=s.zeros(1,350);ell[row]=1
        for k,r in enumerate(rows):ell[r]-=coefficients[k]
        assert ell*D==s.zeros(1,150)
        functional=dict(nonzero_coefficients={str(i):str(a)for i,a in enumerate(ell)if a},pairing=str((ell*Q)[0]))
    result=dict(task_id='E_reduced_nonnormal',source_report=str(report_path),case=record['name'],
        full_D_rank=len(pivots),actual_source_rank=source_rank,normal_count=len(normals),
        all_D_times_directions_zero=True,normal_indices=normal_indices,
        normals=[{f'{i},{j}':str(R[i,j]) for i,j in itertools.combinations(range(5),2) if R[i,j]} for R in normals],
        normal_monomials=list(map(str,monomials)),quadratic_initial_generators=list(map(str,quadrics)),
        zero_dimensional_quadratic_support=bool(gb and gb.is_zero_dimensional),
        grevlex_basis=[str(g.as_expr()) for g in gb.polys] if gb else [],
        literal_Poisson_normal_axes=[i for i,q in enumerate(squares) if q==s.zeros(350,1)],
        one_normal_quadratic_cokernel_functional=functional,
        quadratic_cokernel_functionals=all_functionals,
        full_initial_ideal_not_claimed=True,source_report_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),
        producer_source=Path(__file__).read_text(),seconds=time.monotonic()-start)
    if output is None:return result
    with output.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:result[k] for k in ['case','normal_count','quadratic_initial_generators','zero_dimensional_quadratic_support','literal_Poisson_normal_axes','seconds']}))


def universal_pencils(output):
    """The actual n3 modulus of the p0 source, not another ambient point probe."""
    from sympy.polys.matrices import DomainMatrix
    rpath=ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_COMPONENTS_20261107.json'
    npath=ROOT/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_P0_NORMALS_20261107.json'
    record=json.loads(rpath.read_text())['branches'][0];normal=json.loads(npath.read_text())
    x=s.symbols('h x1:5');lam,t=s.symbols('lam t')
    def matrix(data):
        A=s.zeros(5)
        for key,f in data.items():
            i,j=map(int,key.split(','));A[i,j]=s.sympify(f);A[j,i]=-A[i,j]
        return A
    U=matrix(record['universal_brackets'])
    pt={s.Symbol(p):0 for p in record['full_source_specification']['free_parameters']}
    P=U.subs(pt);T=U.diff(s.Symbol('n3')).subs(pt)
    def projective(A):
        J={ijk:s.expand(sum(A[i,l]*s.diff(A[j,k],x[l])+A[j,l]*s.diff(A[k,i],x[l])+A[k,l]*s.diff(A[i,j],x[l]) for l in range(5)))
           for ijk in itertools.combinations(range(5),3) for i,j,k in [ijk]}
        return [s.expand(sum((-1)**k*x[inds[k]]*J[inds[:k]+inds[k+1:]] for k in range(4)))
                for inds in itertools.combinations(range(5),4)]
    assert all(f==0 for f in projective(P+lam*T))
    normals=list(map(matrix,normal['normals']));coeff=s.symbols('r0:'+str(len(normals)))
    R=sum((a*A for a,A in zip(coeff,normals)),s.zeros(5))
    mixed=[s.expand(f-g-h) for f,g,h in zip(projective(T+R),projective(T),projective(R))]
    equations=[c for f in mixed for c in s.Poly(f,*x).coeffs() if c]
    M,_=s.linear_eq_to_matrix(equations,coeff);kernel=M.nullspace()
    aa=s.symbols('a0:'+str(len(kernel)))
    combo=sum((a*K for a,K in zip(aa,kernel)),s.zeros(len(coeff),1))
    R=R.subs(dict(zip(coeff,combo)),simultaneous=True)
    residual=[s.factor(c) for f in projective(R) for c in s.Poly(f,*x).coeffs() if c]
    residual=list(dict.fromkeys(residual))
    result=dict(task_id='E_reduced_nonnormal',section_parameter='n3',
                full_section_projective_Jacobi_zero=True,
                modulus_compatible_normal_kernel=[list(map(str,k)) for k in kernel],
                compatible_normal={f'{i},{j}':str(s.expand(R[i,j])) for i,j in itertools.combinations(range(5),2) if R[i,j]},
                remaining_Poisson_equations=list(map(str,residual)),
                whole_source_coverage_not_inferred=True,
                producer_source=Path(__file__).read_text())
    with output.open('x') as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:result[k] for k in ['modulus_compatible_normal_kernel','remaining_Poisson_equations']}))


def translation_rank_two_nonplane(output):
    """Two entire exact source branches, then one canonical ambient test each."""
    if output.exists():raise FileExistsError(str(output))
    from sympy.polys.matrices import DomainMatrix
    from classification_three_hyperplanes_quadric_independent import ambient_matrices
    from classification_rank3_mixed_component_review import COLUMNS,jacobi,sparse,minor_certificate
    resource.setrlimit(resource.RLIMIT_CPU,(150,151))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    begun=time.monotonic();h,x,u,v,w=xx=s.symbols('h x u v w');ns=s.symbols('n0:21');r=s.Symbol('r');kap=s.Symbol('kap');t=s.Symbol('t')
    source_path='reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json'
    row=next(z for z in json.loads((ROOT/source_path).read_text())['cases']if z['case']=='j111_rank2')
    env={str(z):z for z in [*xx,*ns]};raw=s.zeros(5)
    for ij,text in row['brackets'].items():
        i,j=map(int,ij.split(','));raw[i,j]=s.sympify(text,locals=env);raw[j,i]=-raw[i,j]
    equations=[s.sympify(f,locals=env)for f in row['remaining_coefficient_equations']]
    def dm(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    def piv(A):return list(dm(A).rref()[1])
    def enc(A):return {f'{i},{j}':str(s.expand(A[i,j]))for i,j in itertools.combinations(range(5),2)if A[i,j]}
    def mat(A):return [[str(v)for v in z]for z in A.tolist()]
    def nz(A):return [[i,j,str(A[i,j])]for i in range(A.rows)for j in range(A.cols)if A[i,j]]
    def strict(P):return [s.expand(sum(P[i,l]*s.diff(P[j,jj],xx[l])+P[j,l]*s.diff(P[jj,i],xx[l])+P[jj,l]*s.diff(P[i,j],xx[l])for l in range(5)))for i,j,jj in itertools.combinations(range(5),3)]
    def push(P,G):return (G*P.subs(dict(zip(xx,G.inv()*s.Matrix(xx))),simultaneous=True)*G.T).applyfunc(s.cancel)
    def from_vector(col):
        A=s.zeros(5)
        for (ij,ex),c in zip(COLUMNS,col):
            if c:A[ij]+=c*s.prod(z**e for z,e in zip(xx,ex))
        return A-A.T
    records=[]
    for branch in [0,1]:
        started=time.monotonic();bindings={ns[11]:s.S.One,ns[18]:s.Integer(branch)};steps=[]
        def eliminate(index,var,allowed_unit):
            expr=s.expand(equations[index].subs(bindings,simultaneous=True));poly=s.Poly(expr,var);assert poly.degree()==1
            coefficient=poly.nth(1);assert s.expand(coefficient/allowed_unit).free_symbols==set() and coefficient!=0
            answer=s.cancel(-poly.nth(0)/coefficient);bindings[var]=answer
            steps.append(dict(original_equation_index=index,pivot=str(var),previously_substituted_equation=str(expr),coefficient=str(coefficient),solution=str(answer),allowed_unit=str(allowed_unit)))
            for rep in range(21):
                updated={j:s.expand(f.subs(bindings,simultaneous=True))for j,f in bindings.items()}
                if updated==bindings:break
                bindings.clear();bindings.update(updated)
        for idx,j in [(2,16),(5,19),(8,20),(9,6),(0,0),(3,2),(1,1),(7,5),(4,3)]:eliminate(idx,ns[j],s.S.One)
        if branch==0:
            partialF=s.expand((-h*raw[2,4]+x*raw[3,4]).subs(bindings,simultaneous=True))
            assert s.expand(partialF.subs(h,0)-ns[17]*x**3)==0
            eliminate(10,ns[12],ns[17]);eliminate(12,ns[13],ns[17]);free=[ns[j]for j in [4,7,8,9,10,14,15,17]]
            point={z:0 for z in free};point.update({ns[8]:-1,ns[14]:-1,ns[17]:1})
        else:
            eliminate(14,ns[13],s.S.One);eliminate(12,ns[12],s.S.One);eliminate(10,ns[17],s.S.One)
            free=[ns[j]for j in [4,7,8,9,10,14,15]];point={z:0 for z in free};point.update({ns[8]:-1,ns[10]:1,ns[14]:-1,ns[15]:1})
        for rep in range(21):bindings={j:s.expand(f.subs(bindings,simultaneous=True))for j,f in bindings.items()}
        assert set(ns)-set(bindings)==set(free)
        assert all(s.expand(f.subs(bindings,simultaneous=True))==0 for f in equations)
        P=raw.subs(bindings,simultaneous=True).applyfunc(s.expand);assert all(f==0 for f in strict(P))
        F=s.expand(-h*P[2,4]+x*P[3,4]);assert s.expand(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3]-x*F)==0
        normalizer=s.eye(5)
        if branch==0:
            normalizer[2,0]=-ns[7];normalizer[2,1]=ns[14]-ns[8];normalizer[2,4]=-ns[9]
            normalizer[3,1]=ns[10]-ns[15];normalizer[3,4]=-ns[4]-2
            normalF=h*h*u-h*x*v+ns[17]*x**3
            kap_values={z:0 for z in free};kap_values.update({ns[8]:-kap,ns[14]:-kap,ns[17]:1})
            fixed=s.eye(5);fixed[3,4]=-2;weights=[0,1,3,2,2];overall=2
        else:
            c0=-ns[7];c1=ns[14]-ns[8];c2=ns[15]-ns[10];c3=bindings[ns[17]]
            normalizer[2,0]=c0;normalizer[2,4]=-ns[9]
            normalizer[3,1]=c3;normalizer[3,4]=-ns[4]-ns[9]
            normalizer[4,0]=c1+c0;normalizer[4,1]=c2+c3
            normalF=(h-x)*(h*u-x*v)+h*x*w
            kap_values={z:0 for z in free};kap_values.update({ns[8]:-kap,ns[10]:kap,ns[14]:-kap,ns[15]:kap})
            fixed=s.eye(5);weights=[0,0,1,1,1];overall=1
        assert normalizer.det()==1
        assert s.expand(F.subs(dict(zip(xx,normalizer.inv()*s.Matrix(xx))),simultaneous=True)-normalF)==0
        normalized_P=push(P,normalizer)
        assert normalized_P[1,2]==x*x and normalized_P[1,3]==h*x and normalized_P[1,4]==0
        model=push(P.subs(kap_values,simultaneous=True),fixed);base=model.subs(kap,0)
        V=s.Matrix([0,0,x,h,0]);Z=s.Matrix([sum(s.diff(base[i,j],xx[i])for i in range(1,5))for j in range(5)])
        assert V[0]==Z[0]==0
        for field in [V,Z]:
            lie=sum((field[j]*base.diff(xx[j])for j in range(5)),s.zeros(5))-field.jacobian(xx)*base-base*field.jacobian(xx).T
            assert lie.applyfunc(s.expand)==s.zeros(5)
        assert (Z.jacobian(xx)*V-V.jacobian(xx)*Z).applyfunc(s.expand)==s.zeros(5,1)
        assert (model-base-kap*(V*Z.T-Z*V.T)).applyfunc(s.expand)==s.zeros(5)
        gt=s.diag(*[t**j for j in weights]);assert (t**(-overall)*push(model,gt)-model.subs(kap,kap*t**overall)).applyfunc(s.cancel)==s.zeros(5)
        pointP=P.subs(point).applyfunc(s.expand);directions=[P.diff(z).subs(point).applyfunc(s.expand)for z in free]
        D,Vert,Orbit,Family=ambient_matrices(pointP,xx,directions)
        matrices=dict(D=D,vertical=Vert,orbit_scale=Orbit,full_source=Family);ranks={name:len(piv(M))for name,M in matrices.items()}
        certificates={name:minor_certificate(M,ranks[name])for name,M in matrices.items()}
        rec=dict(branch='L_h'if branch==0 else'L_h_minus_x',all_passed=True,free_source_parameters=list(map(str,free)),source_units=['n17']if branch==0 else[],
            necessity_triangular_unit_steps=steps,complete_original_n0_to_n20_binding={str(z):str(bindings.get(z,z))for z in ns},all_original_Jacobi_coefficients_zero=True,
            whole_universal_brackets=enc(P),all10_strict_Jacobi_universal_zero=True,actual_cubic=str(F),actual_quintic=str(h*x*F),
            whole_source_cubic_normalizer=mat(normalizer),whole_source_cubic_normalizer_inverse=mat(normalizer.inv()),normal_form_cubic=str(normalF),
            normalizer_preserves_full_cofactor=True,normalizer_determinant='1',whole_source_cubic_normal_form_identity=True,
            all_parameter_absolute_irreducibility_reason='Primitive linear polynomial in u with unit n17 preventing an h common factor'if branch==0 else'Primitive linear polynomial in u,v,w with coefficient gcd1',
            actual_essential_variable_count=4 if branch==0 else 5,canonical_actual_kappa_model=enc(model),canonical_modular_field=list(map(str,Z)),
            canonical_kappa_is_exact_decomposable_Poisson_twist=True,kappa_nonzero_all_one_projective_orbit=True,canonical_torus_weights=weights,canonical_tensor_weight=overall,
            canonical_original_source_values={str(z):str(v)for z,v in point.items()},canonical_point_brackets=enc(pointP),canonical_cubic_normalizer=mat(fixed),
            full_all_source_directions=list(map(enc,directions)),ranks_QQ=ranks,certificates=certificates,all350_D_times_real_source_zero=True,
            full_D_sparse=nz(D),full_F_sparse=nz(Family),projective_actual_lower=ranks['full_source']-25,projective_tangent_upper=125-ranks['D'])
        fc=piv(Family);Find=Family[:,fc];fr=piv(Find.T);Finv=dm(Find[fr,:]).inv().to_Matrix();rr,dc0=dm(D).rref();rr=rr.to_Matrix();dc=list(dc0)
        freecols=[j for j in range(150)if j not in dc];K=s.zeros(150,len(freecols))
        for at,j in enumerate(freecols):
            K[j,at]=1
            for i,col in enumerate(dc):K[col,at]=-rr[i,j]
        projected=K-Find*Finv*K[fr,:];ni=piv(projected);N=projected[:,ni]
        for j in range(N.cols):
            col=N[:,j]*s.ilcm(1,*[s.denom(v)for v in N[:,j]]);g=s.igcd(0,*[abs(int(v))for v in col if v]);N[:,j]=col/g
        assert N[fr,:]==s.zeros(len(fr),N.cols) and D*N==s.zeros(350,N.cols)
        I=s.eye(150)[:,dc];U=I-Find*Finv*I[fr,:];total=Find.row_join(N).row_join(U)
        assert total.shape==(150,150) and len(piv(total))==150
        unit=dm(total).det();assert unit
        rec.update(full_source_independent_columns=fc,true_source_cut_rows=fr,normal_count=N.cols,true_normal_columns=mat(N),
            implicit_standard_columns=dc,true_implicit_columns=mat(U),full150_F_N_U_direct_sum_determinant=str(unit),true_cut_dimension=150-len(fc))
        isolated=N.cols==0;support=dict(status='SMOOTH_SQUEEZE'if isolated else'MULTIPLE_NORMALS_STOP',computed=False)
        if 1<=N.cols<=2:
            normals=[from_vector(N[:,j])for j in range(N.cols)];pairs=list(itertools.combinations_with_replacement(range(N.cols),2))
            sq=[jacobi(sparse(A,xx))for A in normals]
            Q=s.Matrix.hstack(*[sq[i]if i==j else jacobi(sparse(normals[i]+normals[j],xx))-sq[i]-sq[j]for i,j in pairs])
            Dind=D[:,dc];dr=piv(Dind.T);inv=dm(D[dr,dc]).inv().to_Matrix();obstruction=Q-Dind*inv*Q[dr,:]
            qr,qp=dm(obstruction).rref();qr=qr.to_Matrix();variables=s.symbols('t0:'+str(N.cols));mon=[variables[i]*variables[j]for i,j in pairs]
            forms=[s.expand(sum(qr[i,j]*mon[j]for j in range(len(mon))))for i in range(len(qp))]
            functionals=[]
            for rowindex in piv(obstruction.T):
                ell=s.zeros(1,350);ell[rowindex]=1;coeff=D[rowindex,dc]*inv
                for j,rindex in enumerate(dr):ell[rindex]-=coeff[j]
                assert ell*D==s.zeros(1,150)
                functionals.append(dict(nonzero_coefficients={str(i):str(v)for i,v in enumerate(ell)if v},full_quadratic_pairing=list(map(str,list(ell*Q)))))
            support=dict(computed=True,normal_pairs=[list(z)for z in pairs],all350_quadratic_columns_sparse=nz(Q),all350_cokernel_columns_sparse=nz(obstruction),
                implicit_rows=dr,implicit_minor_determinant=str(dm(D[dr,dc]).det()),all_cokernel_functionals=functionals,quadratic_row_space_basis=list(map(str,forms)))
            if N.cols==1:
                isolated=bool(forms);support.update(status='PRIMARY_QUADRATIC'if isolated else'NONPRIMARY_QUADRATIC_STOP',power_identity=None if not forms else dict(power='t0^2',generator=str(forms[0])))
            else:
                degree3=[variables[0]**3,variables[0]**2*variables[1],variables[0]*variables[1]**2,variables[1]**3]
                products=[z*f for f in forms for z in variables];Mac=s.Matrix([[s.Poly(f,*variables).coeff_monomial(mon)for f in products]for mon in degree3])
                mc=piv(Mac)if Mac.cols else[];isolated=len(mc)==4
                support.update(status='PRIMARY_QUADRATIC'if isolated else'NONPRIMARY_QUADRATIC_STOP',degree3_macaulay_matrix=mat(Mac),degree3_macaulay_rank=len(mc))
                if isolated:
                    mi=dm(Mac[:,mc]).inv().to_Matrix();support.update(macaulay_columns=mc,macaulay_unit_determinant=str(dm(Mac[:,mc]).det()),all_cubic_monomial_memberships=mat(mi))
                    assert Mac[:,mc]*mi==s.eye(4)
        rec.update(quadratic_support=support,whole_irreducible_source_is_full_component=isolated,projective_component_dimension=ranks['full_source']-25 if isolated else None,
            full_scheme_smooth=N.cols==0,higher_jets_not_computed=True,component_count_not_changed=True,seconds=time.monotonic()-started)
        records.append(rec)
        checkpoint=dict(task_id='E_reduced_nonnormal',all_passed=False,phase='SOURCE_BRANCH_CHECKPOINT',cases=records,
            producer_source=Path(__file__).read_text(),source_report_sha256=hashlib.sha256((ROOT/source_path).read_bytes()).hexdigest(),seconds=time.monotonic()-begun)
        output.write_text(json.dumps(checkpoint,indent=2)+'\n')
        print(json.dumps(dict(branch=rec['branch'],ranks=ranks,normals=N.cols,status=support['status'],component_dimension=rec['projective_component_dimension'],seconds=rec['seconds'])),flush=True)
    paths=['code/classification_two_hyperplanes_cubic_nilpotent.py','code/classification_three_hyperplanes_quadric_independent.py','code/classification_rank3_mixed_component_review.py',source_path,
        'research/classification/TWO_HYPERPLANES_CUBIC_TRANSLATION_RANK_TWO_NONPLANE_SOURCES.md','research/classification/TWO_HYPERPLANES_TRANSLATION_NONPLANE_CANONICAL_MODELS.md',
        'research/classification/inputs/two_hyperplanes_translation_nonplane_preappend_source_contracts_20261108.json']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,cases=records,exact_nonplane_source_exhaustion_replayed=True,plane_A_B_zero_branch_kept_separate=True,
        source_swap_and_normalization_geometry_requires_pure_note_review=True,new_component_count=0,independent_ambient_and_prior_comparison_required=True,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    output.write_text(json.dumps(result,indent=2)+'\n')


def translation_rank_two_global_graphs(output):
    """Faithful parent recognition at the same two canonical source points."""
    from s61_model import S61,X,pushforward,jacobi_projective,quotient_coefficients
    from classification_degenerate_quadric_q4_nil31_I_graph import darboux
    if output.exists():raise FileExistsError(str(output))
    resource.setrlimit(resource.RLIMIT_CPU,(70,71))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    begun=time.monotonic();path='reports/CLASSIFICATION_TWO_HYPERPLANES_TRANSLATION_NONPLANE_COMPLETION_20261108.json'
    source=json.loads((ROOT/path).read_text());assert source['all_passed']
    old=s.symbols('h x u v w');us=s.symbols('u1:62');model=S61();records=[]
    for row in source['cases']:
        started=time.monotonic();P=s.zeros(5)
        for ij,text in row['canonical_point_brackets'].items():
            i,j=map(int,ij.split(','));P[i,j]=s.sympify(text).subs(dict(zip(old,X)),simultaneous=True);P[j,i]=-P[i,j]
        assert all(f==0 for f in jacobi_projective(P))
        evaluation=[1,1,0,0,0]if row['branch']=='L_h'else[1,1,0,0,1]
        value=P.subs(dict(zip(X,evaluation)))[1:5,1:5]
        pf=s.factor(value[0,1]*value[2,3]-value[0,2]*value[1,3]+value[0,3]*value[1,2]);assert pf==1
        move=s.eye(5)
        for j in range(1,5):move[j,0]=-evaluation[j]
        frame=s.eye(5);frame[1:5,1:5]=darboux(value);pre=frame*move
        normalized=model.normalize(pushforward(P,pre));verification=model.verify_point(normalized,check_equations=True)
        total=normalized['matrix']*pre;inverse=total.inv();assert total*inverse==s.eye(5)
        actual=pushforward(P,total);reconstructed=model.reconstruct(normalized['parameters'])
        assert quotient_coefficients(actual)==quotient_coefficients(reconstructed)
        back=pushforward(reconstructed,inverse);assert quotient_coefficients(P)==quotient_coefficients(back)
        assert all(f==0 for f in jacobi_projective(back))
        results={}
        for name,gpath in [('T39','reports/Q51_UNIVERSAL_SECTION_20261101.json'),('T40','reports/Q52_UNIVERSAL_SECTION_20261101.json')]:
            graph=json.loads((ROOT/gpath).read_text());indices=graph['weight2_u_indices_zero_based']
            residual=[s.cancel(normalized['parameters'][i]-s.sympify(q).subs(dict(zip(us,normalized['parameters'])),simultaneous=True))for i,q in zip(indices,graph['quadratic_section_q'])]
            assert len(residual)==45
            results[name]=dict(canonical_point_in_entire_graph=all(v==0 for v in residual),all45_residuals=list(map(str,residual)),nonzero_count=sum(v!=0 for v in residual))
        rec=dict(branch=row['branch'],canonical_original_brackets=row['canonical_point_brackets'],same_previously_declared_kappa_one_model=True,
            structural_rank4_evaluation=evaluation,evaluation_Pfaffian=str(pf),actual_translation=[[str(v)for v in z]for z in move.tolist()],
            Darboux_frame=[[str(v)for v in z]for z in frame.tolist()],perspective=list(map(str,normalized['perspective'])),
            complete_actual_PGL=[[str(v)for v in z]for z in total.tolist()],complete_actual_inverse_PGL=[[str(v)for v in z]for z in inverse.tolist()],PGL_determinant=str(total.det()),
            all61_parameters=list(map(str,normalized['parameters'])),full150_forward_inverse_reconstructions_modulo_Euler=True,
            all134_original_S61_equations=True,all350_original_projective_Jacobi=True,faithful_verification=verification,
            forward_Euler_difference={f'{i},{j}':str(s.expand((actual-reconstructed)[i,j]))for i,j in itertools.combinations(range(5),2)if actual[i,j]!=reconstructed[i,j]},
            inverse_Euler_difference={f'{i},{j}':str(s.expand((P-back)[i,j]))for i,j in itertools.combinations(range(5),2)if P[i,j]!=back[i,j]},
            global_graph_membership=results,whole_source_covered_if_orbit_closure_bridge_accepted=any(q['canonical_point_in_entire_graph']for q in results.values()),
            seconds=time.monotonic()-started)
        records.append(rec);print(json.dumps(dict(branch=row['branch'],membership={k:v['canonical_point_in_entire_graph']for k,v in results.items()},seconds=rec['seconds'])),flush=True)
    paths=[path,'code/classification_two_hyperplanes_cubic_nilpotent.py','code/s61_model.py','code/classification_degenerate_quadric_q4_nil31_I_graph.py',
        'reports/Q51_UNIVERSAL_SECTION_20261101.json','reports/Q52_UNIVERSAL_SECTION_20261101.json',
        'research/classification/TWO_HYPERPLANES_TRANSLATION_NONPLANE_ORBIT_CLOSURES.md',
        'research/classification/S61_GLOBAL_COMPONENT_EVALUATION_FIBER.md','research/classification/S61_GLOBAL_COMPONENT_FIBER_SEPARATE_REVIEW.md',
        'research/classification/inputs/two_hyperplanes_translation_nonplane_pregraph_source_contracts_20261108.json']
    paths+=['legacy/v02/data/'+p for p in ['symplectic_slice_certificate.json','symplectic_slice_equations.json','symplectic_slice_reduced_equations.json']]
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,cases=records,new_component_count=0,
        scope='Same two structurally fixed kappa1 canonical source tensors only. Positive global graph membership covers the entire corresponding22 orbit closure once its source/flag dimension proof is independently accepted. Negative membership excludes that dense orbit, not every boundary and not other parents. No new anchor, normal-quadratic map, or higher jet.',
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--component-input', type=Path)
    parser.add_argument('--normal-report', type=Path)
    parser.add_argument('--case-index', type=int, default=0)
    parser.add_argument('--universal-pencils', action='store_true')
    parser.add_argument('--translation-nonplane', action='store_true')
    parser.add_argument('--translation-global-graphs', action='store_true')
    args = parser.parse_args()
    if args.translation_global_graphs:
        translation_rank_two_global_graphs(args.output)
    elif args.translation_nonplane:
        translation_rank_two_nonplane(args.output)
    elif args.universal_pencils:
        universal_pencils(args.output)
    elif args.normal_report:
        normal_quadrics(args.normal_report,args.case_index,args.output)
    elif args.component_input:
        components(args.component_input, args.output)
    else:
        main(args.output)
