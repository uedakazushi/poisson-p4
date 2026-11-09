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


def translation_perazzo_semisimple_confluence(output):
    """A literal full-kappa actual curve to the unchanged Perazzo canonical point."""
    from s61_model import X as SX,jacobi_projective
    if output.exists():raise FileExistsError(str(output))
    resource.setrlimit(resource.RLIMIT_CPU,(80,81))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    begun=time.monotonic();H,X,U,V,W=parent=s.symbols('H X U V W');ampA,ampL,period=s.symbols('ampA ampL period')
    h,x,u,v,w=xx=s.symbols('h x u v w');ep,kap=s.symbols('ep kap');P=s.zeros(5)
    entries={(1,2):2*X*U,(1,3):X*V,(1,4):-X*W,(2,3):period*U*V,(2,4):-ampL*H*V-(period+2)*U*W,(3,4):ampA*H*X+period*V*W}
    for ij,val in entries.items():P[ij]=val;P[ij[::-1]]=-val
    def strict(A,variables):return [s.expand(sum(A[i,l]*s.diff(A[j,jj],variables[l])+A[j,l]*s.diff(A[jj,i],variables[l])+A[jj,l]*s.diff(A[i,j],variables[l])for l in range(5)))for i,j,jj in itertools.combinations(range(5),3)]
    def pf(A):return s.expand(A[1,2]*A[3,4]-A[1,3]*A[2,4]+A[1,4]*A[2,3])
    def enc(A):return {f'{i},{j}':str(s.expand(A[i,j]))for i,j in itertools.combinations(range(5),2)if A[i,j]}
    def mat(A):return [[str(s.factor(v))for v in row]for row in A.tolist()]
    assert all(f==0 for f in strict(P,parent))
    parentF=2*ampA*H*X*U+ampL*H*V**2+2*(period+1)*U*V*W
    assert s.expand(pf(P)-X*parentF)==0
    assert s.expand((H*pf(P)).subs({ampA:0,ampL:0,period:1})-4*H*X*U*V*W)==0
    shear=s.eye(5);shear[2,0]=1/ep;shear[3,1]=1/ep;shear[4,0]=1/ep;shear[4,1]=-1/ep
    assert shear.det()==1;inverse_shear=shear.inv();assert shear*inverse_shear==s.eye(5)
    binding={period:-2+2*kap*ep,ampA:1/ep**2-2*kap/ep,ampL:-2/ep+4*kap}
    substituted=P.subs(binding,simultaneous=True).subs(dict(zip(parent,shear*s.Matrix(xx))),simultaneous=True)
    internal=(ep*inverse_shear*substituted*inverse_shear.T).applyfunc(s.cancel)
    assert all(not s.denom(f).free_symbols for f in internal)
    G=s.Matrix([[1,0,0,0,0],[0,-2,0,0,0],[0,0,0,-2,0],[0,0,s.Rational(1,2),0,0],[0,0,s.Rational(1,2),-1,1]])
    assert G.det()==-2;GI=G.inv();assert G*GI==s.eye(5)
    actual=(G*internal.subs(dict(zip(xx,GI*s.Matrix(xx))),simultaneous=True)*G.T).applyfunc(s.cancel)
    assert all(not s.denom(f).free_symbols for f in actual)
    assert all(s.Poly(f,*xx).total_degree()<=2 for f in actual)
    assert all(f==0 for f in strict(actual,xx))
    assert all(f==0 for f in jacobi_projective(actual.subs(dict(zip(xx,SX)),simultaneous=True)))
    source_path='reports/CLASSIFICATION_TWO_HYPERPLANES_TRANSLATION_NONPLANE_COMPLETION_20261108.json'
    source=next(z for z in json.loads((ROOT/source_path).read_text())['cases']if z['branch']=='L_h_minus_x')
    target=s.zeros(5);point=s.zeros(5)
    for data,model in [(source['canonical_actual_kappa_model'],target),(source['canonical_point_brackets'],point)]:
        for ij,text in data.items():
            i,j=map(int,ij.split(','));model[i,j]=s.sympify(text);model[j,i]=-model[i,j]
    endpoint=actual.subs(ep,0).applyfunc(s.expand)
    assert (endpoint-target).applyfunc(s.expand)==s.zeros(5)
    assert (endpoint.subs(kap,1)-point).applyfunc(s.expand)==s.zeros(5)
    F0=(h-x)*(h*u-x*v)+h*x*w
    F1=h*u*(v+w-u)+2*x*v*(w-v);F2=u*v*(2*w-2*v-u)
    expectedD=h*x*(1-2*kap*ep)*(F0+ep*F1+ep**2*F2)
    assert s.expand(h*pf(actual)-expectedD)==0
    assert s.expand(h*pf(endpoint)-h*x*F0)==0
    expected_pull=ep**2*(X*parentF).subs(binding,simultaneous=True).subs(dict(zip(parent,shear*GI*s.Matrix(xx))),simultaneous=True)*G.det()
    assert s.cancel(pf(actual)-expected_pull)==0
    paths=['code/classification_two_hyperplanes_cubic_nilpotent.py','code/s61_model.py',source_path,
        'research/classification/TWO_HYPERPLANES_TRANSLATION_PERAZZO_SS_CONFLUENCE.md','research/classification/TWO_HYPERPLANES_TRANSLATION_NONPLANE_ORBIT_CLOSURES.md',
        'research/classification/TWO_HYPERPLANES_CUBIC_SUPPORT.md','research/classification/PLANE_CUBIC_TORIC_COLLISION_REVIEW.md',
        'research/classification/inputs/two_hyperplanes_translation_pre_ss_confluence_source_contracts_20261108.json']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,branch='L_h_minus_x',new_component_count=0,
        complete_parent_three_parameter_brackets=enc(P),complete_parent_parameter_order=list(map(str,[ampA,ampL,period])),
        parent_all10_strict_Jacobi_zero=True,parent_actual_quintic=str(H*X*parentF),parent_SNC_specialization=dict(ampA=0,ampL=0,period=1,actual_quintic='4*H*X*U*V*W'),
        parent_irreducible_source_component_bridge='The whole irreducible three-parameter parent image lies in an ambient component containing its actual SNC rank4 point; the accepted full SNC-meeting component theorem therefore covers the whole parent by MPS40. No local uniqueness or known-point reverse inference.',
        actual_parent_parameter_binding={str(k):str(v)for k,v in binding.items()},actual_shear_matrix=mat(shear),actual_shear_inverse=mat(inverse_shear),actual_shear_determinant='1',tensor_scale=str(ep),
        complete_pole_free_internal_curve=enc(internal),constant_target_PGL=mat(G),constant_target_inverse_PGL=mat(GI),constant_target_determinant='-2',
        complete_actual_regular_kappa_curve=enc(actual),all_negative_epsilon_powers_cancelled=True,source_kappa_not_inverted=True,
        full_epsilon_kappa_strict10_Jacobi_zero=True,all350_original_projective_Jacobi_zero=True,
        endpoint_full150_equals_same_kappa_canonical_model=True,endpoint_kappa_one_full150_equals_same_previously_declared_point=True,
        complete_endpoint=enc(endpoint),actual_full_quintic=str(expectedD),all_quintic_coefficients_and_actual_Pfaffian_identity=True,endpoint_quintic=str(h*x*F0),
        endpoint_rank4_for_every_kappa=True,whole_nonplane_seven_source_coverage_requires_independent_orbit_bridge=True,
        whole_source_bridge='Once the supplied13+9upper/22orbitlower proof is independently accepted, the entire Lh-minus-x source saturation is the kappa1orbit closure. The literal curve places that same point in MPS40; PGL stability and closedness cover all full source parameters and all normalization boundaries.',
        source_complements=['Lh essential4 is not covered by this result','Plane A=B0 remains separate','All kappa zero boundary is retained directly'],
        no_new_anchor_or_normal_jet=True,input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,whole_kappa_endpoint=True,whole_source_pending_orbit_review=True,seconds=result['seconds'])),flush=True)


def translation_essential_four_semisimple_confluence(output):
    """A second genuine full-kappa curve, to the same essential-four source."""
    from s61_model import X as SX,jacobi_projective
    if output.exists():raise FileExistsError(str(output))
    resource.setrlimit(resource.RLIMIT_CPU,(80,81))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    begun=time.monotonic();H,X,U,V,W=parent=s.symbols('H X U V W');ampA,ampB,period=s.symbols('ampA ampB period')
    h,x,u,v,w=xx=s.symbols('h x u v w');ep,kap=s.symbols('ep kap');P=s.zeros(5)
    for ij,val in {(1,2):3*X*U,(1,3):X*V,(1,4):-X*W,(2,3):period*U*V,(2,4):-ampB*V**2-(period+3)*U*W,(3,4):ampA*H*X+period*V*W}.items():P[ij]=val;P[ij[::-1]]=-val
    def strict(A,variables):return [s.expand(sum(A[i,l]*s.diff(A[j,jj],variables[l])+A[j,l]*s.diff(A[jj,i],variables[l])+A[jj,l]*s.diff(A[i,j],variables[l])for l in range(5)))for i,j,jj in itertools.combinations(range(5),3)]
    def pf(A):return s.expand(A[1,2]*A[3,4]-A[1,3]*A[2,4]+A[1,4]*A[2,3])
    def enc(A):return {f'{i},{j}':str(s.expand(A[i,j]))for i,j in itertools.combinations(range(5),2)if A[i,j]}
    def mat(A):return [[str(s.factor(v))for v in row]for row in A.tolist()]
    def push(A,G):return (G*A.subs(dict(zip(xx,G.inv()*s.Matrix(xx))),simultaneous=True)*G.T).applyfunc(s.cancel)
    assert all(f==0 for f in strict(P,parent))
    parentF=3*ampA*H*X*U+ampB*V**3+3*(period+1)*U*V*W
    assert s.expand(pf(P)-X*parentF)==0
    assert s.expand((H*pf(P)).subs({ampA:0,ampB:0,period:1})-6*H*X*U*V*W)==0
    shear=s.eye(5);shear[2,0]=1/ep;shear[3,1]=1/ep;shear[4,0]=1/ep
    assert shear.det()==1;inverse_shear=shear.inv();assert shear*inverse_shear==s.eye(5)
    binding={period:-3+6*kap*ep,ampA:2/ep**2-6*kap/ep,ampB:-6*ep}
    substituted=P.subs(binding,simultaneous=True).subs(dict(zip(parent,shear*s.Matrix(xx))),simultaneous=True)
    internal=(ep*inverse_shear*substituted*inverse_shear.T).applyfunc(s.cancel)
    assert all(not s.denom(f).free_symbols for f in internal)
    G=s.Matrix([[1,0,0,0,0],[0,1,0,0,0],[0,0,0,1,0],[0,0,0,0,-1],[0,0,s.Rational(1,6),0,s.Rational(1,2)]])
    assert G.det()==-s.Rational(1,6);GI=G.inv();assert G*GI==s.eye(5)
    actual=push(internal,G)
    assert all(not s.denom(f).free_symbols for f in actual)
    assert all(s.Poly(f,*xx).total_degree()<=2 for f in actual)
    assert all(f==0 for f in strict(actual,xx))
    assert all(f==0 for f in jacobi_projective(actual.subs(dict(zip(xx,SX)),simultaneous=True)))
    source_path='reports/CLASSIFICATION_TWO_HYPERPLANES_TRANSLATION_NONPLANE_COMPLETION_20261108.json'
    source=next(z for z in json.loads((ROOT/source_path).read_text())['cases']if z['branch']=='L_h')
    target=s.zeros(5);point=s.zeros(5)
    for data,model in [(source['canonical_actual_kappa_model'],target),(source['canonical_point_brackets'],point)]:
        for ij,text in data.items():
            i,j=map(int,ij.split(','));model[i,j]=s.sympify(text);model[j,i]=-model[i,j]
    endpoint=actual.subs(ep,0).applyfunc(s.expand);assert (endpoint-target).applyfunc(s.expand)==s.zeros(5)
    fixed=s.Matrix([[s.sympify(z)for z in row]for row in source['canonical_cubic_normalizer']])
    assert fixed.det()==1
    assert (push(endpoint.subs(kap,1),fixed.inv())-point).applyfunc(s.expand)==s.zeros(5)
    T=h*u-x*v
    cubic=(1-3*kap*ep)*(h*T+ep*((3*v+6*w)*T-h*u*v)-ep**2*(3*v+6*w)*u*v)+(x+ep*u)**3
    expectedD=h*x*cubic
    assert s.expand(h*pf(actual)-expectedD)==0
    assert s.expand(h*pf(endpoint)-h*x*(h*h*u-h*x*v+x**3))==0
    pulled=ep**2*(X*parentF).subs(binding,simultaneous=True).subs(dict(zip(parent,shear*GI*s.Matrix(xx))),simultaneous=True)*G.det()
    assert s.cancel(pf(actual)-pulled)==0
    paths=['code/classification_two_hyperplanes_cubic_nilpotent.py','code/s61_model.py',source_path,
        'research/classification/TWO_HYPERPLANES_TRANSLATION_ESSENTIAL_FOUR_SS_CONFLUENCE.md','research/classification/TWO_HYPERPLANES_TRANSLATION_NONPLANE_ORBIT_CLOSURES.md',
        'research/classification/TWO_HYPERPLANES_CUBIC_SUPPORT.md','research/classification/PLANE_CUBIC_TORIC_COLLISION_REVIEW.md',
        'research/classification/inputs/two_hyperplanes_translation_pre_essential_four_source_contracts_20261108.json']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,branch='L_h',new_component_count=0,
        complete_parent_three_parameter_brackets=enc(P),complete_parent_parameter_order=list(map(str,[ampA,ampB,period])),parent_all10_strict_Jacobi_zero=True,
        parent_actual_quintic=str(H*X*parentF),parent_SNC_specialization=dict(ampA=0,ampB=0,period=1,actual_quintic='6*H*X*U*V*W'),
        parent_irreducible_source_component_bridge='The entire irreducible three-parameter parent is contained in a component meeting its actual SNC rank4 point; the full SNC-meeting component theorem identifies that component among MPS40. Not a reverse inference from one known point.',
        actual_parent_parameter_binding={str(k):str(v)for k,v in binding.items()},actual_shear_matrix=mat(shear),actual_shear_inverse=mat(inverse_shear),actual_shear_determinant='1',tensor_scale=str(ep),
        complete_pole_free_internal_curve=enc(internal),constant_target_PGL=mat(G),constant_target_inverse_PGL=mat(GI),constant_target_determinant='-1/6',
        complete_actual_regular_kappa_curve=enc(actual),all_negative_epsilon_powers_cancelled=True,source_kappa_not_inverted=True,
        full_epsilon_kappa_strict10_Jacobi_zero=True,all350_original_projective_Jacobi_zero=True,endpoint_full150_equals_same_kappa_canonical_model=True,
        fixed_prior_canonical_PGL=mat(fixed),fixed_prior_canonical_inverse_PGL=mat(fixed.inv()),same_previously_declared_point_recovered_in_all150_coefficients=True,
        complete_endpoint=enc(endpoint),actual_full_quintic=str(expectedD),all_quintic_coefficients_and_actual_Pfaffian_identity=True,endpoint_quintic=str(h*x*(h*h*u-h*x*v+x**3)),
        endpoint_rank4_for_every_kappa=True,whole_nonplane_eight_source_coverage_requires_independent_orbit_bridge=True,
        whole_source_bridge='Once the13+9upper and22orbitlower bridge is independently accepted, the entire Lh source saturation is the kappa1orbit closure. This literal curve places that same prior point, up to its already specified fixed PGL, in MPS40 and therefore covers every source parameter and inverse normalization chart.',
        source_complements=['Plane A=B0 remains separate','Other nilpotent/cofactor sources are not included','All kappa0 boundary is retained directly'],
        no_new_anchor_or_normal_jet=True,input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,whole_kappa_endpoint=True,same_source_point=True,seconds=result['seconds'])),flush=True)


def j21_h_x_semisimple_confluence(output):
    """Whole original six-parameter j21(h,x) source, actual all-s/k SS arc."""
    from s61_model import X as SX, jacobi_projective
    if output.exists():raise FileExistsError(str(output))
    resource.setrlimit(resource.RLIMIT_CPU,(80,81))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    begun=time.monotonic();h,x,u,v,w=xx=s.symbols('h x u v w');ep,kap,amp=s.symbols('ep kap amp')
    H,X,U,V,W=parent=s.symbols('H X U V W');ampA,ampB,period=s.symbols('ampA ampB period')
    def strict(P,variables):return [s.expand(sum(P[i,l]*s.diff(P[j,k],variables[l])+P[j,l]*s.diff(P[k,i],variables[l])+P[k,l]*s.diff(P[i,j],variables[l])for l in range(5)))for i,j,k in itertools.combinations(range(5),3)]
    def pf(P):return s.expand(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3])
    def enc(P):return {f'{i},{j}':str(s.expand(P[i,j]))for i,j in itertools.combinations(range(5),2)if P[i,j]}
    def mat(P):return [[str(s.factor(z))for z in row]for row in P.tolist()]
    def push(P,G):return (G*P.subs(dict(zip(xx,G.inv()*s.Matrix(xx))),simultaneous=True)*G.T).applyfunc(s.cancel)
    def matrix(entries):
        P=s.zeros(5)
        for ij,z in entries.items():
            if isinstance(ij,str):ij=tuple(map(int,ij.split(',')))
            P[ij]=s.sympify(z);P[ij[::-1]]=-P[ij]
        return P
    source_path='reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json'
    source_data_path='code/data/classification_two_hyperplanes_cubic_nilpotent_sources.json'
    source_data=json.loads((ROOT/source_data_path).read_text());spec=next(q for q in source_data['branches']if q['name']=='j21_h_x_irreducible')
    record=json.loads((ROOT/source_path).read_text())['cases'][spec['chart_index']]
    ns=s.symbols('n0:14');binding={s.Symbol(k):s.sympify(z)for k,z in spec['substitution'].items()}
    original=matrix(record['brackets']).subs(binding,simultaneous=True).applyfunc(s.expand)
    assert all(s.expand(s.sympify(z).subs(binding,simultaneous=True))==0 for z in record['remaining_coefficient_equations'])
    assert all(z==0 for z in strict(original,xx))
    normalizer=s.eye(5);normalizer[2,4]=-ns[12];normalizer[2,1]=-ns[12]*ns[8]-ns[7];normalizer[2,0]=-ns[0]-ns[10]*ns[12]
    normalizer[3,1]=-ns[12];normalizer[4,0]=-ns[10];normalizer[4,1]=-ns[8]
    assert normalizer.det()==1;inverse=normalizer.inv();assert normalizer*inverse==s.eye(5)
    target=matrix({(1,2):x*v,(1,3):h*x,(1,4):x*x,(2,3):v*v/2-h*u,(2,4):amp*h*h+kap*x*v-v*w,(3,4):kap*h*x-h*w})
    normalized=push(original,normalizer)
    assert (normalized-target.subs({amp:ns[5],kap:ns[8]},simultaneous=True)).applyfunc(s.expand)==s.zeros(5)
    assert (push(normalized,inverse)-original).applyfunc(s.expand)==s.zeros(5)
    assert all(z==0 for z in strict(target,xx))
    canonicalF=-amp*h**3+x*(v*v/2-h*u);assert s.expand(pf(target)-x*canonicalF)==0
    vf=s.Matrix([0,0,v,h,x]);assert (normalizer*vf-vf.subs(dict(zip(xx,normalizer*s.Matrix(xx))),simultaneous=True)).applyfunc(s.expand)==s.zeros(5,1)
    scale=s.diag(1,1/amp,1,1,1/amp);assert (push(target,scale)-target.subs(amp,1)).applyfunc(s.cancel)==s.zeros(5)
    t=s.Symbol('t');torus=s.diag(1,t**-2,t**2,t,t**-1)
    assert (push(target,torus)/t-target.subs(kap,t*kap)).applyfunc(s.cancel)==s.zeros(5)
    P=matrix({(1,2):3*X*U,(1,3):X*V,(1,4):-X*W,(2,3):period*U*V,(2,4):-ampB*V*V-period*U*W,(3,4):ampA*H*H+period*V*W})
    assert all(z==0 for z in strict(P,parent))
    parentF=3*ampA*H*H*U+ampB*V**3+3*period*U*V*W
    assert s.expand(pf(P)-X*parentF)==0
    assert s.expand((H*pf(P)).subs({ampA:0,ampB:0,period:1})-3*H*X*U*V*W)==0
    shear=s.Matrix([[1,0,0,0,0],[0,1,0,0,0],[0,1/(3*ep),0,0,1],[1/(2*ep**2),0,1,1/(2*ep),0],[-1/ep,0,0,1,0]])
    assert shear.det()==1;shear_inverse=shear.inv();assert shear*shear_inverse==s.eye(5)
    curve_binding={period:1-3*kap*ep,ampA:(1-3*kap*ep)/(2*ep**3),ampB:-8*amp*ep**4}
    substituted=P.subs(curve_binding,simultaneous=True).subs(dict(zip(parent,shear*s.Matrix(xx))),simultaneous=True)
    actual=(ep*shear_inverse*substituted*shear_inverse.T).applyfunc(s.cancel)
    assert all(not s.denom(z).free_symbols for z in actual)
    assert all(s.Poly(z,*xx).total_degree()<=2 for z in actual)
    assert all(z==0 for z in strict(actual,xx))
    assert all(z==0 for z in jacobi_projective(actual.subs(dict(zip(xx,SX)),simultaneous=True)))
    endpoint=actual.subs(ep,0).applyfunc(s.expand);assert (endpoint-target).applyfunc(s.expand)==s.zeros(5)
    pulled=ep**2*(H*X*parentF).subs(curve_binding,simultaneous=True).subs(dict(zip(parent,shear*s.Matrix(xx))),simultaneous=True)
    assert s.cancel(h*pf(actual)-pulled)==0
    assert s.expand(h*pf(endpoint)-h*x*canonicalF)==0
    assert (push(endpoint.subs({amp:ns[5],kap:ns[8]},simultaneous=True),inverse)-original).applyfunc(s.expand)==s.zeros(5)
    paths=['code/classification_two_hyperplanes_cubic_nilpotent.py','code/s61_model.py',source_path,source_data_path,
        'research/classification/TWO_HYPERPLANES_J21_H_X_SS_CONFLUENCE.md','research/classification/TWO_HYPERPLANES_CUBIC_SUPPORT.md',
        'research/classification/PLANE_CUBIC_TORIC_COLLISION_REVIEW.md','research/classification/inputs/two_hyperplanes_pre_j21_hx_ss_source_contracts_20261108.json']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,new_component_count=0,source_branch=spec['name'],original_six_free_parameters=spec['free_parameters'],
        complete_original_binding=spec['substitution'],all_original_remaining_Jacobi_zero=True,complete_original_full_source=enc(original),
        complete_canonical_two_parameter_model=enc(target),canonical_actual_cubic=str(canonicalF),canonical_actual_quintic=str(h*x*canonicalF),
        actual_source_normalizer=mat(normalizer),actual_source_inverse_normalizer=mat(inverse),source_normalizer_determinant='1',source_normalizer_has_no_parameter_denominator=True,
        whole_source_forward_and_inverse_all150_equalities=True,cofactor_preserved_exactly=True,amplitude_scaling=mat(scale),amplitude_scale_requires_only_amp_nonzero=True,
        kappa_nonzero_torus=mat(torus),torus_tensor_scale='1/t',kappa_zero_already_covered_directly=True,
        whole_parent_three_parameter_brackets=enc(P),parent_all10_Jacobi_zero=True,parent_actual_quintic=str(H*X*parentF),
        parent_SNC_specialization={'ampA':0,'ampB':0,'period':1,'actual_quintic':'3*H*X*U*V*W'},
        whole_parent_component_bridge='Irreducible actual A3 parent image is contained in an ambient component containing its SNC point. The full SNC-meeting component theorem identifies that component among MPS40, covering the entire parent.',
        actual_parent_parameter_binding={str(k):str(z)for k,z in curve_binding.items()},actual_coordinate_map=mat(shear),actual_inverse_coordinate_map=mat(shear_inverse),actual_coordinate_determinant='1',tensor_scale='ep',
        complete_pole_free_actual_curve=enc(actual),all_negative_epsilon_powers_cancelled=True,all350_original_projective_Jacobi_zero=True,full_epsilon_amp_kappa_strict10_Jacobi_zero=True,
        all150_endpoint_equals_complete_canonical_model=True,all150_inverse_endpoint_equals_original_whole_six_source=True,
        actual_full_curve_quintic=str(s.factor(h*pf(actual))),full_actual_Pfaffian_matches_parent_pullback=True,complete_endpoint=enc(endpoint),
        amp_kappa_never_inverted_by_actual_curve=True,whole_source_coverage='All original six-source parameters, including every zero magnetic/modulus coefficient, are directly endpoints of genuine MPS-covered parent points after the invertible full-source normalizer. The entire PGL saturation closure is therefore covered; no orbit-dimension or local-uniqueness inference is needed.',
        complements=['amp0 lies outside the irreducible-cubic locus, but the same actual curve retains it','j21(h+x,h) p0, j21(0,h+x), plane residuals and other nilpotent sources remain separate'],
        no_new_anchor_or_normal_jet=True,independent_review_required=True,input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,whole_six_source_direct_actual_curve=True,seconds=result['seconds'])),flush=True)


def j21_zero_hx_modular_confluence(output):
    """Whole u-dependent j21(0,h+x) source via an actual semisimple modular arc."""
    from s61_model import X as SX,jacobi_projective
    if output.exists():raise FileExistsError(str(output))
    resource.setrlimit(resource.RLIMIT_CPU,(80,81));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    begun=time.monotonic();h,x,u,v,w=xx=s.symbols('h x u v w');ep,kap=s.symbols('ep kap');S=h+x;T=S+ep*w;K=-T*u+v*w
    def matrix(entries):
        P=s.zeros(5)
        for ij,z in entries.items():
            if isinstance(ij,str):ij=tuple(map(int,ij.split(',')))
            P[ij]=s.sympify(z);P[ij[::-1]]=-P[ij]
        return P
    def strict(P):return [s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l])for l in range(5)))for i,j,k in itertools.combinations(range(5),3)]
    def pf(P):return s.expand(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3])
    def enc(P):return {f'{i},{j}':str(s.expand(P[i,j]))for i,j in itertools.combinations(range(5),2)if P[i,j]}
    def mat(P):return [[str(s.factor(z))for z in row]for row in P.tolist()]
    def push(P,G):return (G*P.subs(dict(zip(xx,G.inv()*s.Matrix(xx))),simultaneous=True)*G.T).applyfunc(s.cancel)
    def lie(P,L):return (L.jacobian(xx)*P+P*L.jacobian(xx).T-sum((L[i]*P.diff(xx[i])for i in range(5)),s.zeros(5))).applyfunc(s.expand)
    source_path='reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json'
    complete_path='reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_TRANSLATION_COMPONENT_20261107.json'
    spec=json.loads((ROOT/complete_path).read_text())['branches'][0]['full_source_specification']
    record=json.loads((ROOT/source_path).read_text())['cases'][spec['chart_index']]
    ns=s.symbols('n0:19');binding={s.Symbol(k):s.sympify(z)for k,z in spec['substitution'].items()};raw=matrix(record['brackets'])
    original=raw.subs(binding,simultaneous=True).applyfunc(s.expand)
    assert all(s.expand(s.sympify(z).subs(binding,simultaneous=True))==0 for z in record['remaining_coefficient_equations'])
    assert all(z==0 for z in strict(original))
    normalizer=s.eye(5);normalizer[2,0]=-ns[0];normalizer[2,1]=-ns[3];normalizer[4,0]=ns[16]+ns[2];normalizer[4,1]=ns[16]+ns[4];normalizer[4,3]=ns[5]
    inverse=normalizer.inv();assert normalizer.det()==1 and normalizer*inverse==s.eye(5)
    Q=matrix({(1,2):x*v,(1,3):ep*x*v,(1,4):x*T,(2,3):K,(2,4):-K,(3,4):x*v})
    assert all(z==0 for z in strict(Q))
    field=s.Matrix([0,0,v,ep*v,T]);modular=s.Matrix([sum(s.diff(Q[j,i],xx[j])for j in range(1,5))for i in range(5)])
    expected_modular=s.Matrix([0,-2*ep*x,2*v-w-ep*u,ep*v-T,2*h+4*x+2*ep*w])
    assert modular==expected_modular
    assert Q[:,1]==-x*field
    assert lie(Q,field)==s.zeros(5) and lie(Q,modular)==s.zeros(5)
    assert (field.jacobian(xx)*modular-modular.jacobian(xx)*field).applyfunc(s.expand)==s.zeros(5,1)
    actual=(Q+kap*(field*modular.T-modular*field.T)).applyfunc(s.expand)
    assert all(z==0 for z in strict(actual))
    assert all(z==0 for z in jacobi_projective(actual.subs(dict(zip(xx,SX)),simultaneous=True)))
    endpoint=actual.subs(ep,0);normalized=push(original,normalizer)
    assert (normalized-endpoint.subs(kap,ns[16])).applyfunc(s.expand)==s.zeros(5)
    assert (push(normalized,inverse)-original).applyfunc(s.expand)==s.zeros(5)
    assert (push(endpoint.subs(kap,ns[16]),inverse)-original).applyfunc(s.expand)==s.zeros(5)
    baseF=x*v*v+(T+ep*v)*K;c=1+2*kap*ep;F0=-S*S*u+S*v*w+x*v*v
    assert s.expand(pf(Q)-x*baseF)==0
    assert s.expand(pf(actual)-c*x*baseF)==0
    assert s.expand(pf(endpoint)-x*F0)==0
    assert (actual[1,:]-c*x*field.T).applyfunc(s.expand)==s.zeros(1,5)
    semisimple=s.eye(5);semisimple[2,3]=-1/ep;semisimple[4,0]=1/ep;semisimple[4,1]=1/ep
    assert semisimple.det()==1;semisimple_inverse=semisimple.inv();assert semisimple*semisimple_inverse==s.eye(5)
    newfield=(semisimple*field.subs(dict(zip(xx,semisimple_inverse*s.Matrix(xx))),simultaneous=True)).applyfunc(s.cancel)
    assert newfield==s.Matrix([0,0,0,ep*v,ep*w])
    semisimple_tensor=(push(actual,semisimple)/(ep*c)).applyfunc(s.cancel)
    assert semisimple_tensor[1,2]==0 and semisimple_tensor[1,3]==x*v and semisimple_tensor[1,4]==x*w
    assert (push(semisimple_tensor*ep*c,semisimple_inverse)-actual).applyfunc(s.cancel)==s.zeros(5)
    amplitude=s.Symbol('amplitude');scale=s.diag(1,1,amplitude,amplitude,1)
    rawF=s.expand(pf(raw)/x);scaledF=s.expand(amplitude**2*rawF.subs({u:u/amplitude,v:v/amplitude},simultaneous=True))
    assert s.expand(s.diff(scaledF,u)+amplitude*ns[6]*S**2)==0
    paths=['code/classification_two_hyperplanes_cubic_nilpotent.py','code/s61_model.py',source_path,complete_path,
        'research/classification/TWO_HYPERPLANES_J21_ZERO_HX_MODULAR_CONFLUENCE.md','research/classification/TWO_HYPERPLANES_CUBIC_SUPPORT.md',
        'research/classification/PLANE_CUBIC_TORIC_COLLISION_REVIEW.md','research/classification/inputs/two_hyperplanes_pre_j21_zero_hx_source_contracts_20261108.json']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,new_component_count=0,source_branch=spec['name'],complete_original_binding=spec['substitution'],whole_six_free_parameters=spec['free_parameters'],
        complete_original_source=enc(original),all_original_Jacobi_equations_universally_zero=True,
        n6_normalization_scale=mat(scale),n6_coefficient_exactly_rescales_by_amplitude=True,n6_zero_complement_kept_separate=True,
        actual_source_normalizer=mat(normalizer),actual_source_inverse=mat(inverse),normalizer_determinant='1',normalizer_has_no_source_parameter_denominator=True,
        all150_forward_and_inverse_whole_source_equalities=True,complete_canonical_kappa_endpoint=enc(endpoint),endpoint_actual_cubic=str(F0),endpoint_actual_quintic=str(h*x*F0),
        actual_base_curve=enc(Q),full_epsilon_base_strict10_Jacobi_zero=True,actual_base_quintic=str(h*x*baseF),
        Hamiltonian_field=list(map(str,field)),actual_affine_modular_field=list(map(str,modular)),Hamiltonian_and_modular_are_commuting_Poisson_fields=True,
        actual_twist='Q_e+kappa*V_e wedge Z_e',complete_actual_polynomial_kappa_epsilon_curve=enc(actual),all350_original_projective_Jacobi_zero=True,all_actual_strict10_Jacobi_zero=True,
        actual_full_curve_quintic=str(h*c*x*baseF),actual_Pfaffian_equals_literal_c_times_base=True,
        complete_actual_semisimple_coordinate_map=mat(semisimple),complete_actual_semisimple_inverse=mat(semisimple_inverse),semisimple_determinant='1',
        exact_cofactor_scalar=str(ep*c),cofactor_spatial_eigenvalues=['0','1','1'],zero_eigen_translation_is_zero=True,complete_semisimple_tensor=enc(semisimple_tensor),full150_inverse_semisimple_equalities=True,
        source_and_kappa_never_inverted_in_actual_polynomial_curve=True,parameter_complements=['n6=0 not covered by n6 normalization','e=0 and 1+2*kappa*e=0 are excluded only from the semisimple coordinate chart; endpoint directly retained for allkappa'],
        whole_membership_proof='For every fixed kappa the polynomial curve is rank4 with irreducible cubic at e0, hence also for generic e. On e*(1+2kappa*e) nonzero it lies in the accepted full nonzero semisimple diag(1,1,0) source theorem, hence MPS40. Closedness gives allkappa endpoints; the actual no-denominator inverse normalizer and n6 nonzero rescaling cover the entire original source/PGL saturation closure.',
        complements=['j21(0,h+x),n6=0','j21(h+x,h),p0','plane1+1+3 residualthree','other nilpotent sources and totalE'],
        no_new_anchor_or_normal_jet=True,independent_review_required=True,input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,whole_six_source_modular_curve=True,seconds=result['seconds'])),flush=True)


def p0_t06_moving_pair_confluence(output):
    """One whole-modulus actual T06 cubic-Aut frame, complete bounded period class."""
    from sympy.polys.matrices import DomainMatrix
    from s61_model import X as SX,jacobi_projective
    if output.exists():raise FileExistsError(str(output))
    resource.setrlimit(resource.RLIMIT_CPU,(80,81));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    begun=time.monotonic();h,x,u,v,w=xx=s.symbols('h x u v w');yy=s.symbols('x0:5');ep,lam=s.symbols('ep lam');p0,p1,t0,t1=params=s.symbols('p0 p1 t0 t1')
    pairs=list(itertools.combinations(range(5),2));monomials=[xx[i]*xx[j]for i,j in itertools.combinations_with_replacement(range(5),2)]
    def dm(A):return DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    def matrix(entries):
        P=s.zeros(5)
        for ij,z in entries.items():
            if isinstance(ij,str):ij=tuple(map(int,ij.split(',')))
            P[ij]=s.sympify(z);P[ij[::-1]]=-P[ij]
        return P
    def strict(P,variables):return [s.expand(sum(P[i,l]*s.diff(P[j,k],variables[l])+P[j,l]*s.diff(P[k,i],variables[l])+P[k,l]*s.diff(P[i,j],variables[l])for l in range(5)))for i,j,k in itertools.combinations(range(5),3)]
    def canon(P,variables):
        E=s.Matrix(variables);D=s.Matrix([sum(s.diff(P[j,i],variables[j])for j in range(5))for i in range(5)])
        return (P-(E*D.T-D*E.T)/5).applyfunc(s.expand)
    def vec(P):return s.Matrix([s.Poly(P[i,j],*xx).coeff_monomial(m)for i,j in pairs for m in monomials])
    def enc(P):return {f'{i},{j}':str(s.expand(P[i,j]))for i,j in pairs if P[i,j]}
    def mat(P):return [[str(s.factor(z))for z in row]for row in P.tolist()]
    def pf(P,variables):
        out=0
        for i in range(5):
            a,b,c,d=[j for j in range(5)if j!=i]
            out+=(-1)**i*variables[i]*(P[a,b]*P[c,d]-P[a,c]*P[b,d]+P[a,d]*P[b,c])
        return s.expand(out)
    def push(P,G):return (G*P.subs(dict(zip(xx,G.inv()*s.Matrix(xx))),simultaneous=True)*G.T).applyfunc(s.cancel)
    source_path='reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json';data_path='code/data/classification_two_hyperplanes_cubic_nilpotent_sources.json'
    spec=next(q for q in json.loads((ROOT/data_path).read_text())['branches']if q['name']=='j21_hx_h_p0');record=json.loads((ROOT/source_path).read_text())['cases'][spec['chart_index']]
    ns=s.symbols('n0:14');bind={s.Symbol(k):s.sympify(z)for k,z in spec['substitution'].items()};original=matrix(record['brackets']).subs(bind,simultaneous=True).applyfunc(s.expand)
    normalizer=s.eye(5);normalizer[2,0]=ns[10]*ns[3]+ns[5]+ns[10]**2/2;normalizer[2,1]=-ns[2]-ns[3]*ns[8];normalizer[4,0]=-ns[10];normalizer[4,1]=-ns[8]
    assert normalizer.det()==1
    target=matrix({(1,2):x*v,(1,3):x*(h+x),(1,4):h*x,(2,3):h*u-v*v/2+lam*(-h*v+(h+x)*w),(2,4):h*u-w*w/2+lam*h*(w-v),(3,4):h*(v-w)-lam*h*x})
    assert all(z==0 for z in strict(original,xx))
    assert (push(original,normalizer)-target.subs(lam,ns[3])).applyfunc(s.expand)==s.zeros(5)
    assert (push(target.subs(lam,ns[3]),normalizer.inv())-original).applyfunc(s.expand)==s.zeros(5)
    targetC=canon(target,xx);targetvec=vec(targetC);rhs0=targetvec.subs(lam,0);rhs1=targetvec.diff(lam);assert targetvec==rhs0+lam*rhs1
    source_T06='reports/MPS_S61_RAW_CATALOGUE_20261103_v62/T06.json';T06=json.loads((ROOT/source_T06).read_text());P=matrix(T06['homogeneous_bivector'])
    assert all(z==0 for z in jacobi_projective(P))
    F0=(h*(v-w)**2+x*w*w-2*h*x*u)/2;assert s.expand(pf(target,xx)-h*x*F0)==0
    G=s.Matrix([[-1/ep**2,-1/ep**2,-2,-2/ep,0],[0,1,0,0,0],[1,0,0,0,0],[1/ep,0,0,0,1],[0,1/ep,0,1,-1]])
    GI=G.inv();assert G.det()!=0 and G*GI==s.eye(5)
    cubic=yy[0]*yy[1]*yy[2]+yy[2]*yy[4]**2+yy[1]*yy[3]**2
    assert s.cancel(cubic.subs(dict(zip(yy,G*s.Matrix(xx))),simultaneous=True)-2*F0)==0
    parentbasis=[canon(P.diff(z),yy).subs(dict(zip(yy,xx)),simultaneous=True)for z in params]
    parentmatrix=s.Matrix.hstack(*[vec(z)for z in parentbasis]);assert len(dm(parentmatrix).rref()[1])==4
    images=[]
    for z in params:
        pulled=(GI*P.diff(z).subs(dict(zip(yy,G*s.Matrix(xx))),simultaneous=True)*GI.T).applyfunc(s.cancel)
        images.append(canon(pulled,xx))
    bound=8;profiles=list(itertools.product(range(4),range(-bound,bound+1)));rowdata={}
    for j,(i,k)in enumerate(profiles):
        for a,f in enumerate(vec(images[i])):
            for exponent,coefficient in s.Poly(s.expand(ep**bound*f),ep).terms():
                power=exponent[0]-bound+k
                if coefficient and power<=0:rowdata.setdefault((power,a),{})[j]=coefficient
    for a in range(150):rowdata.setdefault((0,a),{})
    keys=sorted(rowdata);A=s.zeros(len(keys),len(profiles));B=s.zeros(len(keys),2)
    for row,key in enumerate(keys):
        for col,z in rowdata[key].items():A[row,col]=z
        if key[0]==0:B[row,0]=rhs0[key[1]];B[row,1]=rhs1[key[1]]
    piv=list(dm(A).rref()[1]);rows=list(dm(A[:,piv].T).rref()[1]);block=A.extract(rows,piv);BI=dm(block).inv().to_Matrix()
    assert block*BI==s.eye(len(piv));assert A[:,piv]*BI*A[rows,:]==A
    Y=s.zeros(len(profiles),2);YP=BI*B[rows,:]
    for j,col in enumerate(piv):Y[col,:]=YP[j,:]
    residual=B-A*Y;gcd=s.Poly(0,lam)
    for row in range(residual.rows):gcd=s.gcd(gcd,s.Poly(residual[row,0]+lam*residual[row,1],lam))
    whole=all(z==0 for z in residual);conditions=[str(residual[i,0]+lam*residual[i,1])for i in range(residual.rows)if any(residual[i,:])]
    functionals=[];functional_gcd=s.Poly(0,lam)
    for row in range(residual.rows):
        if not any(residual[row,:]):continue
        functional=s.zeros(1,A.rows);functional[0,row]=1;coeff=A[row,piv]*BI
        for j,a in enumerate(rows):functional[0,a]-=coeff[j]
        assert functional*A==s.zeros(1,A.cols) and functional*B==residual[row,:]
        functionals.append(dict(row=row,key=keys[row],functional={str(i):str(z)for i,z in enumerate(functional)if z},rhs=list(map(str,residual[row,:]))))
        functional_gcd=s.gcd(functional_gcd,s.Poly(residual[row,0]+lam*residual[row,1],lam))
        if functional_gcd.monic()==gcd.monic():break
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,source_branch=spec['name'],whole_source_parameter_order=spec['free_parameters'],
        source_normalizer=mat(normalizer),source_inverse_normalizer=mat(normalizer.inv()),source_normalizer_has_no_parameter_denominator=True,whole_original_source=enc(original),complete_modulus_model=enc(target),whole_normalizer_all150_forward_inverse=True,
        actual_cubic=str(F0),actual_T06_input=source_T06,source_four_parameter_order=list(map(str,params)),whole_T06_all350_Jacobi_zero=True,
        actual_parent_coordinate_map=mat(G),actual_parent_inverse=mat(GI),actual_parent_determinant=str(G.det()),exact_common_cubic_identity=True,
        canonical_Euler_quotient_definition='P - E wedge div(P)/5',parent_quotient_parameter_rank=4,prior_uniform_pole_bound=8,complete_Laurent_profile_powers=[-8,8],
        period_unknown_order=[[str(params[i]),k]for i,k in profiles],complete_matching_equation_order=[list(z)for z in keys],matching_matrix_shape=list(A.shape),matching_matrix_sparse={f'{i},{j}':str(A[i,j])for i,j in A.todok()},
        matching_rhs_sparse={f'{i},{j}':str(B[i,j])for i,j in B.todok()},matching_rank=len(piv),independent_columns=piv,independent_rows=rows,nonzero_block_inverse=mat(BI),
        full_column_span_identity=True,full_solution_sparse={f'{i},{j}':str(Y[i,j])for i,j in Y.todok()},complete_cokernel_residual_sparse={f'{i},{j}':str(residual[i,j])for i,j in residual.todok()},
        complete_modulus_compatibility_conditions=conditions,compatibility_gcd=str(gcd.monic().as_expr())if gcd else'0',explicit_cokernel_functionals=functionals,
        whole_modulus_actual_T06_confluence_found=whole,negative_scope='If nonzero residual, only this one actual cubic-Aut moving-pair frame class with all meromorphic profiles is excluded at the stated modulus values; not the entire T06 component or other parents.',
        new_component_count=0,no_new_anchor_normal_or_jet=True)
    if whole:
        profiles_actual=[s.expand(sum(ep**k*(Y[j,0]+lam*Y[j,1])for j,(a,k)in enumerate(profiles)if a==i))for i in range(4)]
        actual=sum((profiles_actual[i]*images[i]for i in range(4)),s.zeros(5)).applyfunc(s.cancel)
        assert all(not s.denom(z).free_symbols for z in actual);assert (actual.subs(ep,0)-targetC).applyfunc(s.expand)==s.zeros(5)
        assert all(z==0 for z in jacobi_projective(actual.subs(dict(zip(xx,SX)),simultaneous=True)))
        result.update(actual_parent_parameter_profiles=dict(zip(map(str,params),map(str,profiles_actual))),complete_actual_regular_curve=enc(actual),all_negative_powers_cancelled=True,all150_modulus_endpoint_equalities=True,all350_actual_curve_Jacobi_zero=True,actual_full_curve_quintic=str(s.factor(pf(actual,xx))),endpoint_actual_quintic=str(h*x*F0),independent_review_required=True)
    paths=['code/classification_two_hyperplanes_cubic_nilpotent.py','code/s61_model.py',source_path,data_path,source_T06,'research/classification/TWO_HYPERPLANES_CUBIC_SUPPORT.md','research/classification/TWO_HYPERPLANES_P0_T06_MOVING_MARKED_PAIR.md','research/classification/inputs/two_hyperplanes_pre_p0_t06_source_contracts_20261108.json']
    result.update(input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,whole_modulus_confluence=whole,modulus_condition=result['compatibility_gcd'],matrix_shape=list(A.shape),rank=len(piv),seconds=result['seconds'])),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--component-input', type=Path)
    parser.add_argument('--normal-report', type=Path)
    parser.add_argument('--case-index', type=int, default=0)
    parser.add_argument('--universal-pencils', action='store_true')
    parser.add_argument('--translation-nonplane', action='store_true')
    parser.add_argument('--translation-global-graphs', action='store_true')
    parser.add_argument('--translation-perazzo-ss', action='store_true')
    parser.add_argument('--translation-essential-four-ss', action='store_true')
    parser.add_argument('--j21-h-x-ss', action='store_true')
    parser.add_argument('--j21-zero-hx-modular', action='store_true')
    parser.add_argument('--p0-t06-moving-pair', action='store_true')
    args = parser.parse_args()
    if args.p0_t06_moving_pair:
        p0_t06_moving_pair_confluence(args.output)
    elif args.j21_zero_hx_modular:
        j21_zero_hx_modular_confluence(args.output)
    elif args.j21_h_x_ss:
        j21_h_x_semisimple_confluence(args.output)
    elif args.translation_essential_four_ss:
        translation_essential_four_semisimple_confluence(args.output)
    elif args.translation_perazzo_ss:
        translation_perazzo_semisimple_confluence(args.output)
    elif args.translation_global_graphs:
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
