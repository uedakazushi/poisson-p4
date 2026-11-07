#!/usr/bin/env python3
"""Complete bivectors preserving a smooth projected Veronese surface.

The full280-by160 system imposes q wedge Pi(q)=r(a,b,c)*(dq_a wedge
dq_b wedge dq_c), with every cubic r retained. Only after full Euler24
is checked do we take the projective quotient. No selected tangent ansatz.
"""
from pathlib import Path
import hashlib,itertools,json,resource,signal,time,sys
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,jacobi,minor_certificate
ROOT=Path(__file__).resolve().parents[1]



def dimension():
    """Weak whole-source G bound with finite exact ideal membership lifts."""
    import os,subprocess,tempfile
    from functools import lru_cache
    start=time.monotonic();signal.alarm(115);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    path=ROOT/'reports/CLASSIFICATION_PROJECTED_VERONESE_TANGENT_20261107.json';data=json.loads(path.read_text());zz=s.symbols('z0:'+str(data['complete_projective_tangent_source_dimension']));env={str(v):v for v in zz}
    JC=s.Matrix([[s.Rational(v) for v in row] for row in data['full350_by_parameter_quadric_matrix']]);dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(s.QQ)
    selected=list(dm(JC.T).rref()[1]);gens=[s.Poly(s.sympify(data['all350_projective_Jacobi'][i],locals=env),*zz,domain=s.QQ).clear_denoms()[1].primitive()[1] for i in selected]
    assert len(gens)==54
    script='ring R=0,('+','.join(map(str,zz))+'),dp;\nshort=0;\n'
    script+='ideal I='+','.join(str(g.as_expr()).replace('**','^') for g in gens)+';\nmatrix T; ideal G=liftstd(I,T);\n'
    script+='print("DIM:"+string(dim(G)));\nint i;int j;for(i=1;i<=size(G);i++){print("G:"+string(i)+":"+string(G[i]));}\n'
    script+='for(i=1;i<=nrows(T);i++){for(j=1;j<=ncols(T);j++){if(T[i,j]!=0){print("L:"+string(i)+":"+string(j)+":"+string(T[i,j]));}}}\nquit;\n'
    prefix=ROOT/'work/nonsquare-v62/singular-runtime/prefix';env_run=dict(os.environ)
    env_run['LD_LIBRARY_PATH']='/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu:'+str(prefix/'usr/lib/x86_64-linux-gnu')+':'+env_run.get('LD_LIBRARY_PATH','')
    env_run['SINGULAR_PROCS_DIR']=str(prefix/'usr/lib/x86_64-linux-gnu/singular/MOD')
    env_run['SINGULARPATH']='/workspace/cloud-setup/prefix/usr/share/singular/LIB'
    result={'task_id':'G','input_report':str(path.relative_to(ROOT)),'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'full54_independent_Jacobi_row_indices':selected,'original_ideal_generators':[str(g.as_expr()) for g in gens],'exact_Singular_input':script}
    dest=ROOT/'reports/CLASSIFICATION_PROJECTED_VERONESE_TANGENT_DIMENSION_20261107.json'
    if dest.exists():raise FileExistsError(dest)
    try:
        with tempfile.TemporaryDirectory(prefix='pv-tangent-dim-') as temp:
            task=Path(temp)/'task.sing';task.write_text(script)
            run=subprocess.run([str(prefix/'usr/bin/Singular'),'-q',str(task)],env=env_run,capture_output=True,text=True,timeout=80,check=False)
        assert run.returncode==0 and '?' not in run.stdout,run.stdout[-10000:]+run.stderr
        G={};lifts=[];dim=None
        for line in run.stdout.splitlines():
            if line.startswith('G:'):
                _,j,p=line.split(':',2);G[int(j)-1]=s.Poly(s.sympify(p.replace('^','**'),locals=env),*zz,domain=s.QQ)
            elif line.startswith('L:'):
                _,i,j,p=line.split(':',3);lifts.append((int(i)-1,int(j)-1,s.Poly(s.sympify(p.replace('^','**'),locals=env),*zz,domain=s.QQ)))
            elif line.startswith('DIM:'):dim=int(line.split(':',1)[1])
        assert G and dim is not None
        sums={j:s.Poly(0,*zz,domain=s.QQ) for j in G}
        for i,j,p in lifts:sums[j]+=gens[i]*p
        assert G==sums
        leading=[G[j].monoms(order='grevlex')[0] for j in sorted(G)]
        edges=tuple(sorted(set(sum(1<<i for i,e in enumerate(mon) if e) for mon in leading)))
        @lru_cache(None)
        def cover(edges):
            if not edges:return 0
            edge=min(edges,key=int.bit_count);assert edge
            return 1+min(cover(tuple(e for e in edges if not e&(1<<i))) for i in range(len(zz)) if edge&(1<<i))
        height=cover(edges);bound=len(zz)-height;assert bound==dim
        result.update(all_passed=True,exact_Singular_stdout=run.stdout,Singular_stderr=run.stderr,lifted_polynomials=[str(G[j].as_expr()) for j in sorted(G)],
            exact_ideal_membership_lifts=[{'generator':i,'polynomial':j,'coefficient':str(p.as_expr())} for i,j,p in lifts],all_exact_lift_identities_checked=True,
            leading_monomial_exponents=[list(e) for e in leading],leading_support_minimum_vertex_cover=height,fixed_Poisson_source_dimension_upper=bound,
            smooth_projected_Veroneses_preserving_group_dimension=3,whole_projective_PGL_image_dimension_upper=24+bound-1-3,
            excludes_unknown_G25_to27=24+bound-1-3<25,
            scope='Exact whole fixed-surface Poisson source dimension bound from lifted leading monomials, without rankzero ansatz or fresh ideal-radical assertion. If dimension<=4, its whole PGL saturation cannot be a generic unknown W-dominant component. Global rank2-double-surface reduction is separate geometry; all singular/vertical/non-node complements retained.')
        print('FULL PV TANGENT DIMENSION',bound,'PGL image<=',20+bound,'lift polynomials',len(G),'lifts',len(lifts),flush=True)
    except subprocess.TimeoutExpired as error:
        result.update(all_passed=False,status='finite80second_Singular_limit',partial_stdout=str(error.stdout or '')[-20000:],dimension_bound_proved=None,
            scope='No new dimension theorem obtained. Complete15parameter tangent source and all350Jacobi retained; change mathematical route after at most one additional nonclosing attempt.')
        print('FINITE80SECOND LIMIT; NO DIMENSION PROOF',flush=True)
    result['seconds']=time.monotonic()-start;dest.write_text(json.dumps(result,indent=2)+'\n')

def main():
    start=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    xx=s.symbols('x0:5');a,b,c=s.symbols('a b c');uv=(a,b,c);q=s.Matrix([a*a-b*b,b*b-c*c,a*b,a*c,b*c]);E=s.Matrix(xx)
    triples=list(itertools.combinations(range(5),3));mon=lambda degree:sorted([e for e in itertools.product(range(degree+1),repeat=3) if sum(e)==degree],reverse=True)
    m6,m3=mon(6),mon(3);rows=list(itertools.product(triples,m6));rindex={r:i for i,r in enumerate(rows)};M=s.zeros(280,160)
    for col,(pair,ex) in enumerate(COLUMNS):
        polynomial=s.prod(q[i]**power for i,power in enumerate(ex))
        for tri in triples:
            if set(pair).issubset(tri):
                missing=next(i for i in tri if i not in pair);sign=(-1)**tri.index(missing)
                for e,co in s.Poly(s.expand(sign*q[missing]*polynomial),*uv).terms():M[rindex[tri,e],col]=co
    dq=q.jacobian(uv);tangent3={tri:s.expand(dq.extract(tri,range(3)).det()) for tri in triples}
    for j,e in enumerate(m3):
        monomial=s.prod(v**power for v,power in zip(uv,e))
        for tri in triples:
            for f,co in s.Poly(s.expand(-monomial*tangent3[tri]),*uv).terms():M[rindex[tri,f],150+j]=co
    dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    rank=lambda A:len(dm(A).rref()[1])
    K=dm(M).nullspace().to_Matrix().T;assert M*K==s.zeros(280,K.cols);assert rank(K[:150,:])==K.cols
    def decode(v):
        P=s.zeros(5)
        for coefficient,(ij,ex) in zip(v,COLUMNS):P[ij]+=coefficient*s.prod(x**power for x,power in zip(xx,ex))
        for i,j in PAIRS:P[j,i]=-P[i,j]
        return P.applyfunc(s.expand)
    def vector(P):return s.Matrix([s.Poly(s.expand(P[ij]),*xx).coeff_monomial(ex) for ij,ex in COLUMNS])
    def canonical(P):
        div=s.Matrix([sum(s.diff(P[i,j],xx[i]) for i in range(5)) for j in range(5)])
        return (P-(E*div.T-div*E.T)/5).applyfunc(s.expand)
    V=[]
    for i,j in itertools.product(range(5),repeat=2):
        field=s.zeros(5,1);field[i]=xx[j];V.append(vector(E*field.T-field*E.T))
    V=s.Matrix.hstack(*V);assert rank(V)==24;assert M[:,:150]*V==s.zeros(280,25)
    canonicalK=s.Matrix.hstack(*[vector(canonical(decode(K[:150,i]))) for i in range(K.cols)])
    red,piv=dm(canonicalK.T).rref();B=red.to_Matrix()[:len(piv),:].T
    assert rank(V.row_join(B))==K.cols and B.cols==K.cols-24
    zz=s.symbols('z0:'+str(B.cols));P=decode(B*s.Matrix(zz));print('FULL MATRIX',rank(M),'kernel',K.cols,'projective source',B.cols,flush=True)
    sparse={ij:dict(s.Poly(P[ij],*xx).terms()) for ij in PAIRS if P[ij]!=0};J=jacobi(sparse).applyfunc(s.expand)
    params=[zz[i]*zz[j] for i,j in itertools.combinations_with_replacement(range(len(zz)),2)]
    JC=s.Matrix([[s.Poly(f,*zz).coeff_monomial(m) for m in params] for f in J]);Jrank=rank(JC)
    pf=[]
    for i in range(5):
        j,k,l,m=[t for t in range(5) if t!=i];pf.append((-1)**i*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l]))
    D=s.expand(sum(xx[i]*pf[i] for i in range(5)));Dpoly=s.Poly(D,*xx);Drows=[e for e,co in Dpoly.terms()]
    DC=s.Matrix([[s.Poly(Dpoly.coeff_monomial(e),*zz).coeff_monomial(m) for m in params] for e in Drows]) if Drows else s.zeros(0,len(params))
    contained=rank(JC.col_join(DC))==Jrank;functionals=None
    if contained and DC.rows:
        rowsJ=list(dm(JC.T).rref()[1]);small=JC[rowsJ,:];cols=list(dm(small).rref()[1]);Lsmall=DC[:,cols]*small[:,cols].inv();assert Lsmall*small==DC
        L=s.zeros(DC.rows,350)
        for j,k in enumerate(rowsJ):L[:,k]=Lsmall[:,j]
        assert L*JC==DC;functionals=[{str(j):str(v) for j,v in enumerate(L.row(i)) if v} for i in range(L.rows)]
    print('ALL350 JACOBI rank',Jrank,'actualD row membership',contained,'D monomials',len(Drows),flush=True)
    encode=lambda Q:{str(i)+','+str(j):str(s.expand(Q[i,j])) for i,j in PAIRS if Q[i,j]!=0}
    result={'task_id':'G','all_linear_arithmetic_passed':True,'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),
        'input_sha256':{'code/classification_rank3_mixed_component_review.py':hashlib.sha256((ROOT/'code/classification_rank3_mixed_component_review.py').read_bytes()).hexdigest()},
        'surface_parametrization':list(map(str,q)),'surface_parameters':list(map(str,uv)),'ambient_coordinates':list(map(str,xx)),
        'ambient_coefficient_columns':[{'pair':list(ij),'exponents':list(e)} for ij,e in COLUMNS],
        'restricted_cubic_monomial_exponents':m3,'tangent_three_vector':{str(tri):str(v) for tri,v in tangent3.items()},
        'full280_by160_sparse_matrix':{str(i):{str(j):str(M[i,j]) for j in range(160) if M[i,j]} for i in range(280)},
        'exact_matrix_rank':rank(M),'exact_kernel_dimension':K.cols,'full_kernel_basis_columns':[[str(v) for v in K[:,i]] for i in range(K.cols)],
        'rank_certificate':minor_certificate(M,rank(M)),'full_Euler_vertical_dimension':24,'Euler_is_full_kernel_subspace':True,
        'complete_projective_tangent_source_dimension':B.cols,'canonical_basis':[encode(decode(B[:,i])) for i in range(B.cols)],
        'parameters':list(map(str,zz)),'whole_canonical_source':encode(P),'all350_projective_Jacobi':list(map(str,J)),
        'parameter_quadratic_monomials':list(map(str,params)),'full350_by_parameter_quadric_matrix':[[str(v) for v in row] for row in JC.tolist()],
        'Jacobi_quadric_row_rank':Jrank,'actual_anticanonical_quintic':str(D),'actual_D_coordinate_monomial_exponents':[list(e) for e in Drows],
        'actual_D_coefficient_matrix':[[str(v) for v in row] for row in DC.tolist()],
        'all_D_coefficients_in_full_Jacobi_quadric_row_space':contained,'all_D_from_Jacobi_functionals':functionals,
        'whole_Poisson_source_has_no_rank4':contained,
        'scope':'All global projective bivectors tangent to the fixed smooth projected Veronese, with complete surface bivector cubic, full150 coefficients/Euler24 quotient/all350 Jacobi. If actualD membership is true this excludes all rank4 Poisson tensors preserving such a surface; geometric rank2-double-surface reduction is a separate input. No general irreducible-nonnormal-quintic exhaustion asserted.',
        'seconds':time.monotonic()-start}
    path=ROOT/'reports/CLASSIFICATION_PROJECTED_VERONESE_TANGENT_20261107.json'
    if path.exists():raise FileExistsError(path)
    path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'source_dimension':B.cols,'no_rank4':contained,'seconds':result['seconds']}))

if __name__=='__main__':
    if '--dimension' in sys.argv:dimension()
    else:main()
