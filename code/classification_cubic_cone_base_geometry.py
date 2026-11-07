#!/usr/bin/env python3
"""Split the new cone source into strict base Poisson and a linear vertex fiber."""
from pathlib import Path
from itertools import combinations,combinations_with_replacement
from fractions import Fraction as Q
import hashlib,json,resource,signal,time,sys
from sympy.polys.matrices import DomainMatrix
import sympy as s
from classification_projected_veronese_zero_review import sparse_rref,plus,scale
ROOT=Path(__file__).resolve().parents[1];PAIRS=list(combinations(range(5),2))


def main():
    destination=ROOT/'reports/CLASSIFICATION_CUBIC_CONE_BASE_GEOMETRY_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    path=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_CONORMAL_20261107.json';source=json.loads(path.read_text());case=source['cases'][1]
    B=[[Q(c)for c in row]for row in case['complete_original30_by_reduced_source_matrix']];n=len(B[0]);assert n==18
    vertex=[j for j in range(n)if all(4 in PAIRS[i%10]for i in range(30)if B[i][j])];base=[j for j in range(n)if j not in vertex];assert len(base)==6 and len(vertex)==12
    oldcols=[tuple(v)for v in case['all_reduced_source_quadratic_monomials']];newcols=list(combinations_with_replacement(range(6),2));ci={v:i for i,v in enumerate(newcols)}
    base_rows=[];fiber_rows=[]
    for row in case['all350_reduced_Jacobi_rows']:
        quadratic={};linear={}
        for k,c in row.items():
            pair=oldcols[int(k)];inbase=[j in base for j in pair]
            assert any(inbase)
            if all(inbase):quadratic[ci[tuple(sorted(base.index(j)for j in pair))]]=Q(c)
            else:
                b=next(j for j in pair if j in base);v=next(j for j in pair if j in vertex);linear[base.index(b),vertex.index(v)]=Q(c)
        base_rows.append(quadratic);fiber_rows.append(linear)
    # Different c powers separate pure-base and base/vertex equations.
    assert all(not a or not b for a,b in zip(base_rows,fiber_rows))
    rb,proof,chosen=sparse_rref(base_rows,True);b=s.symbols('b0:6');vv=s.symbols('v0:12')
    polynomials=[sum(s.Rational(c)*b[newcols[j][0]]*b[newcols[j][1]]for j,c in row.items())for row in rb.values()]
    x=s.symbols('x0:4');env={f'x{i}':x[i]for i in range(4)};q=[s.sympify(t,locals=env)for t in case['all_three_ideal_generators']];P=s.zeros(4)
    for i in range(30):
        u,v=PAIRS[i%10]
        if v==4:continue
        P[u,v]+=q[i//10]*sum(s.Rational(B[i][j])*b[k]for k,j in enumerate(base))
    for i,j in combinations(range(4),2):P[j,i]=-P[i,j]
    pf=s.Poly(s.expand(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2]),*x);cert=[];allpf=True
    for exponent,coefficient in pf.terms():
        target={ci[pair]:Q(s.Poly(coefficient,*b).coeff_monomial(b[pair[0]]*b[pair[1]]))for pair in newcols};target={j:c for j,c in target.items()if c};residual=dict(target);w={}
        for p in sorted(rb):
            if p in residual:
                a=residual[p];residual=plus(residual,scale(rb[p],-a));w=plus(w,scale(proof[p],a))
        allpf &=not residual;cert.append({'coordinate_exponents':exponent,'full350_baseJac_functional':w,'residual':residual})
    matrix=[[sum(s.Rational(row.get((k,j),0))*b[k]for k in range(6))for j in range(12)]for row in fiber_rows]
    out={'task_id':'G','input_report':str(path.relative_to(ROOT)),'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'base_six_original18_indices':base,'vertex_twelve_original18_indices':vertex,'all350_base_quadratic_rows':base_rows,'all350_linear_vertex_fiber_matrix':[[str(a)for a in row]for row in matrix],
        'base_quadratic_row_rank':len(rb),'base_rowbasis_polynomials':[str(s.factor(p))for p in polynomials],'base_rowbasis_original350_functionals':proof,
        'all_strict_base_Pfaffian_coefficients_in_baseJac_rowspace':allpf,'base_Pfaffian_certificates':cert,'complete_base_bivector':{str(i)+','+str(j):str(s.expand(P[i,j]))for i,j in combinations(range(4),2)},
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'seconds':time.monotonic()-started,
        'scope':'Exact decomposition of the complete18 abelian-conormal cone source. Base0 gives a12dim rank2 linear space; rank4 requires nonzero strict base Poisson and a vertex vector. No whole cone Gbound is inferred until all nonzero base strata and full linear vertex kernels are bounded.'}
    destination.write_text(json.dumps(out,indent=2,default=str)+'\n')
    print(json.dumps({'base_rank':len(rb),'base_Pf_in_Jac':allpf,'base_polynomials':[str(s.factor(p))for p in polynomials],'seconds':out['seconds']}))


def orbits():
    destination=ROOT/'reports/CLASSIFICATION_CUBIC_CONE_BASE_ORBITS_V2_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    path=ROOT/'reports/CLASSIFICATION_CUBIC_CONE_BASE_GEOMETRY_20261107.json';old=json.loads(path.read_text());x=s.symbols('x0:4');b=s.symbols('b0:6');env={str(v):v for v in (*x,*b)}
    P=s.zeros(4)
    for key,value in old['complete_base_bivector'].items():
        i,j=map(int,key.split(','));P[i,j]=s.sympify(value,locals=env);P[j,i]=-P[i,j]
    coordinate2=[]
    for i,j in combinations_with_replacement(range(4),2):coordinate2.append(x[i]*x[j])
    def vector(M):return s.Matrix([s.Poly(s.expand(M[i,j]),*x).coeff_monomial(mon)for i,j in combinations(range(4),2)for mon in coordinate2])
    coefficient=s.Matrix.hstack(*[vector(P.diff(t))for t in b]);dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    rows=list(dm(coefficient.T).rref()[1]);inverse=coefficient[rows,:].inv();actions={}
    fields={'h':s.Matrix([3*x[0],x[1],-x[2],-3*x[3]]),'e':s.Matrix([0,x[0],2*x[1],3*x[2]]),'f':s.Matrix([3*x[1],2*x[2],x[3],0])}
    for name,X in fields.items():
        images=[]
        for t in b:
            K=P.diff(t);L=s.zeros(4)
            for i,j in combinations(range(4),2):
                L[i,j]=s.expand(sum(X[k]*s.diff(K[i,j],x[k])-K[k,j]*s.diff(X[i],x[k])-K[i,k]*s.diff(X[j],x[k])for k in range(4)));L[j,i]=-L[i,j]
            images.append(vector(L))
        images=s.Matrix.hstack(*images);rho=inverse*images[rows,:];assert coefficient*rho==images;actions[name]=rho
    # A,C,E coordinates on L3; A,B,D on the complementary3space carrying Q2.
    L=s.Matrix([[1,0,0],[0,0,0],[0,1,0],[0,-1,0],[0,0,1],[2,0,0]])
    W=s.Matrix([[1,0,0],[0,1,0],[0,0,0],[0,0,1],[0,6,0],[3,0,0]])
    change=L.row_join(W);assert change.det()!=0;block={}
    for name,rho in actions.items():
        transformed=change.inv()*rho*change;assert transformed[:3,3:]==s.zeros(3)and transformed[3:,:3]==s.zeros(3);block[name]=transformed
    aa,cc,ee=s.symbols('A C E');coords=s.Matrix([aa,cc,ee]);sym=s.symbols('t0:6');mons=[aa*aa,aa*cc,aa*ee,cc*cc,cc*ee,ee*ee];Qtest=sum(u*v for u,v in zip(sym,mons));constraints=[]
    for rho in block.values():
        polynomial=s.expand(sum(s.diff(Qtest,coords[i])*(rho[:3,:3]*coords)[i]for i in range(3)));constraints.extend(s.Poly(polynomial,*coords).coeffs())
    Cmat=s.linear_eq_to_matrix(constraints,sym)[0];null=Cmat.nullspace();assert len(null)==1;invariant=s.factor(sum(a*m for a,m in zip(null[0],mons)))
    for rho in block.values():
        qW=aa*ee-4*cc*cc;assert s.expand(sum(s.diff(qW,coords[i])*(rho[3:,3:]*coords)[i]for i in range(3)))==0
    matrix=s.Matrix([[s.sympify(v,locals={str(t):t for t in b})for v in row]for row in old['all350_linear_vertex_fiber_matrix']])
    assert invariant.subs({aa:0,cc:0,ee:1})!=0 and invariant.subs({aa:0,cc:1,ee:0})==0
    representatives=[('L_semisimple',3,L[:,2]),('L_nilpotent',2,L[:,1]),('Q_nonzero',2,W[:,0])];records=[]
    for name,base_dimension,point in representatives:
        rankmatrix=matrix.subs(dict(zip(b,point)));rank=len(dm(rankmatrix).rref()[1]);kernel=dm(rankmatrix).nullspace().to_Matrix().T;assert rankmatrix*kernel==s.zeros(350,kernel.cols)
        rr=list(dm(rankmatrix.T).rref()[1]);pc=list(dm(rankmatrix).rref()[1]);det=rankmatrix[rr,pc].det();assert det!=0
        total=base_dimension+kernel.cols
        records.append({'name':name,'base_stratum_dimension':base_dimension,'base_representative':list(map(str,point)),
            'full350_by12_vertex_matrix':[[str(v)for v in row]for row in rankmatrix.tolist()],'vertex_matrix_rank':rank,'full_vertex_kernel_dimension':kernel.cols,
            'full_vertex_kernel_basis':[[str(v)for v in row]for row in kernel.tolist()],'rank_minor':{'rows':rr,'columns':pc,'determinant':str(det)},
            'whole_fixed_affine_source_dimension_upper':total,'whole_projective_PGL_image_dimension_upper':24+total-1-8})
        print(name,'base dim',base_dimension,'Vkernel',kernel.cols,'fixed total',total,'PGL',15+total,flush=True)
    result={'task_id':'G','input_report':str(path.relative_to(ROOT)),'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'SL2_linear_fields':{k:list(map(str,v))for k,v in fields.items()},'SL2_action_on_full_base6':{k:[[str(v)for v in row]for row in M.tolist()]for k,M in actions.items()},
        'L3_basis':L.tolist(),'Q2_linear_span_basis':W.tolist(),'SL2_block_actions':{k:[[str(v)for v in row]for row in M.tolist()]for k,M in block.items()},
        'L3_unique_invariant_quadratic':str(invariant),'Q2_invariant_quadratic':'A*D-4*B**2','cases':records,
        'all_nonzero_base_strata_projective_PGL_dimension_at_most24':all(r['whole_projective_PGL_image_dimension_upper']<=24 for r in records),
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'seconds':time.monotonic()-start,
        'supersedes_prior_orbit_scope':'The first ORBITS report mislabeled A=1,C=E=0 as semisimple, although the exact invariant4AC+E^2 makes it nilpotent. Its arithmetic/producerbytes remain immutable and no universal coverage claim may use that report. V2 uses E=1 and explicitly verifies the nonzero invariant before computing its full vertex kernel.',
        'scope':'Allnonzero strictbasePoisson points split into the linearSym2 three-space (semisimple/nilpotent) and the nilcone in anotherSym2. The base0 vertex12 space has rank2 and is irrelevant to irreducible rank4 quintics. Geometric support decomposition/orbit exhaustiveness is proved separately, not inferred from three arbitrary sample points.'}
    destination.write_text(json.dumps(result,indent=2,default=str)+'\n')


if __name__=='__main__':
    if '--orbits'in sys.argv:orbits()
    else:main()
