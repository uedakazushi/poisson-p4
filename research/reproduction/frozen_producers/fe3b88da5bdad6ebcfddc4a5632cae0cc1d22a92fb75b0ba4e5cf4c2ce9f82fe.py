#!/usr/bin/env python3
"""Whole tangent/abelian-conormal sources on both cubic surface orbits.

Retains all150 coefficients before Euler24. The conormal equation is the
intrinsic q_i/x0^2 bracket, with all Euler correction terms, not the rankzero
raw-bracket ansatz. The first finite run uses only linear algebra and one
full-Jacobi rowspace/leading-monomial bound, never a GB search.
"""
import argparse
from functools import lru_cache
from itertools import combinations,combinations_with_replacement
from fractions import Fraction as F
import hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
import classification_projected_veronese_zero_review as a
import classification_cubic_scroll_zero as old
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,minor_certificate

ROOT=Path(__file__).resolve().parents[1]
CASES={
 'smooth_scroll_12':dict(q=['x0*x3-x1*x2','x0*x4-x1*x3','x2*x4-x3**2'],
    variables=['ss','tt','u','v'],mapping=['ss*u','tt*u','ss**2*v','ss*tt*v','tt**2*v'],
    powers=[(1,0,1,0),(0,1,1,0),(2,0,0,1),(1,1,0,1),(0,2,0,1)],group_dimension=6,target_fixed_dimension=7),
 'twisted_cubic_cone_03':dict(q=['x0*x2-x1**2','x0*x3-x1*x2','x1*x3-x2**2'],
    variables=['ss','tt','c'],mapping=['ss**3','ss**2*tt','ss*tt**2','tt**3','c'],
    powers=[(3,0,0),(2,1,0),(1,2,0),(0,3,0),(0,0,1)],group_dimension=8,target_fixed_dimension=9)}


def rank(M):return len(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1])
def lift(poly,powers):
    ans={}
    for e,c in poly.items():
        f=tuple(sum(k*p[j]for k,p in zip(e,powers))for j in range(len(powers[0])))
        ans[f]=ans.get(f,F())+c
    return {e:c for e,c in ans.items()if c}
def sparse_rows(rows):return [{str(j):str(c)for j,c in row.items()}for row in rows]
def encode_matrix(M):return [[str(v)for v in row]for row in M.tolist()]


def linear_source(case):
    X=a.Ring([f'x{i}'for i in range(5)]);x=[X.var(f'x{i}')for i in range(5)]
    q=[X.parse(t)for t in case['q']];dq=[[a.derivative(v,k)for k in range(5)]for v in q]
    assert all(not lift(v,case['powers'])for v in q)
    tangent={};conormal={}
    def insert(target,label,col,value):
        for exponent,c in lift(value,case['powers']).items():
            row=target.setdefault(label+(exponent,),{});row[col]=row.get(col,F())+c
    for col,((u,v),e)in enumerate(COLUMNS):
        m={tuple(e):F(1)}
        for g in range(3):
            contraction=[{}for _ in range(5)]
            contraction[u]=a.times(m,dq[g][v]);contraction[v]=a.scale(a.times(m,dq[g][u]),-1)
            for i,j in PAIRS:
                value=a.plus(a.times(x[i],contraction[j]),a.scale(a.times(x[j],contraction[i]),-1))
                insert(tangent,(g,i,j),col,value)
        for g,h in combinations(range(3),2):
            raw=a.times(m,a.plus(a.times(dq[g][u],dq[h][v]),a.scale(a.times(dq[g][v],dq[h][u]),-1)))
            # {q_g,x0}=m*(dq_g[u]*delta(v,0)-dq_g[v]*delta(u,0)).
            g0=a.times(m,a.plus(a.scale(dq[g][u],int(v==0)),a.scale(dq[g][v],-int(u==0))))
            h0=a.times(m,a.plus(a.scale(dq[h][u],int(v==0)),a.scale(dq[h][v],-int(u==0))))
            intrinsic=a.plus(a.times(x[0],raw),a.scale(a.times(q[h],g0),-2),a.scale(a.times(q[g],h0),2))
            for k in range(5):insert(conormal,(g,h,k),col,a.derivative(intrinsic,k))
    def clean(target):
        labels=sorted(target);rows=[{j:c for j,c in target[l].items()if c}for l in labels]
        return labels,rows
    tl,tr=clean(tangent);cl,cr=clean(conormal)
    rows=tr+cr;M=s.MutableSparseMatrix(len(rows),150,{(i,j):s.Rational(c)for i,row in enumerate(rows)for j,c in row.items()})
    dm=DomainMatrix.from_Matrix(M).convert_to(s.QQ);K=dm.nullspace().to_Matrix().T
    assert M*K==s.zeros(M.rows,K.cols)
    # Full150/Euler24: never replace tangency by coefficientwise vanishing.
    xx=s.Matrix(s.symbols('x0:5'))
    def decode(v):
        P=s.zeros(5)
        for coefficient,(ij,e)in zip(v,COLUMNS):
            if coefficient:P[ij]+=coefficient*s.prod(y**k for y,k in zip(xx,e))
        for i,j in PAIRS:P[j,i]=-P[i,j]
        return P.applyfunc(s.expand)
    def vector(P):return s.Matrix([s.Poly(s.expand(P[ij]),*xx).coeff_monomial(e)for ij,e in COLUMNS])
    def canonical(P):
        divergence=s.Matrix([sum(s.diff(P[i,j],xx[i])for i in range(5))for j in range(5)])
        return (P-(xx*divergence.T-divergence*xx.T)/5).applyfunc(s.expand)
    vertical=[]
    for i in range(5):
        for j in range(5):
            L=s.zeros(5,1);L[i]=xx[j];vertical.append(vector(xx*L.T-L*xx.T))
    V=s.Matrix.hstack(*vertical);assert rank(V)==24 and M*V==s.zeros(M.rows,25)
    CK=s.Matrix.hstack(*(vector(canonical(decode(K[:,j])))for j in range(K.cols)))
    reduced,piv=DomainMatrix.from_Matrix(CK.T).convert_to(s.QQ).rref();B=reduced.to_Matrix()[:len(piv),:].T
    assert rank(V.row_join(B))==K.cols and B.cols==K.cols-24 and M*B==s.zeros(M.rows,B.cols)
    cert=minor_certificate(M,150-K.cols)
    basis_tensors=[decode(B[:,j])for j in range(B.cols)]
    return B,dict(full_tangent_row_labels=tl,full_tangent_rows=sparse_rows(tr),
       full_intrinsic_conormal_firstpartial_labels=cl,full_intrinsic_conormal_rows=sparse_rows(cr),
       combined_linear_shape=list(M.shape),combined_linear_rank=150-K.cols,
       full150_kernel_dimension=K.cols,rank_certificate=cert,full150_kernel_basis=encode_matrix(K),
       tangent_only_rank=rank(M[:len(tr),:]),full_Euler_dimension=24,all_Euler_directions_in_full_kernel=True,
       complete_canonical_projective_source_dimension=B.cols,full150_by_source_basis=encode_matrix(B),
       complete_source_basis=[{str(ij):str(P[ij])for ij in PAIRS if P[ij]!=0}for P in basis_tensors])


def nonlinear(B):
    d=B.cols;cols=list(combinations_with_replacement(range(d),2));ci={e:i for i,e in enumerate(cols)}
    P={}
    for ij in PAIRS:
        P[ij]={(k,e):F(B[j,k])for j,(pair,e)in enumerate(COLUMNS)if pair==ij for k in range(d)if B[j,k]}
        P[ij[1],ij[0]]={t:-v for t,v in P[ij].items()}
    strict={tri:{}for tri in combinations(range(5),3)}
    for tri in strict:
        for i,j,k in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]:
            for ell in range(5):
                for key,value in old.product(P.get((i,ell),{}),P[j,k],ell).items():old.add(strict[tri],key,value)
    rows=[]
    for quad in combinations(range(5),4):
        polynomial={}
        for k,missing in enumerate(quad):
            for (pair,e),value in strict[quad[:k]+quad[k+1:]].items():
                f=list(e);f[missing]+=1;old.add(polynomial,(pair,tuple(f)),(-1)**k*value)
        rows.extend({ci[pair]:value for (pair,f),value in polynomial.items()if f==e}for e in a.monomials(5,4))
    assert len(rows)==350
    D={}
    for missing in range(5):
        i,j,k,ell=[v for v in range(5)if v!=missing]
        for ij,kl,sign in [((i,j),(k,ell),1),((i,k),(j,ell),-1),((i,ell),(j,k),1)]:
            for (pair,e),value in old.product(P[ij],P[kl]).items():
                f=list(e);f[missing]+=1;old.add(D,(pair,tuple(f)),(-1)**missing*sign*value)
    Drows=[{ci[pair]:value for (pair,f),value in D.items()if f==e}for e in a.monomials(5,5)]
    basis,proof,selected=a.sparse_rref(rows,True);Dcert=[];contained=True
    for i,target in enumerate(Drows):
        residual=dict(target);functional={}
        for p in sorted(basis):
            if p in residual:
                c=residual[p];residual=a.plus(residual,a.scale(basis[p],-c));functional=a.plus(functional,a.scale(proof[p],c))
        assert a.plus(a.plus(*(a.scale(rows[i],c)for i,c in functional.items())),residual)==target
        contained&=not residual
        if target:Dcert.append(dict(D_row=i,functional={str(k):str(v)for k,v in functional.items()},residual={str(k):str(v)for k,v in residual.items()}))
    leading=[cols[p]for p in sorted(basis)];edges=tuple(sorted(set((1<<i)|(1<<j)for i,j in leading)))
    @lru_cache(None)
    def cover(edges):
        if not edges:return 0
        edge=min(edges,key=int.bit_count)
        return 1+min(cover(tuple(e for e in edges if not e&(1<<i)))for i in range(d)if edge&(1<<i))
    height=cover(edges)
    return dict(all_source_parameter_quadratic_monomials=cols,full350_Jacobi_rows=sparse_rows(rows),
         full126_actual_D_rows=sparse_rows(Drows),Jacobi_quadric_row_rank=len(basis),
         selected_original350_Jacobi_rows=selected,full_exact_Jacobi_row_basis={str(p):{str(k):str(v)for k,v in row.items()}for p,row in basis.items()},
         full_exact_basis_from_original350_functionals={str(p):{str(k):str(v)for k,v in row.items()}for p,row in proof.items()},
         all_actual_D_rows_in_full_Jacobi_quadric_rowspace=contained,actual_D_rowspace_certificates=Dcert,
         leading_monomial_pairs=leading,leading_support_minimum_vertex_cover=height,
         single_rowspace_initial_ideal_fixed_source_dimension_upper=d-height)


def main(case_name):
    jobpath='research/classification/jobs/cubic_surface_tangent_preflight.json'
    output=ROOT/('reports/CLASSIFICATION_CUBIC_SURFACE_TANGENT_'+case_name.upper()+'_20261107.json')
    if output.exists():raise FileExistsError(output)
    job=json.loads((ROOT/jobpath).read_text())
    for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256']
    resource.setrlimit(resource.RLIMIT_CPU,(115,116));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2));signal.alarm(120)
    start=time.monotonic();case=CASES[case_name];B,linear=linear_source(case)
    result=dict(task_id='G',case=case_name,surface=case,all_linear_arithmetic_passed=True,linear=linear,
      scope='All150 ambient coefficients tangent to the complete cubic surface with intrinsic abelian conormal. Rank2 generic restriction plus singular actual divisor forces this necessary condition; no rankzero ansatz or node/vertical assumption.',
      producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_bytes=Path(__file__).read_text(),exact_inputs_and_hashes=job['exact_inputs_and_hashes'])
    output.write_text(json.dumps(result,indent=2)+'\n')
    print('COMPLETE LINEAR SOURCE',case_name,B.cols,flush=True)
    if B.cols<=26:
        result['nonlinear']=nonlinear(B);bound=result['nonlinear']['single_rowspace_initial_ideal_fixed_source_dimension_upper']
        result['whole_projective_PGL_image_upper']=24+bound-1-case['group_dimension']
        result['rank4_source_empty_by_quadric_rowspace']=result['nonlinear']['all_actual_D_rows_in_full_Jacobi_quadric_rowspace']
        result['weak_G_bound_obtained']=result['whole_projective_PGL_image_upper']<=24 or result['rank4_source_empty_by_quadric_rowspace']
    else:result['weak_G_bound_obtained']=False;result['nonclosing_reason']='Source dimension above26: stop before nonlinear arithmetic; change geometry.'
    result.update(seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result.get(k)for k in ('weak_G_bound_obtained','whole_projective_PGL_image_upper','rank4_source_empty_by_quadric_rowspace','nonclosing_reason','seconds','maximum_RSS_KiB')}),flush=True)


def ranktwo_open():
    """One actual Aut(S)-covered rank2 patch, using only constant-unit solves."""
    output=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_RANK_TWO_OPEN_20261107.json'
    if output.exists():raise FileExistsError(output)
    job=json.loads((ROOT/'research/classification/jobs/cubic_scroll_rank_two_open.json').read_text())
    for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256']
    resource.setrlimit(resource.RLIMIT_CPU,(115,116));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2));signal.alarm(120)
    start=time.monotonic();sourcepath=ROOT/'reports/CLASSIFICATION_CUBIC_SURFACE_TANGENT_SMOOTH_SCROLL_12_20261107.json'
    source=json.loads(sourcepath.read_text());n=source['nonlinear'];t=s.symbols('t0:21');pairs=n['all_source_parameter_quadratic_monomials']
    J=[s.expand(sum(s.Rational(v)*t[pairs[int(j)][0]]*t[pairs[int(j)][1]]for j,v in row.items()))for row in n['full350_Jacobi_rows']]
    sub={t[2]:s.Integer(1)};steps=[]
    for indices,unknowns in [([80],[18]),([26,97],[16,20]),([81],[10]),([57],[4]),([27,99,117],[5,12,13]),([59],[8])]:
        variables=[t[i]for i in unknowns];equations=[s.expand(J[i].subs(sub,simultaneous=True))for i in indices]
        M,rhs=s.linear_eq_to_matrix(equations,variables)
        assert all(not c.free_symbols for c in M)
        det=s.expand(M.det());assert det!=0
        I=M.inv();solution=I*rhs
        update={u:s.expand(v)for u,v in zip(variables,solution)}
        assert all(not any(u in v.free_symbols for u in variables)for v in update.values())
        assert all(s.expand(f.subs(update,simultaneous=True))==0 for f in equations)
        sub.update(update)
        steps.append(dict(original350_rows=indices,solved_variables=unknowns,constant_coefficient_matrix=encode_matrix(M),
                          determinant=str(det),exact_inverse=encode_matrix(I),exact_solution={str(k):str(v)for k,v in update.items()},
                          all_selected_original_equations_identically_zero=True))
    free=[v for v in t if v not in sub];assert len(free)==11
    remaining=[s.Poly(s.expand(f.subs(sub,simultaneous=True)),*free,domain=s.QQ)for f in J]
    allmon=sorted(set().union(*(set(p.monoms())for p in remaining if p)),key=lambda e:(sum(e),e),reverse=True)
    mindex={e:i for i,e in enumerate(allmon)}
    rows=[{mindex[e]:F(c)for e,c in p.terms()if c}for p in remaining]
    basis,proof,selected=a.sparse_rref(rows,True)
    leading=[allmon[p]for p in sorted(basis)];edges=tuple(sorted(set(sum(1<<i for i,p in enumerate(e)if p)for e in leading)))
    @lru_cache(None)
    def cover(edges):
        if not edges:return 0
        if 0 in edges:return len(free)+1
        edge=min(edges,key=int.bit_count)
        return 1+min(cover(tuple(e for e in edges if not e&(1<<i)))for i in range(len(free))if edge&(1<<i))
    height=cover(edges);bound=len(free)-height
    # p=(0,0,1,0,0), chart x2=1, tangent coordinates (k=x0,r=x3).
    B=s.Matrix([[s.Rational(c)for c in row]for row in source['linear']['full150_by_source_basis']])
    spatial=(0,0,2,0,0);row=next(i for i,(ij,e)in enumerate(COLUMNS)if ij==(0,3)and tuple(e)==spatial)
    assert list(B[row,:])==[s.Integer(i==2)for i in range(21)]
    result=dict(task_id='G',scope='All rank2 restrictions on smooth cubic scroll meet this patch via connected actual Aut(S); t2=1 is already tensor projectivization, so normalized target dimension<=6, not7.',
      dense_surface_orbit_point=[0,0,1,0,0],tangent_patch_coordinates=['k=x0/x2','r=x3/x2'],
      restriction_Poisson_at_point='t2*partial_k wedge partial_r',actual_point_evaluation_verified=True,
      tensor_scale_normalization='t2=1',constant_unit_elimination_steps=steps,
      exact_full_source_parameter_substitution={str(k):str(v)for k,v in sub.items()},remaining_parameters=list(map(str,free)),
      full350_remaining_polynomials=[str(p.as_expr())for p in remaining],remaining_monomial_order='graded lexicographic, descending total degree then lex',
      remaining_all_monomial_exponents=allmon,remaining_Jacobi_row_rank=len(basis),
      remaining_exact_row_basis={str(p):{str(k):str(v)for k,v in row.items()}for p,row in basis.items()},
      remaining_exact_basis_original350_functionals={str(p):{str(k):str(v)for k,v in row.items()}for p,row in proof.items()},
      leading_monomial_exponents=leading,leading_minimum_vertex_cover=height,
      normalized_whole_rank2_slice_dimension_upper=bound,whole_projective_PGL_image_upper=18+bound,
      G_non_dominance_obtained=bound<=6,all350_equations_retained=True,no_GB_or_sampling=True,
      producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_bytes=Path(__file__).read_text(),
      exact_inputs_and_hashes=job['exact_inputs_and_hashes'],seconds=time.monotonic()-start,
      maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ('G_non_dominance_obtained','normalized_whole_rank2_slice_dimension_upper','whole_projective_PGL_image_upper','remaining_Jacobi_row_rank','seconds','maximum_RSS_KiB')}),flush=True)



def ranktwo_branch_support():
    """Whole three exact branches; linear support over all base parameters."""
    start=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    jobpath=ROOT/'research/classification/jobs/cubic_scroll_rank_two_branch_support.json'
    output=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_RANK_TWO_BRANCH_SUPPORT_20261107.json'
    if output.exists():raise FileExistsError(output)
    job=json.loads(jobpath.read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    raw=json.loads((ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_RANK_TWO_OPEN_20261107.json').read_text())
    t=s.symbols('t0:21');env={str(v):v for v in t};full=[s.sympify(f,locals=env)for f in raw['full350_remaining_polynomials']]
    product=(t[0]-t[17])*(2*t[0]-t[17])*(2*t[0]-2*t[17]+3*t[19]);assert s.expand(full[73]-s.Rational(5,9)*product)==0
    secondary=[t[i]for i in [1,3,6,7,9,11,14,15]];base=[t[0],t[19]];records=[]
    for name,sub in [('first',{t[17]:t[0]}),('second',{t[17]:2*t[0]}),('third',{t[17]:t[0]+3*t[19]/2})]:
        J=[s.Poly(s.expand(f.subs(sub)),*secondary,domain=s.QQ.poly_ring(*base))for f in full]
        linear=[(i,f)for i,f in enumerate(J)if f.total_degree()<=1 and f.as_expr()!=0]
        A=s.Matrix([[f.coeff_monomial(v)for v in secondary]for _,f in linear]);b=s.Matrix([-f.coeff_monomial(1)for _,f in linear]);dm=DomainMatrix.from_Matrix(A).to_field();_,cc=dm.rref();cc=list(cc)
        _,rr=DomainMatrix.from_Matrix(A[:,cc].T).to_field().rref();rr=list(rr);B=A[rr,cc];det=s.factor(B.det());assert det
        rank=len(cc);assert A.rank()==rank
        record=dict(name=name,exact_branch_substitution={str(k):str(v)for k,v in sub.items()},all350_polynomials=[str(f.as_expr())for f in J],all8_secondary_variables=list(map(str,secondary)),all2_base_variables=list(map(str,base)),linear_original350_rows=[i for i,_ in linear],full_linear_matrix=encode_matrix(A),full_linear_rhs=encode_matrix(b),generic_rank=rank,independent_row_indices=rr,independent_columns=cc,selected_nonzero_minor=str(det),generic_normalized_source_dimension_upper=10-rank,minor_factors=[dict(polynomial=str(f),multiplicity=int(n))for f,n in s.factor_list(det)[1]],generic_augmented_rank=DomainMatrix.from_Matrix(A.row_join(b)).to_field().rank())
        records.append(record)
    result=dict(task_id='G',all_exact_branch_and_linear_identities_passed=True,scope='Complete three linear branch equations from original350 and complete secondary-linear projection. Generic ranks alone do not prove fullsupport dimension; every minor zero component and residual nonlinear equation is retained.',all_G_closed=False,all_three_cases=records,inputs=job['exact_inputs_and_hashes'],producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(cases=[{k:r[k]for k in ['name','generic_rank','selected_nonzero_minor','generic_normalized_source_dimension_upper','generic_augmented_rank']}for r in records],seconds=result['seconds'])),flush=True)

def normal_root_linear_replay():
    """Fresh reconstruction of the unpublished pre-restart linear packet.

    Combines the three lost isolated modes in one reusable mode. This report
    has a fresh identity and fresh measurements; no lost report is fabricated.
    The nine published function ASTs remain unchanged.
    """
    start=time.monotonic();signal.alarm(90);resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,)*2)
    jobpath=ROOT/'research/classification/jobs/cubic_scroll_normal_root_linear_replay.json'
    output=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_NORMAL_ROOT_LINEAR_REPLAY_20261107.json'
    if output.exists():raise FileExistsError(output)
    job=json.loads(jobpath.read_text())
    for binding in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/binding['path']).read_bytes()).hexdigest()==binding['sha256'],binding['path']
    raw=json.loads((ROOT/'reports/CLASSIFICATION_CUBIC_SURFACE_TANGENT_SMOOTH_SCROLL_12_20261107.json').read_text())
    patch=json.loads((ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_RANK_TWO_OPEN_20261107.json').read_text())
    B=s.Matrix([[s.Rational(c)for c in row]for row in raw['linear']['full150_by_source_basis']])
    t=s.symbols('t0:21');r=s.Symbol('r');xx=s.symbols('x0:5');P=s.zeros(5)
    for c,(ij,e)in zip(B*s.Matrix(t),COLUMNS):P[ij]+=c*s.prod(y**k for y,k in zip(xx,e))
    restriction=s.Poly(s.expand((P[1,2]-r*P[0,2]).subs(dict(zip(xx,[1,r,0,0,0])),simultaneous=True)),r)
    assert restriction.degree()<=1
    equations=[s.expand(restriction.coeff_monomial(r**k))for k in range(2)]
    assert equations==[t[18],-t[10]]
    M0,z=s.linear_eq_to_matrix(equations,t);assert z==s.zeros(2,1)
    K0=DomainMatrix.from_Matrix(M0).convert_to(s.QQ).nullspace().to_Matrix().T
    assert M0*K0==s.zeros(2,K0.cols)and K0.cols==19
    env={str(v):v for v in t}
    sub={env[k]:s.sympify(v,locals=env)for k,v in patch['exact_full_source_parameter_substitution'].items()}
    peq=[s.factor(f.subs(sub,simultaneous=True))for f in equations]
    assert s.expand(peq[0]-(2*t[0]-t[17])*(2*t[0]-2*t[17]+3*t[19])/3)==0
    directrix=dict(full_C0_restriction=str(restriction.as_expr()),two_linear_C0_equations=list(map(str,equations)),full_linear_matrix=encode_matrix(M0),rank=2,full_kernel_dimension=19,full_kernel_basis=encode_matrix(K0),complete_ranktwo_patch_C0_equations=list(map(str,peq)),all_patch_substitution={str(k):str(v)for k,v in sub.items()})
    X=a.Ring([f'x{i}'for i in range(5)]);x=[X.var(f'x{i}')for i in range(5)]
    q=X.parse('x2*x4-x3**2');dq=[a.derivative(q,k)for k in range(5)];rootrows={};rawrows={}
    for col,((u,v),e)in enumerate(COLUMNS):
        m={tuple(e):F(1)}
        bq=[a.times(m,a.plus(a.scale(dq[u],int(v==k)),a.scale(dq[v],-int(u==k))))for k in range(5)]
        for k in range(5):
            for exponent,c in lift(bq[k],CASES['smooth_scroll_12']['powers']).items():
                row=rawrows.setdefault((k,exponent),{});row[col]=row.get(col,F())+c
        for k in range(1,5):
            b0k=a.scale(m,int(u==0 and v==k)-int(v==0 and u==k))
            H=a.plus(a.times(x[0],bq[k]),a.scale(a.times(x[k],bq[0]),-1),a.scale(a.times(q,b0k),-2))
            dH=[a.derivative(H,j)for j in range(5)]
            for i,j in PAIRS:
                f=a.plus(a.times(dH[i],dq[j]),a.scale(a.times(dH[j],dq[i]),-1))
                for exponent,c in lift(f,CASES['smooth_scroll_12']['powers']).items():
                    row=rootrows.setdefault((k,i,j,exponent),{});row[col]=row.get(col,F())+c
    matrices=[];records=[]
    for rows in [rootrows,rawrows]:
        labels=sorted(rows);rr=[{j:c for j,c in rows[k].items()if c}for k in labels]
        M=s.MutableSparseMatrix(len(rr),150,{(i,j):s.Rational(c)for i,row in enumerate(rr)for j,c in row.items()})
        A=M*B;matrices.append((M,A))
        records.append(dict(all_row_labels=labels,full150_rows=sparse_rows(rr),full_matrix_shape=list(M.shape),complete21_restricted_matrix=encode_matrix(A)))
    M,A=matrices[0];C=matrices[1][1]
    K=DomainMatrix.from_Matrix(A).convert_to(s.QQ).nullspace().to_Matrix().T
    assert A*K==s.zeros(A.rows,K.cols)and M.shape==(620,150)and K.cols==16
    records[0].update(restricted_rank=5,complete_kernel_dimension=16,complete21_kernel_basis=encode_matrix(K),full150_kernel_basis=encode_matrix(B*K),rank_certificate=minor_certificate(A,5))
    assert matrices[1][0].shape==(84,150)
    total=A.col_join(M0).col_join(C);KT=DomainMatrix.from_Matrix(total).convert_to(s.QQ).nullspace().to_Matrix().T
    assert total*KT==s.zeros(total.rows,KT.cols)and KT.cols==16
    V=[]
    for i in range(5):
        for j in range(5):
            L=s.zeros(5,1);L[i]=xx[j];VP=s.Matrix(xx)*L.T-L*s.Matrix(xx).T
            V.append(s.Matrix([s.Poly(VP[ij],*xx).coeff_monomial(e)for ij,e in COLUMNS]))
    V=s.Matrix.hstack(*V);assert M*V==s.zeros(M.rows,25)
    for source in range(21):
        PS=s.zeros(5)
        for c,(ij,e)in zip(B[:,source],COLUMNS):PS[ij]+=c*s.prod(y**k for y,k in zip(xx,e))
        for i,j in PAIRS:PS[j,i]=-PS[i,j]
        assert all(s.expand(sum(s.diff(PS[i,j],xx[i])for i in range(5)))==0 for j in range(5))
    combined=dict(full21_combined_matrix=encode_matrix(total),shape=list(total.shape),rank=5,complete_kernel_dimension=16,complete21_kernel_basis=encode_matrix(KT),full150_kernel_basis=encode_matrix(B*KT),rank_certificate=minor_certificate(total,5),all21_source_divergences_identically_zero=True,whole_projective_image_upper=33,sufficient_G_bound_from_linear_necessity=False)
    result=dict(task_id='G',replay_provenance='Fresh reconstruction after VM replacement lost three unpublished reports; no original elapsed time or report checksum is claimed. Structural identities match the recorded tool outputs. Published nine functions were preserved.',all_arithmetic_passed=True,directrix=directrix,canonical_rootline=records[0],canonical_raw_contraction=records[1],combined_casimir_necessity=combined,all_Euler24_annihilated_by_intrinsic_rootline=True,all_G_closed=False,scope='Linear diagnostics only. Canonical-root emptiness is instead proved by the separate pure Hartogs/Casimir/divergence note; A25 two types and cubic cone remain open.',inputs=job['exact_inputs_and_hashes'],producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(all_arithmetic_passed=True,directrix=directrix['two_linear_C0_equations'],rootline_shape=records[0]['full_matrix_shape'],rootline_rank=5,combined_dimension=16,linear_G_bound=False,seconds=result['seconds'],maximum_RSS_KiB=result['maximum_RSS_KiB'])),flush=True)



def cone_base_vertex():
    """Whole cone source as base9 plus divergence-free vertex18, no GB."""
    import classification_degenerate_quadric_q4_nil31_I_review as q
    start=time.monotonic();signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    output=ROOT/'reports/CLASSIFICATION_CUBIC_CONE_RANK_TWO_BASE_VERTEX_20261107.json'
    if output.exists():raise FileExistsError(output)
    job=json.loads((ROOT/'research/classification/jobs/cubic_cone_rank_two_base_vertex.json').read_text())
    for b in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    source=json.loads((ROOT/'reports/CLASSIFICATION_CUBIC_SURFACE_TANGENT_TWISTED_CUBIC_CONE_03_20261107.json').read_text())
    B=[[F(v)for v in row]for row in source['linear']['full150_by_source_basis']]
    tensors=[]
    for j in range(27):
        tensor={}
        for row,(ij,e)in enumerate(COLUMNS):
            if B[row][j]:tensor.setdefault(ij,{})[tuple(e)+(0,0,0)]=B[row][j]
        tensors.append(q.antisymmetric(tensor))
    base=[j for j,P in enumerate(tensors)if any(P.get(ij,{})for ij in PAIRS if ij[1]<4)]
    vertex=[j for j in range(27)if j not in base]
    assert base==[0,1,2,3,4,5,6,7,16]and len(vertex)==18
    for j,P in enumerate(tensors):
        for k in range(5):assert not q.add(*(q.der(P.get((i,k),{}),i)for i in range(5)))
        if j in base:
            for i in range(4):
                L=q.add(*(q.der(P.get((k,i),{}),k)for k in range(4)))
                assert P.get((i,4),{})==q.mul(q.atom(1,4),L)
            assert all(not ex[4]for ij in PAIRS if ij[1]<4 for ex in P.get(ij,{}))
        else:
            assert all(not ex[4]for ij in PAIRS for ex in P.get(ij,{}))
            assert not any(q.jacobi(P).values())
    def plus(A,B):return q.antisymmetric({ij:q.add(A.get(ij,{}),B.get(ij,{}))for ij in PAIRS})
    def J(P):return q.jvec(q.jacobi(P))
    def D(P):
        out={}
        for m in range(5):
            i,j,k,l=[v for v in range(5)if v!=m]
            pf=q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{})))
            out=q.add(out,q.scale(q.mul(q.atom(1,m),pf),(-1)**m))
        return out
    base_tensors=[tensors[i]for i in base];v_tensors=[tensors[i]for i in vertex]
    baseJ=[J(P)for P in base_tensors];baseD=[D(P)for P in base_tensors]
    mons=list(combinations_with_replacement(range(9),2));purecols=[];pureD=[]
    for i,j in mons:
        if i==j:purecols.append(baseJ[i]);pureD.append(baseD[i])
        else:
            purecols.append([v-a-b for v,a,b in zip(J(plus(base_tensors[i],base_tensors[j])),baseJ[i],baseJ[j])])
            pureD.append(q.add(D(plus(base_tensors[i],base_tensors[j])),q.scale(baseD[i],-1),q.scale(baseD[j],-1)))
    mixedcols=[];mixedD=[]
    for i,P in enumerate(base_tensors):
        for V in v_tensors:
            mixedcols.append([v-a for v,a in zip(J(plus(P,V)),baseJ[i])])
            mixedD.append(q.add(D(plus(P,V)),q.scale(baseD[i],-1)))
    # All V squares/crosses vanish: independent-c vector fields wedged with d/dc.
    for V,W in combinations(v_tensors,2):assert not any(J(plus(V,W)))
    pure_rows=[{j:c for j,c in enumerate(row)if c}for row in zip(*purecols)]
    mixed_rows=[{j:c for j,c in enumerate(row)if c}for row in zip(*mixedcols)]
    assert len(pure_rows)==len(mixed_rows)==350
    assert all(not(a and b)for a,b in zip(pure_rows,mixed_rows))
    reduced,proof,chosen=a.sparse_rref(pure_rows,True)
    tt=s.symbols('t0:9')
    equations=[s.factor(sum(s.Rational(c)*tt[mons[j][0]]*tt[mons[j][1]]for j,c in row.items()))for row in reduced.values()]
    # D is quadratic in base and affine linear in vertex. Keep all coefficients.
    dmons=a.monomials(5,5)
    def drows(columns):return [{str(j):str(P.get(tuple(e)+(0,0,0),F(0)))for j,P in enumerate(columns)if P.get(tuple(e)+(0,0,0),F(0))}for e in dmons]
    pure_D=drows(pureD);mixed_D=drows(mixedD)
    cdep=[];all_cdep_in_J=True
    for row,e in zip(pure_D,dmons):
        if not e[4]:continue
        residual={int(j):F(v)for j,v in row.items()}
        for pivot in sorted(reduced):
            if pivot in residual:residual=a.plus(residual,a.scale(reduced[pivot],-residual[pivot]))
        all_cdep_in_J &= not residual
        if row:cdep.append(dict(exponent=e,residual={str(j):str(v)for j,v in residual.items()}))
    result=dict(task_id='G',all_arithmetic_passed=True,complete_source_dimension=27,base_source_indices=base,vertex_source_indices=vertex,whole_Euler_canonical_source_identity='Pi=P(x)+(c div(P)+V(x)) wedge d/dc; all V divergencefree, c-independent, whole18',all350_parameter_pure_and_mixed_rows_disjoint=True,all_vertex_quadratics_Jacobi_zero=True,base_quadratic_monomials=mons,base_all350_rows=sparse_rows(pure_rows),base_Jacobi_quadric_row_rank=len(reduced),base_exact_row_basis={str(k):{str(j):str(v)for j,v in row.items()}for k,row in reduced.items()},base_exact_original350_functionals={str(k):{str(j):str(v)for j,v in row.items()}for k,row in proof.items()},base_selected_original350_rows=chosen,base_factored_row_basis=list(map(str,equations)),mixed_parameter_columns=[(i,j)for i in range(9)for j in range(18)],mixed_all350_rows=sparse_rows(mixed_rows),actualD_spatial_monomials=dmons,actualD_pure_base_rows=pure_D,actualD_mixed_rows=mixed_D,all_actualD_cdependent_rows_in_base_Jacobi_quadric_rowspace=bool(all_cdep_in_J),cdependent_actualD_row_remainders=cdep,scope='Whole source structure and one full finite equation presentation only; no nonlinear GB, ranktwo source dimension, orbit exhaustion or global G claimed.',inputs=job['exact_inputs_and_hashes'],producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(all_arithmetic_passed=True,base_quadric_rank=len(reduced),all_actualD_cdependent_rows_in_J=bool(all_cdep_in_J),seconds=result['seconds'],maximum_RSS_KiB=result['maximum_RSS_KiB'],base_equations=list(map(str,equations)))),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--case',choices=CASES,default='smooth_scroll_12');parser.add_argument('--ranktwo-open',action='store_true');parser.add_argument('--ranktwo-branch-support',action='store_true');parser.add_argument('--normal-root-linear-replay',action='store_true');parser.add_argument('--cone-base-vertex',action='store_true');args=parser.parse_args()
    if args.cone_base_vertex:cone_base_vertex()
    elif args.normal_root_linear_replay:normal_root_linear_replay()
    elif args.ranktwo_branch_support:ranktwo_branch_support()
    elif args.ranktwo_open:ranktwo_open()
    else:main(args.case)
