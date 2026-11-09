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


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--case',choices=CASES,default='smooth_scroll_12');args=parser.parse_args();main(args.case)
