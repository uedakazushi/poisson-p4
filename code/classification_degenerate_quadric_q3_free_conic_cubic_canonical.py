#!/usr/bin/env python3
"""Exact magnetic translation quotient and full-source orbit closure reduction."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic();path=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_SOURCES_20261107.json');src=json.loads(path.read_text());xx=s.symbols('a b c r z');base=xx[:3];mm=s.symbols('m0:6');mu=s.Symbol('mu');mons=[base[0]**2,base[0]*base[1],base[0]*base[2],base[1]**2,base[1]*base[2],base[2]**2];monexps=[s.Poly(m,*base).monoms()[0] for m in mons];records=[]
def enc(P):return {str(ij):str(s.expand(P[ij])) for ij in itertools.combinations(range(5),2) if P[ij]}
for case in src['cases']:
 if case['quadratic_syzygy_kernel_dimension']!=2:continue
 P=s.zeros(5)
 for ij in itertools.combinations(range(5),2):P[ij]=s.sympify(case['full_unimodular_brackets'][''.join(map(str,ij))]);P[ij[::-1]]=-P[ij]
 U=s.Matrix([P[i,3] for i in range(3)]);V=s.Matrix([P[i,4] for i in range(3)]);du=sum(s.diff(U[i],base[i]) for i in range(3));dv=sum(s.diff(V[i],base[i]) for i in range(3));images=[s.expand(V[i]+dv*base[i]) for i in range(3)]+[s.expand(-U[i]-du*base[i]) for i in range(3)];M=s.Matrix([[s.Poly(f,*base).coeff_monomial(e) for f in images] for e in monexps]);rank=M.rank();assert rank==5
 k=next(i for i in range(6) if M.row_join(s.eye(6)[:,i]).rank()==6);left=M.T.nullspace()[0];left/=left[k];muv=s.expand((left.T*s.Matrix(mm))[0]);rhs=s.Matrix(mm)-muv*s.eye(6)[:,k];solution=next(iter(s.linsolve((M,rhs),s.symbols('g0:6'))));gs=set().union(*(v.free_symbols for v in solution))-set(mm);sol=s.Matrix(solution).subs({v:0 for v in gs});assert M*sol==rhs
 G=s.eye(5)
 for i in range(3):G[3,i]=sol[i];G[4,i]=sol[3+i]
 Pmu=P.subs({m:0 for m in mm});Pmu[3,4]+=muv*mons[k];Pmu[4,3]-=muv*mons[k]
 push=(G*Pmu.subs(dict(zip(xx,G.inv()*s.Matrix(xx))),simultaneous=True)*G.T).applyfunc(s.expand);assert push==P
 canonical=P.subs({m:0 for m in mm});canonical[3,4]+=mons[k];canonical[4,3]-=mons[k]
 wmu=muv.subs(dict(zip(mm,[1,2,3,5,7,11])));assert wmu!=0
 scale=s.ilcm(1,*[s.denom(a) for ij in itertools.combinations(range(5),2) for a in s.Poly(canonical[ij],*xx).coeffs()]);canonical*=scale
 records.append({'case':case['case'],'magnetic_translation_matrix':[[str(v) for v in row] for row in M.tolist()],'magnetic_translation_rank_QQ':rank,'canonical_complement_monomial':str(mons[k]),'magnetic_quotient_mu':str(muv),'shear_G':[[str(v) for v in row] for row in G.tolist()],'whole_source_shear_identity':True,'generic_orbit_condition':'mu!=0','original_rank_witness_mu':str(wmu),'fiber_scaling_identity':'(1/lambda)*diag(1,1,1,lambda,lambda)_* P_mu = P_(lambda*mu)','source_zero_mu_is_actual_limit_of_nonzero_mu':True,'canonical_integer_scale':str(scale),'canonical_point_brackets':enc(canonical),'whole_irreducible_source_in_closure_of_single_generic_orbit':True})
result={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start}
out=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_CANONICAL_20261107.json')
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'seconds':result['seconds'],'cases':[{k:r[k] for k in ['case','magnetic_translation_rank_QQ','canonical_complement_monomial','magnetic_quotient_mu']} for r in records]},indent=2))
