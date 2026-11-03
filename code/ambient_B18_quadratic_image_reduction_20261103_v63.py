#!/usr/bin/env python3
"""Exact rank2-normalizer reduction of ACTUAL CA/CB intersection at generic B18.

Whole flag fibers and a common regular target base are certified. Actual
intersection is the SUM OF IMAGE KERNELS; no normalizer tensor surrogate.
"""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector,action
from ambient_CB_global_flag_model_20261103_v63 import raw as rawB
ROOT=Path(__file__).resolve().parents[1];h,x,y,z,w=X;ey=s.Matrix([0,0,1,0,0]);PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))
def flags(M):
 L=M*ey;L[2]=0;ell=-M[1,:].T;ell[1]=0;W=[]
 for i,j in PAIRS:
  value=0
  for k in range(2):
   A=s.zeros(2,5);A[0,0]=A[1,1]=1;A[k,:]=-M[k,:];value+=A[:,[i,j]].det()
  W.append(value)
 W[PAIRS.index((0,1))]=0
 return s.Matrix(list(L)+list(ell)+W)
def mj(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cubics=[s.prod((h,x,z,w)[i]for i in m)for m in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in m)];swap=s.eye(5);swap[0,0]=swap[1,1]=0;swap[0,1]=swap[1,0]=1
 def change(P):return(swap*P.subs({h:x,x:h},simultaneous=True)*swap.T).applyfunc(s.expand)
 A=[change(canonical(raw(f,0,0)))for f in[h*h,h*x,x*x]]+[change(canonical(raw(0,1,0)))]+[change(canonical(raw(0,0,H)))for H in cubics]
 B=[change(canonical(rawB(f,0,0)))for f in[h*h,h*x,x*x]]+[change(canonical(rawB(0,h*h,0)))]+[change(canonical(rawB(0,0,H)))for H in cubics];N=[A[0],A[1],*A[3:]]
 D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T;u0=vector(P0);assert canonical(P0)==P0
 O=s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P0))for i,j in itertools.product(range(5),repeat=2)]);oi=list(O.rref()[1]);assert len(oi)==19;OC=O[:,oi];HF=s.Matrix.hstack(*[flags(s.eye(5)[:,i//5]*s.eye(5)[i%5,:])for i in oi])
 # All local flags after the whole intrinsic factor L has been recovered.
 aa,bb,cc,dd,ee=s.symbols('a b c d e');chart=s.eye(5);chart[0,3]=dd;chart[0,4]=ee;chart[1,0]=aa;chart[1,3]=bb;chart[1,4]=cc;chartP=(chart*P0.subs(dict(zip(X,chart.inv()*s.Matrix(X))),simultaneous=True)*chart.T).applyfunc(s.expand);cv=vector(chartP);fiber={}
 for name,basis in [('CA',A),('CB',B),('N',N)]:
  C=s.Matrix.hstack(*map(vector,basis));rr=list(C.T.rref()[1]);eq=(cv-C*C[rr,:].inv()*cv[rr,0]).applyfunc(s.expand);polys=[f for f in eq if f!=0];gb=s.groebner(polys,bb,cc,dd,ee,aa,order='lex');expected=[bb,cc,dd+(2 if name!='CB'else-s.Rational(1,2))*aa,ee,aa*aa];assert list(gb)==expected
  fiber[name]=dict(variables=list(map(str,[bb,cc,dd,ee,aa])),actual_constraints=list(map(str,polys)),entire_groebner_basis=list(map(str,gb)),scheme_length=2)
 # Two exact10-parameter source families; H19 acts outside these core maps.
 p=s.Symbol('p');K=s.symbols('k0:9');extras=[1,2,3,9,10,14,15,18,19];base=s.zeros(20,1);base[0]=-1;base[7]=-1;base[11]=1;base[16]=s.Rational(1,3);records={};Js=[];Ws=[]
 for name,basis,m in [('CA',A,-2),('CB',B,s.Rational(1,2))]:
  M=s.zeros(5);M[1,0]=1;M[0,3]=m;assert M**3==s.zeros(5);C=s.Matrix.hstack(*map(vector,basis));rr=list(C.T.rref()[1]);co=C[rr,:].inv()*vector(action(M,P0))[rr,0];assert C*co==vector(action(M,P0));alpha=base-p*co
  for k,j in zip(K,extras):alpha[j]+=k
  f=alpha[0]*h*h+alpha[1]*h*x+alpha[2]*x*x;H=sum(v*mon for v,mon in zip(alpha[4:,0],cubics));source=change(canonical(raw(f,alpha[3],H)))if name=='CA'else change(canonical(rawB(f,-alpha[3]*f,H)))
  assert source.subs(dict.fromkeys([p,*K],0))==P0;g=s.eye(5)+p*M+p*p*M*M/2;gi=s.eye(5)-p*M+p*p*M*M/2;U=vector((g*source.subs(dict(zip(X,gi*s.Matrix(X))),simultaneous=True)*g.T).applyfunc(s.expand));assert U.diff(p).subs(dict.fromkeys([p,*K],0))==s.zeros(150,1)
  J=OC.row_join(U.jacobian(K).subs(dict.fromkeys([p,*K],0)));assert J.rank()==28;aug=J.row_join(s.zeros(150,1)).col_join(HF.row_join(s.zeros(20,9)).row_join(flags(M)));assert aug.rank()==29
  W=s.Matrix([s.Poly(v,p).nth(2).subs(dict.fromkeys(K,0))for v in U-u0]);Js.append(J);Ws.append(W);records[name]=dict(core_parameters=['p',*map(str,K)],source_fiber_basis_indices=extras,source_kernel_Lie_matrix=mj(M),source_velocity=list(map(str,co)),source_alpha_coefficients=list(map(str,alpha)),core_full150_coordinates=list(map(str,U)),core_group=mj(g),source_etale_augmented_rank=29,entire_source_tangent=mj(aug),linear_J=mj(J),raw_second_coefficient=list(map(str,W)))
 rows=list(Js[0].T.rref()[1]);assert len(rows)==28;dets=[s.factor(J[rows,:].det())for J in Js];assert dets==[-s.Rational(2304,25),-s.Rational(1536,25)];proj=[(W-J*J[rows,:].inv()*W[rows,0]).applyfunc(s.expand)for J,W in zip(Js,Ws)];extra=61;coeffs=[v[extra]for v in proj];assert coeffs==[-4,s.Rational(3,8)];selected=rows+[extra];assert 106 in selected
 # Canonical150 lies in a126-dimensional linear space. Extend the29 actual
 # selected functions to126 raw coordinates for an exact common target T.
 div=s.zeros(25,150)
 for column,((i,j),(m,n))in enumerate(itertools.product(PAIRS,MONS)):
  for ell in [m,n]:
   other=n if ell==m else m
   if ell==i:div[5*j+other,column]+=1
   if ell==j:div[5*i+other,column]-=1
 canonical_basis=s.Matrix.hstack(*div.nullspace());assert canonical_basis.shape==(150,126);order=selected+[i for i in range(150)if i not in selected];piv=list(canonical_basis[order,:].T.rref()[1]);free=[order[i]for i in piv];assert free[:29]==selected and len(free)==126
 reconstruction=canonical_basis*canonical_basis[free,:].inv();assert reconstruction[free,:]==s.eye(126)
 files=[Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py',ROOT/'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63/report.json',ROOT/'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63/report.json']
 result=dict(id='AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63',all_passed=True,field='QQ',generic_B18_point_field=list(map(str,D)),open_orbit_projective_dimension=18,whole_flag_fibers=fiber,orbit_GL_generator_indices=oi,cubic_monomials=list(map(str,cubics)),families=records,common28_linear_target_rows=rows,common28_linear_minors=list(map(str,dets)),common29th_target_row=extra,corrected_quadratic_coefficients=list(map(str,coeffs)),all150_projected_second_coefficients=[list(map(str,v))for v in proj],common_target126_raw_coordinate_indices=free,canonical_reconstruction_from126=mj(reconstruction),common_normalizer_base_cone_dimension=29,common_normalizer_base_projective_dimension=28,both_normalizers_finite_free_rank=2,normalizer_module_basis=['1','p'],projective_scalar_chart='u106=1; this removes one common base parameter',actual_intersection_exact_reduction='T/(ker(phi_CA)+ker(phi_CB)); each phi_i maps97 remaining target generators into its FREE R28-rank2 quadratic normalizer; each kernel is a ONE-variable quadratic elimination, not a normalizer tensor product.',quadratic_coefficients_not_expanded=True,full_local_primary_type_claimed=False,scope='Open B18 orbit in R intersectionN; actualD exact quadratic-image reduction, no completed-ring/nilmodule classification or allboundary extension.',input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in files},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('GENERIC B18 ACTUAL-IMAGE QUADRATIC REDUCTION PASS',result['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
