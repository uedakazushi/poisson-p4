#!/usr/bin/env python3
"""C0 exact S61 first-jet bridge and two honest dominant-slice certificates.

No full134rebuild, Fisher computation, GB, component exhaustion, or secondary jets.
Computes the original150coefficient bridge, both inverse identities, all160Lie
columns via independent finite pullback, and a sparse6 nonzero QQ Jacobian.
"""
import argparse,hashlib,json,resource,signal,sys,time
from pathlib import Path
import sympy as s
from s61_model import S61,BASIS,X
from weight1_structure import J,COORDS,constraints,infinitesimal_action,sp_basis
ROOT=Path(__file__).resolve().parents[1]
def main(output):
 root=ROOT
 start=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(1500000000,1500000000))
 model=S61();w=model.reduced['parameter_weights'];ids=[i for i,x in enumerate(w)if x==1];assert len(ids)==16
 raw=s.zeros(150,16)
 for ri,uid in zip(model.certificate['raw_free_indices'],range(61)):
  if uid in ids:raw[ri,ids.index(uid)]=1
 for rec in model.data['eliminated_raw_variables']:
  for word,num,den in rec['polynomial']:
   if len(word)==1 and word[0] in ids:raw[rec['raw_index'],ids.index(word[0])]+=s.Rational(num,den)
 # Derive first affine jet after Euler dehomogenization, directly from150 coordinates.
 ri={key:i for i,key in enumerate(BASIS)}
 def cf(pair,mon,col):
  sign=1
  if pair[0]>pair[1]:pair=pair[::-1];sign=-1
  return sign*raw[ri[pair,mon],col]
 y=s.Matrix(s.symbols('y1:5'));fields=[];betas=[];Ls=[]
 for col in range(16):
  L=s.zeros(4)
  for i in range(4):
   for j in range(i+1,4):
    for k in range(4):
     value=cf((i+1,j+1),X[0]*X[k+1],col)
     if i==k:value-=cf((0,j+1),X[0]**2,col)
     if j==k:value+=cf((0,i+1),X[0]**2,col)
     L[i,j]+=value*y[k]
    L[j,i]=-L[i,j]
  tau=L[0,1]+L[2,3];a=s.Matrix([s.diff(tau,t)/5 for t in y]);v=J*a
  prim=L-2*(a.dot(y))*J+y*v.T-v*y.T
  assert s.expand(prim[0,1]+prim[2,3])==0
  beta=-J*prim*J;betas.append(s.Matrix([s.expand(beta[i,j]).coeff(y[k])for(i,j),k in COORDS]));Ls.append(L)
 M=s.Matrix.hstack(*betas);D=s.Matrix.vstack(*constraints());assert D*M==s.zeros(8,16)
 B=s.Matrix.hstack(*D.nullspace());rows=list(B.T.rref()[1]);C=B[rows,:].inv()*M[rows,:];assert B*C==M
 Ci=C.inv();assert Ci*C==C*Ci==s.eye(16)
 actions=[]
 for a in sp_basis():
  imgs=s.Matrix.hstack(*(infinitesimal_action(B[:,i],a)for i in range(16)));act=B[rows,:].inv()*imgs[rows,:];assert B*act==imgs;actions.append(act)
 v=s.Matrix(range(1,17));orbit=s.Matrix.hstack(*(a*v for a in actions));S=s.Matrix.hstack(v,*(s.eye(16)[:,i]for i in range(11,16)));T=orbit.row_join(S);det=T.det(method='domain-ge');assert det==8428695240
 # Independent Jacobian check by finite pullback at first order for ten sp generators.
 t=s.symbols('t');ys=list(y)
 for a,act in zip(sp_basis(),actions):
  A=s.eye(4)+t*a
  for col in range(16):
   b=s.zeros(4)
   for val,((i,j),k) in zip(B[:,col],COORDS):b[i,j]+=val*y[k];b[j,i]-=val*y[k]
   shifted=b.subs(dict(zip(ys,A*y)),simultaneous=True);pull=A.T*shifted*A
   vv=s.Matrix([s.expand(pull[i,j]).coeff(t,1).coeff(y[k],1)for(i,j),k in COORDS]);assert vv==B*act[:,col]
 def enc(M):return [[str(x)for x in row]for row in M.tolist()]
 paths=['code/s61_model.py','code/weight1_structure.py','code/primitive_curve_dominant_slice_20261101.py','legacy/v02/data/symplectic_slice_certificate.json','legacy/v02/data/symplectic_slice_equations.json','legacy/v02/data/symplectic_slice_reduced_equations.json','reports/PRIMITIVE_CURVE_DOMINANT_SLICE_20261101.json']
 pin={p:hashlib.sha256((root/p).read_bytes()).hexdigest()for p in paths}
 old=json.loads((root/paths[-1]).read_text());assert s.Matrix(old['dominance_tangent_matrix'])==T and s.Matrix(old['slice_matrix_b'])==S
 out={'all_passed':True,'original_weight1_parameter_indices0':ids,'primitive_basis_coordinate_bridge16':enc(C),'bridge_determinant':str(C.det()),'original_S61_weight1_slice16x6':enc(Ci*S),'primitive_embedding24x16':enc(M),'tangent_determinant':str(det),'ten_actions_all16_columns_independent_first_order_finite_pullback':True,'both_inverse_products':True,'inputs_sha256':pin,'seconds':time.monotonic()-start,'standard_tangent16x16':enc(T),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

 bridge_result=out
 r=ROOT
 p=r/'reports/Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_MINOR_20261104_v63/section6.json';src=json.loads(p.read_text());lam,H,L,ss,qq,eps=s.symbols('lambda H L s0 q0 epsilon');pars=(lam,H,L,ss,qq,eps)
 b=s.Matrix([s.sympify(x.replace('lambda','lam_parameter'),locals={'lam_parameter':lam})for x in src['sparse_beta6']]);D=s.Matrix.vstack(*constraints());B=s.Matrix.hstack(*D.nullspace());rows=list(B.T.rref()[1]);left=B[rows,:].inv();at={lam:1,H:1,L:1,ss:0,qq:0,eps:0};point=B*b.subs(at)
 orbit24=s.Matrix.hstack(*(infinitesimal_action(point,a)for a in sp_basis()));orbit=left*orbit24[rows,:];assert B*orbit==orbit24
 T=orbit.row_join(b.jacobian(pars).subs(at));det=T.det(method='domain-ge');expected=s.Rational(2**33,996451875)*16;assert det==expected
 prepath=r/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63/certificate.json';pre=json.loads(prepath.read_text());bridge=bridge_result;assert s.Matrix(pre['adapter_b16_equals_A_u16'])==s.Matrix(bridge['primitive_basis_coordinate_bridge16'])
 t=s.symbols('t0');b5=b.subs({lam:1,eps:t},simultaneous=True);assert b.subs({eps:lam*t},simultaneous=True)==lam*b5
 model=json.loads((r/'legacy/v02/data/symplectic_slice_reduced_equations.json').read_text());weights=model['parameter_weights'];distribution={}
 for equation in model['equations']:
  ww={sum(weights[i]for i in m)for m,n,d in equation};assert len(ww)==1;w=ww.pop();distribution[w]=distribution.get(w,0)+1
 assert distribution=={4:70,5:64}
 Cinv=s.Matrix(pre['adapter_u16_equals_Ainv_b16']);u=(Cinv*b5).applyfunc(s.expand);stored=s.Matrix([s.sympify(z)for z in pre['source_u16']]);assert u==stored
 paths=['reports/Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_MINOR_20261104_v63/section6.json','reports/Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_MINOR_20261104_v63_REVIEW/report.json','reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63/certificate.json','legacy/v02/data/symplectic_slice_reduced_equations.json']
 sparse_result={'all_passed':True,'Sp4_times_Gm_times_source5_dominates_W16':True,'Sp4_times_source5_alone_not_asserted':True,'witness':{str(k):str(v)for k,v in at.items()},'sparse6_tangent16x16':[[str(x)for x in row]for row in T.tolist()],'determinant':str(det),'expected_factored_at_witness':'2^33*16/996451875','weighted134_distribution':distribution,'lambda_strip_exact_beta6_identity':True,'primitive_original_adapter_matches_prior':True,'u16_input_matches_prior':True,'input_sha256':{p:hashlib.sha256((r/p).read_bytes()).hexdigest()for p in paths},'seconds':time.monotonic()-start}

 bridge_result['sparse5']=sparse_result
 output.parent.mkdir(parents=True,exist_ok=True)
 with output.open('x') as f:json.dump(bridge_result,f,indent=2);f.write('\n')
 print(json.dumps({k:bridge_result[k]for k in ['all_passed','bridge_determinant','tangent_determinant','seconds']}))
 print(json.dumps({'sparse5_determinant':sparse_result['determinant'],'all134_weight_covariance':True}))
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args();main(args.output)
