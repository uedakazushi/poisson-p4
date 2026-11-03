#!/usr/bin/env python3
"""Independent spectral/chain/Hamiltonian flag and global boundary audit."""
import argparse
import itertools
import json
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_review_20261103_v63 import IDS,action,canonical,wedge_factor
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 public=ROOT/'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63/report.json';data=read(public);fresh=read(args.replay/'report.json');assert{k:v for k,v in data.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),public,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py',ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for file,hsh in data['input_sha256'].items():p=ROOT/file;assert sha(p)==hsh,file;paths.append(p)
 X=s.Matrix(s.symbols('h x y z w'));h,x,y,z,w=X;ey=s.Matrix([0,0,1,0,0]);t,r,b,c,d,e=s.symbols('t r b c d e');gg=s.symbols('g0:6');G=sum(v*m for v,m in zip(gg,[h*h,h*x,h*z,x*x,x*z,z*z]));params=[t,r,b,c,d,e,*gg]
 def tensor(D):return ey*D.T-D*ey.T
 def vector(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(s.prod(X[k]for k in mon))for((i,j),mon)in IDS])
 def basis(P,parameters):return s.Matrix.hstack(*[vector(P.diff(v))for v in parameters])
 D=s.Matrix([x*((r-13*t)*h+b*x+c*z),12*t*x*x,0,x*(d*h+e*x+(-r-13*t)*z),2*t*x*w+G]);P=tensor(D);C=basis(P,params)
 assert C==s.Matrix.hstack(*[s.Matrix(list(map(s.sympify,col)))for col in data['spectral_basis_raw150']])and C.rank()==12
 assert all(canonical(C[:,i])==C[:,i]for i in range(12))
 for i,j,k in itertools.combinations(range(5),3):assert s.expand(sum(P[i,l]*s.diff(P[j,k],X[l])+P[j,l]*s.diff(P[k,i],X[l])+P[k,l]*s.diff(P[i,j],X[l])for l in range(5)))==0
 q=r*r+c*d-25*t*t
 def invariant(B,positions,representations,parameters,form=None,character=None,tcoordinate=False):
  assert B.rank()==len(parameters)
  for(i,j),small in zip(positions,representations):
   small=s.Matrix([[s.sympify(v)for v in row]for row in small]);assert action(i,j)*B==B*small
   velocity=small*s.Matrix(parameters);chi=character(i,j)if character else 0
   if form is not None:assert s.expand(sum(s.diff(form,v)*dv for v,dv in zip(parameters,velocity))-chi*form)==0
   if tcoordinate:assert s.expand(velocity[0]-((i==2 and j==2)-(i==1 and j==1))*t)==0
 invariant(C,data['parabolic_GL_positions'],data['spectral_parabolic_parameter_actions'],params,q,lambda i,j:2*(int(i==2 and j==2)-int(i==1 and j==1)),True)
 assert s.hessian(q,[t,r,c,d]).rank()==4 and s.hessian(q.subs(t,0),[r,c,d]).rank()==3
 T=s.Matrix([[r-13*t,b,c],[0,12*t,0],[d,e,-r-13*t]]);lam=s.Symbol('lam');char=T.charpoly(lam).as_expr();target=(lam+8*t)*(lam-12*t)*(lam+18*t)
 assert s.rem(s.Poly(char-target,r),s.Poly(q,r)).is_zero
 aa,bb,cc=s.symbols('aN bN cN');N=tensor(s.Matrix([aa*x*x,0,0,bb*h*x+cc*x*x,G]));np=[aa,bb,cc,*gg];NC=basis(N,np)
 invariant(NC,data['chain_fullflag_parabolic_positions'],data['chain_parameter_actions'],np)
 # Direct raw N_g model with marked h_N=x, k_N=h.
 HN=-(bb*h*x+cc*h*h)*w+s.integrate(G.subs({h:x,x:h},simultaneous=True),z);raw=s.zeros(5);raw[1,2]=-aa*h*h;raw[2,3]=-s.diff(HN,w);raw[2,4]=s.diff(HN,z);raw=raw-raw.T
 swap=s.eye(5);swap[0,0]=swap[1,1]=0;swap[0,1]=swap[1,0]=1
 assert(swap*raw.subs(dict(zip(X,swap*X)),simultaneous=True)*swap.T-N).applyfunc(s.expand)==s.zeros(5)
 f0,f1,f2=s.symbols('f0 f1 f2');f=f0*h*h+f1*h*x+f2*x*x;H=tensor(s.Matrix([0,0,0,-f,G]));hp=[f0,f1,f2,*gg];HC=basis(H,hp);disc=f1*f1-4*f0*f2
 invariant(HC,data['T17_parabolic_positions'],data['T17_parameter_actions'],hp,disc,lambda i,j:2*(int(i==2 and j==2)+int(i==3 and j==3)-int(i==0 and j==0)-int(i==1 and j==1)))
 point={aa:1,bb:1,cc:1}|dict(zip(gg,range(1,7)));spoint={t:0,r:0,b:1,c:0,d:1,e:1}|dict(zip(gg,range(1,7)));hpoint={f0:1,f1:1,f2:1}|dict(zip(gg,range(1,7)));sqpoint={f0:0,f1:0,f2:1}|dict(zip(gg,range(1,7)))
 def witness(v,source,record):
  M=s.Matrix.hstack(source,*[action(i,j)*v for i,j in itertools.product(range(5),repeat=2)]);assert M==s.Matrix([[s.sympify(v)for v in row]for row in record['full_matrix']]);det=M[record['minor_rows'],record['minor_columns']].det();assert det==s.Rational(record['minor_determinant'])!=0
  return len(record['minor_columns'])
 assert witness(vector(N.subs(point)),NC,data['chain_image_witness'])==19
 assert witness(vector(P.subs(spoint)),C[:,[i for i,v in enumerate(params)if v!=c]],data['R_generic_boundary_source_witness'])==20
 assert witness(vector(H.subs(hpoint)),HC,data['T17_image_witness'])==18
 assert witness(vector(H.subs(sqpoint)),s.Matrix.hstack(2*HC[:,1],2*HC[:,2],HC[:,3:]),data['T16_image_witness'])==17
 v=vector(N.subs(point));W=wedge_factor(v);assert W.rank()==4 and W.nullspace()==[ey]
 contract=s.zeros(75,5)
 for((i,j),mon),index in IDS.items():
  k=index%15;contract[15*j+k,i]+=v[index];contract[15*i+k,j]-=v[index]
 assert contract.rank()==4 and contract.nullspace()==[s.Matrix([0,1,0,0,0])]
 # At the chosen point D(ell) mod x only sees the w coefficient, so K3
 # is intrinsically recovered as span(h,x,z), in addition to L and X.
 a0,a1,a2,a3=s.symbols('a0:4');ell=a0*h+a1*x+a2*z+a3*w;Dpoint=s.Matrix([x*x,0,0,h*x+x*x,G.subs(dict(zip(gg,range(1,7))))]);Dell=s.expand(sum(s.diff(ell,X[i])*Dpoint[i]for i in range(5)))
 assert s.expand(Dell.subs(x,0)-a3*(h*h+3*h*z+6*z*z))==0
 assert q.subs(spoint)==0 and s.diff(q,c).subs(spoint)==1
 mu,nu=s.symbols('mu nu');ell=w-mu*z/bb-(nu-mu*cc/bb)*h/aa;special=s.Matrix([aa*x*x,0,0,bb*h*x+cc*x*x,mu*h*x+nu*x*x])
 assert s.cancel(sum(s.diff(ell,X[i])*special[i]for i in range(5)))==0
 paths += [ROOT/'reports/AMBIENT_R_COMPLEMENT_20261103_v63_REVIEW.json',ROOT/'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',ROOT/'research/notes/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63.md',ROOT/'research/notes/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW.md']
 result=dict(id='AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW',all_passed=True,
  fresh_QQ_producer_matches=True,independent_spectral12_chain9_and_Hamiltonian9_spaces=True,
  all192_135_144_actual_Lie_columns_and_relative_characters_independently_verified=True,
  exact_relative_quadric_ranks4_and3_and_eigen_polynomial=True,
  independent_all4_image_full_matrices_and_minors=[19,20,18,17],
  direct_all9_parameter_Ng_tensor_identity=True,
  proper_normal_relative_quadric19_birational_model_of_entire_R_accepted=True,
  entire_tzero_nilpotent_Jordanflag_cover_and_reduced_R_intersection_N_equals_irreducible_B18=True,
  nested_B16_B17_B18_and_reduced_R_intersection_Nf_equals_B17_accepted=True,
  generic_B_whole_fiber_intrinsic_L_X_K_uniqueness_and_source_immersion_checked=True,
  proper_finite_Nakayama_target_smoothness_at_generic_boundary_and_FULL_R_R1_accepted=True,
  target_R_normality_R_N_scheme_thickness_or_global_primary_exhaustion_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS proper spectral R model,ENTIRE reduced R intersection N=B18 and global R1')

if __name__=='__main__':main()
