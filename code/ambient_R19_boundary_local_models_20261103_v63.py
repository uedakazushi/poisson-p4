#!/usr/bin/env python3
"""Exact QQ tangent certificates for generic B17 and a ramified boundary locus.

Proper-fiber and full local-ring arguments are given in the paired proof note.
"""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,vector,action
ROOT=Path(__file__).resolve().parents[1]
h,x,y,z,w=X;ey=s.Matrix([0,0,1,0,0])
def biv(D): return (ey*D.T-D*ey.T).applyfunc(s.expand)
def orbit(P): return s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P)) for i,j in itertools.product(range(5),repeat=2)])
def witness(M):
 cols=list(M.rref()[1]);rows=list(M[:,cols].T.rref()[1]);det=M[rows,cols].det();assert det!=0
 return dict(rank=len(cols),rows=rows,columns=cols,determinant=str(det),matrix=[[str(v) for v in row] for row in M.tolist()])
def main(output):
 output.mkdir(parents=True,exist_ok=False);started=time.monotonic();base=ROOT/'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63/report.json';br=json.loads(base.read_text())
 C=s.Matrix.hstack(*[s.Matrix([s.sympify(v) for v in col]) for col in br['spectral_basis_raw150']]);assert C.rank()==12
 G=h*h+2*h*x+3*h*z+4*x*x+5*x*z+6*z*z
 P=biv(s.Matrix([0,0,0,-h*x,G]));O=orbit(P)
 swap=s.eye(5);swap[0,0]=swap[1,1]=0;swap[0,1]=swap[1,0]=1
 def unvector(col):
  mons=[X[i]*X[j] for i,j in itertools.combinations_with_replacement(range(5),2)];Q=s.zeros(5)
  for k,(i,j) in enumerate(itertools.combinations(range(5),2)):
   Q[i,j]=sum(col[15*k+l]*m for l,m in enumerate(mons));Q[j,i]=-Q[i,j]
  return Q
 def change(col):
  Q=unvector(col);return vector((swap*Q.subs(dict(zip(X,swap*s.Matrix(X))),simultaneous=True)*swap.T).applyfunc(s.expand))
 keep=[i for i in range(12) if i!=3];nilkeep=[i for i in range(12) if i not in (0,3)]
 J1=C[:,keep].row_join(O);J2=s.Matrix.hstack(*[change(C[:,i]) for i in keep]).row_join(O)
 K1=C[:,nilkeep].row_join(O);K2=s.Matrix.hstack(*[change(C[:,i]) for i in nilkeep]).row_join(O)
 R1=witness(J1);R2=witness(J2);RU=witness(J1.row_join(J2));B1=witness(K1);B2=witness(K2);BU=witness(K1.row_join(K2))
 assert [R1['rank'],R2['rank'],RU['rank'],B1['rank'],B2['rank'],BU['rank']]==[20,20,22,19,19,20]
 f0,f1,f2=s.symbols('f0 f1 f2');gg=s.symbols('g0:6');GT=sum(v*m for v,m in zip(gg,[h*h,h*x,h*z,x*x,x*z,z*z]));T=biv(s.Matrix([0,0,0,-f0*h*h-f1*h*x-f2*x*x,GT]));TC=s.Matrix.hstack(*[vector(T.diff(v)) for v in [f0,f1,f2,*gg]]);TJ=TC.row_join(O);TW=witness(TJ);assert TW['rank']==18
 # All surviving coefficients in the one-Casimir chain chart.  After a=b=1,
 # c=0 and removal of hx,x2, G=A*h2+B*hz+E*xz+F*z2.
 A,B,E,F=s.symbols('A B E F');G0=A*h*h+B*h*z+E*x*z+F*z*z;P0=biv(s.Matrix([x*x,0,0,h*x,G0]));J=C[:,keep].row_join(orbit(P0));Jzero=J.subs({A:0,B:0,E:0,F:0});cols=list(Jzero.rref()[1]);rows=list(Jzero[:,cols].T.rref()[1]);V=J[:,cols];det=s.expand(V[rows,:].det());assert det==48
 residual=(J-V*V[rows,:].inv()*J[rows,:]).applyfunc(s.expand)
 nz=[(i,j,str(residual[i,j])) for i in range(residual.rows) for j in range(residual.cols) if residual[i,j]!=0]
 assert nz==[(108,30,'B'),(112,30,'A/6 + E'),(117,30,'F'),(124,30,'-B'),(128,30,'-A/6 - E'),(133,30,'-2*F')]
 special=J.subs({A:1,B:0,E:-s.Rational(1,6),F:0});SP=witness(special);assert SP['rank']==19
 specialP=P0.subs({A:1,B:0,E:-s.Rational(1,6),F:0});CW=witness(orbit(specialP));assert CW['rank']==16
 # Finite normalization of an arbitrary ab!=0 chain; source/target coordinates
 # are kept distinct by the inverse substitution.
 aa,bb,cc=s.symbols('a b c',nonzero=True);g=s.symbols('q0:6');Gc=sum(v*m for v,m in zip(g,[h*h,h*x,h*z,x*x,x*z,z*z]));Pc=biv(s.Matrix([aa*x*x,0,0,bb*h*x+cc*x*x,Gc]));S=s.eye(5);S[0,0]=1/aa;S[3,3]=1/(aa*bb);S[3,0]=-cc/(aa*aa*bb)
 Pc1=(S*Pc.subs(dict(zip(X,S.inv()*s.Matrix(X))),simultaneous=True)*S.T).applyfunc(s.cancel);newG=s.expand(Pc1[2,4]);assert Pc1[2,0]==x*x and Pc1[2,3]==h*x
 mu=s.expand(newG).coeff(h,1).coeff(x,1);nu=s.expand(newG).coeff(x,2);S2=s.eye(5);S2[4,3]=-mu;S2[4,0]=-nu;Pc2=(S2*Pc1.subs(dict(zip(X,S2.inv()*s.Matrix(X))),simultaneous=True)*S2.T).applyfunc(s.expand);assert Pc2[2,4].coeff(h,1).coeff(x,1)==0 and Pc2[2,4].coeff(x,2)==0
 paths=[Path(__file__),base,ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'research/notes/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63.md']
 result=dict(id='AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63',all_passed=True,field='QQ',B17_point_field=[str(v) for v in [0,0,-h*x,G]],R_branch1=R1,R_branch2=R2,R_union=RU,B18_branch1=B1,B18_branch2=B2,B18_union=BU,B17_tangent=TW,
  R_projective_branch_dimensions=[19,19],R_projective_common_tangent_dimension=17,R_projective_tangent_dimension=21,B18_projective_branch_dimensions=[18,18],B18_projective_common_tangent_dimension=17,
  one_Casimir_normalized_chart='D=(x2,0,hx,A*h2+B*hz+E*xz+F*z2)',constant19_minor=dict(rows=rows,columns=cols,determinant=str(det)),entire_residual_entries=nz,
  rank20_criterion='(B,F,A/6+E) != (0,0,0)',exceptional_point_field='D=(x2,0,hx,h2-xz/6)',exceptional_point=SP,exceptional_orbit=CW,chain_normalization_forward=[[str(v) for v in row] for row in S.tolist()],chain_removed_hx_coefficient=str(mu),chain_removed_x2_coefficient=str(nu),chain_final_field=[str(Pc2[2,i]) for i in [0,1,3,4]],
  scope='Generic B17 completed reduced R/B18 images and explicit further ramification; not full Poisson ring or CA/CB scheme intersection on the boundary.',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.monotonic()-started)
 (output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('GENERIC B17 BRANCHES + EXACT RAMIFICATION CHART PASS',result['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
