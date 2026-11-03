#!/usr/bin/env python3
"""Whole finite image algebra at the exceptional C15 orbit of reduced R19.

A polynomial four-source slice generates QQ[B,E,F,p2,p3,Bp,Ep,Fp].
The full 23-dimensional algebra cotangent, together with an etale19 source
chart, proves the completed image statement in the companion proof note.
"""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,vector,action,canonical
ROOT=Path(__file__).resolve().parents[1]
h,x,y,z,w=X;ey=s.Matrix([0,0,1,0,0]);PAIRS=list(itertools.combinations(range(5),2));TRIPLES=list(itertools.combinations(range(5),3))
def biv(D):return (ey*D.T-D*ey.T).applyfunc(s.expand)
def minor(M):
 cols=list(M.rref()[1]);rows=list(M[:,cols].T.rref()[1]);det=M[rows,cols].det();assert det!=0
 return dict(rank=len(cols),rows=rows,columns=cols,determinant=str(det),matrix=[[str(v) for v in row] for row in M.tolist()])
def flags(M):
 # Tangents to Lline, Xline, and wedge3K modulo their own scalar lines.
 L=M*ey;L[2]=0;cx=-M[1,:].T;cx[1]=0;K=[]
 for triple in TRIPLES:
  val=0
  for k in range(3):
   A=s.zeros(3,5)
   for i,row in enumerate([0,1,3]):
    A[i,row]=1
   A[k,:]=-M[[0,1,3][k],:]
   val+=A[:,list(triple)].det()
  K.append(val)
 K[TRIPLES.index((0,1,3))]=0
 return s.Matrix(list(L)+list(cx)+K)
def main(output):
 output.mkdir(parents=True,exist_ok=False);started=time.monotonic();p,B,E,F=s.symbols('p B E F');t=-p/12;r=11*p/12;c=-2*p*p/3
 G=h*h+(B-p/6)*h*z+(E-s.Rational(1,6))*x*z+F*z*z
 D=s.Matrix([x*((r-13*t)*h+x+c*z),12*t*x*x,0,x*(h+(-r-13*t)*z),2*t*x*w+G]);P=biv(D)
 assert canonical(P)==P and s.expand(r*r+c-25*t*t)==0
 M=s.zeros(5);M[3,4]=M[1,0]=1;assert M*M==s.zeros(5)
 g=s.eye(5)+p*M;gi=s.eye(5)-p*M;Q=(g*P.subs(dict(zip(X,gi*s.Matrix(X))),simultaneous=True)*g.T).applyfunc(s.expand);assert canonical(Q)==Q
 # Q is y-independent ey wedge a divergence-zero field, so raw Jacobi vanishes.
 assert all(Q[i,j]==0 for i,j in PAIRS if 2 not in (i,j));assert not Q.has(y)
 u=vector(Q);zero={p:0,B:0,E:0,F:0};u0=u.subs(zero);P0=Q.subs(zero);assert u[106]==1
 # v means target coefficient minus its central value.
 v=u-u0
 # Actual target coefficients recover all eight algebra generators exactly.
 recovered=[v[123]+v[112],v[127],v[132],v[15],-s.Rational(3,2)*v[18],v[108]+v[127]*v[15]+v[15]/6,v[112],v[117]]
 generators=[B,E,F,p*p,p**3,B*p,E*p,F*p];assert all(s.expand(a-b)==0 for a,b in zip(recovered,generators))
 # Every monomial of every target output belongs to the proposed algebra.
 # For an odd p-exponent, either p>=3 or one of B,E,F must occur.
 for value in u-u0:
  for (ip,iB,iE,iF),coef in s.Poly(value,p,B,E,F).terms():
   assert ip%2==0 or ip>=3 or iB+iE+iF>0
 C=s.Matrix([[s.Poly(v,p,B,E,F).coeff_monomial(mon) for mon in generators] for v in u-u0]);assert C.rank()==8
 O=s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P0)) for i,j in itertools.product(range(5),repeat=2)]);OP=O-u0*O[106,:];assert OP.rank()==15
 total=OP.row_join(C);W=minor(total);assert W['rank']==23 and W['determinant']=='-209/162'
 orbit_indices=[i for i in W['columns'] if i<25];assert len(orbit_indices)==15
 H=OP[:,orbit_indices];assert H.row_join(C[:,:3]).rank()==18
 # Full source incidence tangents, including the nonzero vertical p flag.
 HF=s.Matrix.hstack(*[flags(s.eye(5)[:,i//5]*s.eye(5)[i%5,:]) for i in orbit_indices]);pflag=flags(M);assert pflag!=s.zeros(20,1)
 source_tangent=H.row_join(C[:,:3]).row_join(s.zeros(150,1)).col_join(HF.row_join(s.zeros(20,3)).row_join(pflag));SW=minor(source_tangent);assert SW['rank']==19
 aa,bb,ee,ff=s.symbols('a b e f');j=[aa,bb,ee,ff];zz=s.symbols('z0:4');relations=[j[i]*zz[k]-j[k]*zz[i] for i in range(4) for k in range(i+1,4)]+[zz[i]*zz[k]-aa*j[i]*j[k] for i in range(4) for k in range(i,4)]
 sub={aa:p*p,bb:B,ee:E,ff:F,zz[0]:p**3,zz[1]:B*p,zz[2]:E*p,zz[3]:F*p};assert all(s.expand(v.subs(sub))==0 for v in relations)
 paths=[Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_R19_boundary_local_models_20261103_v63.py',ROOT/'reports/AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63/report.json',ROOT/'research/notes/AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63.md']
 result=dict(id='AMBIENT_R19_C15_COMPLETED_IMAGE_20261103_v63',all_passed=True,field='QQ',spectral_slice_parameters=dict(t=str(t),r=str(r),b='1',c=str(c),d='1',e='0',G=str(G)),finite_group=[[str(v) for v in row] for row in g.tolist()],full150_target=list(map(str,u)),actual_generator_formulas=['v123+v112','v127','v132','v15','-3*v18/2','v108+v127*v15+v15/6','v112','v117'],generators=list(map(str,generators)),generator_coefficient_matrix=[[str(v) for v in row] for row in C.tolist()],projective_chart='u106=1',projective_orbit_rank=15,projective_orbit_generator_indices=orbit_indices,entire_algebra_cotangent=W,normalization_source_etale_tangent=SW,vertical_flag_derivative=list(map(str,pflag)),finite_presentation_variables=list(map(str,[aa,bb,ee,ff,*zz])),finite_presentation_relations=list(map(str,relations)),local_ring='C[[u1,...,u15,B,E,F,p²,p³,Bp,Ep,Fp]]',normalization='C[[u1,...,u15,B,E,F,p]]',dimension=19,embedding_dimension=23,depth=16,multiplicity=2,conductor_in_normalization='(B,E,F,p²) * normalization',conductor_in_target='(B,E,F,p²,p³,Bp,Ep,Fp) * target',normalization_quotient='C[[u1,...,u15]] * p',normalization_fiber_length=2,boundary_B18_reduced_completion='C[[u1,...,u15,B,E,F]]',scope='Full completed reduced R image along C15, not full Poisson or actual CA/CB scheme-intersection boundary.',input_sha256={str(v.relative_to(ROOT)):hashlib.sha256(v.read_bytes()).hexdigest() for v in paths},seconds=time.monotonic()-started)
 (output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('C15 FULL FINITE IMAGE + ALGEBRA COTANGENT23 + SOURCE ETALE19 PASS',result['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
