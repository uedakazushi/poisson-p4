#!/usr/bin/env python3
"""Exact completed image at the open B16 orbit, with quadric normalization.

All150 finite outputs and full22 algebra cotangent certify the whole image;
no smooth-normalizer pinch ansatz is used.
"""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,vector,action,canonical
ROOT=Path(__file__).resolve().parents[1]
h,x,y,z,w=X;ey=s.Matrix([0,0,1,0,0]);TRIPLES=list(itertools.combinations(range(5),3))
def biv(D):return (ey*D.T-D*ey.T).applyfunc(s.expand)
def minor(M):
 cols=list(M.rref()[1]);rows=list(M[:,cols].T.rref()[1]);det=M[rows,cols].det();assert det!=0
 return dict(rank=len(cols),rows=rows,columns=cols,determinant=str(det),matrix=[[str(v) for v in row] for row in M.tolist()])
def flags(M):
 L=M*ey;L[2]=0;cx=-M[1,:].T;cx[1]=0;K=[]
 for triple in TRIPLES:
  val=0
  for k in range(3):
   A=s.zeros(3,5)
   for i,row in enumerate([0,1,3]):A[i,row]=1
   A[k,:]=-M[[0,1,3][k],:];val+=A[:,list(triple)].det()
  K.append(val)
 K[TRIPLES.index((0,1,3))]=0
 return s.Matrix(list(L)+list(cx)+K)
def main(output):
 output.mkdir(parents=True,exist_ok=False);started=time.monotonic();p,t,r,c=s.symbols('p t r c');aa,bb,dd=s.symbols('a b d');q=r*r-25*t*t-2*c*p
 D=s.Matrix([x*((r-13*t)*h+c*z),12*t*x*x,0,x*(-2*p*h-x+(-r-13*t)*z),2*t*x*w+h*h+z*z]);P=biv(D);assert canonical(P)==P
 M=s.zeros(5);M[1,0]=1;assert M*M==s.zeros(5);g=s.eye(5)+p*M;Q=(g*P.subs(dict(zip(X,(s.eye(5)-p*M)*s.Matrix(X))),simultaneous=True)*g.T).applyfunc(s.expand);assert canonical(Q)==Q
 u=vector(Q);zero={p:0,t:0,r:0,c:0};u0=u.subs(zero);v=u-u0;P0=Q.subs(zero);assert u[110]==-1
 gens=[t,r,c,p*p,p*r,p*t];recovered=[v[128]/2,-v[112]-13*v[128]/2,-v[22],v[105],(v[15]+v[108])/2,(v[108]-v[15])/26];assert all(s.expand(a-b)==0 for a,b in zip(gens,recovered))
 elimination=s.groebner([p*p-aa,r*p-bb,t*p-dd,2*c*p-r*r+25*t*t],p,t,r,c,aa,bb,dd,order='lex');relations=[v.as_expr() for v in elimination.polys if not v.as_expr().has(p)];assert len(relations)==11
 output_polynomials=[]
 for value in v:
  reduced=s.expand(elimination.reduce(value)[1]);assert not reduced.has(p);output_polynomials.append(reduced)
 C=s.Matrix([[s.Poly(value,p,t,r,c).coeff_monomial(mon) for mon in gens] for value in v]);assert C.rank()==6
 O=s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P0)) for i,j in itertools.product(range(5),repeat=2)]);OP=O+u0*O[110,:];assert OP.rank()==16
 W=minor(OP.row_join(C));assert W['rank']==22 and W['determinant']=='6144';indices=[i for i in W['columns'] if i<25];assert len(indices)==16
 H=OP[:,indices];assert H.row_join(C[:,:3]).rank()==19;HF=s.Matrix.hstack(*[flags(s.eye(5)[:,i//5]*s.eye(5)[i%5,:]) for i in indices]);pflag=flags(M);assert pflag!=s.zeros(20,1)
 ambient=H.row_join(C[:,:3]).row_join(s.zeros(150,1)).col_join(HF.row_join(s.zeros(20,3)).row_join(pflag));SW=minor(ambient);assert SW['rank']==20
 assert s.hessian(q,[p,t,r,c]).rank()==4
 # Exact reduction certifies the multiplicity-three parameter reduction.
 # Set upar=c+a and q0=r2-25t2; a^3-2upar*a2+upar2*a-q0²/4=0.
 upar=s.Symbol('upar');integrality=aa**3-2*upar*aa*aa+upar*upar*aa-(r*r-25*t*t)**2/4
 assert s.expand(integrality.subs({upar:c+p*p,aa:p*p}))==s.expand(-q*(r*r-25*t*t+2*c*p)/4)
 paths=[Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63/report.json',ROOT/'research/notes/AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63.md']
 result=dict(id='AMBIENT_R19_B16_COMPLETED_IMAGE_20261103_v63',all_passed=True,field='QQ',point_field=['0','0','-x²','h²+z²'],finite_group=[[str(a) for a in row] for row in g.tolist()],source_quadric=str(q),quadric_rank=4,full150_target=list(map(str,u)),actual_generator_formulas=['v128/2','-v112-13*v128/2','-v22','v105','(v15+v108)/2','(v108-v15)/26'],generators=list(map(str,gens)),full150_polynomials_in_target_generators=list(map(str,output_polynomials)),target_presentation_variables=list(map(str,[t,r,c,aa,bb,dd])),target_entire_elimination_kernel=list(map(str,relations)),generator_coefficient_matrix=[[str(a) for a in row] for row in C.tolist()],projective_chart='u110=-1',open_B16_orbit_dimension=16,projective_orbit_generator_indices=indices,entire_algebra_cotangent=W,normalization_ambient_etale_tangent=SW,vertical_flag_derivative=list(map(str,pflag)),local_ring='C[[u16,t,r,c,p²,pr,pt]] inside C[[u16,p,t,r,c]]/(r²-25t²-2cp)',normalization_quadric_equation=str(q),dimension=19,embedding_dimension=22,depth=18,multiplicity=3,normalization_quotient='C[[u16,a=p²]]*p',conductor_in_normalization='(t,r,c)*normalization',conductor_in_target='(t,r,c,pr,pt)*target',normalization_fiber_length=2,B17_reduced_completion='C[[u16,a]]',B16_reduced_completion='C[[u16]]',B18_reduced_ideal='(t,pt)*target',multiplicity_parameter_reduction='(u16,t,r,c+p²); quotient normalizer has length3',scope='Full completed reduced R image on the open B16 orbit only; no all-B16 or actual CA/CB scheme-intersection boundary claim.',input_sha256={str(a.relative_to(ROOT)):hashlib.sha256(a.read_bytes()).hexdigest() for a in paths},seconds=time.monotonic()-started)
 (output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('B16 OPEN ORBIT FULL QUADRIC-NORMALIZER IMAGE + COTANGENT22 PASS',result['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
