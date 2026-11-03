#!/usr/bin/env python3
"""Uniform polynomial Darboux map on the actual B18/W curve, including W=0."""
import gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 h,x,y,z,w,lam,r=s.symbols('h x y z w W r');X=s.Matrix([h,x,y,z,w]);P=s.zeros(5)
 P[0,2]=-x*x;P[2,3]=h*x;P[2,4]=h*h+z*z;P[3,4]=lam*x*x;P=P-P.T
 source=ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63/staircase.json.gz'
 data=json.loads(gzip.decompress(source.read_bytes()));assert P==s.Matrix(data['CA']['physical_entire_W_curve'])==s.Matrix(data['CB']['physical_entire_W_curve'])
 assert s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)])==s.zeros(5,1)
 jac=[]
 for i,j,k in itertools.combinations(range(5),3):
  q=s.expand(sum(P[i,a]*s.diff(P[j,k],X[a])+P[j,a]*s.diff(P[k,i],X[a])+P[k,a]*s.diff(P[i,j],X[a])for a in range(5)));assert q==0;jac.append(str(q))
 affine=P.extract([0,2,3,4],[0,2,3,4]).subs(x,1);av=s.Matrix([h,y,z,w])
 Z=z-h*h/2;T=w-h*z*z+s.Rational(2,3)*h**3*z-h**3/3-s.Rational(2,15)*h**5
 forward=s.Matrix([h,y,Z,T]);J=forward.jacobian(av)
 const=s.zeros(4);const[0,1]=-1;const[2,3]=lam;const=const-const.T
 assert(J*affine*J.T-const).applyfunc(s.expand)==s.zeros(4)and J.det()==1
 H,Y,ZZ,TT=s.symbols('H Y Z T');bv=s.Matrix([H,Y,ZZ,TT])
 inverse=s.Matrix([H,Y,ZZ+H*H/2,TT+H*ZZ**2+(1+ZZ)*H**3/3+H**5/20])
 assert(inverse.subs(dict(zip(bv,forward)),simultaneous=True)-av).applyfunc(s.expand)==s.zeros(4,1)
 assert(forward.subs(dict(zip(av,inverse)),simultaneous=True)-bv).applyfunc(s.expand)==s.zeros(4,1)
 F=s.Matrix([h*x**4,x**5,y*x**4,z*x**4-h*h*x**3/2,w*x**4-h*z*z*x*x+s.Rational(2,3)*h**3*z*x-h**3*x*x/3-s.Rational(2,15)*h**5])
 G=s.Matrix([h*x**4,x**5,y*x**4,z*x**4+h*h*x**3/2,w*x**4+h*z*z*x*x+(x+z)*h**3*x/3+h**5/20])
 assert(G.subs(dict(zip(X,F)),simultaneous=True)-x**24*X).applyfunc(s.expand)==s.zeros(5,1)
 assert(F.subs(dict(zip(X,G)),simultaneous=True)-x**24*X).applyfunc(s.expand)==s.zeros(5,1)
 assert s.Poly(F[4],x,h,y,z,w).coeff_monomial(h**5)!=0 and s.Poly(G[4],x,h,y,z,w).coeff_monomial(h**5)!=0
 diag=s.diag(1,1,r,1,1);changed=(diag*P.subs({y:y/r},simultaneous=True)*diag.T/r).applyfunc(s.expand)
 assert(changed-P.subs(lam,lam/r)).applyfunc(s.expand)==s.zeros(5)
 # At W=0, each leaf is the y-cone over a basepoint-free rational quintic.
 S,A,c,d=s.symbols('S A c d');curve=s.Matrix([A*S**4,S**5,c*S**5+A*A*S**3/2,d*S**5+c*c*A*S**4+(1+c)*A**3*S*S/3+A**5/20])
 assert all(s.Poly(f,S,A).total_degree()==5 for f in curve)
 assert curve.subs(S,0)==s.Matrix([0,0,0,A**5/20])and curve[1]==S**5
 assert s.cancel(curve[0]/curve[1])==A/S
 payload=dict(physical_P=rows(P),all10_strict_Jacobi_polynomials=jac,affine_coordinate_order=list(map(str,av)),affine_P=rows(affine),forward=list(map(str,forward)),inverse=list(map(str,inverse)),constant_target=rows(const),homogeneous_degree5_forward=list(map(str,F)),homogeneous_degree5_inverse=list(map(str,G)),both_compositions_common_factor='x^24',all_nonzero_W_orbit_matrix='diag(1,1,r,1,1)',orbit_tensor_scalar='1/r',orbit_parameter_change='W -> W/r',W0_Casimirs=[str(Z),str(T)],W0_leaf_projective_curve=list(map(str,curve)),leaf_curve_parameter_recovery='h/x=A/S',leaf_curve_degree_exact=5,leaf_surface='projective cone over the rational quintic, vertex [0:0:1:0:0]')
 (out/'certificate.json').write_text(json.dumps(payload,indent=2)+'\n')
 note=ROOT/'research/notes/AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63.md';paths=[Path(__file__),note,source]
 report=dict(id='AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63',all_passed=True,actual_entire_W_family_bound_to_CA_CB_source=True,full_canonical_divergence_and_all10_strict_Jacobi_zero=True,uniform_polynomial_Darboux_and_inverse_all_W_including_zero=True,parameter_W_not_inverted_in_Darboux_map=True,exact_Cremona_bidegree=[5,5],all_nonzero_W_single_projective_PGL_orbit=True,W0_two_global_affine_polynomial_Casimirs=True,W0_each_affine_leaf_is_A2=True,W0_each_leaf_closure_rational_cone_degree=5,nonzero_W_affine_open_single_symplectic_leaf_A4=True,all_W_local_completion_types_deduced_without_generic_W_input=False,global_all_B18_special_strata_or_component_classification_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main(ROOT/'reports/AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63')
