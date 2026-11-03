#!/usr/bin/env python3
"""A fixed cubic polynomial Darboux map for every physical Qa, including a=0."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 a=s.Symbol('a');x,y,z,w=s.symbols('x y z w');u=s.Matrix([x,y,z,w])
 P=s.Matrix([[0,0,0,a],[0,0,-s.Rational(2,9),4*x/3],[0,s.Rational(2,9),0,2*y/3+x*x],[-a,-4*x/3,-2*y/3-x*x,0]])
 shift=6*x*z-s.Rational(3,2)*y*y-s.Rational(9,2)*x*x*y
 forward=s.Matrix([x,y,z,w+shift]);inverse=s.Matrix([x,y,z,w-shift])
 assert (forward.subs(dict(zip(u,inverse)),simultaneous=True)-u).applyfunc(s.expand)==s.zeros(4,1)
 assert (inverse.subs(dict(zip(u,forward)),simultaneous=True)-u).applyfunc(s.expand)==s.zeros(4,1)
 jac=forward.jacobian(u);assert jac.det()==1
 constant=s.Matrix([[0,0,0,a],[0,0,-s.Rational(2,9),0],[0,s.Rational(2,9),0,0],[-a,0,0,0]])
 assert (jac*P*jac.T-constant).applyfunc(s.expand)==s.zeros(4)
 assert (inverse.jacobian(u)*constant*inverse.jacobian(u).T-P.subs(dict(zip(u,inverse)),simultaneous=True)).applyfunc(s.expand)==s.zeros(4)
 h=s.Symbol('h');v=s.Matrix([h,x,y,z,w]);hh=s.symbols('H X Y Z W');V=s.Matrix(hh)
 F=s.Matrix([h**3,h*h*x,h*h*y,h*h*z,h*h*w+6*h*x*z-s.Rational(3,2)*h*y*y-s.Rational(9,2)*x*x*y])
 H,X,Y,Z,W=hh
 G=s.Matrix([H**3,H*H*X,H*H*Y,H*H*Z,H*H*W-6*H*X*Z+s.Rational(3,2)*H*Y*Y+s.Rational(9,2)*X*X*Y])
 assert (F.subs(dict(zip(v,G)),simultaneous=True)-H**8*V).applyfunc(s.expand)==s.zeros(5,1)
 assert (G.subs(dict(zip(V,F)),simultaneous=True)-h**8*v).applyfunc(s.expand)==s.zeros(5,1)
 assert s.gcd_list(list(F))==s.gcd_list(list(G))==1
 Ph=s.zeros(5);Ph[1:,1:]=h*h*P.subs({x:x/h,y:y/h,z:z/h,w:w/h},simultaneous=True)
 Ph=Ph.applyfunc(s.cancel);Ch=s.zeros(5);Ch[1:,1:]=H*H*constant
 assert (F.jacobian(v)*Ph*F.jacobian(v).T-Ch.subs(dict(zip(V,F)),simultaneous=True)).applyfunc(s.expand)==s.zeros(5)
 assert (G.jacobian(V)*Ch*G.jacobian(V).T-Ph.subs(dict(zip(v,G)),simultaneous=True)).applyfunc(s.expand)==s.zeros(5)
 cert=dict(parameter_scope='QQ[a], no localization',affine_coordinates=list(map(str,u)),forward=list(map(str,forward)),inverse=list(map(str,inverse)),Jacobian_determinant='1',transformed_poisson=[[str(q)for q in row]for row in constant.tolist()],homogeneous_forward=list(map(str,F)),homogeneous_inverse=list(map(str,G)),forward_after_inverse_scalar='H^8',inverse_after_forward_scalar='h^8',both_coordinate_gcds='1',a0_polynomial_Casimirs=['x',str(w+shift)])
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes/Q_QA_UNIFORM_DARBOUX_20261103_v63.md',ROOT/'research/notes/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63.md']
 result=dict(id='Q_QA_UNIFORM_DARBOUX_20261103_v63',all_passed=True,
  parameter_scope='All finite complex a including0, uniformly over QQ[a]',
  fixed_parameter_independent_polynomial_automorphism=True,Jacobian_determinant=1,
  affine_Poisson_normal_form='a*d_x wedge d_W - (2/9)*d_y wedge d_z',
  homogeneous_raw5x5_Poisson_pushforward_identities=True,
  fixed_projective_Cremona_bidegree=[3,3],
  a_nonzero_affine_open_is_one_symplectic_leaf=True,
  a0_polynomial_Casimirs_and_all_affine_leaves_A2=True,
  PGL_equivalence_or_component_exhaustion_or_hyperplane_boundary_classification_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('UNIFORM ALL-Qa FIXED CUBIC DARBOUX/(3,3) CREMONA PASS',result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
