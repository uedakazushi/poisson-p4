#!/usr/bin/env python3
"""All-magnetic node/two-node-tangents source as an actual MPS degeneration.

The geometric source Jac(f)+lambda W*dW wedge Euler_ABC contains a rank-four
SNC toric point. A collision of its two independent hyperplanes with the two
node tangents has polynomial coefficients after tensor scaling. Integrating
the magnetic quadratic in B preserves every one of its six parameters.
"""
import hashlib,itertools,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(4),2))

def jac(P,variables):
 return [s.expand(sum(P[i,k]*s.diff(P[j,l],variables[k])+P[j,k]*s.diff(P[l,i],variables[k])+P[l,k]*s.diff(P[i,j],variables[k]) for k in range(4))) for i,j,l in itertools.combinations(range(4),3)]

def matrix(entries):
 P=s.zeros(4)
 for (i,j),a in zip(PAIRS,entries):P[i,j]=s.expand(a);P[j,i]=-P[i,j]
 return P

def encode(P):return {str(i+1)+','+str(j+1):str(s.factor(P[i,j])) for i,j in PAIRS if P[i,j]!=0}

def main():
 start=time.monotonic();h,x,y,z,w,e=s.symbols('h x y z w epsilon');A,B,C,W=s.symbols('A B C W')
 r=s.symbols('r0:6');lam=s.symbols('lambda');g=s.symbols('g0:4')
 G=3*r[0]*C*C*B+3*r[1]*C*A*B+s.Rational(3,2)*r[2]*C*B*B+3*r[3]*A*A*B+s.Rational(3,2)*r[4]*A*B*B+r[5]*B**3
 f0=-lam*A*B*C+sum(g[i]*A**(3-i)*C**i for i in range(4));f=f0+e*G
 original=matrix([s.diff(f,C),-s.diff(f,B),-lam*W*A,s.diff(f,A),-lam*W*B,-lam*W*C])
 assert jac(original,(A,B,C,W))==[0]*4
 pf=s.expand(original[0,1]*original[2,3]-original[0,2]*original[1,3]+original[0,3]*original[1,2])
 assert s.expand(pf+3*lam*W*f)==0
 # H is Casimir. On the affine H=h slice the inverse coordinate formulas are
 # x=A, y=B, z=(C-H)/e, w=(W-A)/e. Multiply the tensor by e.
 sub={A:x,B:y,C:h+e*z,W:x+e*w}
 family=matrix([e*s.diff(f,C),-s.diff(f,B),-lam*W*A,s.diff(f,A),-lam*W*B+s.diff(f,C),-lam*h*w-lam*e*z*w-s.diff(G,B)])
 family=family.applyfunc(lambda a:s.expand(a.subs(sub,simultaneous=True)))
 assert jac(family,(x,y,z,w))==[0]*4
 # Check the raw rational pushforward against the displayed polynomial tensor.
 T=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,1/e,0],[-1/e,0,0,1/e]])
 rational=e*T*original*T.T
 assert all(s.cancel((rational[i,j].subs(sub,simultaneous=True))-family[i,j])==0 for i,j in PAIRS)
 pf_family=s.factor(family[0,1]*family[2,3]-family[0,2]*family[1,3]+family[0,3]*family[1,2])
 expected_pf=-3*lam*(x+e*w)*f.subs(sub,simultaneous=True)
 assert s.expand(pf_family-expected_pf)==0
 assert all(s.Poly(a,h,x,y,z,w).total_degree()==2 for a in family if a!=0)
 # The standard model at f=ABC, lambda=-1 is a genuine SNC point: all five
 # projective factors H,W,A,B,C are independent and its Pfaffian is 3ABCW.
 toric=matrix([A*B,-B*C,W*A,B*C,W*B,W*C])
 # Build this independently from the universal formula to avoid a sign typo.
 toric=original.subs({e:0,lam:-1,**dict(zip(g,(0,0,0,0)))}).subs({})
 toric=matrix([A*B,-A*C,W*A,B*C,W*B,W*C])
 assert jac(toric,(A,B,C,W))==[0]*4
 toric_pf=s.factor(toric[0,1]*toric[2,3]-toric[0,2]*toric[1,3]+toric[0,3]*toric[1,2])
 assert toric_pf==3*A*B*C*W
 # Match every coefficient of the existing normalized node source, with all
 # six magnetic parameters retained. f0=ABC-(A-C)^3/8 is its exact cubic.
 node=family.subs({lam:-1,g[0]:-s.Rational(1,8),g[1]:s.Rational(3,8),g[2]:-s.Rational(3,8),g[3]:s.Rational(1,8)})
 Z,V=s.symbols('Z V');fiber=s.Matrix([[-s.Rational(8,3),0],[s.Rational(1,8),s.Rational(1,8)]])
 change=s.eye(4);change[2:4,2:4]=fiber
 transformed=(24*change*node*change.T).subs({z:-3*Z/8,w:8*V+3*Z/8},simultaneous=True).applyfunc(s.expand)
 source_path=ROOT/'reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json';packet=json.loads(source_path.read_text());case=next(c for c in packet['cases']if c['case']=='node_two_node_tangents')
 env={'h':h,'x':x,'y':y,'z':Z,'w':V,**{str(a):a for a in r}}
 target=s.zeros(4)
 for key,value in case['full_six_magnetic_source'].items():
  i,j=map(int,key.strip('()').split(','));target[i-1,j-1]=s.sympify(value,locals=env);target[j-1,i-1]=-target[i-1,j-1]
 assert transformed.subs(e,0)==target
 node_pf=s.factor(transformed[0,1]*transformed[2,3]-transformed[0,2]*transformed[1,3]+transformed[0,3]*transformed[1,2])
 q=x*h*y-(x-h)**3/8
 assert s.expand(node_pf.subs(e,0)+576*x*q)==0
 here=Path(__file__);out=dict(task_id='E_reduced_nonnormal',status='WHOLE_NODE_TWO_NODE_TANGENTS_SOURCE_IN_MPS_CLOSED_UNION',all_passed=True,
  coefficient_ring='QQ[lambda,g0,g1,g2,g3,r0,...,r5,epsilon,h,x,y,z,w]',
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [here,source_path]},producer_bytes=here.read_text(),
  model_f0=str(f0),magnetic_primitive_G=str(G),model_bivector=encode(original),
  model_whole_strict_Jacobi_zero=True,model_affine_pfaffian=str(-3*lam*W*f),
  collision_old_coordinates=['H=h','W=x+epsilon*w','A=x','B=y','C=h+epsilon*z'],collision_determinant='-epsilon**2 up to the stated order',tensor_scale='epsilon',
  full_polynomial_collision_bivector=encode(family),rational_pushforward_equals_polynomial_family=True,
  full_family_strict_Jacobi_zero=True,all_coefficients_homogeneous_quadratic=True,
  actual_affine_pfaffian=str(expected_pf),actual_projective_divisor='H*(x+epsilon*w)*f(x,y,h+epsilon*z)',
  epsilon_zero_fiber_B_coefficient='-lambda*h*w-G_B(x,y,h)',
  toric_model_brackets=encode(toric),toric_model_pfaffian=str(toric_pf),toric_projective_divisor='3*H*W*A*B*C',toric_five_coordinate_factors_independent=True,
  node_f0=str(A*B*C-(A-C)**3/8),node_lambda='-1',fiber_coordinate_change=['Z=-8*z/3','V=(z+w)/8'],fiber_determinant=str(fiber.det()),final_tensor_scale='24',
  exact_whole_existing_six_parameter_source_identity=True,full_node_family_brackets=encode(transformed),node_epsilon_zero_pfaffian=str(-576*x*q),
  inverted_elements=['epsilon only in explaining the PGL map for nonzero epsilon; the final family is polynomial everywhere'],
  complements='epsilon=0 is the exact entire six-magnetic target. No r coefficient is inverted and no magnetic boundary is removed. lambda=-1 is fixed for the node source. All generic epsilon members belong to the accepted model closed union regardless of whether their cubic is smooth.',
  geometric_dependency='The entire irreducible model source lambda!=0,f any ternary homogeneous cubic contains the SNC toric point f=ABC. Its image closure lies in a component meeting SNC, hence in the MPS40 list. This finite closed union remains invariant under PGL and nonzero tensor scaling. The polynomial collision and closedness include epsilon=0.',
  limitation='Exact producer identities are checked here. The SNC-to-known-components theorem and the irreducible-image/closed-union argument require independent mathematical review. No new component count or coverage of the other three plane sources is claimed.',seconds=time.monotonic()-start)
 dest=ROOT/'reports/CLASSIFICATION_PLANE_CUBIC_TORIC_COLLISION_20261107.json'
 if dest.exists():raise FileExistsError(str(dest))
 dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status=out['status'],whole_source_identity=True,seconds=out['seconds'])))
if __name__=='__main__':main()
