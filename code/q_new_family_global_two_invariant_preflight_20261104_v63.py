#!/usr/bin/env python3
"""Intrinsic H-foliation/chord invariant; no point or finite-equivalence search."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_GLOBAL_TWO_INVARIANT_PREFLIGHT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 refs={
  'geometry':'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json',
  'target_review':'reports/Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63_REVIEW/report.json',
  'source_review':'reports/Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW/report.json',
  'compression_review':'reports/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63_REVIEW/report.json',
  'inverse':'reports/Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63/report.json',
  'inverse_review':'reports/Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63_REVIEW/report.json',
 }
 data={k:json.loads((ROOT/p).read_text())for k,p in refs.items()}
 assert all(d['all_passed']for d in data.values())
 u,v,x,y,z,w,h=s.symbols('u v x y z w h');n0,n1,n2=s.symbols('n0 n1 n2');coords=(x,y,z,w)
 P=s.Matrix(data['geometry']['affine_P']).subs({n0:u,n1:1,n2:v})
 def degree2(f):return sum(c*s.prod(a**e for a,e in zip(coords,m))for m,c in s.Poly(f,*coords).terms()if sum(m)==2)
 lead=P.applyfunc(degree2);F=x+v*y
 V=s.Matrix([64*w*F,-u*v*y*y/2-u*x*y-w*w+x*x/64-64*w*y,0,x*F])
 ez=s.Matrix([0,0,1,0]);assert(lead-(ez*V.T-V*ez.T)).applyfunc(s.expand)==s.zeros(4)
 Q=s.expand(-x*x+64*w*w+64*u*x*y+32*u*v*y*y)
 original_Q=s.sympify(data['geometry']['homogeneous_quadric']).subs({n0:u,n1:1,n2:v,h:0})
 assert s.expand(Q-72*original_Q)==0
 def deriv(f):return s.expand(sum(V[i]*s.diff(f,a)for i,a in enumerate(coords)))
 lp,lm=x+8*w,x-8*w
 assert s.expand(deriv(lp)-8*F*lp)==0
 assert s.expand(deriv(lm)+8*F*lm)==0
 assert s.expand(deriv(Q)+u*F*Q)==0
 assert s.expand(deriv(lp*lm))==0
 assert s.expand(-u*F+s.Rational(1,16)*u*(8*F+8*F))==0
 assert s.simplify(2+(u/16-1)+(-u/16-1))==0
 qcoords=(x,y,w);QM=s.hessian(Q,qcoords)/2
 assert s.factor(QM.det())==-2048*u*(32*u+v)
 assert Q.subs({x:0,y:1,w:0})==32*u*v
 li=s.Matrix([1,0,8]);lj=s.Matrix([1,0,-8]);Qi=QM.inv()
 aa=s.factor((li.T*Qi*li)[0]);bb=s.factor((lj.T*Qi*lj)[0]);cc=s.factor((li.T*Qi*lj)[0])
 J=s.factor(256*cc*cc/(aa*bb));assert s.cancel(J-(v/u+16)**2)==0
 A,b,j=s.symbols('A b J');bprime=-32*A-b
 assert s.expand((bprime+16*A)**2-(b+16*A)**2)==0
 poly=b*b+32*A*b+(256-j)*A*A
 assert s.discriminant(poly,b)==4*A*A*j
 assert s.cancel(poly.subs({A:u*u,b:u*v,j:J},simultaneous=True))==0
 cert=dict(
  family_parameters=['u','v'],intrinsic_H='unique multiplicity3 hyperplane in actual h^3 Q',
  intrinsic_vertex='unique vertex of Q restricted to H; [x:y:z:w]=[0:0:1:0]',
  leading_full_tensor=[[str(c)for c in row]for row in lead.tolist()],
  induced_P2_coordinates=['x','y','w'],vector_field=[str(V[i])for i in(0,1,3)],
  conic_polynomial=str(Q),conic_matrix=[[str(c)for c in row]for row in QM.tolist()],
  conic_determinant=str(s.factor(QM.det())),conic_at_line_intersection='32*u*v',
  invariant_lines=[str(lp),str(lm)],cofactors=['8*(x+v*y)','-8*(x+v*y)','-u*(x+v*y)'],
  closed_log_form='dlog(Q/(lplus*lminus))+(u/16)*dlog(lplus/lminus)',
  residues=dict(Q='1',lplus='u/16-1',lminus='-u/16-1'),
  no_other_geometric_generic_invariant_line_proof='restriction to any other line gives integer orders a+(u/16)b=0; generic u implies a=b=0; lplus/lminus constant forces line through O; Q(O)!=0 contradicts Q/(lplus*lminus) constant',
  dual_conic_values=dict(lplus_self=str(aa),lminus_self=str(bb),cross=str(cc)),
  intrinsic_chord_invariant='256*(lplus^T Q^-1 lminus)^2/((lplus^T Q^-1 lplus)*(lminus^T Q^-1 lminus))',
  chord_restriction=str(J),
  first_global_invariant='A=1600*r_firstjet=u^2, independent of evaluation point and Poisson scalar on dense actual fullsource',
  target_field_sandwich='k(A,J) subset k(X_NEW25)^PGL5 subset k(A,b), A=u^2,b=u*v',
  degree2_equation=str(poly),discriminant='4*A^2*J',
  endpoints_both_pure_rational_dimension2=True,
  undetermined_generic_extra_equivalence='sigma:(u,v)->(u,-32*u-v), besides known simultaneous sign',
  no_sigma_existence_or_nonexistence_proved=True,
  no_special_parameter_exhaustion=True)
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),*[ROOT/p for p in refs.values()]]
 report=dict(id=ID,all_passed=True,actual_full_H_leading_tensor_logarithmic_and_chord_identities_verified=True,
  global_u_squared_invariance_proof_and_geometric_generic_line_uniqueness_recorded=True,
  actual_generic_fullPGL_moduli_field_rational2_by_degree2_sandwich=True,
  WHICH_endpoint_and_sigma_classification_proved=False,
  no_new_point_prime_high_invariant_GB_or_full_equivalence_solve=True,
  independent_conceptual_review_pending=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
  evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','actual_generic_fullPGL_moduli_field_rational2_by_degree2_sandwich','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
