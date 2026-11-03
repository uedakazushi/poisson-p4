#!/usr/bin/env python3
"""Uniform a!=0 Poisson Cremona, rational leaves, and actual LPB/Log sources."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_CA_SELF_PAIR_CREMONA_20261103_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def sign(seq):return (-1)**sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 parent=ROOT/'reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63';pr=load(parent/'report.json');assert pr['all_passed']
 for n,digest in pr['input_sha256'].items():assert sha(ROOT/n)==digest,n
 for n,digest in pr['evidence_sha256'].items():assert sha(parent/n)==digest,n
 source=load(parent/'certificate.json.gz')
 h,x,y,z,w,a=s.symbols('h x y z w a');X=s.Matrix([h,x,y,z,w]);H,Z,W=s.symbols('H Z W')
 V=s.Matrix([x*x+12*a*h*x/5,-3*a*x*x/5,0,h*x+2*a*x*z/5,h*h+z*z-8*a*w*x/5]);ey=s.Matrix([0,0,1,0,0]);P=ey*V.T-V*ey.T
 pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
 def vec(P,variables=X):return s.Matrix([s.Poly(P[i,j],*variables).coeff_monomial(variables[r]*variables[t]) for i,j in pairs for r,t in mons])
 def canonical(P,variables=X):
  d=s.Matrix([sum(s.diff(P[i,j],variables[i]) for i in range(5)) for j in range(5)])
  return (P-(variables*d.T-d*variables.T)/5).applyfunc(s.expand)
 actual=s.Matrix([sum(s.Rational(c)*a**power for power,c in q) for q in source['first_source_physical150']])
 assert vec(P)==actual and canonical(P)==P
 D=s.Matrix([(V[0]-h*V[1]).subs(x,1),(V[3]-z*V[1]).subs(x,1),(V[4]-w*V[1]).subs(x,1)])
 assert D==s.Matrix([1+3*a*h,h+a*z,h*h+z*z-a*w])
 Hexpr=h+1/(3*a);Zexpr=z-h/(2*a)-1/(2*a*a)
 hv=H-1/(3*a);zv=Z+H/(2*a)+1/(3*a*a)
 target=s.Poly(s.expand(hv*hv+zv*zv),H,Z)
 Q=sum(c*H**i*Z**j/(a*(3*i+j+1)) for (i,j),c in target.terms())
 Q=s.expand(Q);assert s.expand(a*(3*H*s.diff(Q,H)+Z*s.diff(Q,Z)+Q)-hv*hv-zv*zv)==0
 expected=(4*a*a+1)*H*H/(28*a**3)+H*Z/(5*a*a)+Z*Z/(3*a)+(1-2*a*a)*H/(12*a**4)+Z/(3*a**3)+(a*a+1)/(9*a**5)
 assert s.cancel(Q-expected)==0
 Wexpr=w-Q.subs({H:Hexpr,Z:Zexpr},simultaneous=True)
 deriv=lambda f:s.expand(sum(D[i]*s.diff(f,v) for i,v in enumerate([h,z,w])))
 assert s.cancel(deriv(Hexpr)-3*a*Hexpr)==s.cancel(deriv(Zexpr)-a*Zexpr)==s.cancel(deriv(Wexpr)+a*Wexpr)==0
 aff=s.Matrix([Hexpr,y,Zexpr,Wexpr]);oldaff=s.Matrix([h,y,z,w]);ja=aff.jacobian(oldaff);assert ja.det()==1
 vaff=s.Matrix([D[0],0,D[1],D[2]]);eaff=s.Matrix([0,1,0,0]);paff=eaff*vaff.T-vaff*eaff.T
 newv=s.Matrix([3*a*Hexpr,0,a*Zexpr,-a*Wexpr]);assert (ja*paff*ja.T-(eaff*newv.T-newv*eaff.T)).applyfunc(s.cancel)==s.zeros(4)
 A=h+x/(3*a);B=z-h/(2*a)-x/(2*a*a)
 Qhom=lambda U,T,S:s.expand(S*S*Q.subs({H:U/S,Z:T/S},simultaneous=True))
 C=s.expand(x*w-Qhom(A,B,x));phi=s.Matrix([A,x,y,B,C/x]);jac=phi.jacobian(X)
 assert s.cancel(jac.det()-1)==0 and (jac*X-phi).applyfunc(s.cancel)==s.zeros(5,1)
 HH,XX,YY,ZZ,WW=s.symbols('HH XX YY ZZ WW');new=s.Matrix([HH,XX,YY,ZZ,WW]);newey=s.Matrix([0,0,1,0,0])
 diagV=s.Matrix([3*a*XX*HH,0,0,a*XX*ZZ,-a*XX*WW]);diagP=canonical(newey*diagV.T-diagV*newey.T,new)
 assert (jac*P*jac.T-diagP.subs(dict(zip(new,phi)),simultaneous=True)).applyfunc(s.cancel)==s.zeros(5)
 forward=s.Matrix([A*x,x*x,y*x,B*x,C])
 backward=s.Matrix([XX*(HH-XX/(3*a)),XX*XX,XX*YY,XX*(ZZ+HH/(2*a)+XX/(3*a*a)),XX*WW+Qhom(HH,ZZ,XX)])
 assert (backward.subs(dict(zip(new,forward)),simultaneous=True)-x**3*X).applyfunc(s.expand)==s.zeros(5,1)
 assert (forward.subs(dict(zip(X,backward)),simultaneous=True)-XX**3*new).applyfunc(s.expand)==s.zeros(5,1)
 # Casimirs and a full birational coordinate system on the dense Z!=0 open.
 f=Hexpr/Zexpr**3;g=Zexpr*Wexpr
 assert s.cancel(deriv(f))==s.cancel(deriv(g))==0
 c1,c2,t,Y=s.symbols('c1 c2 t Y');inverse={h:c1*t**3-1/(3*a),y:Y,z:t+c1*t**3/(2*a)+1/(3*a*a),w:c2/t+Q.subs({H:c1*t**3,Z:t},simultaneous=True)}
 assert s.cancel(f.subs(inverse,simultaneous=True)-c1)==s.cancel(g.subs(inverse,simultaneous=True)-c2)==0
 assert s.cancel(Zexpr.subs(inverse,simultaneous=True)-t)==0
 # Quotient curve; y is the cone coordinate. Degree follows from basepoint
 # freeness at0 and infinity and inverse t=B/x, not from a sampled implicit GB.
 curve=[s.expand(t*inverse[h]),t,s.expand(t*inverse[z]),s.expand(c2+t*Q.subs({H:c1*t**3,Z:t},simultaneous=True))]
 assert [s.degree(p,t) for p in curve]==[4,1,4,7]
 qHH=s.factor(s.Poly(Q,H,Z).coeff_monomial(H*H));qHZ=s.factor(s.Poly(Q,H,Z).coeff_monomial(H*Z))
 assert qHH==(4*a*a+1)/(28*a**3) and qHZ==1/(5*a*a)
 assert s.cancel(s.Poly(curve[3],t).coeff_monomial(t**7)-c1*c1*qHH)==0
 assert s.Poly(curve[3],t).coeff_monomial(t**6)==0
 assert s.Poly(curve[3],t).coeff_monomial(t**5)==c1*qHZ and s.Poly(curve[3],t).coeff_monomial(1)==c2
 # Explicit existing Log(1,1,1,2) source, in the exact old builder convention.
 factors=[A,B,x,C];r=s.Matrix([1,-3,2,0]);rr=s.Matrix([0,1,-3,1]);degrees=s.Matrix([1,1,1,2]);assert r.dot(degrees)==rr.dot(degrees)==0
 residue=r*rr.T-rr*r.T;pa,pb,pc=s.Rational(1,4),-s.Rational(3,4),s.Rational(7,4)
 assert [residue[i,j] for i,j in itertools.combinations(range(4),2)]==[4*pa,4*pb,-2*(pa+pb),4*pc,2*(pa-pc),2*(pb+pc)]
 omega={};dual={}
 for i,j in pairs:
  om=0
  for k,l in itertools.combinations(range(4),2):
   om+=residue[k,l]*s.prod(factors[v] for v in range(4) if v not in [k,l])*(s.diff(factors[k],X[i])*s.diff(factors[l],X[j])-s.diff(factors[k],X[j])*s.diff(factors[l],X[i]))
  omega[i,j]=s.expand(om);comp=[u for u in range(5) if u not in [i,j]];u,v,q=comp
  dual[i,j]=sign(comp+[i,j])*(X[u]*P[v,q]-X[v]*P[u,q]+X[q]*P[u,v])
  assert s.cancel(omega[i,j]+dual[i,j]/a)==0
 dform={}
 for i,j,k in itertools.combinations(range(5),3):dform[i,j,k]=s.expand(s.diff(omega[j,k],X[i])-s.diff(omega[i,k],X[j])+s.diff(omega[i,j],X[k]))
 logP=s.zeros(5)
 for i,j in pairs:
  comp=tuple(v for v in range(5) if v not in [i,j]);logP[i,j]=sign([i,j,*comp])*dform[comp];logP[j,i]=-logP[i,j]
 assert (logP+5*P/a).applyfunc(s.cancel)==s.zeros(5)
 cofactors=[s.Rational(12,5)*a*x,s.Rational(2,5)*a*x,-s.Rational(3,5)*a*x,-s.Rational(11,5)*a*x]
 for F,K in zip(factors,cofactors):assert s.cancel(sum(V[i]*s.diff(F,X[i]) for i in range(5))-K*F)==0
 # Existing LPB source: permute constant vector e_y to legacy center e4.
 U=s.Matrix(s.symbols('u0:5'));perm=[0,1,4,2,3];sub={X[i]:U[perm[i]] for i in range(5)};lpb=s.zeros(5)
 for i,j in pairs:lpb[perm[i],perm[j]]=P[i,j].subs(sub,simultaneous=True);lpb[perm[j],perm[i]]=-lpb[perm[i],perm[j]]
 vv=s.Matrix([V[i].subs(sub,simultaneous=True) for i in [0,1,3,4]]+[0]);ee=s.Matrix([0,0,0,0,1])
 assert lpb==ee*vv.T-vv*ee.T and canonical(lpb,U)==lpb
 oldpath=ROOT/'legacy/v14/data/six_rank2_components.json';old=load(oldpath);cases=old if isinstance(old,list) else old.get('cases',old.get('components'))
 oldcases={name:next(c for c in cases if c['name']==name) for name in ['Log1112','LPB']}
 for name,dimension in [('Log1112',28),('LPB',39)]:
  c=oldcases[name];assert c['projective_component_dimension']==dimension and c['ambient_affine_tangent_dimension']==c['family_affine_differential_rank']==dimension+1
 # Powers of a only: no hidden exceptional factor is divided in this model.
 for value in [Q,A,B,C,*phi,*backward]:
  denominator=s.denom(s.cancel(value));assert denominator.free_symbols<={a,x}
  assert len(s.Poly(denominator,a,x).terms())==1
 files=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',parent/'report.json',parent/'certificate.json.gz',oldpath,
        ROOT/'legacy/v14/code/rank2_lower_bounds/lpb_family.py',ROOT/'legacy/v14/code/rank2_lower_bounds/verify_rank2_components.py']
 payload=dict(Q=str(Q),affine_coordinates=list(map(str,aff)),
  affine_inverse=['H-1/(3*a)','y','Z+H/(2*a)+1/(3*a^2)',str(W+Q)],
  generic_leaf_parameterization={str(v):str(q) for v,q in inverse.items()},
  homogeneous_linear_A=str(A),homogeneous_linear_B=str(B),homogeneous_quadratic_C=str(C),
  forward_quadratic=list(map(str,forward)),inverse_quadratic=list(map(str,backward)),
  rational_Casimirs=[str(f),str(g)],generic_quotient_curve=list(map(str,curve)),
  degree7_coefficient=str(c1*c1*qHH),degree5_coefficient=str(c1*qHZ),
  Log1112_divisors=list(map(str,factors)),Log1112_residues=[list(map(int,r)),list(map(int,rr))],legacy_log_parameters=['1/4','-3/4','7/4'],
  all10_log_form_coefficients={str(k):str(v) for k,v in omega.items()},all150_log_tensor_coefficients=list(map(str,vec(logP))),
  LPB_center_permutation=perm,LPB_vector_field=list(map(str,vv)),existing_component_certificates=oldcases)
 (out/'certificate.json').write_text(json.dumps(payload,indent=2)+'\n')
 result=dict(id=STEM,all_passed=True,field='QQ(a), all a!=0',all150_bound_to_actual_polynomial_self_pair=True,
  uniform_quadratic_polynomial_affine_automorphism=True,only_parameter_a_inverted=True,
  homogeneous_quadratic_Cremona_both_compositions_verified=True,full_actual_Poisson_pushforward_verified=True,
  diagonal_affine_Poisson='d_y wedge(3aH*d_H+aZ*d_Z-aW*d_W)',rational_center='C(H/Z^3,ZW)',
  generic_leaf_dense_open_is_C_times_Gm=True,generic_projective_leaf_closure_cone_degree=7,
  exceptional_nonzero_parameter_degree_condition='4a^2+1=0',exceptional_generic_leaf_closure_degree=5,
  degree_requires_nonzero_generic_Casimirs_c1_c2=True,
  exact_existing_LPB39_and_Log1112_28_membership_for_all_nonzero_a=True,
  entire_polynomial_line_including_a0_in_both_closed_components=True,
  Log1112_actual_tensor_factor='-5/a',all10_form_and150_tensor_identities_exact=True,
  new_component_count=0,CB_or_other_component_membership_or_full_conductor_claimed=False,
  old_large_component_certificates_reused_not_rerun=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in files},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS',STEM,result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();main(args.output.resolve())
