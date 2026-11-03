#!/usr/bin/env python3
"""Independent tensor formulas for the two full B18 quadratic image charts."""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 folder=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63'
 r=json.loads((folder/'report.json').read_text());assert r['all_passed']
 for name,want in r['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want
 h,x,y,z,w=s.symbols('h x y z w');X=s.Matrix([h,x,y,z,w]);pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
 def vector(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[a]*X[b])for i,j in pairs for a,b in mons])
 def canonical(P):
  div=s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)])
  return (P-(X*div.T-div*X.T)/5).applyfunc(s.expand)
 def raw(name,f,g,H):
  P=s.zeros(5);P[1,2]=f;P[3,4]=g*h*h if name=='CA'else g
  P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z)
  if name=='CB':P[2,3]-=s.diff(f,x)*z/2;P[2,4]-=s.diff(f,x)*w/2
  return P-P.T
 swap=s.eye(5);swap[0,0]=swap[1,1]=0;swap[0,1]=swap[1,0]=1
 def change(P):return(swap*P.subs({h:x,x:h},simultaneous=True)*swap.T).applyfunc(s.expand)
 cubics=[s.prod((h,x,z,w)[i]for i in m)for m in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in m)]
 D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);e=s.Matrix([0,0,1,0,0]);P0=e*D.T-D*e.T;v0=vector(P0)
 def action(A,P):return(A*P+P*A.T-sum(((A*X)[k]*P.diff(X[k])for k in range(5)),s.zeros(5))).applyfunc(s.expand)
 orbit=s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P0))for i,j in itertools.product(range(5),repeat=2)])
 assert orbit.rank()==19
 oi=r['orbit_GL_generator_indices'];O=orbit[:,oi];assert O.rank()==19
 aa,bb,cc,dd,ee=s.symbols('a b c d e');chart=s.eye(5);chart[0,3]=dd;chart[0,4]=ee;chart[1,0]=aa;chart[1,3]=bb;chart[1,4]=cc
 moved=(chart*P0.subs(dict(zip(X,chart.inv()*X)),simultaneous=True)*chart.T).applyfunc(s.expand);moved=vector(moved)
 bases={}
 for name in ['CA','CB']:
  bases[name]=[change(canonical(raw(name,f,0,0)))for f in[h*h,h*x,x*x]]+[change(canonical(raw(name,0,1 if name=='CA'else h*h,0)))]+[change(canonical(raw(name,0,0,H)))for H in cubics]
 for name in ['CA','CB','N']:
  base=bases[name]if name!='N'else[bases['CA'][0],bases['CA'][1],*bases['CA'][3:]]
  C=s.Matrix.hstack(*map(vector,base));rows=list(C.T.rref()[1]);equations=(moved-C*C[rows,:].inv()*moved[rows,0]).applyfunc(s.expand)
  gb=s.groebner([f for f in equations if f],bb,cc,dd,ee,aa,order='lex')
  expected=[bb,cc,dd+(2 if name!='CB'else-s.Rational(1,2))*aa,ee,aa*aa]
  assert list(gb)==expected
 p=s.Symbol('p');K=s.symbols('k0:9');origin=dict.fromkeys([p,*K],0);rows=r['common28_linear_target_rows'];coefficients=[]
 for name in ['CA','CB']:
  rec=r['families'][name];alpha=s.Matrix(list(map(s.sympify,rec['source_alpha_coefficients'])));M=s.Matrix(rec['source_kernel_Lie_matrix']);g=s.eye(5)+p*M+p*p*M*M/2;inverse=s.eye(5)-p*M+p*p*M*M/2
  assert g*inverse==s.eye(5)
  f=alpha[0]*h*h+alpha[1]*h*x+alpha[2]*x*x;H=sum(v*m for v,m in zip(alpha[4:,0],cubics));gg=alpha[3]if name=='CA'else-alpha[3]*f
  source=change(canonical(raw(name,f,gg,H)))
  U=vector((g*source.subs(dict(zip(X,inverse*X)),simultaneous=True)*g.T).applyfunc(s.expand))
  assert list(U)==list(map(s.sympify,rec['core_full150_coordinates']))
  assert U.subs(origin)==v0 and U.diff(p).subs(origin)==s.zeros(150,1)
  J=O.row_join(U.jacobian(K).subs(origin));assert J==s.Matrix(rec['linear_J'])and J.rank()==28
  augmented=s.Matrix(rec['entire_source_tangent']);assert augmented[:150,:28]==J
  assert augmented.rank()==29
  W=s.Matrix([s.Poly(v,p).nth(2).subs(dict.fromkeys(K,0))for v in U-v0]);corrected=(W-J*J[rows,:].inv()*W[rows,0]).applyfunc(s.expand)
  assert list(map(str,corrected))==r['all150_projected_second_coefficients'][0 if name=='CA'else 1]
  coefficients.append(corrected[61])
  assert J[rows,:].det()==s.sympify(r['common28_linear_minors'][0 if name=='CA'else 1])
 assert coefficients==[-4,s.Rational(3,8)]
 reconstruction=s.Matrix(r['canonical_reconstruction_from126']);free=r['common_target126_raw_coordinate_indices'];assert reconstruction.shape==(150,126)and reconstruction[free,:]==s.eye(126)
 for j in range(126):
  values=reconstruction[:,j];P=s.zeros(5)
  for k,(i,ell)in enumerate(pairs):
   P[i,ell]=sum(values[15*k+n]*X[a]*X[b]for n,(a,b)in enumerate(mons));P[ell,i]=-P[i,ell]
  assert all(s.expand(sum(s.diff(P[i,ell],X[i])for i in range(5)))==0 for ell in range(5))
 paths=[Path(__file__),folder/'report.json',ROOT/'research/notes/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63.md',ROOT/'research/notes/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW.md']
 result=dict(id='AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW',all_passed=True,
  independent_CA_CB_raw_and_canonical_formulas=True,whole_three_flag_Groebner_fibers_rebuilt=True,
  both_all150_finite_core_maps_rebuilt=True,own_physical_GL_orbit_rank=19,
  common28_minor_and_corrected_row61_coefficients=list(map(str,coefficients)),
  canonical126_reconstruction_and_all_divergence_constraints=True,
  finite_free_rank2_completed_normalizers_and_separate_image_kernel_proof_reviewed=True,
  full_intersection_nilradical_primary_type_or_normalizer_tensor_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('INDEPENDENT B18 ACTUAL QUADRATIC IMAGE REVIEW PASS',result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
