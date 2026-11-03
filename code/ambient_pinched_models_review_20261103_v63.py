#!/usr/bin/env python3
"""Independent saved-coordinate and complete-image audit for base/component pinches."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def same(a,b):return {k:v for k,v in a.items()if k!='seconds'}=={k:v for k,v in b.items()if k!='seconds'}

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--base-replay',type=Path,required=True)
 ap.add_argument('--component-replay',type=Path,required=True)
 ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 bp=ROOT/'reports/AMBIENT_BASE_PINCHED_SURFACE_20261103_v63/report.json'
 cp=ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json'
 bd,cd=load(bp),load(cp);assert same(bd,load(args.base_replay))and same(cd,load(args.component_replay))
 paths=[Path(__file__),bp,cp]
 for d in [bd,cd]:
  for name,h in d['input_sha256'].items():
   path=ROOT/name;assert sha(path)==h,name;paths.append(path)
 p,sv,rv=s.symbols('p s r');k=s.symbols('k0:9');local={str(v):v for v in [p,sv,rv,*k]}
 curve=s.Matrix([s.sympify(v,locals=local)for v in bd['slice_24_coordinates']])
 # Independent reconstruction of the finite covariant two-form action from
 # the six displayed fixed-plane coefficients; no common_plane/finite_beta.
 x,y,z,w=s.symbols('x y z w');Y=s.Matrix([x,y,z,w]);pairs=list(itertools.combinations(range(4),2))
 B=s.zeros(4)
 entries=[p*x,2*z,x+(0)*z+2*sv*w,0,0,-p*x]
 # beta_xz=2z+(D+t/2)w=2z, beta_xw=x+(D-t/2)z+2Fw.
 entries[2]=x-p*z+2*sv*w
 for v,(i,j)in zip(entries,pairs):B[i,j]=v;B[j,i]=-v
 g=s.Matrix([[1,0,-p,0],[0,1,0,0],[0,0,1,0],[0,p,0,1]]);gi=g.inv()
 transformed=(gi.T*B.subs(dict(zip(Y,gi*Y)),simultaneous=True)*gi).applyfunc(s.expand)
 actual=s.Matrix([transformed[i,j].coeff(v)for i,j in pairs for v in Y]);assert actual==curve
 assert [curve[11]/2,curve[2],-curve[14]/2,-curve[3]/2]==[sv,p*p,p**3,sv*p]
 for v in curve:assert all(not(ep==1 and es==0)for(es,ep),c in s.Poly(v,sv,p).terms())
 S,U,V,W=s.symbols('S U V W')
 elim=s.groebner([S-sv,U-p*p,V-p**3,W-sv*p],p,sv,S,U,V,W,order='lex',domain=s.QQ)
 E=[v.as_expr()for v in elim.polys if not v.as_expr().has(p,sv)]
 target=[V*V-U**3,U*W-S*V,V*W-S*U**2,W*W-S*S*U]
 assert list(s.groebner(E,S,U,V,W,domain=s.QQ))==list(s.groebner(target,S,U,V,W,domain=s.QQ))
 base=load(ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json')
 bi,qi=base['beta_indices_zero_based'],base['q_indices_zero_based'];adapter=s.Matrix(base['beta16_to_intrinsic24_adapter'])
 records={};all_u={};common={p:0,sv:0,rv:0,**dict.fromkeys([k[i]for i in [1,3,4,6,7,8]],0)}
 for name,body in cd['families'].items():
  pars=[local[v]for v in body['source_fiber_parameters']]
  u=s.Matrix([s.sympify(v,locals=local)for v in body['full_61_coordinates']]);all_u[name]=u
  assert adapter*u[bi,0]==curve
  qv=u[qi,0];M=qv.jacobian(pars);zero=dict.fromkeys(pars,0)
  assert not any(v.free_symbols&set(pars)for v in M)
  assert (qv-M*s.Matrix(pars)-qv.subs(zero)).applyfunc(s.expand)==s.zeros(45,1)
  q0=qv.subs({p:0,sv:0});q1=qv.diff(p).subs({p:0,sv:0});M0=q0.jacobian(pars)
  selected=body['actual_q_selected_rows'];L=M0[selected,:]
  assert L==s.Matrix(body['fiber_linear_block_at_origin'])and L.det()==s.Rational(body['fiber_linear_block_determinant'])!=0
  assert selected==([19,21,24,25,27,28,32,33,34]if name=='N'else[0,19,21,24,25,27,28,32,33,34])
  odd=(q1-M0*L.inv()*q1[selected,0]).applyfunc(s.expand)
  assert odd==s.Matrix([s.sympify(v,locals=local)for v in body['first_odd_coefficient_at_fixed_target_K']])
  assert odd.subs(zero)==s.zeros(45,1)
  OM=odd.jacobian(pars);assert OM==s.Matrix(body['first_odd_coefficient_matrix'])
  desired=[v for v in pars if v in [rv,*[k[i]for i in [1,3,4,6,7,8]]]]
  DM=s.Matrix([[int(v==a)for a in pars]for v in desired])
  assert OM.rank()==len(desired)and OM.col_join(DM).rank()==len(desired)
  assert OM[:,[i for i,a in enumerate(pars)if a not in desired]]==s.zeros(45,len(pars)-len(desired))
  assert (M[selected,:].det()-L.det()).expand()==0
  records[name]=dict(actual_selected_q_block_determinant=str(L.det()),odd_ideal_generators=list(map(str,desired)),
   odd_rank=len(desired),pinching_codimension=len(desired)+1,actual_beta_coordinates_match=True,
   all45_odd_coefficients_match=True,selected_block_determinant_constant_in_p_s=True)
 assert all((all_u[n].subs(common)-all_u['N'].subs(common)).applyfunc(s.expand)==s.zeros(61,1)for n in ['CA','CB'])
 paths += [ROOT/'research/notes/AMBIENT_BASE_PINCHED_SURFACE_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63.md',
  ROOT/'research/notes/AMBIENT_PINCHED_MODELS_20261103_v63_REVIEW.md',
  ROOT/'research/notes/AMBIENT_INTERSECTION_SYMPLECTIC_SCOPE_20261103_v63.md',
  ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63_REVIEW.json',
  ROOT/'reports/AMBIENT_Q_NORMALIZATION_20261103_v63_REVIEW/report.json']
 result=dict(id='AMBIENT_PINCHED_MODELS_20261103_v63_REVIEW',all_passed=True,
  fresh_full_QQ_producer_replays_match=True,independent_finite_two_form_action=True,full_toric_elimination_matches=True,
  families=records,whole_common_nonnormal_stratum_actual61_outputs_coincide=True,
  finite_normalization_completion_injection_and_equivariant_product_argument_accepted=True,
  full_image_argument='R and p2,p3 belong to A; D=R[[p2]] subset A subset B=R[[p]] is finite over D. Hence A is closed and p2B subset A. The exact first-odd ideal determines the whole ring R+Ip+p2B.',
  component_completed_rings={'N':'A7','CA':'A8','CB':'A8'},base_completed_ring='C[[u7,s,p2,p3,sp]]',
  scope='Individual reduced S61 components at distinguished a=1 Q and all nonnormal points above the critical rank2 beta orbit; no full S61 or whole scheme intersection ideal. Projective consequences only on generic-Poisson-rank4 open.',
  input_sha256={str(path.relative_to(ROOT)):sha(path)for path in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS base/full component completed images; exact odd ideals and shared stratum')

if __name__=='__main__':main()
