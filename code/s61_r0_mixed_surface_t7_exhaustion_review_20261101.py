#!/usr/bin/env python3
"""Parent coefficient/proof review for complete T7 surface-orbit removal."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
def review():
 start=time.time();path=ROOT/'reports/S61_R0_MIXED_SURFACE_T7_EXHAUSTION_20261101.json';r=json.loads(path.read_text());assert r['all_passed'];oldpath=ROOT/'reports/S61_R0_MIXED_SURFACE_T7_FORMAL_20261101.json';old=json.loads(oldpath.read_text())['universal'];m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];b=s.symbols('b:11');H=s.Matrix(r['beta_kernel']).applyfunc(s.Rational);_,F=specialize(m,dict(zip(one,H*s.Matrix(b))),True);poly,*v=ring(b,s.QQ);field=s.QQ.frac_field(*b);fixed=r['chart_fixed_q_index'];chart=r['chart_columns'];J=s.Matrix(old['jacobian']).applyfunc(s.sympify);left=s.Matrix([r['left_kernel_row']]).applyfunc(s.sympify)
 assert (left.to_DM().convert_to(field)*J.to_DM().convert_to(field)).is_zero_matrix
 minor=J[r['minor_rows'],r['minor_columns']].to_DM().convert_to(field).det();assert minor==field.from_sympy(s.sympify(r['jacobian_minor']))
 den=[None,v[8]**2*v[9]**2,v[8]**4*v[9]**5,v[8]**6*v[9]**8];ld=v[8]**3
 def clear(a,d):
  nu,de=s.fraction(s.sympify(a));return (poly.from_expr(nu)*d).exquo(poly.from_expr(de))
 W={}
 for degree,key,col in [(1,'tangent',2),(2,'quadratic_correction',5),(3,'cubic_correction',9)]:
  src=[clear(row[col],den[degree]) for row in old[key]];q=[poly.zero]*46
  for i,j in enumerate(chart):q[j]=src[i]
  q[45]=src[44];assert q[45]==0;W[degree]=q
 # Construct each raw equation's fourth coefficient before the left projection.
 # This reverses the contraction order used by the producer.
 fourth=[]
 for equation in F:
  total=poly.zero
  for mon,c in equation.items():
   ij=list(mon)+[45]*(2-len(mon));i,j=ij
   total+=poly.from_expr(c)*(W[1][i]*W[3][j]+W[2][i]*W[2][j]+W[3][i]*W[1][j])
  fourth.append(total)
 numerator=sum((clear(lv,ld)*q4 for lv,q4 in zip(left,fourth)),poly.zero)
 expected=s.sympify(r['quartic_coefficient']);assert numerator==clear(expected,ld*den[2]**2)
 print('independent raw fourth coefficient PASS',time.time()-start,flush=True)
 z,bu,be=primitive_beta(m);tensor=s.Matrix([be[i,j].diff(z[k]) for (i,j),k in COORDS]);B=tensor.jacobian(bu);br=list(B.T.to_DM().rref()[1]);inv=B[br,:].inv();act=s.Matrix.hstack(*[inv*infinitesimal_action(tensor,A)[br,:] for A in sp_basis()]);act=act.subs(dict(zip(bu,H*s.Matrix(b))),simultaneous=True)
 results=[]
 for c in r['exceptions']:
  var=s.Symbol(c['elimination_variable']);value=s.sympify(c['elimination_value']);other=[x for x in b if x!=var];sub={var:value};bb=(H*s.Matrix(b)).subs(sub,simultaneous=True);D=bb.jacobian(other).row_join(act.subs(sub,simultaneous=True));saved=s.Matrix(c['full_source_derivative']).applyfunc(s.sympify)
  assert (D.to_DM().convert_to(field)-saved.to_DM().convert_to(field)).is_zero_matrix
  K=s.Matrix(c['left_kernel']).applyfunc(s.sympify);assert K.rank()==16-c['beta_image_dimension'];assert (K.to_DM().convert_to(field)*D.to_DM().convert_to(field)).is_zero_matrix
  vals={s.Symbol(k):s.Rational(v) for k,v in c['sample_values'].items()};M=D.subs(vals);det=M[c['sample_minor_rows'],c['sample_minor_columns']].det();assert det==s.Rational(c['sample_minor']) and det!=0
  f=s.sympify(c['factor']);assert s.degree(f,var)==1 and s.gcd(s.diff(f,var),f.subs(var,0))==1
  results.append(dict(factor=c['factor'],independent_beta_image_dimension=c['beta_image_dimension']))
 paths=[Path(__file__),path,oldpath,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'research/notes/S61_R0_MIXED_SURFACE_T7_EXHAUSTION_20261101.md']
 return dict(id='S61_R0_MIXED_SURFACE_T7_EXHAUSTION_REVIEW_20261101',all_passed=True,independent_fourth_coefficient=True,exception_checks=results,proof_audit=['Universal full quadratic/cubic vanishings and t-jets are supplied by the independently reviewed T7 formal certificate.','The nonzero normal source minor gives an actual smooth fixed-beta surface with tangent v2=0 on the whole stated open.','Normal quartic nonzero gives transverse upper4; same-surface generic lower4 gives finite freeR2 rank4 and t=0.','The irreducible factors b8,b9,P all have proper full group beta images, covering every bad factor.','Hence the entire T7 orbit over sufficiently general beta is removed, without including T5 or binary/flag endpoints.'],limitations=['Saved universal 3-jets and the row-103 left projection are verified/used as certificates, not independently re-solved.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=review()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({k:r[k] for k in ['all_passed','seconds']}))
