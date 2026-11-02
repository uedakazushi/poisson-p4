"""Independent proof/source review with small saved-certificate arithmetic."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_HYPERPLANE_BFREE_HORIZONTAL_REVIEW_20261102'
NOTE='research/notes/'+STEM+'.md'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def verify():
 start=time.monotonic();primary='reports/S61_BINARY_HYPERPLANE_BFREE_HORIZONTAL_20261102.json';q=json.loads((ROOT/primary).read_text());assert q['all_passed']and not q['resumed_run'];hashes={}
 for path,h in q['input_sha256'].items():assert sha(path)==h,path;hashes[path]=h
 for st in ['QUADRATIC_CODIMENSION','CUBIC_CODIMENSION','GENERIC_RANKS']:
  pp='reports/S61_BINARY_HYPERPLANE_'+st+'_20261102.json';z=json.loads((ROOT/pp).read_text());assert z['all_passed'];hashes[pp]=sha(pp)
  for key in ['input_sha256','primary_input_sha256']:
   for path,h in z.get(key,{}).items():assert sha(path)==h,path;hashes[path]=h
 assert q['leading_constant_quadratic_rank']==132 and q['leading_cubic_rank']==31
 assert q['first_obstruction_rank']==41 and q['second_obstruction_rank']==18
 assert q['kernel_dimensions']==[106-31,75-41,34-18]==[75,34,16]
 assert q['universal_maximal_minor_valuation_lower']==sum(q['kernel_dimensions'])==125
 for a,n,r in zip(q['kernel_audits'],[75,34,16],[31,41,18]):assert a['all_passed']and a['kernel_columns']==n and a['rank_lower']==r and a['free_anchor_diagonal']and all(a['denominator_point_values'])
 mod=q['modular_certificate'];p=mod['prime'];assert mod['all_passed']and p==32003
 assert mod['generic_degree_bound_from_weights']==mod['selected_row_weight_sum']+mod['column_weight_sum']==384
 assert mod['conditional_QQ_d_degree_upper']==384-125==mod['determinant_degree']==259
 coeff=mod['determinant_coefficients_ascending'];assert len(coeff)==260 and coeff[-1]==6872
 values=mod['independent_determinant_evaluations']
 if isinstance(values,dict):values=[[int(k),int(v)]for k,v in values.items()]
 for pair in values:
  if isinstance(pair,dict):x=pair.get('point',pair.get('d'));v=pair.get('value',pair.get('determinant'))
  else:x,v=pair
  y=0
  for c in coeff[::-1]:y=(y*x+c)%p
  assert y==v
 det=1
 for x in mod['module_basis_diagonal']:assert x;det=det*x%p
 assert len(mod['module_basis_diagonal'])==238 and det==mod['module_basis_determinant']
 assert mod['module_rank']==238 and mod['module_degree']==0 and mod['row_module_uses_only_unit_constant_divisions']
 assert mod['rational_coefficients_p_integral']and mod['all328x238_coefficient_specializations_match']and mod['all321_generic_row_bounds_rechecked']
 h=q['heavy_square_certificate'];assert h['all_passed'];a,e,d=s.symbols('a e d');K=s.QQ.frac_field(a,e)
 cf=[s.Poly(s.sympify(z),d,domain=K)for z in h['raw_heavy_square_coefficients']]
 raw=sum((s.Poly(s.sympify(v),d,domain=K)*cf[int(i)]for i,v in h['bezout_raw'].items()),s.Poly(0,d,domain=K));assert raw==s.Poly(1,d,domain=K)
 ud=s.sympify(h['implicit_common_denominator'])
 actual=sum((s.Poly(s.sympify(v),d,domain=K)*s.Poly(cf[int(i)].as_expr()/ud,d,domain=K)for i,v in h['bezout_actual'].items()),s.Poly(0,d,domain=K));assert actual==s.Poly(1,d,domain=K)
 assert len([f for f in cf if f])==21
 low=q['Bzero_low_ranks'];assert low['all_passed']and[(v['rank'],v['minor'])for v in low['matrices']]==[(6,9925),(57,24692)]
 aa,ee,dd=map(s.Rational,low['point_rational']);B=50*aa*dd+18*dd*ee**2-9*dd*ee-2*dd+60*ee-40;assert B==0
 F=500*aa+90*ee**2+3*ee-24;D=50*aa+9*ee**2+3*ee-6;C=100*aa+18*ee**2-9*ee-2;assert aa*F*D*C
 assert q['weighted_Hilbert_function']==[1,7,23,34]and sum(q['weighted_Hilbert_function'])==q['weighted_algebra_free_rank_over_Kd']==65
 for path in [str(Path(__file__).relative_to(ROOT)),NOTE,primary,'research/notes/S61_BINARY_HYPERPLANE_BFREE_HORIZONTAL_20261102.md']:hashes[path]=sha(path)
 return dict(id=STEM,all_passed=True,accepted_scope='K=QQ(a,e); nuisance-zero weighted transverse algebra free65 over K[d], t^2=0; no vertical/full-core/ambient conclusion',proof_source_review=True,primary_fresh_rerun=False,raw135_second_implementation=False,large_polynomial_kernels_second_implementation=False,small_saved_arithmetic=dict(raw_and_actual_H2_Bezout_independently_multiplied=True,finite_field_det_evaluations=True,constant_triangular_diagonal_product=True,Bzero_point_known_open=True,rank_and_degree_chain=True),primary_report_sha256=sha(primary),primary_hashes_checked=len(q['input_sha256']),input_sha256=hashes,proof_note=NOTE,seconds=time.monotonic()-start)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
 with a.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps({'all_passed':True,'seconds':r['seconds'],'hashes':len(r['input_sha256'])}))
