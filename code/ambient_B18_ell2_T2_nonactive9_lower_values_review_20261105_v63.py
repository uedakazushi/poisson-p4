#!/usr/bin/env python3
"""Independent nine fixed first-product postcheck, no search or source solve."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_NONACTIVE9_LOWER_VALUES_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and r['code_sha256']==sha(ROOT/'code/ambient_B18_ell2_T2_nonactive9_lower_values_v2_20261105_v63.py')
 assert r['note_sha256']==sha(ROOT/'research/notes'/f'{BASE}.md')
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 own=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63_REVIEW/report.json');assert own['all_passed'] and sha(ROOT/'reports/AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63_REVIEW/report.json')==r['source_peer_sha256']
 ep=f'reports/{BASE}/actual-four-lower-by-nine-source-values-and-fixed-rank.json';d=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 s=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63/exact171-decomposition-162-zero-9-lower-first-source-input.json')
 live=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_LIVE_BOUND_20261105_v63/report.json');assert live['all_passed'] and live['proven_whole_live_bound']==7728<11408
 Hr,He=list(map(sp,s['Hstar1_r_selected99_columns'])),list(map(sp,s['Hstar1_ell_selected99_columns']))
 cols=[];D=int(s['common_integer_lower_D']);bound=s['all36_lower_positive_D_U']
 # Independent contraction order: contract each H column with chi first,
 # then the complete b1 vector. Preserve both ordinary mixed products.
 for j,rec in enumerate(s['remaining9_sources']):
  br,be=sp(rec['b1_r_selected99']),sp(rec['b1_ell_selected99']);col={}
  for a,chi in enumerate([(0,1,4),(11,0,270)]):
   hr=[sum((Q(c)*v.get(k,Q(0)) for k,c in enumerate(chi)),Q(0)) for v in Hr]
   he=[sum((Q(c)*v.get(k,Q(0)) for k,c in enumerate(chi)),Q(0)) for v in He]
   rr=-sum((hr[i]*x for i,x in br.items()),Q(0))
   re=-sum((hr[i]*x for i,x in be.items()),Q(0))-sum((he[i]*x for i,x in br.items()),Q(0))
   for row,x in [(2*a,rr),(2*a+1,re)]:
    Dj,Uj=map(int,bound[4*j+row]);assert Dj%x.denominator==0 and abs(x.numerator)*(Dj//x.denominator)<=Uj
    assert not x
  cols.append(col)
 assert cols==list(map(sp,d['actual_lower4_by9_columns']))==[{}]*9
 assert list(map(sp,d['integer_columns']))==[{}]*9 and int(d['common_integer_denominator'])==D==5628
 assert d['known_lower_quadratic_relation_rank']==0 and not d['selected_original9_columns'] and not d['normalized_basis']
 assert not d['fixed_rank_receipts'] and d['fixed_original_minor_rows']==[] and d['fixed_integer_minor_determinant']=='1'
 assert list(map(sp,d['all9_original_basis_factorizations']))==[{}]*9
 assert list(map(sp,d['kernel_source9_combinations']))==[{i:Q(1)} for i in range(9)]
 assert d['upper_components_not_zero_filled'] and not d['new_source_physical_gauge_or_implicit_coefficients']
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);e=out/'independent36-first-products-and-rank0-scope.json'
 e.write_text(json.dumps({'all36_exact_lower_values_zero':True,'fixed_rank':0,'kernel':'identity9','prior_whole_live':7728,'upper18_components_not_computed':True,'no_HF2_equality_total_or_torsion_claim':True,'proof':'Source peer proves the lower image exactly -chi(Hstar1*b1), including Schur signs and coherent normalization modulo Q2. Contracting H first independently reproduces all36 zeros; the rank0 factorization and I9 kernel are then literal. The nine upperQ2 components remain unknown.'},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'all36_lower_values_zero':True,'four_lower_image_rank':0,'upper_components_or_fullHF2_not_computed':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{e.name:sha(e)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'lower_rank':0}))
if __name__=='__main__':main()
