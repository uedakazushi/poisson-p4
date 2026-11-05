#!/usr/bin/env python3
"""Independent truncated free-polynomial replay of the formal high3 fold."""
import json
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_UPPER563_HIGH3_FORMAL_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 pins={}
 def read(p):pins[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and pins[rp]['sha256']=='9a938f415c5f80700d35a33fc55e5853b05314b9f5cae60bec28ee76a87e5f53'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/two-formal-direct-old-highestthird-sensitivity-rows.json';e=read(ep);assert pins[ep]['sha256']==r['evidence_sha256'][Path(ep).name]
 gb='reports/AMBIENT_B18_ELL2_T2_UPPER563_HIGH3_INPUT_GATE_20261105_v63';g=read(gb+'/report.json');d=read(gb+'/direct-old-highestthird-origin-dual-source-and-positive-input.json');assert g['all_passed'] and sha(ROOT/gb/'report.json') in json.dumps(read(gb+'_REVIEW/report.json'))
 maximum=0
 def ck(x):
  nonlocal maximum
  maximum=max(maximum,abs(x.numerator).bit_length(),x.denominator.bit_length());assert maximum<=57;return x
 def add(a,b,c=Q(1)):
  z=dict(a)
  for m,x in b.items():
   v=ck(z.get(m,Q(0))+ck(c*x))
   if v:z[m]=v
   else:z.pop(m,None)
  return z
 def mul(a,b):
  z={}
  for ma,x in a.items():
   for mb,y in b.items():
    m=ma+mb
    if len(m)<2:z=add(z,{m:ck(x*y)})
  return z
 poly={};originbound={i:tuple(map(int,v)) for i,v in d['positive_node_origin_envelopes']};bound={i:{j:tuple(map(int,v)) for j,v in row} for i,row in d['positive_node_high3_sensitivity_envelopes']}
 for i,n in d['typed_source_program']:
  if n[0]=='actual_C0_scalar_function':v={(n[1],):Q(1)}
  elif n[0]=='c':v={():Q(n[1],n[2])} if n[1] else {}
  elif n[0]=='sum':
   v={}
   for j,a,b in n[1]:v=add(v,poly[j],Q(a,b))
  else:
   assert n[0]=='product';v={():Q(n[1],n[2])}
   for j in n[3]:v=mul(v,poly[j])
  poly[i]=v
  for mon,x in v.items():
   D,U=originbound[i] if not mon else bound[i][mon[0]];assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
 target=[{},{}]
 for term in d['original_factor_terms']:
  v=poly[term['root']] if term['root_kind']=='typed' else {(term['root'],):Q(1)};target[term['channel']]=add(target[term['channel']],v,Q(term['coefficient_factors'][0])*Q(term['coefficient_factors'][1]))
 outrows=[]
 for row in target:
  assert not row.get((),Q(0));outrows.append({m[0]:x for m,x in row.items() if m})
 assert outrows==list(map(sp,e['two_direct_OLD_highestthird_symbol_coefficients'])) and all(set(row)=={327603,327627,505938} for row in outrows)
 owned=2*count(d)+45*26+2*26+600+1000;assert owned<7618
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ve=out/'independent-free-polynomial-linear-truncation-high3-proof.json';ve.write_text(json.dumps({'two_formal_direct_old_high3_rows':[[[i,str(x)] for i,x in sorted(row.items())] for row in outrows],'all45_free_affine_polynomials_and_positive_bounds_checked':True,'remaining_symbol_ids':[327603,327627,505938],'all_symbols_unassigned':True,'maximum_known_coefficient_bits':maximum,'conservative_owned':owned,'scope':'OnlydirectOLDGammaXc3highestthirdsymbolcoefficients. No scalar third, full upper relation, source-order necessity or native/carry cancellation conclusion.'},separators=(',',':'))+'\n');rr={'id':STEM,'all_passed':True,'primary_sha256':pins[rp]['sha256'],'formal_six_coefficients_verified':True,'no_actual_scalar_or_upper_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ve.name:sha(ve)}};(out/'report.json').write_text(json.dumps(rr,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
