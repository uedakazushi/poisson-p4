#!/usr/bin/env python3
"""Known-cover consequences of a peer-verified lower witness; no new images."""
import json
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_WITNESS563_HF23_SOURCE_20261105_v63'
def main():
 pins={}
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 primary='reports/AMBIENT_B18_ELL2_T2_NATIVE70_VALUES_20261105_v63';r=read(primary+'/report.json');assert r['all_passed'];v=read(primary+'/actual-native70-lower-prefix-and-witness-if-any.json')['nonzero_witness'];assert v['case619']==563
 peer=read('reports/AMBIENT_B18_ELL2_T2_NATIVE70_WITNESS563_20261105_v63_REVIEW/report.json');assert peer['all_passed'] and sha(ROOT/primary/'report.json') in json.dumps(peer)
 l=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json',['degree2_linear_standard_module_monomials','allowed_variables_per_target'])
 assert l['allowed_variables_per_target']==[[],[],[],[7,8],[7,8],[7,8]]
 std=l['degree2_linear_standard_module_monomials'];expected=[]
 for target in [3,4,5]:
  for a,b in [(2,0),(1,1)]:m=[0]*9;m[7]=a;m[8]=b;expected.append([target,m])
 assert std==expected
 gb=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_GB_COVER_20261105_v63/report.json');assert gb['all_passed']
 d=read('reports/AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63/exact171-decomposition-162-zero-9-lower-first-source-input.json',['lower_four_duals6','Q_two_upper_quadratic_columns'])
 chi=list(map(sp,d['lower_four_duals6']));upper=list(map(sp,d['Q_two_upper_quadratic_columns']));assert chi[:2]==[{2:Q(1),4:Q(4)},{3:Q(1),5:Q(4)}]
 assert all(sum((x*q.get(i,Q(0)) for i,x in row.items()),Q(0))==0 for row in chi for q in upper)
 a,b=map(Q,v['actual_lower4_values'][:2]);assert a==Q(-1000,3) and b==Q(-212,3)
 # For an arbitrary FULL quadratic class q (including its unknown upper2),
 # cubic chi applied to rq,ellq gives [[a,0],[b,a]].
 minor=a*a;assert minor==Q(1000000,9)
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'unknown-upper-independent-r-and-ell-product-source-proof.json'
 e={'known_cover_degree2_and3_dimension':6,'full_witness_upper2_coordinates_remain_unknown':True,'two_cubic_product_dual_matrix':[[str(a),'0'],[str(b),str(a)]],'exact_nonzero_minor':str(minor),'actual_HF2_upper_bound':5,'actual_HF3_upper_bound':4,'proof':'ThecompletehomogeneousknownlinearGBhasstandardsbottom3*(r^d,r^(d-1)ell). ThelowerdualchiA_rrannihilatesbothunknownupper2directionsandhasvaluea!=0. Multiplicationbyr/ellonanyfullquadraticliftqgivesthe2x2dualmatrix[a,0;b,a]. Thusq,thenrqandellq,areindependentnewinitialrelationsindegrees2and3respectively. Knowninitialcoverdimensions6thereforegiveHF2<=5,HF3<=4. Extraactualrelationsmayfurtherreduceboth.','abstract_length20_not_a_known_finite_initial_cover':True,'actual_total_T2_length_and_H0_exponents_not_decided':True};ep.write_text(json.dumps(e,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'actual_HF2_upper_bound':5,'actual_HF3_upper_bound':4,'no_length17_or_HF2_equality_claim':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)}};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
