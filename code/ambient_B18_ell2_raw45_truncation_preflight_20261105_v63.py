#!/usr/bin/env python3
"""Pure ell² source truncation theorem/input inventory; no rank or values."""
import argparse
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_RAW45_TRUNCATION_PREFLIGHT_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'source':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
      'source_peer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json',
      'product_peer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_CURVE_BASE_20261104_v63_REVIEW/report.json',
      'DVR':'reports/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63/DVR20-raw-cokernel26-and-comparison-preflight.json',
      'DVR_peer':'reports/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW/report.json'}
    data={k:json.loads((ROOT/p).read_text()) for k,p in paths.items()}
    for k in ('source_peer','product_peer','DVR_peer'):assert data[k]['all_passed']
    s=data['source'];plane={}
    for k in ('beta2_before_ell','beta3_before_ell'):
        plane[k]=[[[m[7:],c] for m,c in row if not any(m[:7])] for row in s[k]]
    assert sum(map(len,plane['beta2_before_ell']))==9
    assert sum(map(len,plane['beta3_before_ell']))==4
    assert len(s['gA_plane3'])==len(s['gB_own_plane3'])==4
    assert s['YB_in_YA_linear'][7:]==[
      ['0']*7+['7/12','5/6'],['0']*7+['5/12','1/6']]
    d=data['DVR'];assert d['both_actual_CI_curve_modules_are_DVR_free_rank3']
    assert d['product_curve_ring_DVR_type']=={'free_rank':6,'D_mod_ell_squared_copies':7,'D_mod_ell_copies':7,'special_fiber_dimension':20}
    source={'beta_graph_plane_through3':plane,'gA_plane3':s['gA_plane3'],
      'gB_own_plane3':s['gB_own_plane3'],'YB1_bottom_plane2_by2':[row[7:] for row in s['YB_in_YA_linear'][7:]],
      'ell_substitution':'sA=r+ell/5, ell²=0',
      'z_dependent_beta3_change':'unit I+M2 removes these terms; changes graph only in degree>=4',
      'same_source_physical_rows':s['actual_physical_rows9']}
    proof={
      'scope':'pure exact truncation theorem, not a computed ell² quotient',
      'union_D_basis6':['(1,1)','(r,r)','(r²,0)','(0,r²)','(ell²,0)','(ell*r,0)'],
      'matching_over_A3':'a0-b0 in ell²D; a1-b1 in ellD; a2,b2 arbitrary',
      'Weierstrass_coefficient_valuation_lower_bounds':[1,2,3],
      'union_m4_bound_including_Tor':True,
      'product_m6_bound':True,
      'exact_truncation':'(IA*IB,ell²)=(IA_le3*IB_le3,ell²)',
      'proof_of_equality':'apply J*m4S=0 and m6S=0 in both actual and truncated rings; generator changes are O(m4)',
      'module_generators8':['1']+[f'z{i}' for i in range(7)],
      'module_relation_families':['21 beta_j zi-beta_i zj','7 gA(zi+beta_i)','7 gB zi'],
      'zi_action':'zi*1=zi; zi*zj=-beta_j*zi',
      'z_closure':'pair and gA rows map to zero; gB rows map to beta multiples of gB rows',
      'ambient_monomials':{'plane':11,'z_sectors':7*9,'total':74},
      'row_bound':21*5+14*3,'dense_slots_upper_bound':74*(21*5+14*3),
      'conservative_input_nonzeros_bound':1470,
      'predicted_actual_product_dimension_from_accepted_DVR_type':33,
      'predicted_raw_cokernel_dimension':45,
      'rank_source_coefficient_or_IFT_job_run':False,
      'physical12_action_and_honest_T2_still_separate_gates':True}
    rejected={'rejected_statement':'U/ell² injects into RA/ell² directsum RB/ell² because U is ell-regular',
      'reason':'Tor1_D(A3,D/ell²)=A3 survives; ell-regularity does not remove this boundary',
      'repair':'explicit six-element fiber-product D basis proves the m4 bound including the attachment',
      'numerical_job_before_repair':False}
    evidence={}
    for name,obj in [('literal-existing-cubic-source-input.json',source),('union-six-basis-truncation-and-module-preflight.json',proof),('rejected-tensor-injection-shortcut.json',rejected)]:
        p=out/name;p.write_text(json.dumps(obj,separators=(',',':'))+'\n');evidence[name]=sha(p)
    report={'id':STEM,'all_passed':True,'scope':'source theorem/preflight only; numerical quotient not run',
      'conditional_module_shape_bound':[147,74],'conditional_input_nonzeros_bound':1470,
      'input_sha256':{p:sha(ROOT/p) for p in paths.values()},'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':evidence}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps({'all_passed':True,'shape_bound':[147,74],'numerical_quotient_run':False}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
