#!/usr/bin/env python3
"""Source-bound ell² torsion criterion; no coefficient or rank evaluation."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_TORSION_LENGTH_PREFLIGHT_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'length11':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json',
      'DVR':'reports/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63/DVR20-raw-cokernel26-and-comparison-preflight.json',
      'DVR_peer':'reports/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW/report.json',
      'rank4':'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json',
      'complete_source':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/report.json',
      'raw26':'reports/AMBIENT_B18_RAW_COKER26_CLOSED_ACTIONS_20261104_v63/report.json'}
    data={k:json.loads((ROOT/p).read_text()) for k,p in paths.items()}
    for k in ['length11','DVR_peer','rank4','complete_source','raw26']:assert data[k]['all_passed']
    t=data['length11'];assert (t['actual_T_length'],t['actual_mu_D_K'],t['actual_mu_D_torsion_K'])==(11,7,3)
    assert data['rank4']['rank_D_K']==4
    source=data['DVR']
    def length(typ,n):return typ['free_rank']*n+typ['D_mod_ell_squared_copies']*min(n,2)+typ['D_mod_ell_copies']
    product=[length(source['product_curve_ring_DVR_type'],n) for n in [1,2]]
    raw=[length(source['raw_cokernel_DVR_type'],n) for n in [1,2]]
    assert product==[20,33] and raw==[26,45]
    assert len(source['all12_physical_difference_actions_on_raw_cokernel'])==12
    evidence={
      'DVR':'D=k[[ell]], ell=5(sA-rA), not rA',
      'actual_K_DVR_rank':4,'actual_C_DVR_rank':4,'actual_K_torsion_cyclic_summands':3,
      'all_n_exact_sequence':'0 -> K/ell^nK -> T_n -> C[ell^n] -> 0',
      'actual_T_n_definition':'ker(V mod ell^n)/image(X mod ell^n), same allorder original source',
      'image_V_ell_regular_reason':'image V is a submodule of the free R-module R^3',
      'exact_difference':'d_n=len(T_n)-len(C/ell^nC)=sum_{i=1}^3 min(n,e_i)',
      'd1':11-8,'d2_range':[3,6],
      'positive_criterion':'d2=3 iff every e_i=1 iff length H0(D0)=3',
      'if_d2_gt3':'d2-3 counts precisely the exponents at least2; total length remains open',
      'general_count':'d_n-d_(n-1)=number of e_i>=n',
      'product_curve_ell_quotient_dimensions':product,
      'raw_cokernel_ell_quotient_dimensions':raw,
      'conditional_C2_matrix_shape':[45,12*45],
      'conditional_C2_dense_QQ_slots':45*12*45,
      'C2_source_gate':'actual ell² raw45 quotient basis and coordinate/physical actions, not the old raw26',
      'T2_source_gate':'honest ell² kernel/source-lift or complete compatible resolution before finite-base tensor',
      'unsafe_kernel_counterexample':{'R':'k[[ell,x]]','V':'x:R->R','X':'0',
        'true_T_n':0,'after_x_squared_quotient':'ker(x on R/(ell^n,x²))=x*(R/(ell^n,x²)) !=0'},
      'all_coefficient_rank_GB_implicit_jobs_run':False,
      'torsion_exponents_or_H0_length_computed':False}
    ep=out/'ell2-exact-criterion-and-source-gates.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    report={'id':STEM,'all_passed':True,'scope':'source-only exact criterion; ell² case not evaluated',
      'actual_H0_length_computed':False,'conditional_C2_QQ_slots':24300,
      'input_sha256':{p:sha(ROOT/p) for p in paths.values()},
      'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{ep.name:sha(ep)}}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps({'all_passed':True,'conditional_C2_QQ_slots':24300,'actual_H0_length_computed':False}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
