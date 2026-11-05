#!/usr/bin/env python3
"""Verify explicit source-coordinate obstructions, no adjoint solve."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_GRAPH_INVARIANT_ADJOINT6_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','cert':f'reports/{A}/actual-two-graph-invariant-six-adjoint-certificates.json',
'source':'reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json',
'source_peer':'reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def sp(p):return {i:Q(c) for i,c in p if Q(c)}
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and D['source_peer']['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 S=D['source'];C=D['cert'];cols=list(map(sp,S['six_rr_rell_source_adjoint_columns_in234_domain']));G1=list(map(sp,S['actual_G1linear3_by26_columns']));f=list(map(sp,S['two_invariant_functionals_on_all234_F1_coefficients']))
 assert len(cols)==6 and all(not p for p in cols) and not G1[12]
 assert sorted({k%9 for p in G1 for k in p})==list(range(7))
 for i,p in enumerate(C['results']):
  d=sp(p['restricted_support_obstruction234']);assert d=={115:Q(1)} and not p['rr_rell_representation_exists']
  assert all(not sum((c*d.get(k,Q(0)) for k,c in col.items()),Q(0)) for col in cols)
  assert sum((c*d.get(k,Q(0)) for k,c in f[i].items()),Q(0))==Q(p['nonzero_value'])!=0
 assert C['fixed_known_rank']==0 and C['fixed_receipts']==[]
 ep=o/'independent-dual234-obstructions-and-fullquadratic-limitation.json';dump(o,ep.name,{'inputs':inputs,'all_six_adjoint_columns_literal_zero':True,'both_nonzero_invariant_dual_pairings_verified':True,'same_obstruction_source_coefficient':'F1 source12 rA, scalar slot115','stronger_fullquadratic_obstruction':'G1linear sourcecolumn12 is literallyzero. Thus anF1 supportedon source12*rA is killedbytheENTIRE degree2G1linearmap, notonlyitsrr/rell outputs. Both desired invariants pairnontrivially withthissameinput. Neither canbe reconstructedfromanylinearfunctionalofall132 quadraticoutputs.',
'actual_M1_graph_or_invariant_values_not_evaluated':True,'no_zero_value_or_higher_jet_claim':True})
 report={'id':STEM,'all_passed':True,'both_six_adjoint_nonexistence_certificates_verified':True,'both_also_fail_fullquadratic_G1linear_output_reconstruction':True,'actual_M1_or_graph_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json')}))
if __name__=='__main__':main()
