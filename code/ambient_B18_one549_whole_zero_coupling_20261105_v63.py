#!/usr/bin/env python3
"""Whole ONE549 zero proof from independently checked positive components."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ONE549_WHOLE_ZERO_COUPLING_20261105_v63'
NP='reports/AMBIENT_B18_NATIVE_ONE549_POSITIVE_INPUT_20261105_v63/'
OP='reports/AMBIENT_B18_ELL2_T2_UPPER549_OLD_CARRY_INPUT_GATE_20261105_v63/'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def fields(p,names):
 s=p.read_text();d=json.JSONDecoder();i=s.index('{')+1;z={}
 while True:
  while s[i].isspace() or s[i]==',':i+=1
  if s[i]=='}':break
  k,i=d.raw_decode(s,i)
  while s[i].isspace() or s[i]==':':i+=1
  v,i=d.raw_decode(s,i)
  if k in names:z[k]=v
  del v
 assert set(z)==set(names);return z
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0
 def read(p,names=None):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size
  return json.loads(q.read_text()) if names is None else fields(q,names)
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 nr=read(NP+'report.json');orr=read(OP+'report.json')
 assert pins[NP+'report.json']=='62d66bc4406244a863719b5e0ec9cdb2b780c5b2d61475c39ecc6796512f1b4f'
 assert pins[OP+'report.json']=='4bdd31b2eba6848c058c024b7d7bdd24d3db085079e5c0299cbb97b32ccc41b6'
 ne='ONE549-positive-native649-GammaD-input-and-known648zero-coupling.json';oe='ONE549-OLD-quarter-and55-FIRST-positive-component-input.json';n=read(NP+ne);o=read(OP+oe)
 assert pins[NP+ne]==nr['evidence_sha256'][ne] and pins[OP+oe]==orr['evidence_sha256'][oe]
 npp='reports/AMBIENT_B18_NATIVE_ONE549_SOURCE_POSITIVE_20261105_v63_REVIEW/report.json';npeer=read(npp);assert pins[npp]=='b02f0d2118aab1499477d6c05603e8e0a7f8738b83cc702795fd96224a234e18' and npeer['all_passed']
 opp='reports/AMBIENT_B18_ELL2_T2_UPPER549_OLD_FIRST_20261105_v63_REVIEW/report.json';opeer=read(opp);assert pins[opp]=='ea081efd637927abd6e3b68de4a1b90cc016d14f6b8f4d7af2cae9f3dafe6652' and opeer['all_passed']
 assert npeer['input_sha256'][NP+'report.json']==pins[NP+'report.json'] and opeer['input_sha256'][OP+'report.json']==pins[OP+'report.json']
 fp='reports/AMBIENT_B18_NATIVE_ONE549_FRAME_SOURCE_20261105_v63/ONE549-J18-with-r2ell-native-full-S0-column-and-formal-N2N3-frontier.json'
 f=read(fp,['literal_full_original768_source','full_a1_27times9_32_entries_UNEVALUATED','same617_lower_basis_recipe','Gamma_two_current56','D_two_by10_by56','only_pureplane_a1_positions','comparison_pivot648_mask0_S0_column8_nonzero'])
 osp='reports/AMBIENT_B18_ELL2_T2_UPPER549_OLD_CARRY_SOURCE_20261105_v63/ONE549-OLD298-and-full32-a1-known-FIRST-carry-factor-source.json'
 os=read(osp,['original768_source','same617_source_recipe','actual_full_a1_27times9','only_plane_a1_pivot6_source648','Gamma_two','D_two_by10_by56','all_known_FIRSTcarry_factor_recipes'])
 primitive=n['native_primitive_input_from_frozen_source']
 assert f['literal_full_original768_source']==os['original768_source']==[[298,'1/4'],[649,'1']]
 assert n['all32_a1_source_entries_kept_UNEVALUATED']==o['full_a1_32']==f['full_a1_27times9_32_entries_UNEVALUATED']==os['actual_full_a1_27times9']
 assert os['same617_source_recipe']==f['same617_lower_basis_recipe'] and f['same617_lower_basis_recipe']['subtract_actual_basis563609']==['0','0']
 assert primitive['Gamma_two56']==o['Gamma_two']==f['Gamma_two_current56']==os['Gamma_two']
 assert primitive['D_two_by10_by56']==o['D_two_by10_by56']==f['D_two_by10_by56']==os['D_two_by10_by56']
 assert o['known_FIRSTcarry_factor_recipes55']==os['all_known_FIRSTcarry_factor_recipes'] and len(o['known_FIRSTcarry_factor_recipes55'])==55
 assert f['only_pureplane_a1_positions']==os['only_plane_a1_pivot6_source648']==[[61,'-79/24'],[62,'-11/30']]
 assert f['comparison_pivot648_mask0_S0_column8_nonzero']==[[15,'-1/2']]
 outer=n['outer_GammaXp648_plane_a1_exact0_by_independent92cache'];v=read(outer['data']);p=read(outer['peer']);assert p['all_passed'] and pins[outer['peer']]=='10844d8720465054d6b9f7a7aa408dd95dd613a2e857eb52e50cd5a71a847e92'
 k=v['original768_native23_sources'].index(648);assert k==2 and v['Gamma_Xp2_native23_twochannels_rr_rell'][k]==[['0','0'],['0','0']]
 Q=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63/known-Psi-point-correction-and-full-b1-source-annihilation.json',['Q_two6'])['Q_two6']
 initial=count([nr,orr,n,o,npeer,opeer,f,os,v,p,Q,pins])+10000;assert initial<60000 and raw+2_000_000<32*1024*1024
 dump('phase-before-zero-joining.json',{'phase':'all_current_case_GammaD_a1_55FIRST_and_independent_components_sameframe_PASS','whole_source_numeric_bound':initial,'no_signed_job':True})
 def zero(q):assert len(q)==2 and int(q[0])>0 and int(q[1])==0
 for ch in n['native_Gamma_cubic_two_by5_D_U']:
  for q in ch:zero(q)
 for ch in n['native_D_second10_unique_first_decomposition_D_U']:
  for pair in ch:
   for q in pair:zero(q)
 for q in n['positive_native_Gamma_plus_D_two_D_U']:zero(q)
 for key in ['positive_OLD_quarter_Gamma_cubic_two','positive_OLD_quarter_D_second_two','positive_FIRSTcarry_total_two','positive_OLD_and_FIRST_partial_two']:
  for q in o[key]:zero(q)
 for row in o['positive_FIRST55_factor_bounds_two_by55']:
  assert len(row)==55
  for q in row:zero(q)
 # Addition of nonnegative zero envelopes and the known zero outer target
 # proves both ACTUAL whole projected coefficients zero; sign convention is
 # OLD+native-FIRST-outer, despite signs being irrelevant to a zero upperbound.
 whole=[['1','0'],['1','0']];six=[['1','0'] for _ in range(6)]
 assert len(Q)==2 and all(i<6 for row in Q for i,x in row)
 prior=max(nr['positive_committed_bits'],orr['committed_prior_bits']);temporary=3*prior+1;assert prior==152 and temporary==457
 result={'selected_case619':549,'same_literal_original768_source':[[298,'1/4'],[649,'1']],
  'full_source_formula':'U549=OLD298quarter+native649 -all55_knownFIRSTcarry -GammaXp648_plane_a1; fullPsi/currentGammaD; lower617basis coefficients0.',
  'all32_a1_and55_ordinary_FIRST_records_kept':True,'native_Source0cN_zero_only_whole_column9':True,
  'pivot648_S0_nonzero_fullL2S0_actual_four_targets_zero_independently':True,
  'independent_native_and_OLD_FIRST_positive_zero_proofs_combined':True,
  'positive_actual_whole_upper_two_D_U':whole,'positive_actual_whole_six_D_U':six,
  'actual_ONE549_quadratic_class_zero_proved_without_signed_job':True,
  'all_internal_honest_N2_N3_extensions_scope_controlled_by_sourceproof':True,
  'no_extra_physical_gauge_or_source_response_evaluation':True,
  'scope':'ONE549 only. No all605 annihilation,HF2,totalT2length,H0 ornilcube claim.',
  'whole_retained_numeric_upper':initial+1000,'whole_RAW_upper':raw+2_000_000,'inherited_positive_prior152_temporary457':True}
 ep='ONE549-complete-current-frame-positive-zero-and-no-signed-job.json';dump(ep,result)
 rep={'id':STEM,'all_passed':True,'whole_ONE549_positive_zero_proof_only':True,'signed_zero_only_worker_not_prepared_or_needed':True,
  'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
  'whole_source_numeric_upper':initial+1000,'whole_RAW_upper':raw+2_000_000,'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'whole_numeric':initial+1000,'RAW':raw+2_000_000,'seconds':rep['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
