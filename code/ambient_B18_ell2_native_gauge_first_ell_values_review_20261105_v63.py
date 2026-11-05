#!/usr/bin/env python3
"""Postverify provided30 gauge derivatives; no solve or source evaluation."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
from ambient_B18_ell2_native_gauge_first_ell_availability_review_20261105_v63 import sparse,mv
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import add,enc
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_20261105_v63';STEM=A+'_REVIEW'
G='AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_PREFLIGHT_20261105_v63'
P={
 'primary':f'reports/{A}/report.json','value':f'reports/{A}/actual30-firstell-and-native-S1ell.json',
 'gate':f'reports/{G}/actual-first-ell-input-and-old8-source-crosschecks.json',
 'gate_peer':f'reports/{G}_REVIEW/report.json',
 'S':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'PROVIDED30_CANONICAL_RESPONSE_POSTCHECK','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and D['gate_peer']['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 E=D['value'];GATE=D['gate'];phi=list(map(Q,GATE['actual_phi_first_ell']));gamma=list(map(Q,E['actual_gauge30_first_ell']))
 assert E['actual_phi_first_ell']==GATE['actual_phi_first_ell'] and len(gamma)==30
 JP=list(map(sparse,GATE['Jg_phi_sparse_rows']));Jgamma=list(map(sparse,GATE['Jg_gamma_sparse_rows']));target=list(map(sparse,GATE['G32_nu2_target_Jacobian_sparse_rows']))
 rhs=mv(JP,phi);assert list(map(str,rhs))==E['RHS30']
 assert all(a+b==0 for a,b in zip(mv(Jgamma,gamma),rhs))
 response=mv(target,phi+gamma);assert response[32:]==[0,0]
 g=[response[i*4:i*4+4] for i in range(8)]
 assert [list(map(str,p)) for p in g]==E['actual_G8_C4_first_ell'] and E['actual_nuA_nuB_first_ell']==['0','0']
 native=D['S']['native_S12_by32'];assert native==GATE['native_S12_by32']
 out=[{} for _ in range(128)]
 for row,j,k,sign in native['koszul28_in_first8_rows']:
  for mask,c in enumerate(g[k]):
   if not c:continue
   for source_mask in range(4):
    if mask & source_mask:continue # Both original monics and ell derivatives vanish here.
    idx=4*j+source_mask;out[idx]=add(out[idx],{4*row+(mask|source_mask):Q(sign)*c})
 assert list(map(enc,out))==E['actual_native_S1_ell_48_by128_sparse_columns']
 assert all(not p for p in out[112:])
 assert E['no_native_L1_or_N1_evaluated'] and E['no_individual_actual_native_X1_claim']
 ep=o/'independent-provided30-G8-and-nativeS1ell-proof.json'
 dump(o,ep.name,{'inputs':inputs,'all30_provided_derivatives_satisfy_canonical_constraints':True,'uniqueness_premise':'same fixedJ30 unit from frozen independent availability peer; no inverse or unknown solve repeated','all34_actual_target_responses_checked':True,'actual_monic_ell_response_zero':True,'all48_by128_intrinsic_native_S1ell_entries_reconstructed':True,'intrinsic_native_firstell_nonzeros':sum(map(len,out)),'free4_source_column_derivatives_zero':True,'effective_use_only':'L0*S1ell*cNative on originalsourcecombinations withS0*cNative=0; no individual actualnativeX1','new_physical_source_or_higherjet_or_N1_or_T2_value_computed':False})
 report={'id':STEM,'all_passed':True,'provided30_firstell_gauge_G8_and_intrinsic_S1ell_accepted':True,'sourceframe_and_effective_only_scope_accepted':True,'unknown_solve_or_physical_coefficient_or_actual_T2_value_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
