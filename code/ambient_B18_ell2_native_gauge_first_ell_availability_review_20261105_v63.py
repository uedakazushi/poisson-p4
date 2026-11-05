#!/usr/bin/env python3
"""Availability/source postcheck only: no new first-ell gauge value."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_PREFLIGHT_20261105_v63';STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json','gate':f'reports/{A}/actual-first-ell-input-and-old8-source-crosschecks.json',
 'J':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63/phase1-actual-J135-blocks-and-target-Jacobian.json',
 'J_peer':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63_REVIEW/report.json',
 'phi':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase1-actual105-B56-cache-composition.json',
 'values':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase2-all105-equations-and-primitive56-source-columns.json',
 'values_peer':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63_REVIEW/report.json',
 'oldJ':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase2-fixed30-Jacobian.json',
 'gamma':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase3-gauge30-firsts.json',
 'G8':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
 'nu':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase1-source-firsts.json',
 'native_peer':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63_REVIEW/report.json',
 'S':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def sparse(p):return {int(i):Q(c) for i,c in p if Q(c)}
def mv(A,v):return [sum((c*v[j] for j,c in row.items()),Q(0)) for row in A]
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'OLD8_POSTCHECK_AND_NEWELL_AVAILABILITY_ONLY','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert all(D[k]['all_passed'] for k in ['primary','J_peer','values_peer','native_peer'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 pin=D['J_peer']['input_sha256']['J'];assert pin==inputs['J']
 J=D['J'];E=D['gate'];assert J['source_dependent_monic_carry_included'] and J['all_base_G8_and_gauge_constraints_checked']
 assert J['physical105_order']==D['phi']['unknown105_order']==E['physical105_order']
 assert J['canonical_gauge30_order']==E['canonical_gauge30_order']==['even0..14','odd0..14']
 Aphi=list(map(sparse,J['Jg_phi30_by105']));Agamma=list(map(sparse,J['Jg_gamma30']));Ai=list(map(sparse,J['Jg_gamma_inverse']));target=list(map(sparse,J['properG32_and_nu2_target_Jacobian34_by135']))
 for key,values in [('Jg_phi_sparse_rows',J['Jg_phi30_by105']),('Jg_gamma_sparse_rows',J['Jg_gamma30']),('Jg_gamma_inverse_sparse_rows',J['Jg_gamma_inverse']),('G32_nu2_target_Jacobian_sparse_rows',J['properG32_and_nu2_target_Jacobian34_by135'])]:assert E[key]==values
 fixed=D['oldJ'];assert Q(fixed['det_J30'])==Q(4,9)
 expected=[]
 for i in range(30):
  row={}
  for j in range(30):
   c=Q(fixed['A0'][i%15][j%15]) if (i//15)==(j//15) else Q(fixed['A1_eta_delta'][i-15][j]) if i>=15 and j<15 else Q(0)
   if c:row[j]=c
  expected.append(row)
 assert Agamma==expected and Ai==[{j:Q(c) for j,c in enumerate(row) if Q(c)} for row in fixed['J30_inverse']]
 for j in range(30):
  u=[Q(i==j) for i in range(30)];assert mv(Agamma,mv(Ai,u))==u and mv(Ai,mv(Agamma,u))==u
 oldgamma=D['gamma']['even15_first9']+D['gamma']['odd15_first9'];origin=[Q(p[0]) for p in oldgamma]
 assert origin[:15]==[0]*15 and origin[15:]==list(map(Q,D['S']['group30_to120']['group30_odd_origin']))
 basis=[tuple(m) for m in D['phi']['basis56']];ix={m:i for i,m in enumerate(basis)};phi=list(map(sparse,D['phi']['all105_cached_values_B56']))
 first=[]
 for k in range(9):
  m=tuple(int(i==k) for i in range(9));first.append([p.get(ix[m],Q(0)) for p in phi])
 assert list(map(str,first[8]))==E['actual_phi_first_ell']
 checks=[]
 for k in range(8):
  g=[Q(row[k+1]) for row in oldgamma];rhs=mv(Aphi,first[k]);assert all(a+b==0 for a,b in zip(mv(Agamma,g),rhs))
  response=mv(target,first[k]+g);expectedG=[Q(D['G8']['G8'][i][4*(k+1)+mask]) for i in range(8) for mask in range(4)]
  assert response[:32]==expectedG and response[32:]==[Q(row[k+1]) for row in D['nu']['nu_first9']]
  record=E['all8_old_direction_crosschecks'][k];assert record=={'direction':k,'all30constraints_zero':True,'gauge30':list(map(str,g)),'G32_nu2':list(map(str,response))}
  checks.append({'direction':k,'all30provided_gauge_constraints_and34target_responses_verified':True})
 # These monic coefficients were already evaluated in the source cache.
 # Checking them here computes no new gauge coefficient.
 assert all(not any(j>=105 for j in row) for row in target[32:])
 known_nu=[sparse(D['values'][name]).get(ix[(0,)*8+(1,)],Q(0)) for name in ['nuA','nuB']]
 assert mv(target[32:],[*first[8],*([Q(0)]*30)])==known_nu==[0,0]
 assert E['new_ell_gauge_and_G8_values_evaluated'] is False and E['physical_implicit_higherjet_or_N1_queries'] is False
 assert E['native_S12_by32']==D['S']['native_S12_by32']
 assert E['whole_slot_bound']==sum(E['slot_breakdown'].values())==12051<20000
 assert E['limits']=={'seconds':60,'bytes':2_000_000_000,'artifact_bytes':32*1024*1024,'slots':20000,'coefficient_bits':384}
 ep=o/'independent-same-branch-Jacobian-old8-and-firstell-source-gate.json'
 dump(o,ep.name,{'inputs':inputs,'same_actual_physical105_and_even15_odd15_gauge_orders':True,'saved_actual_odd_center_and_J30_det4over9_and_bothinverseproducts_verified':True,'all8_old_provided_gauge_and_G32_nu2_responses_verified':checks,'phiell_genuine_common105_cache_and_sA_equals_rA_plus_ell_over5':True,'source_trace_functions_not_replaced_by_F9s':True,'known_monic_firstell_zero_independently_checked':True,'new30_gammaell_and_G8ell_values_not_evaluated':True,'unique_new_response_formula':'gammaell=-Jgamma_inverse*Jphi*phiell; targetell=Jtarget*(phiell,gammaell)','whole_slot_bound':12051,'conditional_numeric_worker_scope':'30gaugecoefficients and derivedG32/nu/S1ell only; noProperAD/physical105solve/N1','effective_native_use_requires_S0cNative_zero':True})
 report={'id':STEM,'all_passed':True,'same_canonical_branch_and_existing_Jacobian_firstell_availability_accepted':True,'all8_old_source_gauge_target_and_monic_crosschecks_accepted':True,'new_firstell_gauge_or_target_values_evaluated':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
