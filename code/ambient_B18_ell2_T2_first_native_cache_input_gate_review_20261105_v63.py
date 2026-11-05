#!/usr/bin/env python3
"""Independent cached native C16 input gate; no N or S2 reconstruction."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time,resource,signal
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,mv
from ambient_B18_ell2_T2_first_X2_source_fold_review_20261105_v63 import jetmul
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','gate':f'reports/{A}/actual-native-Rell4-cache-source-input-and-adapter-gate.json',
'raw':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
'raw_peer':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json',
'ms_report':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/report.json',
'ms_peer':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW/report.json',
'gamma':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/phase2-canonical30-exact-values.json',
'proper':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/phase3-actual-kappa9-and-normalized-nilframe.json',
'nu':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/phase1-source150-readback-fingerprints.json',
'phi':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',
'phi_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
'B56':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase1-actual105-B56-cache-composition.json',
'firstell':'reports/AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_20261105_v63/actual30-firstell-and-native-S1ell.json',
'firstell_peer':'reports/AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_20261105_v63_REVIEW/report.json',
'old':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
'old_peer':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json',
'physical':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json'}
MON=[(0,0),(1,0),(0,1),(1,1)]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def dense(p):return [str(p.get(m,Q(0))) for m in MON]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'NATIVE_CACHED_INPUT_ONLY','inputs':inputs});D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert D['primary']['all_passed'] and all(D[k]['all_passed'] for k in ['raw_peer','ms_report','ms_peer','phi_peer','firstell_peer','old_peer'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 E=D['gate'];R=D['raw'];images=[]
 for sector,r,e in R['actual33_basis']:images.append({(r,e):Q(1)} if sector==0 and r<=1 and e<=1 else {})
 def project(v):
  out={}
  for i,c in v.items():
   for k,w in images[i].items():out[k]=out.get(k,Q(0))+c*w
  return {k:c for k,c in out.items() if c}
 actions=[list(map(sp,M)) for M in R['coordinate_actions9_sparse_columns']]
 for i,(sector,r,e) in enumerate(R['actual33_basis']):
  for j in range(33):
   v={j:Q(1)}
   for k in ([sector-1] if sector else [])+[7]*r+[8]*e:v=mv(actions[k],v)
   assert project(v)==jetmul(images[i],images[j])
 for k,M in enumerate(actions):
  var={} if k<7 else {MON[k-6]:Q(1)}
  for j,p in enumerate(M):assert project(p)==jetmul(var,images[j])
 assert list(map(dense,images))==E['raw33_to_Rell4_columns']
 def scalar(v):return project({i:Q(c) for i,c in enumerate(v) if Q(c)})
 def trace(p):
  vals=[Q(0)]*16
  for i,c in enumerate(map(Q,p)):
   base,mask=divmod(i,4)
   for m,w in images[base].items():vals[4*MON.index(m)+mask]+=w*c
  return list(map(str,vals))
 even=list(map(scalar,D['gamma']['actual_even15_B33']));odd=list(map(scalar,D['gamma']['actual_odd15_B33']));actual=list(map(trace,D['proper']['actual14_C132']));kappa=list(map(trace,D['proper']['kappa9_C132']));nu=list(map(scalar,D['nu']['nu_A_B33']))
 assert list(map(dense,even))==E['actual_even15_Rell4'] and list(map(dense,odd))==E['actual_odd15_Rell4'] and actual==E['actual14_C16'] and kappa==E['actual_kappa9_C16'] and list(map(dense,nu))==E['actual_nu_A_B_Rell4']
 assert [p.get((0,1),Q(0)) for p in even+odd]==list(map(Q,D['firstell']['actual_gauge30_first_ell']))
 assert all(set(p)<= {(1,1)} for p in nu) and [p.get((0,1),Q(0)) for p in nu]==list(map(Q,D['firstell']['actual_nuA_nuB_first_ell']))
 assert [v[:8] for v in actual]==[list(map(str,map(Q,v[:8]))) for v in D['old']['actual14_C_B44']] and [v[:8] for v in kappa]==[list(map(str,map(Q,v[:8]))) for v in D['old']['kappa9_C_B44']]
 phi={i:project(sp(p)) for i,p in zip(D['phi']['joint105_unknown_order'],D['phi']['all105_unknown_final_values'])};bm=list(map(tuple,D['B56']['basis56']))
 for i,p in zip(D['B56']['unknown105_order'],D['B56']['all105_cached_values_B56']):
  q={}
  for j,c in sp(p).items():
   mon=bm[j]
   if not any(mon[:7]) and (mon[7],mon[8]) in MON:q[mon[7],mon[8]]=c
  assert q==phi[i]
 assert [[i,dense(v)] for i,v in sorted(phi.items())]==E['actual105_Rell4']
 physical=[]
 for p in D['physical']['physical12_coefficients']:
  parts=[project(sp(p[k])) for k in ['delta','bA','bB']];v=[Q(0)]*16
  for m in MON:
   for mask,(q,sign) in enumerate(zip(parts,[1,1,-1])):v[4*MON.index(m)+mask]=sign*q.get(m,Q(0))
  physical.append(list(map(str,v)))
 assert physical==E['actual_physicalF12_C16'] and [p['physical_row'] for p in D['physical']['physical12_coefficients']]==[18,19,28,29,60,67,72,63,64,68,73,74]
 assert sum(E['slot_breakdown'].values())==65484==E['whole_slot_bound']<150000
 ep=o/'independent-native-C16-source-frame-and-cache-inputs.json';dump(o,ep.name,{'inputs':inputs,'all1089_raw33_ring_products_and297_action_intertwiners':True,'actual105_B56_and_raw33_projection_overlap_exact':True,'actual30_canonical_gamma_ell_matches_independent_firstell':True,'all14_and_kappa9_constant_and_rfirst_match_accepted_B11':True,'all12_physical_tensor_rows_and_signs_verified':True,'actual_nu_A_B_rell':list(map(dense,nu)),'whole_slot_bound':65484,'source_scope':'Same actual canonical even15+odd15*tauA/auxetaA branch, physical original105 retained; no F9 trace substitution. All normalizer finite telescope references are hashbound by primarymanifest.','required_future_value_checks':E['worker_required_checks'],'actual_N1ell_intrinsicS2rell_native_values_HF1_not_computed':True})
 report={'id':STEM,'all_passed':True,'native_C16_cache_source_gate_accepted':True,'actual_N1ell_or_S2_or_HF1_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()
