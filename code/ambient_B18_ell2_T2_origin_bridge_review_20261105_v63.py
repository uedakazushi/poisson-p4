#!/usr/bin/env python3
"""Existing-origin bridge to the new minimal26 frame; no X1 query."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv,mp,eye
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_ORIGIN_BRIDGE_20261105_v63_REVIEW'
P={
 'origin':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'origin_peer':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63_REVIEW/report.json',
 'adapter':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/phase2-native-adapter.json',
 'S':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'frame':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json',
 'frame_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json',
 'G2_peer':'reports/AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def mu_transpose(mat):
 out=[]
 for k in range(12):
  for b in range(4):
   p={}
   for i in range(12):
    for a,c in enumerate(mat[k][i]):
     if Q(c) and not(a&b):p=add(p,{4*i+(a|b):Q(c)})
   out.append(p)
 return out
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'EXISTING27_UNITS_AND741_ORIGIN_BRIDGE','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert all(D[k]['all_passed'] for k in P if k.endswith('_peer'))
 O=D['origin'];assert O['C0_native_algebra']=='QQ[tauA,tauB]/(tauA^2,tauB^2)'
 X=[{i:Q(O['full_X0_56_by768'][i][j]) for i in range(56) if Q(O['full_X0_56_by768'][i][j])} for j in range(768)]
 free=O['free_kernel_rows32'];assert len(free)==32
 basis={};reprs={};selected=[];combos=[]
 for j,r in enumerate(O['all768_fixed_elimination_receipts_and_constant_combinations']):
  assert r['source_column']==j;v={i:X[j][k] for i,k in enumerate(free) if k in X[j]};w={j:Q(1)}
  for p,c in r['subtract']:
   assert v.get(p)==Q(c);v=add(v,scale(basis[p],-Q(c)));w=add(w,scale(reprs[p],-Q(c)))
  if r.get('zero'):
   assert not v and not mv(X,w);combos.append((j,w))
  else:
   p=r['pivot'];c=Q(r['divide']);assert c and min(v)==p and v[p]==c and p not in basis
   v=scale(v,1/c);w=scale(w,1/c);basis[p]=v;reprs[p]=w;selected.append(j)
   assert enc(v)==r['normalized']
   actual=mv(X,w);assert {i:actual.get(k,Q(0)) for i,k in enumerate(free) if actual.get(k,Q(0))}==v
  assert enc(w)==r['constant_original_source_combination']
 assert len(selected)==27 and len(combos)==741 and selected==O['selected_independent_columns']
 assert [j for j,c in combos]==[j for j in range(768) if j not in selected]
 assert all(c[j]==1 and all(k==j or k in selected for k in c) for j,c in combos)
 assert sum(j<640 for j in selected)==4 and sum(j>=640 for j in selected)==23
 M=list(map(sp,D['frame']['actual_Mnew0_columns']));assert len(M)==56 and all(not mv(M,x) for x in X)
 J=D['frame']['same_selected26_V_columns'];C=[i for i in range(56) if i not in J];assert len(C)==30
 # Unit M_J identifies kerM0 with its30 complementary coordinates, so
 # the same27 independent source columns remain units in the free part.
 assert Q(D['frame']['new_selected26_minor_determinant'])!=0
 G=[list(map(Q,p)) for p in D['adapter']['G8_C0']];S=D['S']['native_S12_by32'];S0=[{} for _ in range(128)]
 for i,j,k,sign in S['koszul28_in_first8_rows']:
  for mask,c in enumerate(G[k]):
   for b in range(4):
    if c and not(mask&b):S0[4*j+b]=add(S0[4*j+b],{4*i+(mask|b):Q(sign)*c})
 for i,j,c in S['free4_zero_source_columns']:
  for b in range(4):S0[4*j+b]=add(S0[4*j+b],{4*i+b:Q(c)})
 L=mu_transpose(D['adapter']['N']);Li=mu_transpose(D['adapter']['N_inverse'])
 assert mp(L,Li)==mp(Li,L)==eye(48)
 native=mp(L,S0);assert native==[{i:c for i,c in X[j+640].items() if i<48} for j in range(128)]
 assert all(not any(i>=48 for i in X[j]) for j in range(640,768))
 native_checks=[]
 for j,c in combos:
  old={k:a for k,a in c.items() if k<640};nc={k-640:a for k,a in c.items() if k>=640}
  assert not mv(X,old) and not mv(S0,nc)
  native_checks.append({'nonpivot_source':j,'constant_original768_combination':enc(c)})
 ep=o/'independent27units-741source-combinations-and-nativeS0zero.json'
 dump(o,ep.name,{'inputs':inputs,'same_frozen_origin_no_new_point':True,'all768_fixed_origin_operations_replayed':True,'selected27_source_columns':selected,'all741_original_constant_zero_source_combinations':native_checks,'old4_and_native23_parts_separately_zero_for_all741':True,'intrinsicS0_rebuilt_and_L0_both_inverse_products_verified':True,'new_minimal26_M0_annihilates_all768_X0':True,'new_free30_coordinate_projection_preserves_the27_units':True,'no_X1_or_T2_mu_or_rank_search_computed':True})
 report={'id':STEM,'all_passed':True,'all27_origin_units_741_zero_combinations_and_nativeS0_separately_accepted':True,'new26_frame_point_bridge_accepted':True,'new_X1_or_T2_mu_or_source_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
