#!/usr/bin/env python3
"""Literal raw33->plane9 homomorphism; reduce missing HEAD29 plane4 to one."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_HEAD29_CACHE_AVAILABILITY_20261105_v63'
PATHS={
 'head29':'reports/AMBIENT_B18_C_CONNECTING_HEAD29_SOURCE_20261105_v63/coherent-head29-high-operator-and-three-selected-plane4-source-gate.json',
 'head29_report':'reports/AMBIENT_B18_C_CONNECTING_HEAD29_SOURCE_20261105_v63/report.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'raw33_peer':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'E7':'reports/AMBIENT_B18_C_GENUINE_T7_CACHED_VALUES_20261105_v63/four-cached-coefficients-actual-E7-and-genuine-four-target-coordinates.json',
 'E7_peer':'reports/AMBIENT_B18_C_GENUINE_T7_CACHED_VALUES_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 st=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={'raw33':['coordinate_order9','actual33_basis','coordinate_actions9_sparse_columns','actual_source_equality'],
  'C2':['physical12_coefficients'],'E7':['all4_actual_cached_coefficients']}
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 assert d['head29_report']['all_passed']
 b=d['raw33']['actual33_basis'];assert len(b)==33
 plane=[(0,0),(1,0),(0,1),(2,0),(1,1),(3,0),(2,1),(4,0),(3,1)];index={m:i for i,m in enumerate(plane)}
 P=[({index[(r,e)]:Q(1)}if z==0 else{})for z,r,e in b]
 actions=d['raw33']['coordinate_actions9_sparse_columns'];assert len(actions)==9
 def proj(x):
  v={}
  for j,q in x:
   for k,a in P[j].items():v[k]=v.get(k,Q(0))+Q(q)*a
  return{k:q for k,q in v.items()if q}
 def plane_action(i,col):
  if i<7:return{}
  v={}
  for j,q in col.items():
   a,e=plane[j];m=(a+int(i==7),e+int(i==8));k=index.get(m)
   if k is not None:v[k]=v.get(k,Q(0))+q
  return{k:q for k,q in v.items()if q}
 dig=hashlib.sha256();checks=0
 for i in range(9):
  for j in range(33):
   assert proj(actions[i][j])==plane_action(i,P[j]),(i,j)
   checks+=1;dig.update(json.dumps([i,j,sorted(proj(actions[i][j]).items())],default=str,separators=(',',':')).encode())
 assert checks==297 and P[0]=={0:Q(1)}
 assert[P[j]for j,m in enumerate(b)if m[0]==0]==[{i:Q(1)}for i in range(9)]
 rec=next(r for r in d['C2']['physical12_coefficients']if r['physical_row']==63);delta=dict(rec['delta'])
 ir4=b.index([0,4,0]);ir3e=b.index([0,3,1]);assert(ir4,ir3e)==(31,32)
 known={('CACHED_delta63_plane4','r4'):Q(delta[ir4]),('CACHED_delta63_plane4','r3ell'):Q(delta[ir3e])}
 assert list(known.values())==[Q(-1236125,1728),Q(-91325,288)]
 replacement={('NEW_delta63_plane4','r4'):('CACHED_delta63_plane4','r4'),('NEW_delta63_plane4','r3ell'):('CACHED_delta63_plane4','r3ell')}
 polys=[]
 for p in d['head29']['four_genuine_and_canonical_high_polynomials_UNPAIRED']:
  mapped={}
  for mon,q in p:
   key=tuple(sorted(replacement.get(tuple(label),tuple(label))for label in mon));mapped[key]=mapped.get(key,Q(0))+Q(q)
  polys.append([[[list(label)for label in mon],str(q)]for mon,q in sorted(mapped.items())if q])
 missing={tuple(label)for p in polys for mon,q in p for label in mon if label[0].startswith('NEW_')}
 assert missing=={('NEW_delta63_plane4','r2ell2'),('NEW_CACHED_bB63_ell3',)}
 sensitivity=[]
 for p in polys:
  terms=[]
  for mon,q in p:
   if ['NEW_delta63_plane4','r2ell2']in mon:
    assert mon==[['NEW_delta63_plane4','r2ell2']];terms.append(Q(q))
  sensitivity.append(sum(terms,Q(0)))
 assert sensitivity==[Q(1),Q(0),Q(5,3),Q(-121,48)]
 data={
  'scope':'Literal existing-cache availability refinement ONLY. Head29sourcepolynomialsremainUNPAIRED; no selectedfourthIFTfunctional/cacheell³valuecomputed.',
  'raw33_to_actual_plane9':'z7→0,r→r,ell→ell, QQ[r,ell]/(ell²,(r,ell)^5)',
  'plane9_basis':[list(m)for m in plane],
  'actual_projection9_by33_sparse_columns':[[[i,str(q)]for i,q in col.items()]for col in P],
  'all297_coordinate_basis_intertwiners_verified':True,'all297_intertwiner_hash':dig.hexdigest(),
  'algebra_homomorphism_proof':'Raw33iscycliconunitwiththeaccepted9actualcoordinateactions. P1=1andP*x_i=plane_x_i*P onALL33basisvectors forALL9generators impliesP(f)=f(z=0,r,ell) forANYformalfunction(bothfinite). ThisissameR0homomorphism, notanassumedplainmonomialNF. Surjectivebecauseall9planeunitwordsareliteralraw33basis.',
  'ideal_source_proof':'ActualIA=(z7,gAorder3),IB=(z7+beta2+...,gBorder3). Atz=0theirproducthasplaneorder≥5 (gA*beta2) or6; henceitvanishesmodplane m5. Ell²cutisretained. The297actualactionchecksindependentlycertifythisideal/sourceframe.',
  'two_actual_existing_delta63_fourth_inputs':[[list(k),str(q)]for k,q in known.items()],
  'literal_cached_input_locations':[['C2','physical_row63','delta',31,'r4'],['C2','physical_row63','delta',32,'r3ell']],
  'actual_delta63_C2_source_roots':rec['roots'],
  'refined_four_head29_polynomials_UNPAIRED':polys,
  'remaining_callback_labels':[list(k)for k in sorted(missing)],
  'remaining_one_selected_plane4':'delta63_r²ell² only; raw33 hasell²=0anddoesNOTdeterminethiscoefficient.',
  'remaining_one_existing_cache_composition':'bB63_ell³ fromgivenjoint105third/root38648; previousB6r²0/ell³0jobdoesNOTdeterminit.',
  'coherent_delta63_r2ell2_sensitivity_A_B_rg_k2':list(map(str,sensitivity)),
  'selected_IFT4_next_SOURCE_gate':'Onlythecoherentdelta63_r²ell²functional needsafreshselectedIFT4sourceprofile. Samefixedρ=gradδ63 J105^-1 andRHS4fromgivenplaneΦ≤3. Actualρ,closure,height,ownership,RAW,preparedcode/custody/peer/ROOTONEGOremainOPEN. No fullplane4/source4jet isneededbythisrefinement.',
  'prior_frozen_head29_source_kept':'Its3unreadplane4variableswereaconservativeprofile; twoarenowboundtoliteralexistinginputslotswithoutchangingthatfrozenpacket.',
  'no_head29_response_or_original_H0_nil_claim':True,
  'selective_source_receipts':receipts,
 }
 owned=S.count(d)+S.count(data)+12000;assert owned<100000;data['whole_owned_numeric_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'literal-raw33-plane9-map-and-one-selected-fourth-head29-gate.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_RAW33_PLANE9_ALL297_TWO_CACHED_FOURTH_ONE_MISSING',
  'two_head29_plane4_inputs_already_cached':True,'remaining_selected_plane4_callbacks':1,'remaining_existing_cache03_compositions':1,
  'head29_values_or_newsource4_computed':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},'whole_numeric_upper':owned,'seconds':time.monotonic()-st}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
