#!/usr/bin/env python3
"""Literal accepted-alpha support/typed dependency proof only; no contraction."""
from pathlib import Path
import itertools,json,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_ALPHA_ACTIVE_SUPPORT_SOURCE_20261105_v63'
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 a='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_20261105_v63/';r=rd(a+'report.json');assert r['all_passed'] and pins[a+'report.json']=='7b58e5bfaa3107eaa672535eee86c44f79d300bd6af2d6b1968bce3d98f81dff'
 p='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_20261105_v63_REVIEW/report.json';assert rd(p)['all_passed'] and pins[p]=='2b9a5ee58eda108ba19df265442cc226194afde86d84012af2a242d081184b5d'
 d=rd(a+'actual-two-full918-alpha-covectors-and836-G3-identities.json');support=[[i for i,x in row] for row in d['alpha_two_by918_sparse_rows']];assert list(map(len,support))==[34,52]
 b='reports/AMBIENT_B18_ELL2_T2_FULL918_RIGHTSOLVE_SOURCE_20261105_v63/';rr=rd(b+'report.json');assert rr['all_passed'] and pins[b+'report.json']=='8741d10b4cc1253c086e7f3c89163cc3f87caac2374a62f17402dcbcdf2fb494'
 e=rd(b+'actual-factored918-rightsolve-source-circuit-and-all-generator-identities.json')
 for stem,rep in [(a,r),(b,rr)]:
  for p,h in rep['evidence_sha256'].items():assert sha(ROOT/stem/p)==h
 va=[0,1,4,5,7,8];m9=e['quadratic_pair_order44'];m6=[x for x in itertools.combinations_with_replacement(range(6),2) if x!=(5,5)];allowed={44*t0+m9.index(sorted((va[i],va[j]))) for t0 in range(17) for i,j in m6};assert len(allowed)==340
 seed=[]
 for row in support:
  assert all(i//9<42 and i%9 in va for i in row);nodes={e['full_source_outputs918'][i] for i in row};seed.append(sorted(nodes))
  assert all(e['full_source_outputs918'][i]==e['active_source_outputs252'][6*(i//9)+va.index(i%9)] for i in row)
 program=e['scalar_linear_nodes'];closures=[];inputsets=[]
 for s in seed:
  seen=set();todo=set(s)
  while todo:
   i=todo.pop()
   if i in seen:continue
   seen.add(i);n=program[i]
   if n[0]=='scale':kids=[n[2]]
   elif n[0]=='sum':kids=n[1]
   elif n[0]=='unit28_component':kids=e['active_unit_block28']['input_nodes28']
   else:assert n[0] in ('constant','input');kids=[]
   for k in kids:assert k<i;todo.add(k)
  inputs=sorted({program[i][1] for i in seen if program[i][0]=='input'});assert set(inputs)<=allowed
  closures.append(sorted(seen));inputsets.append(inputs)
 proof={'all86_alpha_slots_onlyactive42_and_activevariables':True,'alpha_support_slots34_52':support,'seed_node_ids':seed,'all_reachable_scalar_nodes':closures,'all_possible_rho_input_slots_by_literal_DAG_dependency':inputsets,'all340_allowed_top17_activequadratic_output_slots':sorted(allowed),'nonactive_tensor_sources_never_seeded':True,'rho_outside340_structurally_zero':'ActualalphaRonlyseedsactive252Routputs. Thoseoutputsdependonlyonactive340topquadratics,includingthefixed28unitblock. No306/102/54/24tensoroutputcontributes. Reverseallunitblockinputsaretakenwithoutassumingsparseinverse.','no_signed_alpha_times_R_or_rho_value_evaluated':True,'future_scope':'Freshtranspose/28unitheight/wholejobGOandall918rhoFidentitycheckstillrequired. ThisdoesnotproveanyRHSsourcefunctionaldescendsraw33oranyphysicalhighercallbackavailability.'}
 ep='literal-alpha-active-support-and-pruned-rightsolve-dependency.json';(o/ep).write_text(json.dumps(proof,separators=(',',':'))+'\n');owned=count([r,d,support,rr,e,allowed,seed,closures,inputsets,proof,pins])+6000;assert owned<60000
 report={'id':STEM,'all_passed':True,'SOURCE_metadata_only_no_signed_contraction':True,'alpha_supports':[34,52],'reachable_input_counts':list(map(len,inputsets)),'structural_rho_support_within_top_active340':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};(o/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'inputcounts':list(map(len,inputsets))}))
if __name__=='__main__':main()
