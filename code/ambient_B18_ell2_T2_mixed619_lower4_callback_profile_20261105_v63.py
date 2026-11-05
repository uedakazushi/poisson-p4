#!/usr/bin/env python3
"""Literal callback/support profile only; no fold, gradient or callback value."""
import collections,json,resource,signal,time
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_CALLBACK_PROFILE_20261105_v63'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):
  inputs[p]={'path':p,'sha256':sha(ROOT/p)}
  return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 prefix='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_X2_SOURCE_20261105_v63';r=read(prefix+'/report.json');assert r['all_passed']
 gp=prefix+'/actual-four-Gamma-X2-source-rows-and-lazy619-recipe.json';g=read(gp);assert sha(ROOT/gp)==r['evidence_sha256'][Path(gp).name]
 Gamma=[[(i//8,i%8,x) for i,x in row] for row in g['Gamma4_on_source56_times8']]
 c=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations']
 m=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['old640_source_columns','old640_typed_closure_nodes'])
 B=[dict(v) for v in m['old640_source_columns']];ops=dict(m['old640_typed_closure_nodes'])
 assert len(ops)==1302 and all(i<48 for row in Gamma for i,p,x in row)
 inputnumeric=count([g,c,m]);prior=inputnumeric+3*8*1302+8000;assert prior<150000
 oldcount=[0]*4;nativecount=[0]*4;paircounts=collections.Counter();sourcecount=collections.Counter();native_rows=set()
 # Stream counters, never expand the four-by619 coefficient tensors.
 for rec in c:
  assert rec['full56_Xgamma1_equals_Xp0a1']
  for source,coef in rec['original768_combination']:
   for target,row in enumerate(Gamma):
    for i,p,a in row:
     if source>=640:nativecount[target]+=1;native_rows.add((i,p));continue
     root=B[source].get(i)
     if root is None or ops[root][0]=='c':continue
     oldcount[target]+=1;paircounts[(root,p)]+=1;sourcecount[root]+=1
 closure=set()
 def visit(k):
  if k in closure:return
  closure.add(k);n=ops[k]
  if n[0]=='sum':children=[j for j,a,b in n[1]]
  elif n[0]=='product':children=n[3]
  else:children=[]
  for j in children:visit(j)
 for root in sourcecount:visit(root)
 leaves=sorted({ops[k][1] for k in closure if ops[k][0]=='actual_C0_scalar_function'})
 nonlinear=sorted(k for k in closure if ops[k][0]=='product' and len(ops[k][3])>1)
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'literal-oldB-native-X2-callback-and-firstcarry-profile.json'
 ev={'Gamma_support4':[len(v) for v in Gamma],'all_Gamma_rows_less48_no_direct_P8_second':True,
 'old_uncombined_term_counts4':oldcount,'native_uncombined_term_counts4':nativecount,
 'old_root_pair_counts':[[list(k),n] for k,n in sorted(paircounts.items())],'old_distinct_roots':sorted(sourcecount),'reachable_typed_node_ids':sorted(closure),
 'reachable_actual_C0_scalar_labels':leaves,'reachable_nontrivial_product_nodes':nonlinear,'native_output_row_pair_slots':sorted(native_rows),
 'source_prior_numeric':prior,'actual_profile_owned':count([g,c,m,paircounts,closure,leaves])+8000,
 'firstcarry_scope':'RetainactualXp1*a1forall619. Original8directionsarecached; individualnativepivotellrequiresL0Se+LeS0fromacceptedC16,notonlyeffectiveL0Se. No firstcarry numericalvalueiscomputedhere.',
 'opaque_scope':'EverydirectC0secondandnonlinearproductdescendantremainsformal. Thisprofiledoesnotinfercancellation,descentorcacheavailabilityfromnames.',
 'native_scope':'WeightednativeX2requiresseparateS0cN0beforeL2drop; retainL0S2+LrS_e+LeS_r/mixedanalogues inactualframe. No native secondcallbackevaluated.',
 'no_signed_fold_or_gradient_or_X2_values':True}
 ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'source_only':True,'old_distinct_roots':len(sourcecount),'old_distinct_root_pairs':len(paircounts),'reachable_nodes':len(closure),'actual_C0_labels':len(leaves),'nonlinear_products':len(nonlinear),'no_signed_fold_gradient_or_values':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'oldterms':oldcount,'nativeterms':nativecount,'roots':len(sourcecount),'pairs':len(paircounts),'nodes':len(closure),'leaves':len(leaves),'nonlinear':len(nonlinear),'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
