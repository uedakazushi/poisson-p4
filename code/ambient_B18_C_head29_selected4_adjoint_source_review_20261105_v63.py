#!/usr/bin/env python3
"""Independent forward sparse first-jet proof of grad(delta63)=rho=0."""
import collections,hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_HEAD29_SELECTED4_ADJOINT_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PRE=0
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def check(a,b,m=False):
 global PRE
 if a and b:
  na,da=abs(a.numerator).bit_length(),a.denominator.bit_length();nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
  h=max(na+nb,da+db)if m else max(na+db,nb+da)+1;PRE=max(PRE,h);assert h<=2048
 return a*b if m else a+b
def plus(a,b):
 z=dict(a)
 for k,q in b.items():
  v=check(z.get(k,Q(0)),q)
  if v:z[k]=v
  else:z.pop(k,None)
 return z
def scale(a,q):return{k:check(v,q,True)for k,v in a.items()if v and q}
def mm(A,B):
 out=[]
 for row in A:
  z={}
  for k,q in row.items():z=plus(z,scale(B[k],q))
  out.append(z)
 return out
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None,parse=True):
  nonlocal raw
  p=pathlib.Path(p);k=str(p.relative_to(ROOT));s=sha(p)
  if h:assert h==s,(k,h,s)
  if k not in pins:raw+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text())if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json','0a5fffe39baf1f7cb40b3f603ebfe39007ea09a3ad849dcfee670e61f5a58e64')
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 pin(ROOT/'code/ambient_B18_C_head29_selected4_adjoint_source_20261105_v63.py',pr['code_sha256'],False)
 en=next(iter(pr['evidence_sha256']));ep=ROOT/'reports'/PRI/en;pin(ep,pr['evidence_sha256'][en],False)
 completion=pin(ROOT/'reports'/PRI/'whole-source-counter-completion-no-math-rerun.json','d7f00043aab9d7859529f8d9a99a193003f07fb8daf2e026e4a4c620cc769343')
 assert completion['all_passed'] and completion['primary_report_sha256']==sha(ROOT/'reports'/PRI/'report.json')
 assert sum(completion['whole_owned_components'].values())==completion['whole_owned_numeric_upper']==149030<150000
 prior=pin(ROOT/'reports/AMBIENT_B18_C_CONNECTING_HEAD29_SOURCE_20261105_v63_REVIEW/report.json','295b19f1a80e7fa512c11218501f32123d8987db2b18b9bf79b84cd527c0b4b7');assert prior['all_passed']
 original=ROOT/next(p for p in pins if p.endswith('compiled-original164-DAG-and-provided-pure-ell-inputs.json'))
 joint=ROOT/next(p for p in pins if p.endswith('joint-degree3-certificate.json'))
 records=field(original,'operator_records_originalID_op_args');provided=field(original,'provided105_originalID_pure_ell_coefficients');eq=field(original,'equation105_order')
 omap={r[0]:r for r in records};todo=[20144,40922];need=set()
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i);r=omap[i]
  if r[1]not in['c','input']:todo.extend(r[2:])
 assert len(records)==18584 and len(need)==4387
 phase0=count([records,provided,eq,pins,pr,completion])+len(omap)+len(need)+12000
 assert phase0<150000
 selected=[omap[i]for i in sorted(need)];del records,omap
 U=field(joint,'all105_unknown_original_node_order');F=field(joint,'all105_equation_original_node_order')
 jets=field(joint,'all105_unknown_jets');given={i:Q(v[0])for i,v in jets};del jets
 assert len(U)==len(F)==105 and F==eq and field(joint,'CA_F9_targets')==['0']*7+['r','s']
 assert given=={i:Q(v[0])for i,v in provided};ui={i:j for j,i in enumerate(U)};assert set(given)==set(U)
 uses=collections.Counter([20144,40922])
 for r in selected:
  if r[1]not in['c','input']:uses.update(r[2:])
 # At most128 live nodes, each with point plus105 formal linear coefficients.
 # We reserve two slots per derivative key/value even before cancellations.
 static=count([selected,provided,eq,U,F,given,ui,pins,pr,completion,prior])+len(selected)+len(need)+2*len(uses)+16000
 forward_positive=static+128*(1+2*105)+4*(1+2*105)
 assert forward_positive<150000 and raw+1000000<32<<20
 save(out/'positive-before-independent-forward-first-jets.json',{'phase0_numeric_upper':phase0,'forward_numeric_upper':forward_positive,'live_nodes_prior':128,'gradient_dimension':105,'hard_bits':2048,'no_fourth_coefficient_or_RHS_evaluation':True,'original_unused_records_released_before_AD':True})
 cache={};peaklive=0;peakslots=0;units=[];roots={};finger=hashlib.sha256()
 for rec in selected:
  i,op=rec[:2]
  if op=='input':val=(given[i],{ui[i]:Q(1)})
  elif op=='c':
   PRE0=max(abs(rec[2]).bit_length(),abs(rec[3]).bit_length());assert PRE0<=2048;val=(Q(rec[2],rec[3]),{})
  elif op=='a':
   a,da=cache[rec[2]];b,db=cache[rec[3]];val=(check(a,b),plus(da,db))
  elif op=='n':
   a,da=cache[rec[2]];val=(-a,scale(da,Q(-1)))
  elif op=='m':
   a,da=cache[rec[2]];b,db=cache[rec[3]];val=(check(a,b,True),plus(scale(da,b),scale(db,a)))
  elif op=='i':
   a,da=cache[rec[2]];assert a==1;val=(Q(1),scale(da,Q(-1)));units.append([i,rec[2],'1'])
  else:raise AssertionError(op)
  cache[i]=val;peaklive=max(peaklive,len(cache));peakslots=max(peakslots,sum(1+2*len(v[1])for v in cache.values()))
  assert len(cache)<=128 and static+peakslots+4*(1+2*105)<forward_positive
  if i in[20144,40922]:roots[i]=(val[0],dict(val[1]));uses[i]-=1
  if op not in['c','input']:
   for j in rec[2:]:
    uses[j]-=1
    if not uses[j]:cache.pop(j,None)
  if not uses[i]:cache.pop(i,None)
  finger.update(json.dumps([i,str(val[0]),[[k,str(q)]for k,q in sorted(val[1].items())]],separators=(',',':')).encode())
 grad=plus(roots[20144][1],scale(roots[40922][1],Q(-1)));assert not grad and roots[20144][0]==roots[40922][0]
 assert not field(ep,'origin_gradient_g_wrt105_sparse') and not field(ep,'fixed_rho_equation_order_sparse')
 assert units==field(ep,'origin_unit_inverse_checks')
 # Release original records, jets and last-use controls before parsing either
 # stored105 inverse. No dense inverse is computed; both fixed products checked.
 del selected,provided,need,uses,cache,roots
 matrices=[];matrix_peak=0
 for key in['fixed105_Jacobian','fixed105_Jacobian_inverse']:
  dense=field(joint,key);assert len(dense)==105 and all(len(r)==105 for r in dense)
  sp=[{j:Q(q)for j,q in enumerate(row)if Q(q)}for row in dense]
  matrix_peak=max(matrix_peak,count([dense,sp,matrices,U,F,given,ui,pins,pr,completion])+20000)
  assert matrix_peak<150000;matrices.append(sp);del dense
 J,Inv=matrices;identity=[{i:Q(1)}for i in range(105)]
 assert mm(Inv,J)==identity and mm(J,Inv)==identity
 assert mm([grad],Inv)==[{}] and mm([{}],J)==[{}]
 assert [20144,40922]==field(ep,'actual_g_original_root_pair')==field(ep,'future_functional_roots')
 assert field(ep,'future_stream_lastuse_nodes')==128
 assert len(field(ep,'future_needed_original_records'))==4387
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'independent_forward_sparse_firstjet_algebra_NOT_author_reverseAD':True,'all4387_original_scalar_nodes_replayed_at_given105_origin':True,
  'original_delta63_roots':[20144,40922],'same105_origin_equation_and_CA_F9_order_verified':True,
  'all105_origin_gradient_components_zero':True,'rho_grad_times_stored_inverse_zero':True,'both105_stored_inverse_products_verified':True,
  'selected_fourth_dependence_proof':'For any two implicit lifts with identical Phi<=3, their difference begins in degree4. Taylor expansion of the same regular rational scalar g gives (g(Phi)-g(Phi_trunc))_4=(D_Phi g)(Phi0)*Phi4. Every one of105 gradient components is exactly zero, so delta63_r2ell2 equals its lower-input-only value. This does not assert truePhi4=0 and does not evaluate that lower-input value.',
  'no_RHS4_or_selectedfourth_or_bB63ell3_or_head29_value_evaluated':True,
  'source_record_alias_and_completion_scope':'Producer main123655 was not whole; immutable-record alias convention plus separately charged controllers/views/temporaries gives149030. Primary bytes remain unchanged. Reviewer prunes unused original records before forward AD and releases AD buffers before stored matrix checks.',
  'producer_whole_countercompletion149030_passed':True,'review_phase0_upper':phase0,'review_forward_prior':forward_positive,'review_matrix_prior':matrix_peak,
  'actual_firstjet_live_nodes':peaklive,'actual_firstjet_numeric_map_slots':peakslots,'actual_preoperation_bits':PRE,'RAW':raw,
  'forward_firstjet_fingerprint':finger.hexdigest(),'higher_response_CODE_budget_custody_ROOT_GO_not_supplied':True,
  'no_originalH0_Kbar_exponent_nil_claim':True}
 evidence=out/'independent-forward105-firstjets-zero-adjoint-and-no-fourth-IFT-dependence.json';save(evidence,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_only':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{evidence.name:sha(evidence)},'all105_gradient_and_rho_zero':True,'selected_fourth_response_computed':False,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'forward_prior':forward_positive,'prebits':PRE,'seconds':report['seconds']}))
if __name__=='__main__':main()
