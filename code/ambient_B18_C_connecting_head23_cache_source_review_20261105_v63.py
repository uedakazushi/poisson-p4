#!/usr/bin/env python3
"""Independent single-root compilation/majorants/unrun worker audit."""
import ast,hashlib,json,math,pathlib,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_CONNECTING_HEAD23_CACHE_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
MON=[(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)];IDX={m:i for i,m in enumerate(MON)}
HEIGHT=0
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def mark(x):
 global HEIGHT
 h=max(abs(x[0].numerator).bit_length(),x[0].denominator.bit_length()) if isinstance(x[0],Q) else (x[0]*x[1]).bit_length()+1
 HEIGHT=max(HEIGHT,h);assert HEIGHT<=103;return x
def bound(q):return mark(((abs(q.numerator)+q.denominator-1)//q.denominator,q.denominator))
def add(a,b):
 if isinstance(a[0],Q) and isinstance(b[0],Q):return mark((a[0]+b[0],1))
 if isinstance(a[0],Q):a=bound(a[0])
 if isinstance(b[0],Q):b=bound(b[0])
 return mark((a[0]+b[0],math.lcm(a[1],b[1])))
def mul(a,b):
 if isinstance(a[0],Q) and isinstance(b[0],Q):return mark((a[0]*b[0],1))
 if isinstance(a[0],Q):a=bound(a[0])
 if isinstance(b[0],Q):b=bound(b[0])
 if not a[0] or not b[0]:return (0,1)
 return mark((a[0]*b[0],a[1]*b[1]))
def zero():return [(Q(0),1)]+[(0,1)]*5
def conv(a,b):
 v=zero()
 for i,x in enumerate(MON):
  for j,y in enumerate(MON):
   k=IDX.get((x[0]+y[0],x[1]+y[1]))
   if k is not None:v[k]=add(v[k],mul(a[i],b[j]))
 return v
def worker_count(x):
 if isinstance(x,(bool,int,float,Q)):return 1
 if isinstance(x,str):return int(bool(__import__('re').fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple,set)):return sum(worker_count(v) for v in x)
 if isinstance(x,dict):return sum(worker_count(k)+worker_count(v) for k,v in x.items())
 return 0
def main():
 global HEIGHT
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None,parse=True):
  nonlocal RAW
  p=pathlib.Path(p);v=sha(p)
  if h:assert v==h,(str(p),v,h)
  k=str(p.relative_to(ROOT))
  if k not in pins:RAW+=p.stat().st_size
  pins[k]=v
  return json.loads(p.read_text()) if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json','e9f6a9af2398b26121b8ace947b73e0693c3069022e3a8aaf2aebe61b0abbaad')
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 pin(ROOT/'code/ambient_B18_C_connecting_head23_cache_source_20261105_v63.py',pr['code_sha256'],False)
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'],False)
 C=pin(ROOT/'reports'/PRI/'compiled-single-bB63-root-and-provided-joint105-cache.json',pr['evidence_sha256']['compiled-single-bB63-root-and-provided-joint105-cache.json'])
 P=pin(ROOT/'reports'/PRI/'single-root-head23-coherent-positive-height-and-whole-budget.json',pr['evidence_sha256']['single-root-head23-coherent-positive-height-and-whole-budget.json'])
 headpath=ROOT/'reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63'
 hr=pin(headpath/'report.json','05bc2794abc9f1d9bb7fe671516fcfc41fc8818ae0fa5d59a5b0d1ff6ee394e4')
 H=pin(headpath/'full-coherent-head23-two-socle-cache-polynomials-and-low-two-row-source.json','25cb781012de5c2d77bd9a0bf262202790458e759224810eed623aaabef28529')
 hp=pin(ROOT/'reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63_REVIEW/report.json','25cb1e164eb54b85f6835a51781df9cfa0f5adae5279193e22c1a48e34d2b281');assert hp['all_passed']
 # Original immutable closure; no author compiler imports.
 orig=ROOT/next(p for p in pr['input_sha256'] if p.endswith('compiled-original164-DAG-and-provided-pure-ell-inputs.json'))
 records=field(orig,'operator_records_originalID_op_args');index={r[0]:r for r in records};wanted=set();todo=[38648]
 while todo:
  i=todo.pop()
  if i in wanted:continue
  wanted.add(i);r=index[i]
  if r[1] not in ['c','input']:assert all(j<i for j in r[2:]);todo.extend(r[2:])
 selected=[index[i] for i in sorted(wanted)]
 assert len(selected)==2194 and selected==C['operator_records_originalID_op_args']
 source_peak1=count(records)+len(index)+len(wanted)+count([C,P,H,pr,hr,hp])+15000
 assert source_peak1<150000
 del records,index,wanted
 joint=ROOT/next(p for p in pr['input_sha256'] if p.endswith('compiled-six-even-roots-and-provided-joint105-plane3-cache.json'))
 given=field(joint,'provided105_originalID_r_s_coefficients10')
 assert given==C['provided105_originalID_r_s_coefficients10'] and len(given)==105
 assert C['B6_monomials']==[list(m) for m in MON]
 assert field(joint,'provided_monomials10')==[[0,0],[1,0],[0,1],[2,0],[1,1],[0,2],[3,0],[2,1],[1,2],[0,3]]
 assert field(joint,'all630_known105_lower_overlap_checked')
 # Expected lowerroot labels are exactprior readbacks, not guesses.
 old=ROOT/next(p for p in pr['input_sha256'] if p.endswith('actual-source105-ell2-raw33-values.json'))
 priorroots=field(old,'all_requested164_root_values');prior_peak=count(priorroots)
 oldroot=dict(dict(priorroots)[38648]);del priorroots
 ell=ROOT/next(p for p in pr['input_sha256'] if p.endswith('all164-cached-pure-ell-readbacks-and-physical12-second-values.json'))
 rootell=dict(field(ell,'all164_actual_pure_ell_root_values3'))[38648]
 lower=[oldroot.get(0,'0'),oldroot.get(8,'0'),rootell[1],oldroot.get(25,'0'),rootell[2]]
 assert lower==C['expected_root_lower_B6_order_1_r_ell_rell_ell2'] and Q(lower[0])==Q(rootell[0])
 # Pure positive-majorant propagation only. No signed positiveorder root.
 bounds={};rawparse=0
 for i,arr in given:
  q=list(map(Q,arr));rawparse=max(rawparse,max(max(abs(x.numerator).bit_length(),x.denominator.bit_length()) for x in q))
  v=zero();v[0]=(q[0],1);v[1]=add(bound(q[1]),bound(q[2]));v[2]=mul((Q(1,5),1),bound(q[2]))
  v[3]=mul((Q(1,5),1),add(bound(q[4]),mul((Q(2),1),bound(q[5]))))
  v[4]=mul((Q(1,25),1),bound(q[5]));v[5]=mul((Q(1,25),1),add(bound(q[8]),mul((Q(3),1),bound(q[9]))));bounds[i]=v
 uses=Counter()
 for r in selected:
  if r[1] not in ['c','input']:uses.update(r[2:])
 uses[38648]+=1;live={};peak=0;hashes=hashlib.sha256();units=[]
 for r in selected:
  i,op=r[:2]
  if op=='input':v=bounds[i]
  elif op=='c':assert r[3];v=zero();v[0]=(Q(r[2],r[3]),1)
  elif op=='a':v=[add(a,b) for a,b in zip(live[r[2]],live[r[3]])]
  elif op=='n':v=list(live[r[2]]);v[0]=(-v[0][0],1)
  elif op=='m':v=conv(live[r[2]],live[r[3]])
  else:
   assert op=='i';x=live[r[2]];assert x[0][0]==1;v=zero();v[0]=(Q(1),1)
   for j,m in enumerate(MON[1:],1):
    for k,n in enumerate(MON[1:],1):
     t=IDX.get((m[0]-n[0],m[1]-n[1]))
     if t is not None and t<j:v[j]=add(v[j],mul(x[k],v[t]))
   conv(x,v);conv(v,x);units.append([i,r[2],'1'])
  live[i]=v;peak=max(peak,len(live));hashes.update((json.dumps([i,[[str(a),str(b)] for a,b in v]],separators=(',',':'))+'\n').encode())
  if op not in ['input','c']:
   for j in r[2:]:
    uses[j]-=1
    if not uses[j]:del live[j]
  if not uses[i]:del live[i]
 assert peak==120 and units==[[21179,21178,'1']]
 assert hashes.hexdigest()==P['all_node_positive_majorant_sha256']
 rootbound=live[38648][5];assert [str(x) for x in rootbound]==P['root_lambda_majorant_U_D']
 assert HEIGHT==P['root_stored_and_diagnostic_height_majorant']==41 and rawparse<=41
 leaves={tuple(k):bound(Q(q)) for k,q,locator in H['all_input_values_already_cached']};leaves['NEW_bB63_lambda12',]=rootbound
 polynomials=H['head1_genuine_high_A_B_UNPAIRED']+H['head1_canonical_socle_coordinates_rg_k2_UNPAIRED']+[p for pair in H['five_head_low_rows_rH_ellH_UNPAIRED'] for p in pair]
 outputs=[]
 for poly in polynomials:
  value=(0,1)
  for labels,weight in poly:
   term=(Q(weight),1)
   for label in labels:term=mul(term,leaves[tuple(label)])
   value=add(value,term)
  outputs.append(value)
 witness=mul(outputs[5],outputs[3]);assert HEIGHT==103 and 2*HEIGHT+1==207
 assert [[str(a),str(b)] for a,b in outputs]==P['head_output_polynomial_majorants_U_D']
 assert [str(a) for a in witness]==P['simpletarget_minor_majorant_U_D']
 # Exact code/custody snapshot and strict before-operation behavior.
 custodypath=ROOT/'reports/AMBIENT_B18_C_CONNECTING_HEAD23_VALUE_CUSTODY_20261105_v63'
 custody=pin(custodypath/'prepared-worker-code-source-custody.json','c54d5fa005d9f17da6e4e5c35929ebe980075f27387e0fe7f01c7b94ca9034ed')
 cr=pin(custodypath/'report.json');assert cr['all_passed'] and cr['custody_sha256']==pins[str((custodypath/'prepared-worker-code-source-custody.json').relative_to(ROOT))]
 for p,h in custody['input_sha256'].items():pin(ROOT/p,h,False)
 worker=ROOT/'code/ambient_B18_C_connecting_head23_cached_values_20261105_v63.py';pin(worker,'ee08ceaeef6fa355f770bb3a098464e8d2ed56c484d72ae0988c7d04a3d31665',False)
 assert custody['prepared_worker_sha256']==pins[str(worker.relative_to(ROOT))]
 code=worker.read_text();tree=ast.parse(code);functions={n.name:ast.get_source_segment(code,n) for n in tree.body if isinstance(n,ast.FunctionDef)}
 for k,expr in [('add','a+b'),('mul','a*b')]:assert functions[k].index('before(')<functions[k].index('return check('+expr+')')
 assert "h>PRELIMIT or h>2048" in functions['before']
 assert 'before(max(abs(a).bit_length(),b.bit_length()))' in functions['parse']
 assert "PRELIMIT=83" in functions['main'] and "PRELIMIT=207" in functions['main']
 assert "both_products_identity" in functions['main'] and 'conv(x,v)' in functions['inv'] and 'conv(v,x)' in functions['inv']
 assert "save('phase1-single-root-existing-cache-readback.json'" in functions['main']
 assert functions['main'].index("save('phase1-single-root-existing-cache-readback.json'")<functions['main'].index('del given,cache,uses,conversion,result')<functions['main'].index('PRELIMIT=207')
 assert "sha(Path(__file__))in peer['input_sha256'].values()" in functions['main']
 assert "assert owned_peak<=40000" in functions['main'] and 'RAW<=2000000' in functions['main']
 assert "STOP_PRESERVE_NO_AUTOMATIC_RETRY" in functions['main']
 # Postphase references all14 literal source expressions exactly; no term
 # may be pruned by its expectedvalue, and all outputs/factors are saved.
 assert 'for mon,w in poly:' in functions['main'] and 'terms.append([mon,w,term])' in functions['main']
 assert "witness=mul(low[0][1],k2)" in functions['main']
 assert custody['whole_positive_owned_upper']==36102 and custody['static_owned_cap']==40000
 assert custody['whole_fresh_RAW_positive_upper']==683553 and custody['fresh_RAW_cap']==2000000
 assert custody['guard_inverse_fixture_result'].startswith('PURE_GUARD_AND_B6_INVERSE_FIXTURES_PASS')
 assert 'usage:' in custody['help_result'] and custody['actual_value_folder_absent']
 assert not (ROOT/'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63').exists()
 # Originalchecker's own positive-onlySOURCE allocation, separatefrom
 # the MUCHsmaller future worker's input/lastuse/copy layout.
 own=max(source_peak1,count([pr,hr,hp,C,P,H,selected,given,bounds,live,uses,leaves,outputs,custody,cr,pins])+prior_peak+25000)
 assert own<150000 and RAW+1000000<32<<20
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'complete_original2194_root38648_closure_verified':True,'all_given105x10_source_values_exact_previous_joint_table':True,
  'existing630_lower_provenance_and5_expectedroot_lower_sources_pinned':True,
  'all_node_positive_majorant_stream_exact':True,'root_H41_pre83_andpost_H103_pre207':True,
  'all_14_output_and_minor_positive_majorants_exact':True,'lastuse_peak120_andoneunitcenter1':True,
  'strict_worker_parse_add_mul_conv_inverse_diagnostics_prechecks':True,
  'phase1_saved_before_postphase_cache_release_no_auto_retry':True,
  'fresh_coherent_source_math_peer_and_cache_code_peer_both_required':True,
  'worker_only_readies_original_cached_composition_not_new_physics':'root38648 isphysicalbB63, SAMEoriginalDAG/providedcommon105 branch ands=r+ell/5,ordinaryr*ell2. Onlyitsonepreexistingthirdcoefficient andthe14sourceexpressions/minor may be evaluated underparentONEGO.',
  'no_actual_root_or_edge_coefficient_evaluated':True,
  'future_worker_whole_positive_owned':36102,'static':40000,'hard':60000,
  'future_RAW_upper':683553,'future_RAW_hard':2000000,
  'own_SOURCE_positive_only_owned':own,'own_SOURCE_RAW':RAW,
  'parent_conditional_ONE_readiness':True}
 ep=out/'independent-single-root-closure-majorants-unrun-code-and-whole-custody.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'original_cache_code_and_whole_custody_PASS':True,'parent_conditional_ONE_readiness':True,
  'actual_cache_root_or_head_values_computed':False,'seconds':time.monotonic()-started}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'own_SOURCE':own,'seconds':report['seconds']}))
if __name__=='__main__':main()
