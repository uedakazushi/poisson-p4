#!/usr/bin/env python3
"""Fresh HB pivot03 original-source/majorant/code/custody gate; no root values."""
import ast,hashlib,json,pathlib,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
import ambient_B18_C_full_eta4_cache_source_review_20261105_v63 as M
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
SOURCE='6cbc0f67ea1b88f6d5e14d06176f38b58d876b34afb5823a45c73bc78a4023d0'
WORKER='code/ambient_B18_C_full_HB_pivot03_values_20261105_v63.py'
WORKERSHA='54ee6475ecb38426280117d0a730964924dbee1b4b1b28ad671bcce974e94100'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;assert out.exists() and not (out/'report.json').exists();pins={};rawbytes=0
 def pin(p,h=None):
  nonlocal rawbytes
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h,(str(p),s,h)
  k=str(p.relative_to(ROOT))
  if k not in pins:rawbytes+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text()) if p.suffix=='.json' else p.read_text()
 mathproof=pin(out/'independent-full-HB-formal-six-coordinate-source-proof.json',
  '6408bb4690c53e7ae8d9b73ca481acaa206cd7a18c6c2fe7d97ac539fb5e19c3')
 for p,h in mathproof['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 pin(ROOT/'code/ambient_B18_C_full_HB_pivot03_formal_review_20261105_v63.py',mathproof['code_sha256'])
 previous=json.loads((ROOT/'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_SOURCE_20261105_v63_REVIEW/report.json').read_text())
 pin(ROOT/'code/ambient_B18_C_full_eta4_cache_source_review_20261105_v63.py',previous['code_sha256'])
 s=pin(ROOT/'reports'/PRI/'report.json',SOURCE);assert s['all_passed']
 c=pin(ROOT/'reports'/PRI/'compiled-six-pivot-even-roots-and-pure-ell4-inputs.json',s['evidence_sha256']['compiled-six-pivot-even-roots-and-pure-ell4-inputs.json'])
 prof=pin(ROOT/'reports'/PRI/'full-HB-pivot03-frame-and-whole-operational-profile.json',s['evidence_sha256']['full-HB-pivot03-frame-and-whole-operational-profile.json'])
 formal=pin(ROOT/'reports'/PRI/'full-HB6-by2-affine-source-operator-and-old27-leaves.json',s['evidence_sha256']['full-HB6-by2-affine-source-operator-and-old27-leaves.json'])
 rec=c['operator_records_originalID_op_args'];meta=c['root_pairs_physical3']
 roots=[r[k] for r in meta for k in ['aA_node','aB_node']]
 assert roots==[20114,40892,20138,40916,20153,40931]
 original=pin(ROOT/'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63/compiled-original164-DAG-and-provided-pure-ell-inputs.json',
  'cec7274d71c0f361e1d580e6001396c4374f8ffe70359665bd9789882f440399')
 omap={r[0]:r for r in original['operator_records_originalID_op_args']};wanted=set();stack=list(roots)
 while stack:
  i=stack.pop()
  if i in wanted:continue
  wanted.add(i);r=omap[i]
  if r[1]not in ['c','input']:stack.extend(r[2:])
 assert sorted(wanted)==[r[0] for r in rec] and all(omap[r[0]]==r for r in rec) and len(rec)==7095
 prior=original['provided105_originalID_pure_ell_coefficients'];eqs=original['equation105_order']
 own_source_peak=M.nc(original)+M.nc(c)+M.nc(prof)+M.nc(formal)+M.nc(mathproof)+len(omap)+len(wanted)+12000
 assert own_source_peak<150000
 del original,omap,wanted,stack
 joint=pin(ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
  '1909547a8d80ad59905518282c8d266d2c92f791c2c9c120b01abc3cd0860d44')
 assert joint['all105_equations_zero_through_degree3'] and joint['all105_equation_original_node_order']==eqs
 old=pin(ROOT/'reports/AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63/compiled-six-even-roots-and-provided-joint105-plane3-cache.json',
  '5a10b3d66333cc85d59446278f228cd867909eac47a6e0c6db4e8918670ccaa8')
 assert joint['all105_unknown_jets']==old['provided105_originalID_r_s_coefficients10']
 assert rec==old['operator_records_originalID_op_args'] and meta==old['root_pairs_physical3']
 lower=dict(prior);provided=[]
 for i,v in joint['all105_unknown_jets']:
  w=[Q(v[k])/Q(5**j) for j,k in enumerate([0,2,5,9])]
  assert w[:3]==list(map(Q,lower[i]));provided.append([i,list(map(str,w))])
 assert len(provided)==105 and provided==c['provided105_originalID_pure_ell4'] and prior==c['provided105_originalID_pure_ell3_prior']
 ell=pin(ROOT/'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
  'f213f27f92b77d734f864dbc48f49d964e47980038c3139a958608ad3647b408')
 expected=dict(ell['all164_actual_pure_ell_root_values3'])
 assert c['expected6_roots_prior_pure_ell3']==[[i,expected[i]] for i in roots]
 own_source_peak=max(own_source_peak,M.nc(c)+M.nc(prof)+M.nc(formal)+M.nc(joint)+M.nc(old)+M.nc(ell)+20000)
 assert own_source_peak<150000
 del joint,old,ell,lower
 # Absolute majorants only, using the independent closed inverse envelopes
 # and reversed convolution from our already-frozen own SOURCE checker.
 M.height=0
 phi={i:[(Q(v[0]),1)]+[M.mark((M.ceiling(Q(q)),Q(q).denominator)) for q in v[1:]] for i,v in provided}
 uses=Counter()
 for r in rec:
  if r[1]not in['c','input']:uses.update(r[2:])
 for i in roots:uses[i]+=1
 values={};peak=0;digest=hashlib.sha256();inv=[];constants=[]
 for r in rec:
  i,op=r[:2]
  if op=='input':a=phi[i]
  elif op=='c':a=M.zseries();a[0]=M.mark((Q(r[2],r[3]),1));constants.append([i,r[2],r[3]])
  elif op=='a':a=[M.sum_bound(x,y) for x,y in zip(values[r[2]],values[r[3]])]
  elif op=='n':a=list(values[r[2]]);a[0]=M.mark((-a[0][0],1))
  elif op=='m':a=M.product_series(values[r[2]],values[r[3]])
  elif op=='i':
   a=M.inverse_bounds(values[r[2]]);M.product_series(a,values[r[2]]);M.product_series(values[r[2]],a);inv.append([i,r[2],'1'])
  else:raise AssertionError(op)
  digest.update((M.enc([i,[[str(x),str(y)] for x,y in a]])+'\n').encode())
  values[i]=a;peak=max(peak,len(values))
  if op not in ['c','input']:
   for j in r[2:]:
    uses[j]-=1
    if not uses[j]:del values[j]
  if not uses[i]:del values[i]
 deltas=[]
 for row in meta:
  a,b=values[row['aA_node']],values[row['aB_node']]
  v=[M.mark((a[0][0]-b[0][0],1))]+[M.sum_bound(x,y) for x,y in zip(a[1:],b[1:])]
  deltas.append([row['physical_row'],[[str(x),str(y)] for x,y in v]])
 assert peak==162 and inv==c['unit_inverse_centers']
 assert digest.hexdigest()==prof['every_node_majorant_sha256'] and constants==prof['literal_all_n_d_constants']
 assert deltas==prof['root_delta_final_majorants'];rootheight=M.height;assert rootheight==35
 leafbounds={tuple(k):(M.ceiling(Q(q)),Q(q).denominator) for k,q in formal['known27_old_cached_leaf_values']}
 for p,vals in deltas:leafbounds[('mu03',p,-1)]=tuple(map(int,vals[3]))
 entrybounds=[]
 for row in formal['full_HB_6_by2_affine_operators']:
  rowbounds=[]
  for entry in row:
   v=(0,1)
   for k,w in reversed(entry):v=M.sum_bound(v,M.product_bound((Q(w),1),leafbounds[tuple(k)]))
   rowbounds.append(v)
  entrybounds.append(rowbounds)
 for row in entrybounds:M.sum_bound(row[1],M.product_bound((Q(10),1),row[0]))
 assert [[[str(x),str(y)] for x,y in row] for row in entrybounds]==prof['full_HB_entry_majorants']
 assert M.height==65 and 2*rootheight+1==prof['root_phase_before_operation_bound']==71
 assert 12*M.height+3==prof['post_affine_rank_inverse_operational_bound']==783<2048
 # No mathematical retry in the producer custody path correction.
 custody=pin(ROOT/'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_VALUE_CUSTODY_20261105_v63/report.json',
  '4b751cee5358efc9e945993834f294c9bba1e2ea03a9b1b330e343195047abad')
 for p,h in custody['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 for p,h in custody['evidence_sha256'].items():cw=pin(ROOT/'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_VALUE_CUSTODY_20261105_v63'/p,h)
 history=pin(ROOT/'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_CUSTODY_PATH_HISTORY_20261105_v63/path-only-failure-and-lossless-correction.json',
  '5e1c8711cb58e0bead9eef76c9e660b4b87408defedfe1e0a04bf1daa6b595f2')
 assert history['status']=='V1_CUSTODY_FILE_NOT_FOUND_BEFORE_ANY_INPUT_PARSE_OR_MATH'
 assert sha(ROOT/history['snapshot'])==history['snapshot_sha256']
 worker=pin(ROOT/WORKER,WORKERSHA);tree=ast.parse(worker)
 fs={n.name:ast.get_source_segment(worker,n) for n in tree.body if isinstance(n,ast.FunctionDef)}
 assert 'parse(str(r[2])+\'/\'+str(r[3]))' in fs['main']
 for name,expr in [('add','a+b'),('times','a*b')]:assert fs[name].index('check(')<fs[name].index('return note('+expr+')')
 assert 'add(v[n],neg(times(x[k],v[n-k])))' in fs['inverse']
 assert 'conv(cache[r[2]],v)==one==conv(v,cache[r[2]])' in fs['main']
 assert "assert MAX_PRECHECK<=profile['root_phase_before_operation_bound']" in fs['main']
 assert "phase1-six-cached-pivot-pure-ell4-values.json" in fs['main']
 assert fs['main'].index('rp.write_text')<fs['main'].index('leaves={')
 assert "POST_ROOT_ONLY_STOP_CHECKPOINT.json" in fs['main'] and 'no_root_rerun' in fs['main']
 assert "MAX_PRECHECK<=profile['post_affine_rank_inverse_operational_bound']" in fs['main']
 assert 'check(max(abs(det.numerator).bit_length(),det.denominator.bit_length()))' in fs['main']
 assert 'invdet=note(Q(det.denominator,det.numerator))' in fs['main']
 assert 'for A,B in[(block,inv),(inv,block)]' in fs['main']
 assert "assert all(row[1]==times(parse('-10'),row[0])for row in matrix)" in fs['main']
 future=M.nc(s)+M.nc(c)+M.nc(prof)+M.nc(formal)+5000+2*len(rec)+525+4*162+16+14000
 assert future==58222==cw['owned_with5000_future_peer_numeric_slots'] and future<75000
 assert prof['source_owned_numeric_peak_bound']<150000 and prof['future_owned_numeric_upper']==61001
 assert not (ROOT/'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_VALUES_20261105_v63').exists()
 assert rawbytes<2<<20
 proof={'formal_math_proof_sha256':sha(out/'independent-full-HB-formal-six-coordinate-source-proof.json'),
  'new7095_original_ops_and_exact_closure_verified':True,
  'same105_namespace315_input_and18_old_root_overlap':True,
  'fresh_closed_inverse_majorants_and_reverse_convolution_all_nodes':digest.hexdigest(),
  'lastuse_peak':peak,'root_height_bound':rootheight,'root_preop_bound':71,
  'positive_affine_entry_height_bound':65,'all_postfold_rank_inverse_preop_bound':783,
  'whole_own_source_upper':own_source_peak,'whole_value_upper_with5000_peer':future,
  'unrun_worker_literal_guard_lifetime_phase1_and_poststop_pass':True,
  'custody_V1_path_only_failure_preserved_before_any_math':True,
  'own_raw33_missing_ell2_reader_failure_preserved_frame_binding_corrected':True,
  'no_mu_or_HB_value_pairing':True,
  'parent_conditional_one_cache_job_readiness':True,
  'rank2_consequence':'HF4=4,b2,e3or4; no exact e',
  'rank1_scope':'known two-source submodule only; HF4 still4or5',
  'hard_source_slots':150000,'hard_value_slots':100000,'static_value_slots':75000,
  'hard_bits':2048,'hard_RAW_bytes':32<<20,'hard_seconds':60,'hard_memory_bytes':2<<30}
 ep=out/'independent-full-HB-original-cache-whole-guard-and-unrun-code-gate.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
  'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep),'independent-full-HB-formal-six-coordinate-source-proof.json':sha(out/'independent-full-HB-formal-six-coordinate-source-proof.json')},
  'full_source_operator_code_custody_and_whole_budget_PASS':True,
  'requested_mu_or_HB_values_evaluated':False,'parent_conditional_ONE_job_ready':True,
  'elapsed_seconds':time.monotonic()-start}
 assert M.nc(report)<5000
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
