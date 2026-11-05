#!/usr/bin/env python3
"""Unrun full tauAB worker/input/positive-height and whole-life gate review."""
import ast,hashlib,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
import ambient_B18_C_full_eta4_cache_source_review_20261105_v63 as helper
import ambient_B18_C4_HB_residual_source_review_20261105_v63 as poly
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
SOURCE='2debebb789b978069fc3dd5863cb5ff5f834366add91716c34065209f3dd852c'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;assert out.exists() and not (out/'report.json').exists();pins={};rawbytes=0
 def pin(p,h=None):
  nonlocal rawbytes
  p=pathlib.Path(p);v=sha(p)
  if h:assert v==h,(str(p),v,h)
  k=str(p.relative_to(ROOT))
  if k not in pins:rawbytes+=p.stat().st_size
  pins[k]=v
  return json.loads(p.read_text()) if p.suffix=='.json' else p.read_text()
 mathproof=pin(out/'independent-full-tauAB-scalar-symbol-polynomial-proof.json',
  '8cc36264674b8bb2f33541a5a2cd9f6457d89f798cffeed2699793edfa4a1480')
 for p,h in mathproof['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 pin(ROOT/'code/ambient_B18_C_full_tauAB12_formal_review_20261105_v63.py',mathproof['code_sha256'])
 s=pin(ROOT/'reports'/PRI/'report.json',SOURCE)
 c=pin(ROOT/'reports'/PRI/'full-tauAB12-quartic-scalar-circuit-and-source-proof.json',s['evidence_sha256']['full-tauAB12-quartic-scalar-circuit-and-source-proof.json'])
 bounds=pin(ROOT/'reports'/PRI/'all-scalar-majorants-no-value-evaluation.json',s['evidence_sha256']['all-scalar-majorants-no-value-evaluation.json'])
 gr=pin(ROOT/'reports/AMBIENT_B18_C_FULL_TAUAB12_WHOLE_GATE_20261105_v63/report.json',
   'fabe4a9b1328d38888a24cd18b090538a39b3ed0f4090e98e3c4ede43130587f')
 for p,h in gr['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 gate=pin(ROOT/'reports/AMBIENT_B18_C_FULL_TAUAB12_WHOLE_GATE_20261105_v63/whole-source-values-rank-and-C4-positive-profile.json',
  '3283a21c3896b0c6bb5cffeebee7220db97ac0521895f600d6abf8ad72afa292')
 worker=pin(ROOT/'code/ambient_B18_C_full_tauAB12_values_20261105_v63.py',
  '2e1897055596d5abe4fcdb3c2d15ab515540185910a91aeae97f513b6dbd8312')
 iface=pin(ROOT/'code/ambient_B18_C_full_tauAB12_input_interface_20261105_v63.py',
  'cd5f54310fa100f91936d8faa15ee3125b702f79668a3975cdeef2dd3b0c0acd')
 # Literal39 reader reproduced one source file at a time; no SLP output.
 groups={}
 for row in c['provided_inputs']:groups.setdefault(row['source_path'],[]).append(row)
 checks=[];peak=literal_raw=0
 for path,records in sorted(groups.items()):
  doc=pin(ROOT/path,s['input_sha256'][path]);peak=max(peak,helper.nc(doc));literal_raw+=(ROOT/path).stat().st_size
  for row in records:
   kind,*args=row['key']
   if kind=='delta':actual=dict(next(r['delta'] for r in doc['physical12_coefficients'] if r['physical_row']==args[0])).get(args[1],'0')
   elif kind=='nu_sparse33':actual=dict(next(v for i,v in doc['all_requested164_root_values'] if i==args[0])).get(args[1],'0')
   elif kind=='nu_ell2':actual=doc['actual_nuA_B_pure_ell3'][args[0]][2]
   elif kind=='lambda':actual=next(r['delta_r_ell2'] for r in doc['all3_lambda_readbacks'] if r['physical_row']==args[0])
   elif kind=='mu03':actual=next(r['delta_ell3'] for r in doc['all3_pivot_mu03_readbacks'] if r['physical_row']==args[0])
   else:raise AssertionError(kind)
   assert Q(actual)==Q(row['value']);checks.append([row['node'],row['key'],str(Q(actual)),path,row['locator']])
  del doc
 assert len(checks)==39 and len(groups)==5 and peak==3777 and literal_raw==83374
 assert checks==gate['all39_literal_inputs_match_pinned_source_locators']
 # Independently recompute every positive bound, never signed responses.
 b=[];H=Op=0;byid={r['node']:r for r in c['provided_inputs']}
 for i,node in enumerate(c['nodes']):
  if node[0] in ['c','input']:
   q=Q(node[1],node[2]) if node[0]=='c' else Q(byid[i]['value'])
   U=(abs(q.numerator)+q.denominator-1)//q.denominator;D=q.denominator;operation=max(abs(q.numerator).bit_length(),D.bit_length())
  else:
   a,d=[b[j] for j in node[1:]]
   if node[0]=='a':U,D=a[0]+d[0],math.lcm(a[1],d[1]);operation=2*max(a[2],d[2])+1
   elif node[0]=='m':U,D=a[0]*d[0],a[1]*d[1];operation=a[2]+d[2]
   else:raise AssertionError(node[0])
  h=max((U*D).bit_length(),D.bit_length());b.append((U,D,h));H=max(H,h);Op=max(Op,operation)
 assert [[str(u),str(d),h] for u,d,h in b]==bounds['bounds_U_D_height']
 assert H==86 and Op==171 and 18*H+5==1553<2048
 assert gate['circuit_before_operation_bits']==171 and gate['rank_and_factorization_before_operation_bits']==1553
 # Full4-source Cramer factors require the larger18H+5 bound:
 # det/numerators4H+1, ratios8H+2, basisproduct9H+2,
 # two-term rawadd18H+5. No magnitude lower bound on det is used.
 functions={n.name:ast.get_source_segment(worker,n) for n in ast.parse(worker).body if isinstance(n,ast.FunctionDef)}
 for name,expr in [('add','a+b'),('times','a*b')]:assert functions[name].index('check(')<functions[name].index('return note('+expr+')')
 assert 'return times(a,note(Q(b.denominator,b.numerator)))' in functions['divide']
 assert functions['divide'].index('check(')<functions['divide'].index('Q(b.denominator,b.numerator)')
 maincode=functions['main']
 assert 'iface.read39(source,c,parse)' in maincode and '954' in maincode
 assert 'for h in[SOURCE_SHA,sha(ROOT/GATE_REPORT),sha(Path(__file__))]' in maincode
 assert 'abs(v.numerator)<=U*v.denominator and D%v.denominator==0' in maincode
 assert 'check(U.bit_length()+v.denominator.bit_length())' in maincode
 assert 'for X in[v]+norm:' in maincode
 assert 'for A,B in[(block,inv),(inv,block)]' in maincode
 assert 'all(add(times(v[a],ab[0]),times(W[a],ab[1]))==X[a]for a in range(6))' in maincode
 assert "cut_positive=[(j,k)for j in range(3)for k in[1,3,4]" in maincode
 assert 'MAX_PRECHECK<=gate[\'rank_and_factorization_before_operation_bits\']' in maincode
 assert 'del doc' in iface and "source['input_sha256'][path]" in iface
 assert "record['delta']" in iface and "r['physical_row']==physical" in iface
 # A new W has arbitrary first coordinate. Recheck the generalized C4
 # minor, retaining all lower symbols instead of reusing HB1_0=6.
 v=list(map(Q,c['known_full_HB_vector']));minorchecks=[]
 for k in [1,3,4]:
  x,y,p,q,u,w=[poly.var(j) for j in range(6)]
  A=[[poly.const(v[0]),x,{},{}],[poly.const(v[k]),y,{},{}],[p,q,poly.const(1),{}],[u,w,{},poly.const(1)]]
  dp={0:poly.const(1)}
  for r in range(4):
   nxt={}
   for mask,z in dp.items():
    for col in range(4):
     if mask>>col&1:continue
     sign=(-1)**sum(mask>>j&1 for j in range(col+1,4));m=mask|1<<col
     nxt[m]=poly.add(nxt.get(m,{}),poly.scale(poly.mul(z,A[r][col]),Q(sign)))
   dp=nxt
  expected=poly.add(poly.scale(y,v[0]),poly.scale(x,-v[k]));assert dp[15]==expected
  minorchecks.append([k,[[list(m),str(q)] for m,q in sorted(expected.items())]])
 own=C4=pin(ROOT/'reports/AMBIENT_B18_C4_HB_RESIDUAL_SOURCE_20261105_v63_REVIEW/report.json',
  'ccfcb7f5326a7c47243119c5e1454b6a5cf98a55ade68cc2a625a39677c7c75f')
 assert C4['all_passed']
 future=helper.nc(s)+helper.nc(c)+helper.nc(bounds)+helper.nc(gr)+helper.nc(gate)+5000+954+39+peak+14000
 assert future<45000<60000 and gate['value_static_numeric']==45000
 assert rawbytes<2<<20 and gate['existing_RAW_bytes_including_all5_literal_input_files']+310000<1000000
 assert not (ROOT/'reports/AMBIENT_B18_C_FULL_TAUAB12_VALUES_20261105_v63').exists()
 proof={'independent_full_symbol_polynomial_proof_sha256':sha(out/'independent-full-tauAB-scalar-symbol-polynomial-proof.json'),
  'all39_inputs_five_grouped_literal_source_files_and_hashes_verified':True,
  'positive_every954_bound_reconstructed_H86_op171':True,
  'fixed_rank_inverse_and_all4_Cramer_factorback_preop_bound':1553,
  'no_determinant_magnitude_lower_bound_or_unknown_inverse':True,
  'unrun_worker_all_guards_interface_lifetimes_and_stop_checkpoint_pass':True,
  'generalized_C4_four_source_minor_checks':minorchecks,
  'whole_future_owned_with5000_peer':future,'source_formal_whole_upper':113222,
  'response_values_not_evaluated':True,'parent_conditional_ONE_ready':True,
  'scope':'known HB+tauAB12 rank; nonzero residual givesHF4=4,b2,e3or4; nonzero cut-row1/3/4 givese3; allzero subsetnegative only',
  'hard_static_value':45000,'hard_value_numeric':60000,'hard_global_numeric':150000,
  'hard_bits':2048,'hard_RAW':32<<20,'hard_seconds':60,'hard_memory_bytes':2<<30}
 ep=out/'independent-whole-tauAB-worker-interface-height-and-C4-gate.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
  'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep),'independent-full-tauAB-scalar-symbol-polynomial-proof.json':sha(out/'independent-full-tauAB-scalar-symbol-polynomial-proof.json')},
  'full_formal_source_code_interface_whole_budget_PASS':True,'response_values_not_evaluated':True,
  'parent_conditional_ONE_cache_job_ready':True,'elapsed_seconds':time.monotonic()-start}
 assert helper.nc(report)<5000
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
