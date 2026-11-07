#!/usr/bin/env python3
"""Independent final20 SOURCE/code/lifetime/custody; no response computation."""
import argparse,ast,hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER20_WHOLE_V4_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;cbits=0
 def read(p):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size;assert raw<32*1024*1024;return json.loads(q.read_text())
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 def qck(q):
  nonlocal cbits
  cbits=max(cbits,abs(q.numerator).bit_length(),q.denominator.bit_length());assert cbits<=32;return q
 gp='reports/AMBIENT_B18_ELL2_T2_UPPER20_WHOLE_INPUT_GATE_V4_20261105_v63/'
 gr=read(gp+'report.json');assert gr['all_passed'] and pins[gp+'report.json']=='cbdbe69c3a701ab3809d6ac5220adb3abc5e0eef04b3ba37c0f7c373b3e73de4'
 d=read(gr['value_input_path']);assert pins[gr['value_input_path']]==gr['value_input_sha256']=='dd8c9d761d6954cee49150166e11fa76c50071702f6450857c1523decf7e9342'
 cp='reports/AMBIENT_B18_ELL2_T2_UPPER20_VALUES_V4_PREPARED_CUSTODY_20261105_v63/report.json';cu=read(cp);assert cu['all_passed'] and pins[cp]=='49dacd5da40887594c86a1beb480991fd5fe8bd2910969c27f4d51de3d444d9c'
 op='reports/AMBIENT_B18_ELL2_T2_UPPER20_PRUNED_SOURCE_20261105_v63/';orr=read(op+'report.json');os=read(op+'pruned20-exceptional-OLD-FIRST-outer-source-and605-coverage.json')
 np='reports/AMBIENT_B18_NATIVE20_PRUNED_SOURCE_20261105_v63/';nr=read(np+'report.json');ns=read(np+'pruned-native20-formal16-toactual15-and-explicit10zero-kernel-source.json')
 o=d['OLD'];n=d['NATIVE'];assert o==os and o['exceptional_cases20']==ns['exceptional_cases20']
 cases=o['exceptional_cases20'];assert len(cases)==20
 for path,h in gr['input_sha256'].items():assert sha(ROOT/path)==h;pins[path]=h
 assert d['committed_prior']==354 and d['temporary_prior']==1063 and d['hard_guards']==[384,1153]
 expectedprimitive={k:v for k,v in ns['native_primitive_known_cache'].items() if k not in ['Gamma_two56','D_two_by10_by56']};assert n['primitive']==expectedprimitive
 assert ns['native_primitive_known_cache']['Gamma_two56']==o['Gamma_two'] and ns['native_primitive_known_cache']['D_two_by10_by56']==o['D_two_by10_by56']
 assert n['raw_program']==ns['raw_auxiliary_program'] and n['restrictions']==ns['actual_kernel_restrictions'] and n['raw_point_S0_columns48']==ns['raw_point_S0_columns48']
 defs={r['kernel105']:r['original_native_source_combination'] for r in ns['actual_kernel_restrictions']+ns['explicit_zero_kernel_certificates']}
 fines={r['kernel105']:r['two_channels_Gamma_D_sum_DU'] for r in ns['actual_kernel_fine_Gamma_D_sum_DU']+ns['explicit_zero_kernel_certificates']}
 assert {r['kernel105']:r['two_channels_Gamma_D_sum_DU'] for r in n['kernel_fine']}==fines
 zids=[r['kernel105'] for r in ns['explicit_zero_kernel_certificates']];assert n['source_zero_kernel_ids']==zids and len(zids)==10
 assert all(int(q[1])==0 for k in zids for row in fines[k] for q in row)
 coords={r['case619']:r for r in ns['all20_constant_native_coordinate_recipes']};assert dict(n['case_coordinates'])=={j:r['full_kernel105_coordinate_pairs'] for j,r in coords.items()}
 assert sum(len(r['kernel105_coordinate_pairs']) for r in coords.values())==28
 for j,r in coords.items():assert sorted(r['kernel105_coordinate_pairs']+r['explicit_proven_zero_kernel_terms'])==r['full_kernel105_coordinate_pairs']
 recipes=dict(o['corrected20_source_recipes']);assert recipes[617]==['-20','0'] and all(v==['0','0'] for j,v in recipes.items() if j!=617)
 c0=dict(o['original_source_allc0_raw21']);assert set(c0)==set(cases)|{563} and o['required_raw_cases21']==sorted(c0)
 before=count([gr,d,cu,orr,os,nr,ns,defs,fines,coords,c0,pins])+20000;assert before<150000
 dump('before-constant-source-and-lifetime-checks.json',{'no_native_or20_responses':True,'constant_source_prior32_temp97':True,'higher354_temp1063_inherited':True,'source_owned_before':before})
 def put(v,i,q):
  z=qck(v.get(i,F(0))+q)
  if z:v[i]=z
  else:v.pop(i,None)
 for j in cases:
  lhs={i-640:F(q) for i,q in c0[j] if i>=640}
  for base,q in zip([563,609],recipes[j]):
   if F(q):
    for i,a in c0[base]:
     if i>=640:put(lhs,i-640,qck(-F(q)*F(a)))
  rhs={}
  for k,q in coords[j]['full_kernel105_coordinate_pairs']:
   for i,a in defs[k]:put(rhs,i,qck(F(q)*F(a)))
  assert lhs==rhs
 coverage=o['all605_coverage'];zero=set(coverage['SOURCE_positive_zero584']);assert len(zero)==584 and coverage['independently_actual_zero554']==554
 assert len(zero|{554}|set(cases))==605 and not(zero&set(cases)) and 554 not in zero|set(cases)
 opaque={579,580,593,594,597,598,605,606,611,613,615,616};assert zero|{554}|set(cases)==set(range(619))-opaque-{563,609}
 # Boolean SOURCE shapes only: no signed jets or convolutions are computed.
 basis=[tuple(v) for v in o['ordinary_B21_basis']];bi={v:i for i,v in enumerate(basis)};tab=[[bi.get(tuple(x+y for x,y in zip(a,b))) for b in basis] for a in basis]
 def prod(a,b):return {tab[i][j] for i in a for j in b if tab[i][j] is not None}
 leaves={i:{k for k,q in enumerate(v) if F(q)} for i,mode,v in o['boundary_cache']};C=dict(leaves);T={}
 for i,r in o['literal_C0_program']:
  if r[0]=='c':v={0} if r[1] else set()
  elif r[0]=='n':v=C[r[1]]
  elif r[0]=='a':v=C[r[1]]|C[r[2]]
  else:assert r[0]=='m';v=prod(C[r[1]],C[r[2]])
  C[i]=v
 for i,r in o['typed_program']:
  if r[0]=='actual_C0_scalar_function':v=leaves[r[1]]
  elif r[0]=='sum':v=set().union(*(T[k] for k,a,b in r[1] if a)) if r[1] else set()
  else:
   assert r[0]=='product';v={0} if r[1] else set()
   for k in r[3]:v=prod(v,T[k])
  T[i]=v
 L0={j:{i for i,q in v if F(q)} for j,v in o['native_L0_selected_columns']};Le={j:{i for i,q in v if F(q)} for j,v in o['native_Le_selected_columns']};first={}
 for rec in o['FIRST_column_recipes']:
  r=rec['first_column_recipe']
  if r['kind']=='literal56':v={i for i,q in r['column'] if F(q)}
  else:
   assert r['kind']=='native_ell_SUM';v=set()
   for label,mat in [('Se',L0),('S0',Le)]:
    for j,q in r[label]:
     if F(q):v|=mat[j]
  first[rec['pivot27'],rec['variable9']]=v
 lc=sum(map(len,leaves.values()));cc=sum(len(C[i]) for i,r in o['literal_C0_program']);tc=sum(map(len,T.values()));fc=sum(map(len,first.values()));assert (cc,tc,fc)==(44,223,481)
 factors=sum((slot%9==a)+(slot%9==b and a!=b) for j,v in o['original_full_a1_raw21'] for slot,q in v for a,b in o['pair_order10']);assert factors==871
 ledger=d['derived_workspace_ledger'];derived=sum(ledger.values());assert derived==16549
 # Exact tuple-key storage needs two integersperFIRSTkey. The sourcebook's
 # 200point allowance has178spare records (actual16keys+3pointpairs=22).
 # Rebucket62 ofthose spare records intoFIRST, preserving the WHOLEbound.
 first_actual=2*fc+2*len(first);native_actual=16*5+count({i:dict(v) for i,v in n['raw_point_S0_columns48']})+25*3
 assert first_actual==1086 and native_actual==177
 assert ledger['actual_FIRST62_sparse']+ledger['native_aux_point_kernel_responses']==1379>=first_actual+native_actual==1263
 rebucket=dict(ledger);rebucket['actual_FIRST62_sparse']+=len(first);rebucket['native_aux_point_kernel_responses']-=len(first)
 assert sum(rebucket.values())==derived and rebucket['native_aux_point_kernel_responses']>=native_actual
 n0=d['input_numeric_before_finalmetadata'];knownNEW=d['NEWstate_numeric_upper_for_two_extra_codec_copies'];guard=d['guard_extra_codec_reserve']
 assert knownNEW==7826 and guard==2*knownNEW+1000==16652
 assert d['shared_codec_extra_numeric_upper']==max(guard,2*d['largest_success_codec_numeric_upper']+1000,2*d['rank_codec_numeric_upper']+1000,2000)
 assert d['whole_numeric_cap']==2*(n0+1000)+derived+guard==58883<60000 and count(d)<=n0+1000
 assert d['whole_RAW_cap']==2_000_000+58883*(2*math.ceil(354*math.log10(2))+50)==17545112<32*1024*1024
 E={(name,t):tuple(map(int,b)) for name,t,b in d['positive_global_envelopes']};U0,U1=E['upper',0][1],E['upper',1][1];minorU=2*U0*U1;compareU=2*max(U0,U1)*minorU
 assert minorU==int(d['rank_integer_minor_bound']) and compareU==int(d['rank_TEMP_comparison_bound']) and minorU.bit_length()<=334 and compareU.bit_length()<=502
 worker='code/ambient_B18_ell2_T2_upper20_values_v4_20261105_v63.py';pins[worker]=sha(ROOT/worker);assert pins[worker]==cu['worker_sha256']=='cdb2418923e09c9d12330e02d446fdf0251e41eef7bf79263f36041e494851ba'
 source=(ROOT/worker).read_text();allowed={'argparse','gc','hashlib','json','resource','signal','time','traceback','fractions','pathlib'}
 for node in ast.walk(ast.parse(source)):
  if isinstance(node,ast.Import):assert all(a.name in allowed for a in node.names)
  elif isinstance(node,ast.ImportFrom):assert node.module in allowed
 assert source.index('traceback.clear_frames(tb)')<source.index("dump(out/'FAILURE_PREFIX.json'") and 'error.__traceback__=None;error.__context__=None;error.__cause__=None' in source
 assert 'signal.alarm(55)' in source and 'signal.alarm(5)' in source and 'gc.collect()' in source
 assert "work['current_env_check']" in source and "work['last_completed_binary_result']" in source
 assert cu['source_gate_sha256']==pins[gp+'report.json'] and cu['whole_numeric_cap']==58883 and cu['whole_RAW_cap']==17545112
 for path,h in cu['worker_and_source_sha256'].items():assert sha(ROOT/path)==h;pins[path]=h
 result=subprocess.run([sys.executable,str(ROOT/worker),'--help'],capture_output=True,timeout=10,cwd=ROOT);assert result.returncode==0 and not result.stderr and hashlib.sha256(result.stdout).hexdigest()==cu['safe_help_stdout_sha256']
 owned=count([gr,d,cu,orr,os,nr,ns,defs,fines,coords,c0,recipes,coverage,zero,basis,bi,tab,leaves,C,T,L0,Le,first,ledger,rebucket,E,pins])+15000;assert owned<150000
 ep='independent-final20-source-frame-lifetime-owner-code-custody-proof.json';ev={'SOURCE_CODE_only_no20responses':True,'native16FORMAL_to15actualkernels_and10explicitzero_certificates':True,
  'all28active_plus_explicitzero_case_map_terms_verified':True,'all20_full128constant_source_reconstructions_exact':True,'case617_nativecorrected_once_OLD_FIRST_outer_raw617plus20raw563':True,
  'all584SOURCEzero_plus554actualzero_plus20_complete605_coverage':True,'OLD_physical_vs_only79_arbitrary_extension_scope_preserved':True,
  'allBoolean_C0_44_typed223_FIRST481_shapes_and871factors':True,'full_native_ell_SUM_L0Se_LeS0_and_native92_fullL2S0_preserved':True,
  'independent_joint_FIRST_NATIVE_record_book':{'FIRST_tuple_key_records':first_actual,'native_aux_point_kernel_records':native_actual,'joint_allocated':1379,'joint_used_upper':1263,'rebucketed_FIRST':rebucket['actual_FIRST62_sparse'],'rebucketed_native':rebucket['native_aux_point_kernel_responses']},
  'unchanged_whole_numeric58883_hard60000_andRAW17545112':True,'exception_frame_context_cause_cycle_release_before_GUARD':True,'success_guard_codecs_disjoint_after_release':True,
  'external_ROOT_timeout60_REQUIRED_for_absolute_deadline':'55s calculation+cleanup+remainingcheckpoint window; ROOT runner adds external60s withoutcodechange.',
  'committed354_TEMP1063_andrank334_502_verified':True,'allsource_old_failedbudget_andV2V3withdrawnbytes_hashpreserved':True,'exactV4worker_AST_stdlib_safehelp_andcustody':True,
  'no_new_case_phi_gauge_NATIVE_or_higher_value':True,'source_constant_bits':cbits,'reviewer_owned':owned}
 dump(ep,ev);report={'id':STEM,'all_passed':True,'fresh_SOURCE_CODE_whole_custody_PASS':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep),'before-constant-source-and-lifetime-checks.json':sha(out/'before-constant-source-and-lifetime-checks.json')},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'reviewer_owned_numeric':owned,'reviewer_RAW_upper':raw+sum(p.stat().st_size for p in out.iterdir()),'seconds':time.monotonic()-start}
 dump('report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owned,'source_constant_bits':cbits}))
if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--output',required=True,type=Path);a=pa.parse_args();main(a.output)
