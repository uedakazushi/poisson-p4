#!/usr/bin/env python3
"""Independent whole554 source/assembly/height/custody; no response arithmetic."""
import argparse,ast,hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_UPPER554_WHOLE_GATE_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def leaves(v):
 if isinstance(v,(list,tuple)):
  for x in v:yield from leaves(x)
 else:yield F(v)
def main(out,nativepeer,custody):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;metric=0
 def read(p):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size;assert raw<32*1024*1024;return json.loads(q.read_text())
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 def ck(n):
  nonlocal metric
  metric=max(metric,abs(n).bit_length());assert metric<=384 and 3*metric+1<=1153;return n
 def env(q):q=F(q);return(q.denominator,abs(q.numerator)) if q else(1,0)
 def add(a,b):
  d=ck(math.lcm(a[0],b[0]));return d,ck(a[1]*(d//a[0])+b[1]*(d//b[0]))
 def mul(a,b):return (ck(a[0]*b[0]),ck(a[1]*b[1])) if a[1] and b[1] else(1,0)
 def dec(q):return int(q[0]),int(q[1])
 def total(v):
  z=(1,0)
  for q in v:z=add(z,dec(q))
  return z
 wp='reports/AMBIENT_B18_ELL2_T2_UPPER554_WHOLE_INPUT_GATE_20261105_v63/'
 wr=read(wp+'report.json');assert pins[wp+'report.json']=='e250de02066f1df5abdbc512c4a268a9e3f83ca1722e0f34ba97986c01d9cdfa'
 ep='ONE554-whole-positive-input-fine-assembly-and-native339-height.json';w=read(wp+ep);assert pins[wp+ep]==wr['evidence_sha256'][ep]
 op='reports/AMBIENT_B18_ELL2_T2_UPPER554_OLD_FIRST_20261105_v63_REVIEW/report.json';oldpeer=read(op);assert oldpeer['all_passed'] and pins[op]=='130f7b209b26c56cd68d8c0a2666a93028d7a627dab88b5fd5d0041bfa4c524c'
 np=read(nativepeer);assert np['all_passed'] and '554' in np['id']
 assert np['input_sha256']['reports/AMBIENT_B18_NATIVE_ONE554_POSITIVE_INPUT_20261105_v63/report.json']=='e9dbeecb12be941d83a97b701ce7b723c6db06200bdc698d8233300d5b6e11bd'
 n=w['NATIVE_outer_component'];o=w['OLD_FIRST_component'];p=n['full_primitive_input_J18_N_G_nu_GammaD']
 assert n['selected_case619']==o['selected_case619']==w['selected_case619']==554
 originals={'OLD_FIRST_component':('reports/AMBIENT_B18_ELL2_T2_UPPER554_OLD_CARRY_INPUT_GATE_20261105_v63/ONE554-OLD-two-columns-and57-FIRST-positive-component-input.json',o),
  'NATIVE_outer_component':('reports/AMBIENT_B18_NATIVE_ONE554_POSITIVE_INPUT_20261105_v63/ONE554-positive-four-native-and-nonzero644-outer-input.json',n)}
 for key,(path,x) in originals.items():assert read(path)==x
 assert p['Gamma_two56']==o['Gamma_two'] and p['D_two_by10_by56']==o['D_two_by10_by56']
 assert n['all29_a1_kept_UNEVALUATED']==o['full_a1_29'] and n['plane_a1_pivot5_original644']==[[52,'-175/3'],[53,'-377/60']]
 c0=sorted(o['OLD_source2']+[[r['original768'],r['constant_source_coefficient']] for r in n['parameterized_four_native_source_records']]);assert w['literal_source_all6']==c0
 # Literal value certificate for lower4=0, distinct from the617 source recipe.
 lp='reports/AMBIENT_B18_ELL2_T2_NATIVE70_VALUES_20261105_v63/actual-native70-lower-prefix-and-witness-if-any.json';lower=read(lp);rec=next(v for v in lower['computed_prefix'] if v['case619']==554);assert rec['actual_lower4_values']==['0']*4
 lpr='reports/AMBIENT_B18_ELL2_T2_NATIVE70_PREFIX14_20261105_v63_REVIEW/report.json';lr=read(lpr);assert lr['all_passed'] and pins[lpr]=='868ba1c858dcc03946180b2cbf952437cb8f0f7323bc1dc3e48363b43911bb3f' and lr['input_sha256'][lp]==pins[lp]
 rp='reports/AMBIENT_B18_NATIVE_PIVOT23_SECOND_VALUES_20261105_v63/actual92-Gamma-nativepivot23-plane2-and322-full-product-pieces.json';rv=read(rp);i=rv['original768_native23_sources'].index(644);assert i==1
 assert n['actual644_GammaXp_rr_rell_twochannels_UNPAIRED']==rv['Gamma_Xp2_native23_twochannels_rr_rell'][i]
 assert any(F(q) for ch in n['actual644_GammaXp_rr_rell_twochannels_UNPAIRED'] for q in ch)
 # Only envelope products, never the signed outer pairing.
 plane=[F('-377/60'),F('-175/3')]
 for t in range(2):
  for k in range(2):assert mul(env(n['actual644_GammaXp_rr_rell_twochannels_UNPAIRED'][t][k]),env(plane[k]))==dec(n['positive_outer644_two_by2'][t][k])
 basis=[tuple(v) for v in p['basis18']];bi={v:i for i,v in enumerate(basis)};assert len(basis)==18 and p['N_known_positions']==[0,1,2,3,4,5,6,15,16]
 assert {e[1] for r in n['parameterized_four_native_source_records'] for e in r['intrinsic_source_entries']}=={0,1,2,3}
 assert all(r['full_point_S0_column48']==[] for r in n['parameterized_four_native_source_records'])
 for a in p['nu_two_J18']:
  for b in p['nu_two_J18']:assert all(tuple(x+y for x,y in zip(basis[i],basis[j])) not in bi for i,u in enumerate(a) if F(u) for j,v in enumerate(b) if F(v))
 h=w['native_common_denominator_height'];A=[q for row in p['Gamma_two56'] for i,q in row]+[q for ch in p['D_two_by10_by56'] for row in ch for q in row]
 cls={'known_partial_N9':p['N_known12_by12_by9_by4'],'canonical_G8':p['G8_C72'],'monics':p['nu_two_J18'],'Gamma_and_D_scalar_coefficients':A,'native_c0_constants':[r['constant_source_coefficient'] for r in n['parameterized_four_native_source_records']]}
 D=1;U=2**28
 for name,v in cls.items():
  xs=list(leaves(v));d=1
  for q in xs:d=ck(math.lcm(d,q.denominator))
  u=max([d]+[ck(abs(q.numerator)*(d//q.denominator)) for q in xs]);r=h['scalar_operand_classes'][name];assert str(d)==r['D'] and str(u)==r['U_including_scalar_one'];D=ck(D*d);U=ck(U*u)
 assert 4*2*12*16*9*18*19*11==h['explicit_native_direct_paths_upper']<2**28
 assert str(D)==h['common_D'] and str(U)==h['common_absolute_scaled_U_with_path2pow28'] and max(D.bit_length(),U.bit_length())==339
 metric=0;ng=[];nd=[];outer=[];upper=[]
 for t in range(2):
  ng.append(total(n['positive_Gamma_weighted_cubic_two_by5'][t]));nd.append(total([v for x in n['positive_D_weighted_second_two_by10_unique_first_decomp'][t] for v in x]));outer.append(total(n['positive_outer644_two_by2'][t]))
  # Different ordering: pair nativeΓ withOLDΓ, nativeD withOLD D, and
  # FIRST withouter; combinepositivegroupsinreversechannelindependentorder.
  a=add(ng[t],dec(o['positive_OLD_Gamma_cubic_two'][t]));b=add(nd[t],dec(o['positive_OLD_D_second_two'][t]));c=add(outer[t],dec(o['positive_FIRSTcarry_total_two'][t]));upper.append(add(add(c,b),a))
 for key,v in [('positive_fine_native_Gamma_two',ng),('positive_fine_native_D_two',nd),('positive_outer_total_two',outer),('positive_WHOLE_upper_two',upper)]:assert [dec(q) for q in w[key]]==v
 six=[];Q=w['fixed_Q_embedding_two6']
 for i in range(6):six.append(add(mul(env(dict(Q[1]).get(i,'0')),upper[1]),mul(env(dict(Q[0]).get(i,'0')),upper[0])))
 assert [dec(q) for q in w['positive_WHOLE_full_six']]==six and metric<=214
 # Literal ownership and semantic dimensions, no producer count import.
 before={k:v for k,v in w.items() if k not in {'input_numeric_before_finalmeta','metadata_reserve_per_input_copy','explicit_derived_workspace','success_live_upper','guard_live_upper','future_whole_numeric_cap','future_whole_RAW_cap','guard_prefix_rule'}}
 n0=count(before);assert n0==w['input_numeric_before_finalmeta']==13805
 ws=w['explicit_derived_workspace'];assert ws['native_knownN_required4rows_by_all12targets_9base4mask']==4*12*9*4
 assert ws['native_Sc4targets18base4mask_sparse_index_value']==4*18*4*2
 dest={i//4 for ch in p['D_two_by10_by56'] for row in ch for i,q in enumerate(row) if q!='0' and i<48}|{i//4 for row in p['Gamma_two56'] for i,q in row if i<48};assert len(dest)==7
 assert ws['native_NSc_targets_by18base4mask_sparse_index_value']==len(dest)*18*4*2
 fs=sum(len(v[2]) for v in o['positive_Xp_FIRST_columns28']);assert ws['point_FIRST_cache_possible_support_index_value']==2*fs+4*28
 derived=sum(ws.values());success=2*(n0+1000)+derived+4000;guard=success+12000;assert success==43515==w['success_live_upper'] and guard==55515==w['guard_live_upper']
 assert count(w)<=n0+1000 and w['future_whole_numeric_cap']==guard and w['future_whole_RAW_cap']==2_000_000+guard*(2*math.ceil(339*math.log10(2))+42)==15767720
 # Worker hash,stdlib-onlyAST, safehelp and code-specific control proof are
 # inspectedbythereviewerBEFOREthiscertificateisrun.
 cr=read(custody);assert pins[custody]=='cc3ecd488be3cc6e2f345117621b2f1d457464eef804a08263e49edce4902924' and cr['all_passed'] and cr['source_gate_sha256']==pins[wp+'report.json']
 worker='code/ambient_B18_ell2_T2_upper554_values_v2_20261105_v63.py';source=(ROOT/worker).read_text();pins[worker]=sha(ROOT/worker);assert pins[worker]==cr['worker_sha256']=='f1ccb1cf6fafdd032d27dd48c77d8634b482911b6dd81e8064f6e7ea825cd930'
 assert cr['case554_lower_zero_literal_record']==rec and cr['accepted_lower_zero_peer_sha256']==pins[lpr]
 assert cr['guard_live_encoded_extra_two_copies_upper']==2*4080==8160<cr['guard_extra_reserve']==12000
 assert cr['whole_numeric_cap']==55515 and cr['committed_prior_bits']==339 and cr['temporary_prior_bits']==1018
 # An explicit NEWstate-only guard byte book is smaller than the job RAWcap.
 guard_raw=cr['RAW_prior_read_and_hash_files_upper']+4080*(2*math.ceil(339*math.log10(2))+50)+100000
 assert guard_raw<15767720
 allowed={'argparse','hashlib','json','math','resource','signal','sys','time','fractions','pathlib'}
 for node in ast.walk(ast.parse(source)):
  if isinstance(node,ast.Import):assert all(a.name in allowed for a in node.names)
  if isinstance(node,ast.ImportFrom):assert node.module in allowed
 run=subprocess.run([sys.executable,str(ROOT/worker),'--help'],capture_output=True,check=False,timeout=10,cwd=ROOT);assert run.returncode==0 and not run.stderr
 assert hashlib.sha256(run.stdout).hexdigest()==cr['safe_help_stdout_sha256']
 for label,hsh in cr['worker_and_source_sha256'].items():
  q=ROOT/label;assert sha(q)==hsh;pins[label]=hsh
 for path,hsh in wr['input_sha256'].items():assert sha(ROOT/path)==hsh;pins[path]=hsh
 owner=count([wr,w,oldpeer,np,n,o,p,lower,rec,lr,rv,h,A,cls,upper,six,Q,ws,cr,pins])+15000;assert owner<150000
 assert raw+15767720<32*1024*1024
 ev={'original6c0_and29a1_currentGammaD_andouter644_exact':True,'accepted554_lower4_actualzero_bound_to_prefix14_peer':True,'full6_Qupper_valid_only_after_lowerzero':True,
  'OLD_FIRST_own_dense_positive_peer':pins[op],'native68_fresh_numeric_peer':pins[nativepeer],
  'direct_C72_one_each_inputclass_path_bound339_temp1018':True,'all_double_nu_histories_supportzero_before_products':True,
  'fine_envelope_positive_assembly_and_Q_in_different_order_verified':True,'own_fine_assembly_bits':metric,
  'whole_inputs_and_success_guard_dimensions_verified':True,'success_upper':success,'guard_upper':guard,'future_RAW':15767720,
  'literal_worker_custody_stdlibAST_safehelp_PASS':True,'V1_preserved_V2_only_guard_cursor_timeout_changes':True,'guard_extra_two_encoded_copies8160_within12000':True,'explicit_guard_RAW_upper':guard_raw,'manual_algorithm_frame_order_and_preproduct_support_pruning_checked':True,
  'no_signed_native_OLD_FIRST_outer_upper_response_computed':True,'onlyONE554_not605_HF_orlength':True,'reviewer_owned':owner}
 e='independent-whole554-source-height-ownership-and-worker-custody-checks.json';dump(e,ev)
 rep={'id':STEM,'all_passed':True,'source_code_positive_whole_only':True,'input_sha256':pins,'evidence_sha256':{e:sha(out/e)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'owned_numeric':owner,'whole_RAW_upper':raw+15767720,'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owner,'assembly_bits':metric}))
if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--output',required=True,type=Path);pa.add_argument('--native-peer',required=True);pa.add_argument('--custody',required=True);a=pa.parse_args();main(a.output,a.native_peer,a.custody)
