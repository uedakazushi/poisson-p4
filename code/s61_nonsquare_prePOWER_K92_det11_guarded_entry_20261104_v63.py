#!/usr/bin/env python3
"""ONE guardeddet11 attempt; restores exact prefix, never promotes incomplete sum."""
import argparse,hashlib,importlib.util,json,math,resource,shlex,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def tc(x):
 if isinstance(x,list):
  if len(x)==2 and isinstance(x[0],list)and len(x[0])==3 and all(isinstance(i,int)for i in x[0])and isinstance(x[1],str):return 1
  return sum(tc(v)for v in x)
 if isinstance(x,dict):return sum(tc(v)for v in x.values())
 return 0
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={};cache={};retained=peak=0;events=[];dump(out/'phase.json',dict(phase='prepareONEguardedattempt; no source arithmetic yet'))
 def source(rel):
  nonlocal retained,peak
  if rel not in cache:
   p=ROOT/rel;pins[rel]=sha(p);x=read(p);full=tc(x);peak=max(peak,retained+full);assert peak+15000<=150000
   if 'CUBIC_GENUINE14_H5_PREFLIGHT' in rel:x={k:x[k]for k in['raw_compact14x5','known_row_units_rho14']}
   elif 'H22_PURE3_FIXED_MINOR' in rel:x={'node_values':x['node_values']}
   elif 'RAW8_BOTTOM24_SOURCE' in rel:x={'bottom4x6':x['bottom4x6']}
   elif 'CUBIC_ZW25_PAIR_DUAL_PREFLIGHT' in rel:x={'Lframe_dual6x2':x['Lframe_dual6x2']}
   elif 'RAW_H_GENERIC_RANK' in rel:x={'nonzero_selected_minor_Delta':x['nonzero_selected_minor_Delta'],'cleared_upper_row_coordinates7xr':[[],x['cleared_upper_row_coordinates7xr'][1],[],x['cleared_upper_row_coordinates7xr'][3],x['cleared_upper_row_coordinates7xr'][4]]}
   kept=tc(x);events.append(dict(path=rel,full_terms=full,kept_terms=kept,parse_peak=retained+full));retained+=kept;cache[rel]=x
  return cache[rel]
 base='reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63';gate=source(base+'/report.json');meta=source(base+'/source-circuit-metadata.json');norm=source(base+'/normalization.json');assert pins[base+'/report.json']=='9a2a547edf9c0fae2d4776463b28990b4867e320eb5cbb9ebda9676fac0b86bf'and gate['guarded_attempt_gate_passes']
 sb='reports/S61_NONSQUARE_PREPOWER_K92_H3_SIGN_SCHEDULE_20261104_v63';assert source(sb+'/report.json')['all_passed'];schedule=source(sb+'/schedule.json')['ten_scalar_terms'];assert len(schedule)==10
 assert source('reports/S61_NONSQUARE_PREPOWER_RAWCROSS18_ENTRY_20261104_v63_REVIEW/report.json')['all_passed']
 q=source('reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json');cp=source('reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');pd=source('reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json');hc=source('reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json');rk=source('reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63/certificate.json');N=meta['nodes'];roots=meta['matrix11x8_input_nodes'];assert len(N)==575
 need=set();todo=[i for r in roots for i in r]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i);todo+=N[i]['args']
 def atom(n):
  r=n['reference']
  if 'constant'in r:return [[[0,0,0],r['constant']]]
  if 'named_FULLU_unit'in r:return [[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
  if r['path']=='accepted_Q2L_source_atom':
   assert not r.get('known_B_denominator_power',0);a=r['definition'];return cp['matrix25x32'][a[1]][a[2]][0]if a[0]=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1][0]if a[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]][0]
  x=source(r['path'])
  if 'old_DAG_root'in r:return next(v['Laurent_polynomial']for v in x['node_values']if v['old_DAG_root']==r['old_DAG_root'])
  if 'h9_target'in r:return x['source_h9'][x['original_target_labels'].index(r['h9_target'])]
  x=x[r['key']]
  for j in r.get('index',[]):x=x[j]
  return x
 def normalized(poly,unit=Q(1),ashift=0,expected=None):
  p=[([m[0]+ashift,m[1],m[2]],Q(v)/unit)for m,v in poly if Q(v)];assert 3*len(p)<=15000;D=math.lcm(*(v.denominator for m,v in p))if p else 1;nums=[v.numerator*(D//v.denominator)for m,v in p];G=__import__('functools').reduce(math.gcd,map(abs,nums),0);c=Q(G,D)if G else Q(1)
  if expected is not None:assert c==Q(expected)
  ints=[(m,v//G)for(m,_),v in zip(p,nums)]if G else[];assert all(abs(v).bit_length()<=1024 for m,v in ints);assert abs(c.numerator).bit_length()<=4096 and c.denominator.bit_length()<=4096;return c,ints
 extras=[hc['raw_compact14x5'][r][j]for r in[10,11,8,6,12,9]for j in[1,3,2]]+[rk['cleared_upper_row_coordinates7xr'][3][0],rk['cleared_upper_row_coordinates7xr'][4][0]]+rk['cleared_upper_row_coordinates7xr'][1]+[rk['nonzero_selected_minor_Delta']];assert len(extras)==24
 ip=out/'native-input.txt';bindings=[]
 with ip.open('w')as f:
  f.write('575 24\n')
  for i,n in enumerate(N):
   if i not in need:f.write('x 0 1\n');continue
   upper=n['bound']['support']if n['bound']else 0;c=n['content_lower']
   if n['op']=='atom':
    qc,p=normalized(atom(n),Q(n['QQ_unit']),n['a_shift'],c);f.write(f'p {upper} {qc}\n{len(p)}\n')
    for m,v in p:f.write(' '.join(map(str,m))+f' {v}\n')
    bindings.append(dict(node=i,reference=n['reference'],QQ_unit=n['QQ_unit'],a_shift=n['a_shift'],content=str(qc),terms=len(p)))
   elif n['op']=='known_a_shift':f.write(f"k {upper} {c} {n['args'][0]} {n['a_exponent_shift']}\n")
   else:f.write(dict(add='a',sub='s',mul='m')[n['op']]+f' {upper} {c} '+' '.join(map(str,n['args']))+'\n')
  for poly in extras:
   c,p=normalized(poly);f.write(f'{c} {len(p)}\n')
   for m,v in p:f.write(' '.join(map(str,m))+f' {v}\n')
  f.write(str(len(need))+' '+' '.join(map(str,sorted(need)))+'\n')
  for r in roots:f.write(' '.join(map(str,r))+'\n')
  f.write('18 '+' '.join(map(str,range(18)))+'\n6 18 19 20 21 22 23\n4\n')
  for rr,p in zip([[0,2,4],[0,2,5],[1,2,4],[1,2,5]],norm['cofactor_content_Newton_proofs']):f.write(' '.join(map(str,rr))+' '+p['exact_content_from_Gauss']+'\n')
  f.write(norm['content_lower_divisor']+' '+str(norm['known_output_a_power'])+'\n')
 assert ip.stat().st_size<=3*1024*1024
 dump(out/'source-binding.json',dict(source_gate=base+'/report.json',source_gate_sha256=pins[base+'/report.json'],actual_source575_nodes_needed=sorted(need),actual_pure_matrix_roots=roots,input_atoms=bindings,extra_H3_rows=[10,11,8,6,12,9],extra_H3_columns=[1,3,2],extra_lambda_order=['415','442','505','419','428','Delta'],H3_literal_rows=[[0,2,4],[0,2,5],[1,2,4],[1,2,5]],normative_sign_schedule=schedule,
  NO_double_sign_on_literal_H3=True,g_is_kept_factored=True,source_rows_are_current15_with_all_Q_corrections=True,all10_terms_and_all_tail_effects_retained=True))
 dump(out/'preparation-resource.json',dict(load_events=events,explicit_source_input_peak=peak,with15000reserve=peak+15000,input_bytes=ip.stat().st_size,input_byte_cap=3*1024*1024,numeric_cellcap=135000,work15000=True,QQcontent_bitcap=4096,coefficient_bitcap=1024,prefix_RAW_upper=26980736))
 dump(out/'phase.json',dict(phase='sourceinputfrozen beforecompile/fixture; scalarNOTstarted'));cache.clear();del q,cp,pd,hc,rk,N,extras,bindings,meta,norm
 cpp=Path(__file__).with_suffix('.cpp');parser=ROOT/'code/s61_nonsquare_prePOWER_K92_restore_sparse_state_20261104_v63.py';note=ROOT/f'research/notes/{ID}.md';exe=out/'native-worker'
 cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+shlex.quote(str(cpp))+' -o '+shlex.quote(str(exe))+' $(pkg-config --cflags --libs flint)';run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,55-(time.monotonic()-began)));(out/'compile.log').write_text(run.stdout)
 if run.returncode:
  (out/'compile-failed-source.cpp').write_bytes(cpp.read_bytes());dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,arithmetic_STARTED=False,reason='compile-onlyfailure, sourcebytes/logpreserved',code_sha256=sha(cpp)));signal.alarm(0);return
 cmd='source /workspace/cloud-setup/activate.sh && exec '+shlex.quote(str(exe))+' --fixture '+shlex.quote(str(out));fixture=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=4);(out/'fixture.log').write_text(fixture.stdout);assert fixture.returncode==0
 spec=importlib.util.spec_from_file_location('restore',parser);restore=importlib.util.module_from_spec(spec);spec.loader.exec_module(restore);fixtureState=restore.read_state(out/'fixture.bin');restore.check_fixture(fixtureState);dump(out/'fixture-restore.json',dict(all_passed=True,version=restore.VERSION,fixture_sha256=sha(out/'fixture.bin'),product_and_precommit_cell_bit_guards=True,parser_sha256=sha(parser)))
 dump(out/'phase.json',dict(phase='compile/fixtures/restorerPASS beforeONE actualscalar'));assert sum(p.stat().st_size for p in out.iterdir())<2*1024*1024+3*1024*1024
 cmd='source /workspace/cloud-setup/activate.sh && exec '+shlex.quote(str(exe))+' '+shlex.quote(str(ip))+' '+shlex.quote(str(out))+' '+sha(ip);run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,58-(time.monotonic()-began)));(out/'native.log').write_text(run.stdout)
 resultpath=out/'native.json'if(out/'native.json').exists()else out/'frontier.json';nr=read(resultpath)if resultpath.exists()else dict(all_passed=False,reason='no nativecompletionrecord',returncode=run.returncode);state=restore.read_state(out/'state.bin');assert state['header']['input_sha256']==sha(ip);assert state['total_terms']<=135000
 full=bool(nr.get('full_scalar_complete'));assert full==state['header']['complete'];stateHash=sha(out/'state.bin');statecheck=dict(all_passed=True,version=restore.VERSION,state_sha256=stateHash,complete=full,active_polynomials=len(state['polynomials']),exact_saved_terms=state['total_terms'],has_uncommitted_update=state['pending']is not None,header=state['header'],restored_but_arithmetic_NOT_resumed=True);dump(out/'restored-state-check.json',statecheck)
 for p in[Path(__file__),cpp,parser,note]:pins[str(p.relative_to(ROOT))]=sha(p)
 total=sum(p.stat().st_size for p in out.iterdir());assert total<=26980736 and total<32*1024*1024
 dump(out/'phase.json',dict(phase='FINAL completeactualdet11'if full else'FINAL exactprefixonly; noautoretry; NOdetnonzero-proof'))
 final=dict(id=ID,all_passed=full,frozen=True,source_scalar_attempt_once=True,safety_and_restore_checks_passed=True,full_scalar_complete=full,actual_det11_nonzero=nr.get('output_nonzero')if full else None,output_terms=nr.get('output_terms')if full else None,completed_H3_terms=nr.get('completed_terms'),guard_reason=nr.get('reason'),coefficient_updates=nr.get('coefficient_updates'),peak_owned_cells=nr.get('peak_owned_cells'),peak_owned_plus_work=nr.get('peak_owned_plus_work'),all_active_state_restored=True,state_version=restore.VERSION,state_sha256=stateHash,
  exact_restore_command=f'.venv/bin/python code/{parser.name} --check {out.relative_to(ROOT)}/state.bin',restore_does_NOT_resume_arithmetic=True,no_prefix_to_fullscalar_inference=True,no_other_entry_retry_or_global57_claim=True,total_RAW_bytes_including_binary_before_report=total,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-began)
 dump(out/'report.json',final);signal.alarm(0);print(json.dumps({k:v for k,v in final.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
