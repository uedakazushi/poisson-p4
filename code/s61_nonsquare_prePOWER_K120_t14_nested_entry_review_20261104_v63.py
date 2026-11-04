#!/usr/bin/env python3
"""Fixed source + ordinary symbolic nmod N7, no producer imports."""
import hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE='S61_NONSQUARE_PREPOWER_K120_T14_NESTED_ENTRY_20261104_v63';GATE='S61_NONSQUARE_PREPOWER_K120_T14_NESTED_SOURCE_PREFLIGHT_20261104_v63';ID=BASE+'_REVIEW';OUT=ROOT/'reports'/ID;P=32003
BP=[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={};cache={};denchecks=0
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;k=str(p.relative_to(ROOT));pins[k]=sha(p)
  if k not in cache:cache[k]=json.loads(p.read_text())
  return cache[k]
 def reduce(c):
  nonlocal denchecks
  c=Q(c);assert c.denominator%P;denchecks+=1;return (c.numerator%P)*pow(c.denominator,-1,P)%P
 r=load(BASE);assert pins[f'reports/{BASE}/report.json']=='180d03a78234db32a00548dedbebca37e0dbfd4be42a84b48d5ff00af9a9508f'
 for k,h in r['input_sha256'].items():assert sha(ROOT/k)==h
 for k,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/BASE/k)==h
 assert load(GATE+'_REVIEW')['all_passed'] and pins[f'reports/{GATE}_REVIEW/report.json']=='2c48439b98402610c26b5d9c22c3e34b5b06b36bfd27bacca901841c70eda2bd'
 plan=load(GATE,'nested-plan.json');old=load(*plan['external_scalar_nodes']['path'].split('/')[1:]);N=old['nodes']+plan['appended_scalar_nodes'];q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63','certificate.json');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63','certificate.json');hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json')
 def source(n):
  x=n['reference']
  if 'constant'in x:return [[[0,0,0],x['constant']]]
  if 'named_FULLU_unit'in x:return BP
  if x.get('derived_existing_coefficient_difference'):
   z={tuple(m):Q(c)for m,c in hc['raw_compact14x5'][11][1]}
   for m,c in hc['raw_compact14x5'][9][1]:z[tuple(m)]=z.get(tuple(m),Q(0))-Q(c)
   return [[list(m),str(c)]for m,c in z.items()if c]
  if x['path']=='accepted_Q2L_source_atom':
   a=x['definition'];assert a==q['scalar_DAG']['nodes'][x['old_DAG_root']]
   num,den=cp['matrix25x32'][a[1]][a[2]]if a[0]=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1]if a[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
   assert len(den)==1 and den[0][0][1:]==[0,0]and Q(n['QQ_unit'])==Q(den[0][1])and n['a_shift']==-den[0][0][0];return num
  z=load(*x['path'].split('/')[1:])
  if 'old_DAG_root'in x:return next(v['Laurent_polynomial']for v in z['node_values']if v['old_DAG_root']==x['old_DAG_root'])
  if 'h9_target'in x:return z['source_h9'][z['original_target_labels'].index(x['h9_target'])]
  z=z[x['key']]
  for j in x.get('index',[]):z=z[j]
  return z
 needed=set()
 for e in plan['outer21']:
  needed.update(e['weighted2_source_nodes'])
  for t in e['inner10']:needed.update(t['source_nodes'])
 atoms={i:source(N[i])for i in needed if N[i]['op']=='atom'};assert sum(len(x)for x in atoms.values())==11990
 # Before the one arithmetic worker: exact source/input normalization and p-units.
 own=OUT/'independent-input.txt';audit=[]
 with own.open('w')as f:
  f.write(f'{len(N)}\n')
  for i,n in enumerate(N):
   if i not in needed:f.write(f'{i} x\n');continue
   upper=n['bound']['support']if n['bound']else 0;assert int(n['QQ_denominator'])%P
   if n['op']=='atom':
    raw=atoms[i];s=min((m[0]+n['a_shift']for m,c in raw if Q(c)),default=0);unit=Q(n['QQ_unit']);assert unit.numerator%P and unit.denominator%P
    vals=[([m[0]+n['a_shift']-s,m[1],m[2]],reduce(Q(c)/unit))for m,c in raw if Q(c)];f.write(f'{i} p {upper} {s} {len(vals)}\n')
    for m,c in vals:f.write(' '.join(map(str,m+[c]))+'\n')
    audit.append(dict(node=i,reference=n['reference'],unit=str(unit),shift=s,terms=len(vals)))
   elif n['op']=='known_a_shift':f.write(f"{i} k {upper} {n['args'][0]} {n['a_exponent_shift']}\n")
   else:f.write(f"{i} "+dict(add='a',sub='s',mul='m')[n['op']]+' '+str(upper)+' '+' '.join(map(str,n['args']))+'\n')
  b=plan['output_bound'];f.write(' '.join(map(str,[v for z in b['box']for v in z]+b['degree']))+'\n21\n')
  for e in plan['outer21']:
   need=e['weighted2_source_nodes'];f.write(f"{e['weighted2_root']} {e['outer_sign']} {len(need)} "+' '.join(map(str,need))+'\n10\n')
   for t in e['inner10']:
    need=t['source_nodes'];f.write(f"{t['left2minor']} {t['right3minor']} {t['sign']} {len(need)} "+' '.join(map(str,need))+'\n')
  expected=load(BASE,'native.json')['normalized_N7'];assert int(expected['global_denominator'])%P;f.write(expected['global_denominator']+' '+str(len(expected['integer_coefficients']))+'\n')
  for m,c in expected['integer_coefficients']:f.write(' '.join(map(str,m+[c]))+'\n')
 # Native ASCII and saved QQ denominator side independently source-bound.
 ip=ROOT/'reports'/BASE/'native-input.txt';tokens=iter(ip.read_text().split());assert (next(tokens),next(tokens))==('11990','11523')
 for i,n in enumerate(N):
  k=next(tokens)
  if i not in needed:assert k=='x';continue
  upper,D=next(tokens),next(tokens);assert int(upper)==(n['bound']['support']if n['bound']else 0)and D==n['QQ_denominator']and int(D)%P
  if k=='p':
   sh,u,count=int(next(tokens)),next(tokens),int(next(tokens));raw=atoms[i];amin=min((m[0]for m,c in raw),default=0);assert sh==n['a_shift']+amin and Q(u)==Q(n['QQ_unit'])and count==len(raw)
   payload={(int(next(tokens)),int(next(tokens)),int(next(tokens))):Q(next(tokens))for _ in range(count)};assert payload=={(m[0]-amin,m[1],m[2]):Q(c)for m,c in raw}
  elif k=='k':assert int(next(tokens))==n['args'][0]and int(next(tokens))==n['a_exponent_shift']
  else:assert dict(a='add',s='sub',m='mul')[k]==n['op']and [int(next(tokens)),int(next(tokens))]==n['args']
 # The following plan fields are separately bound by the completed source-only peer.
 maxTop=max(e['top5_dense_slots']for e in plan['outer21']);live=78338+11990+39549+maxTop+1000;assert live==146711 and live<150000
 dump(OUT/'preflight.json',dict(primary_report=pins[f'reports/{BASE}/report.json'],all_denominator_checks=denchecks,all_source_atom_normalizations=True,all11523_native_source_records_bound=True,max_live_upper=live,method='ordinary symbolic nmod source and inner products; fused outer only',no_parameter_point_or_newprime=True))
 dump(OUT/'phase.json',dict(phase='preflight complete; before one independent symbolic worker'))
 dump(OUT/'source-binding.json',dict(atoms=audit,phase_plan_sha256=pins[f'reports/{GATE}/nested-plan.json'],source_peer=pins[f'reports/{GATE}_REVIEW/report.json']))
 cache.clear();del atoms,N,q,cp,pd,hc,old,plan,expected,tokens
 cpp=Path(__file__).with_suffix('.cpp');exe=OUT/'independent-worker';cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)';s=subprocess.run(['bash','-lc',cmd],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=max(1,57-(time.monotonic()-start)));(OUT/'compile.log').write_text(s.stdout);assert s.returncode==0
 cmd='source /workspace/cloud-setup/activate.sh && exec '+str(exe)+' '+str(own)+' '+str(OUT/'native.json')+' '+str(OUT/'progress.log');s=subprocess.run(['bash','-lc',cmd],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=max(1,58-(time.monotonic()-start)));(OUT/'native.log').write_text(s.stdout);assert s.returncode==0
 nr=json.loads((OUT/'native.json').read_text());assert nr['all_passed']and nr['all14986_stored_QQ_coefficients_reduce_to_same_polynomial']
 for p in[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 files={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name not in['report.json','independent-worker']};assert sum(p.stat().st_size for p in OUT.iterdir()if p.is_file()and p.name!='independent-worker')<32*1024*1024
 dump(OUT/'phase.json',dict(phase='FINAL symbolic polynomial NONZERO, no full QQ coefficient identity'));files['phase.json']=sha(OUT/'phase.json')
 report=dict(id=ID,all_passed=True,frozen=True,source_bound_symbolic_modp_nonzero=True,full_QQ_coefficient_identity_NOT_claimed=True,all14986_saved_coefficients_reduce_correctly=True,actual_N7_NONZERO_char0=True,input_sha256=pins,evidence_sha256=files,seconds=time.monotonic()-start,scope='ONE K120/target14; no other source, point, prime, or full rank claim')
 dump(OUT/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
def finalize_existing():
 pins={};primary=json.loads((ROOT/'reports'/BASE/'report.json').read_text())
 for rel,h in primary['input_sha256'].items():assert sha(ROOT/rel)==h;pins[rel]=h
 pins[f'reports/{BASE}/report.json']=sha(ROOT/'reports'/BASE/'report.json')
 for rel in [f'reports/{GATE}_REVIEW/report.json',f'reports/{BASE}/native.json']:
  pins[rel]=sha(ROOT/rel)
 for path in [Path(__file__),Path(__file__).with_suffix('.cpp'),ROOT/f'research/notes/{ID}.md']:pins[str(path.relative_to(ROOT))]=sha(path)
 nr=json.loads((OUT/'native.json').read_text());assert nr['all_passed'] and nr['all14986_stored_QQ_coefficients_reduce_to_same_polynomial']
 dump(OUT/'phase.json',dict(phase='FINAL source-bound symbolic polynomial NONZERO; full QQ coefficient identity not claimed'))
 dump(OUT/'report.json',dict(id=ID,all_passed=True,frozen=True,source_bound_symbolic_modp_nonzero=True,full_QQ_coefficient_identity_NOT_claimed=True,all14986_saved_coefficients_reduce_correctly=True,actual_N7_NONZERO_char0=True,arithmetic_seconds=nr['seconds'],peak_visible_slots=nr['peak_visible_slots_with1000reserve'],compile_only_failure_preserved=True,arithmetic_worker_runs=1,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name not in['report.json','independent-worker']},scope='ONE K120/target14; no other source, point, prime, or full rank claim'))
if __name__=='__main__':
 import sys
 if '--finalize-existing' in sys.argv:finalize_existing()
 else:main()
