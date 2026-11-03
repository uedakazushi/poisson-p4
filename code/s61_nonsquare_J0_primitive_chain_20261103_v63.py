#!/usr/bin/env python3
"""Bounded J0 primitive-span continuation, exporting after EVERY exact pivot."""
import argparse,hashlib,json,math,os,resource,shlex,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def job(command,log,cpu,wall):
 def bounds():resource.setrlimit(resource.RLIMIT_AS,(1800000000,1800000000));resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu+1))
 started=time.monotonic()
 with log.open('w')as stream:
  p=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT,preexec_fn=bounds,start_new_session=True)
  try:p.wait(timeout=wall);status='completed'if p.returncode==0 else'refused'if p.returncode==2 else'exit_'+str(p.returncode)
  except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait();status='wall_timeout'
 return {'status':status,'returncode':p.returncode,'seconds':time.monotonic()-started,'cpu_limit':cpu,'wall_limit':wall}
def profile(path):
 with path.open()as f:
  n,m,k=map(int,f.readline().split());assert len(f.readline().split())==k;ids=list(map(int,f.readline().split()));assert len(ids)==n
  rows=[]
  for i in range(n):
   counts=[];degree=0
   for _ in range(m):
    z=f.readline().split();nt=int(z[0]);assert len(z)==1+4*nt;counts.append(nt)
    for j in range(nt):
     ex=list(map(int,z[2+4*j:5+4*j]));assert ex[2]==0;degree=max(degree,sum(ex))
   rows.append({'row':i,'origin':ids[i],'terms':sum(counts),'max_degree':degree,'nonzero_entries':sum(c>0 for c in counts)})
  assert not f.read().strip()
 return [n,m,k],sorted(rows,key=lambda r:(r['terms'],r['max_degree'],r['row']))
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',type=Path,required=True);ap.add_argument('--expected-sha256',required=True)
 ap.add_argument('--output',type=Path,required=True);ap.add_argument('--budget',type=int,default=300);ap.add_argument('--per-row-cap',type=int,default=60)
 ap.add_argument('--priority-rows',default='',help='Initial zero-based row priorities; candidate ordering only, not an invariant-submodule assertion.')
 args=ap.parse_args();assert 0<args.budget<=300;args.output=args.output.resolve();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic();source=args.source.resolve();assert sha(source)==args.expected_sha256
 units=ROOT/'work/nonsquare-J0-pid-v63/units.txt';searchsrc=ROOT/'reports/S61_NONSQUARE_J0_PRIMITIVE_UNIT_20261103_v63/primitive-search.cpp';pivotsrc=ROOT/'code/polynomial_one_witness_pivot_20261103_v62.cpp'
 flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
 search=args.output/'search';pivot=args.output/'pivot'
 for cpp,exe in [(searchsrc,search),(pivotsrc,pivot)]:subprocess.run(['c++','-std=c++17','-O3',str(cpp),*flags,'-lgmp','-o',str(exe)],check=True)
 stages=[];state=source;initial_shape,_=profile(state);total=0
 preferred={int(x)for x in args.priority_rows.split(',')if x};assert all(0<=i<initial_shape[0]for i in preferred)
 priority=[i in preferred for i in range(initial_shape[0])]
 while args.budget-(time.monotonic()-started)>5:
  number=len(stages)+1;folder=args.output/f'stage-{number}';folder.mkdir();shape,rows=profile(state)
  rows.sort(key=lambda r:(not priority[r['row']],r['terms'],r['max_degree'],r['row']))
  if shape[0]==0:break
  stage={'source':str(state.relative_to(ROOT)),'source_sha256':sha(state),'source_shape':shape[:2],'profiles':rows,'attempts':[]};found=None
  for row in rows:
   remaining=args.budget-(time.monotonic()-started)
   if remaining<5:break
   witness=folder/f'row{row["row"]}.witness';log=folder/f'row{row["row"]}.log';cpu=max(1,min(args.per_row_cap,math.floor(remaining)-1));wall=min(remaining,cpu+2)
   attempt=job([str(search.resolve()),str(state),str(units),'1',str(row['row']),str(witness.resolve()),'8'],log,cpu,wall)
   text=log.read_text();attempt.update(profile=row,original_unit_identity_verified='EXACT_KNOWN_UNIT_DENOMINATOR_LIFT_PASS'in text,witness_exists=witness.exists());stage['attempts'].append(attempt)
   print('STAGE',number,'ROW',row['row'],attempt['status'],'WITNESS',attempt['original_unit_identity_verified'],'SECONDS',time.monotonic()-started,flush=True)
   if attempt['status']=='completed'and attempt['original_unit_identity_verified']and witness.exists():found=(row,witness);break
  remaining=args.budget-(time.monotonic()-started)
  if found and remaining>3:
   row,witness=found;out=folder/'residual.txt';cpu=max(1,math.floor(remaining)-1)
   applied=job([str(pivot.resolve()),str(state),str(out.resolve()),str(units),'1',str(witness.resolve()),str(row['row'])],folder/'pivot.log',cpu,remaining)
   stage['pivot']=applied
   if applied['status']=='completed':
    newshape,_=profile(out);assert newshape==[shape[0]-1,shape[1]-1,1]
    assert 'PASS 1'in(folder/'pivot.log').read_text();stage.update(output_shape=newshape[:2],output_sha256=sha(out),output_bytes=out.stat().st_size,output_path=str(out.relative_to(ROOT)),witness_path=str(witness.relative_to(ROOT)),witness_sha256=sha(witness))
    with out.open()as f:f.readline();f.readline();origins=list(map(int,f.readline().split()))
    priority=[priority[i]for i in origins]
    total+=1;state=out;print('EXACT_ONE_UNIT_EXPORTED',total,'SHAPE',newshape[:2],'SECONDS',time.monotonic()-started,flush=True)
  stages.append(stage);(folder/'receipt.json').write_text(json.dumps(stage,indent=2)+'\n')
  if 'output_shape'not in stage:break
  (args.output/'latest-confirmed.json').write_text(json.dumps({'additional_exact_unit_pivots':total,'shape':stage['output_shape'],'state':str(state.relative_to(ROOT)),'sha256':sha(state)},indent=2)+'\n')
 report={'id':'S61_NONSQUARE_J0_PRIMITIVE_CHAIN_20261103_V63','base_ring':'QQ[a,e,1/(a*F*Bs*C*(3e+2))]','J_not_inverted':True,
  'initial_shape':initial_shape[:2],'additional_exact_unit_pivots':total,'last_confirmed_state':str(state.relative_to(ROOT)),'last_confirmed_sha256':sha(state),
  'last_confirmed_shape':profile(state)[0][:2],'stages':stages,'export_after_every_exact_unit':True,'old_field_gcd_searches_repeated':False,'traced_Smith_repeated':False,
  'full_coker_unit_proved':profile(state)[0][0]==0,'J0_defect_support_finite_proved':profile(state)[0][0]==0,
  'initial_priority_rows':sorted(preferred),'priority_is_heuristic_only':True,
  'budget_seconds':args.budget,'per_row_cap_seconds':args.per_row_cap,'AS_limit_bytes':1800000000,'seconds':time.monotonic()-started,
  'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in [Path(__file__),source,units,searchsrc,pivotsrc]}}
 (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['additional_exact_unit_pivots','last_confirmed_shape','seconds']},indent=2))
if __name__=='__main__':main()
