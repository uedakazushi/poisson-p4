#!/usr/bin/env python3
"""Resume saved J0 Bezout witnesses, then use only the remaining time budget.

No first-stage gcd searches or traced Smith computations are repeated.
Exact replay states are exported after selected saved steps.  All output is
new; the frozen first checkpoint and its input files are read-only.
"""
import argparse,hashlib,json,math,os,re,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WRAPPER=ROOT/'work/nonsquare-v62/singular-runtime/run-singular-optimized-v6.sh'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def run_job(command,logfile,cpu,wall):
 def limits():resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu+1))
 started=time.monotonic()
 with logfile.open('w')as log:
  process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits,start_new_session=True)
  try:process.wait(timeout=wall);status='completed'if process.returncode==0 else'exit_'+str(process.returncode)
  except subprocess.TimeoutExpired:
   os.killpg(process.pid,signal.SIGKILL);process.wait();status='wall_timeout'
 return {'status':status,'seconds':time.monotonic()-started,'cpu_limit':cpu,'wall_limit':wall,'returncode':process.returncode}

def export_code(path,step):
 rows=36-step
 return f'''if(remainingRows!={rows}){{ERROR("saved step row count failed");}}
print("SAVED_STEP_{step}_VERIFIED_ROWS_{rows}");
write("{path}","ring R=(0,e),(a),dp;matrix M[{rows}][70];");
for(int exportI=1;exportI<={rows};exportI++){{for(int exportJ=1;exportJ<=70;exportJ++){{if(Mcur[exportI,exportJ]!=0){{write("{path}","M["+string(exportI)+","+string(exportJ)+"]="+string(Mcur[exportI,exportJ])+";");}}}}}}
write("{path}","poly KnownUnit=a*(500*a+90*e^2+3*e-24)*(300*a*e+54*e^3+9*e^2-24*e-4)*(100*a+18*e^2-9*e-2)*(3*e+2);");
print("SAVED_STATE_{step}_EXPORT_COMPLETE");
'''

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
 ap.add_argument('--checkpoint',type=Path,default=ROOT/'work/nonsquare-J0-pid-v63/QQe-a-rowgcd-verified')
 ap.add_argument('--budget',type=int,default=600);ap.add_argument('--replay-segment-cap',type=int,default=300)
 args=ap.parse_args();assert 0<args.budget<=600 and 0<args.replay_segment_cap<=600
 args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic();cp=args.checkpoint.resolve()
 proof=cp/'proof.sing';library=cp/'replay-library.sing';matrix=cp.parent/'QQe-a-smith/matrix.sing'
 text=proof.read_text();matches=list(re.finditer(r'(?m)^matrix T(\d+)\[70\]\[1\];\n',text));assert len(matches)==12
 segments={}
 for i,match in enumerate(matches):
  n=int(match[1]);end=matches[i+1].start()if i+1<len(matches)else len(text)
  body=text[match.start():end];assert f'list step{n}=ReplayStep('in body;segments[n]=body
 stage=0;best=0;states={};runs=[]
 while best<12:
  remaining=args.budget-(time.monotonic()-started)
  if remaining<5:break
  stage+=1;folder=args.output/f'replay-{stage}';folder.mkdir();input_matrix=matrix if best==0 else states[best]
  script=f'execute(read("{input_matrix}"));\nexecute(read("{library}"));matrix Mcur=M;int remainingRows={36-best};\n'
  targets={}
  for n in range(best+1,13):
   script+=segments[n]
   script+=f'print("SAVED_COFACTOR_STEP_{n}_APPLIED_ROWS_"+string(remainingRows));\n'
   if n in [4,8,11,12]:
    path=(folder/f'state-after-{n}.sing').resolve();targets[n]=path;script+=export_code(path,n)
  script+='print("EXACT_SAVED_12_STEP_PREFIX_REPLAY_PASS");quit;\n';runfile=folder/'run.sing';runfile.write_text(script)
  cpu=max(1,min(args.replay_segment_cap,int(remaining)-1));wall=min(remaining,cpu+5)
  result=run_job([str(WRAPPER),'--cpus=0','-q',str(runfile.resolve())],folder/'solver.log',cpu,wall)
  log=(folder/'solver.log').read_text();error=bool(re.search(r'(?m)^\s*\?',log));result['arithmetic_error']=error
  result['start_after_step']=best;result['applied_steps']=[int(m)for m in re.findall(r'SAVED_COFACTOR_STEP_(\d+)_APPLIED',log)]
  completed=[n for n,path in targets.items()if f'SAVED_STATE_{n}_EXPORT_COMPLETE'in log and path.exists()]
  if not error:
   for n in completed:states[n]=targets[n]
  result['exported_steps']=completed;runs.append(result)
  (folder/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
  previous=best
  if states:best=max(states)
  print('SAVED_REPLAY_STAGE',stage,'STATUS',result['status'],'BEST_EXPORT',best,'SECONDS',time.monotonic()-started,flush=True)
  if error or best<=previous:break
 continuation=None
 if best==12:
  remaining=args.budget-(time.monotonic()-started)
  if remaining>=5:
   out=args.output/'continued-row-gcd';cpu=max(1,math.floor(remaining)-2)
   command=['python',str(ROOT/'code/s61_nonsquare_J0_pid_rowgcd_20261103_v63.py'),'--matrix-input',str(states[12]),'--output',str(out.resolve()),'--timeout',str(cpu)]
   # The child producer itself enforces the same AS/CPU cap on Singular.
   continuation=run_job(command,args.output/'continued-row-gcd.stdout',max(1,int(remaining)),remaining)
   if (out/'receipt.json').exists():continuation['mathematical_receipt']=json.loads((out/'receipt.json').read_text())
 report={'id':'S61_NONSQUARE_J0_PID_RESUME_20261103_V63','budget_seconds':args.budget,'AS_limit_bytes':2200000000,
  'saved_prefix_searches_repeated':False,'traced_Smith_repeated':False,'replay_runs':runs,'best_exact_export_step':best,
  'states':{str(n):{'path':str(path.relative_to(ROOT)),'bytes':path.stat().st_size,'sha256':sha(path)}for n,path in states.items()},
  'saved_12_step_replay_confirmed':best==12,'continuation':continuation,'seconds':time.monotonic()-started,
  'input_sha256':{str(f.resolve().relative_to(ROOT)):sha(f)for f in [Path(__file__),proof,library,matrix,ROOT/'code/s61_nonsquare_J0_pid_rowgcd_20261103_v63.py']}}
 (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
