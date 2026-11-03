#!/usr/bin/env python3
"""Bounded new primitive-entry unit-span search on saved J0 bivariate23x57.

This uses the frozenV9 method, not a repeated rational-function gcd search.
Per-entry unit removal is used only to DISCOVER a witness.  Native code then
clears all known-unit denominators and verifies the ORIGINAL polynomial row.
"""
import argparse,hashlib,json,math,os,resource,shlex,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(command,log,cpu,wall,limit=1800000000):
 def bounds():resource.setrlimit(resource.RLIMIT_AS,(limit,limit));resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu+1))
 start=time.monotonic()
 with log.open('w')as f:
  proc=subprocess.Popen(command,stdout=f,stderr=subprocess.STDOUT,preexec_fn=bounds,start_new_session=True)
  try:proc.wait(timeout=wall);status='completed'if proc.returncode==0 else'refused'if proc.returncode==2 else'exit_'+str(proc.returncode)
  except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait();status='wall_timeout'
 return {'status':status,'returncode':proc.returncode,'seconds':time.monotonic()-start,'cpu_limit':cpu,'wall_limit':wall,'AS_limit_bytes':limit}

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--budget',default=300,type=int);ap.add_argument('--max-rows',default=4,type=int);args=ap.parse_args();assert 0<args.budget<=300;args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 source=ROOT/'work/nonsquare-J0-pid-v63/J0-residual.txt';units=ROOT/'work/nonsquare-J0-pid-v63/units.txt'
 frozen=ROOT/'reports/S61_NONSQUARE_J0_PID_20261103_v63/report.json';old=json.loads(frozen.read_text());assert sha(source)==old['omitted_easy_to_regenerate_outputs'][str(source.relative_to(ROOT))]['sha256']
 profiles=[]
 with source.open()as f:
  assert list(map(int,f.readline().split()))==[23,57,13];assert len(f.readline().split())==13;origins=list(map(int,f.readline().split()));assert len(origins)==23
  for row in range(23):
   counts=[];maxdeg=0
   for _ in range(57):
    tokens=f.readline().split();nt=int(tokens[0]);assert len(tokens)==1+4*nt;counts.append(nt)
    for j in range(nt):
     ex=list(map(int,tokens[2+4*j:5+4*j]));assert ex[2]==0;maxdeg=max(maxdeg,sum(ex))
   profiles.append({'row':row,'previous_origin':origins[row],'terms':sum(counts),'nonzero_entries':sum(n>0 for n in counts),'max_degree':maxdeg})
 selected=sorted(profiles,key=lambda r:(r['terms'],r['max_degree'],r['row']))[:args.max_rows]
 original=ROOT/'work/nonsquare-v62/primitive-unit-span-v9/primitive-unit-span-optimized.cpp';snapshot=args.output/'primitive-search.cpp';snapshot.write_bytes(original.read_bytes())
 flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
 exe=args.output/'search';subprocess.run(['c++','-std=c++17','-O3',str(snapshot),*flags,'-lgmp','-o',str(exe)],check=True)
 attempts=[];found=None
 for profile in selected:
  remaining=args.budget-(time.monotonic()-start)
  if remaining<3:break
  row=profile['row'];witness=args.output/f'row{row}.witness';log=args.output/f'row{row}.log';cpu=max(1,math.floor(remaining)-1)
  result=run([str(exe.resolve()),str(source),str(units),'1',str(row),str(witness.resolve()),'8'],log,cpu,remaining)
  text=log.read_text();result.update(profile=profile,witness_file_exists=witness.exists(),original_polynomial_identity_pass='EXACT_KNOWN_UNIT_DENOMINATOR_LIFT_PASS'in text)
  attempts.append(result);print('PRIMITIVE_ROW',row,result['status'],'WITNESS',result['original_polynomial_identity_pass'],'SECONDS',time.monotonic()-start,flush=True)
  if result['status']=='completed'and result['original_polynomial_identity_pass']and witness.exists():found={'row':row,'witness':witness,'search_result':result};break
  if result['status']not in ['refused']:break
 pivot=None
 if found:
  remaining=args.budget-(time.monotonic()-start)
  if remaining>=3:
   native=ROOT/'code/polynomial_unit_span_pivots_20261103_v62.cpp';pivotexe=args.output/'pivot';subprocess.run(['c++','-std=c++17','-O3',str(native),*flags,'-lgmp','-o',str(pivotexe)],check=True)
   result=args.output/'residual.txt';remaining=args.budget-(time.monotonic()-start);cpu=max(1,math.floor(remaining)-1)
   pivot=run([str(pivotexe.resolve()),str(source),str(result.resolve()),str(units),'1','colclean',str(found['witness'].resolve()),str(found['row'])],args.output/'pivot.log',cpu,remaining)
   if pivot['status']=='completed':
    with result.open()as f:shape=list(map(int,f.readline().split()))
    assert shape[2]>=1;pivot.update(residual_shape=shape[:2],additional_known_unit_pivots=shape[2],residual_sha256=sha(result))
 inputs=[Path(__file__),source,units,original,frozen,ROOT/'code/polynomial_unit_span_pivots_20261103_v62.cpp']
 report={'id':'S61_NONSQUARE_J0_PRIMITIVE_SPAN_20261103_V63','locus':'J=0,a*F*Bs*C*(3e+2)!=0','base_ring':'QQ[a,e,1/(a*F*Bs*C*(3e+2))]',
  'only_declared_nonconstant_units_used':True,'J_not_inverted':True,'old_field_gcd_searches_repeated':False,'traced_Smith_repeated':False,
  'starting_shape':[23,57],'profiles':profiles,'candidate_rows':selected,'search_attempts':attempts,
  'new_original_row_unit_witness_verified':found is not None,'new_witness_row':found['row']if found else None,'pivot_replay':pivot,
  'full_coker_unit_or_allpoint_rank57_claimed':bool(pivot and pivot.get('residual_shape',[1])[0]==0),
  'budget_seconds':args.budget,'seconds':time.monotonic()-start,
  'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in inputs}}
 (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['new_original_row_unit_witness_verified','pivot_replay','seconds']},indent=2))
if __name__=='__main__':main()
