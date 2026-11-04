#!/usr/bin/env python3
"""ONE combinedbudget fourremainingclosedlinear Qcomponent tests; no othercomponent."""
import hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_Q_LINEAR_BATCH_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);cpp=ROOT/'code/s61_nonsquare_raw_H_B_Q_first_component_unit_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-Q-linear-batch';worker=ROOT/'code/s61_nonsquare_raw_H_B_Q_linear_component_worker_20261104_v63.py'
 with(out/'compile.log').open('w')as f:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=50)
 exesha=sha(exe);results=[]
 for idx in[1,2,3,4]:
  child=out/('component-'+str(idx));log=out/('component-'+str(idx)+'.log')
  with log.open('w')as f:subprocess.run([str(ROOT/'.venv/bin/python'),'-u',str(worker),'--factor-index',str(idx),'--native-executable',str(exe),'--output',str(child)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
  r=json.loads((child/'report.json').read_text());results.append({'factor_index':idx,'dimension_N':r['finite_QQ_algebra_dimension_N'],'matrix_shape':r['matrix_shape'],'matrix_nonzero':r['matrix_nonzero'],'cofactor_witness':r['selected_Q_linear_component_FULLU_empty_proved'],'child_report_sha256':sha(child/'report.json')});assert sum(x['matrix_nonzero']for x in results)<150000
 exe.unlink();allunit=all(x['cofactor_witness']for x in results);paths=[Path(__file__),worker,cpp,ROOT/'research/notes'/(ID+'.md')];report={'id':ID,'results':results,'all_four_remaining_linear_components_FULLU_empty_proved':allunit,'all_five_linear_Q_components_closed_including_previous':allunit,'untested_factor_degrees':[6,6,9,11],'no_auto_other_component_fiber_point_prime_GB_resultant_jobs':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},'compiled_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.rglob('*')if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':main()
