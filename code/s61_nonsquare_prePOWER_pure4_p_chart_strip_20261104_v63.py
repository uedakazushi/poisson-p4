#!/usr/bin/env python3
"""ONE at-most-two exact tests of existing H5 chart factor p, not factorization."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_PURE4_P_CHART_STRIP_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);assert not(out/'report.json').exists();phase=read(out/'phase.json')
 bp=ROOT/'reports/S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63';assert sha(bp/'report.json')==phase['source_R_named_strip_report_SHA256'];d=read(bp/'native.json');assert d['all_passed'];R=d['remaining_polynomial'];hp=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63/certificate.json';p=read(hp)['p'];assert len(R)==434 and len(p)==21 and max(m[1]for m,c in R)==16 and max(m[1]for m,c in p)==4
 forecast=dict(input_terms=455,source_R_degree20=True,p_total_degree6=True,p_e_degree4=True,two_tests_are_explicit_budget_not_maximality_bound=True,initial_wrong_p_e6_assumption_preserved_before_native_job=True,conservative_total_live_terms_upper=100000,strict150000_dynamic_guard=True,no_otherfactor_or_minor=True)
 (out/'preflight.json').write_text(json.dumps(forecast,indent=2)+'\n');text=['455']
 for a in[R,p]:text.append(str(len(a)));text.extend(' '.join(map(str,m))+' '+c for m,c in a)
 ip=out/'native-input.txt';ip.write_text('\n'.join(text)+'\n');cpp=Path(__file__).with_suffix('.cpp');exe=out/'p-chart-native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,53-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:run=subprocess.run([str(exe),str(ip),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,56-(time.monotonic()-start)))
 n=read(out/'native.json')if(out/'native.json').exists()else dict(all_passed=False,returncode=run.returncode);ok=run.returncode==0 and n['all_passed']
 if ok:assert n['literal_R_equals_p_power_times_remaining_verified']and n['peak_live_terms_with_input_reserve']<=150000
 paths=[Path(__file__),cpp,ROOT/'code/s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md',bp/'report.json',bp/'native.json',hp]
 result=dict(id=ID,all_passed=ok,frozen=True,independent_review_pending=True,exact_chart_p_power=n.get('exact_chart_p_power'),p_power_maximality_proved=n.get('p_power_maximality_proved'),remaining_unknown_terms=n.get('remaining_terms'),peak_live_terms_with_input_reserve=n.get('peak_live_terms_with_input_reserve'),no_global_unit_full9_full57_or_complement_removal_claim=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
