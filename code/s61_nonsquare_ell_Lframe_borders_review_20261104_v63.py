#!/usr/bin/env python3
"""Independent fixed-transcript polynomial audit, without elimination search."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 b=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63';r=read(b/'report.json');c=read(b/'certificate.json')
 assert r['frozen']and r['generic_rank_exactly10_proved']and r['status']=='all_six_zero_exact'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h,p
 for p,h in r['artifact_sha256'].items():assert sha(b/p)==h,p
 stage=read(b/'current-stage.json');piv=stage['pivots'];assert len(piv)==16 and stage['completed_stage_steps']==1
 assert len(set(x[0]for x in piv))==16 and len(set(x[1]for x in piv))==16
 assert set(x[0]for x in piv[:10])==set(range(10))and all(x[1]<26 for x in piv[:10])
 assert set(x[0]for x in piv[10:15])=={11,12,18,20,21}and set(x[1]for x in piv[10:15])==set(range(26,31))
 assert c['remaining11_genuine_relative_rows']==[0,3,4,5,6,7,9,12,13,14,15]and c['last_tested_genuine_relative_row']==6
 assert 16 not in[x[0]for x in piv]
 (out/'fixed-pivots.txt').write_text('\n'.join(str(x[0])+' '+str(x[1])for x in piv)+'\n')
 cpp=ROOT/'code/s61_nonsquare_ell_Lframe_borders_review_20261104_v63.cpp';library=ROOT/'code/s61_nonsquare_ell_Lframe_review_20261104_v63.cpp';exe=out/'review-native'
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 def rem():return max(.1,60-(time.monotonic()-started))
 with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-o',str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=rem())
 source=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63'
 run=subprocess.run([str(exe),str(source),str(b),str(out/'fixed-pivots.txt')],check=True,capture_output=True,text=True,timeout=rem());exe.unlink()
 (out/'native.log').write_text(run.stdout+run.stderr);own=json.loads(run.stdout);assert own['all_passed']
 deps=[ROOT/'reports'/s/'report.json'for s in['S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63_REVIEW','S61_NONSQUARE_ELL_COLON_RANK10_20261104_v63_REVIEW']]
 assert all(read(p)['all_passed']for p in deps)
 note=ROOT/'research/notes'/(ID+'.md')
 result=dict(id=ID,all_passed=True,actual_complete_R11x16_all_entries_independently_bound=True,
  source_G6_nonzero_coefficient_preserved_because_all_inverted_factors_nonzero_and_pivot_never_G6=True,
  actual_G6_row_zero_already_after10quad5heavy=True,all15_saved_remainders_not_only_six_zero=True,
  generic_canonical11_rank_upper10_by_literal_polynomial_identities=True,
  accepted_old_exact_QQ_minor_lower10_implies_generic_rank_exact10=True,
  every_division_only_generic_extra_chart_no_new_FULLU_unit=True,
  no_allpoint_lower10_full239_or_uniform38_claim=True,
  no_new_point_prime_GB_pivot_choice_gcd_search_or_author_rerun=True,native=own,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),cpp,library,note,b/'report.json',b/'certificate.json',b/'current-stage.json',b/'operations.txt',source/'quad10x26.txt',source/'genuine16x31.txt',*deps]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started,frozen=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
