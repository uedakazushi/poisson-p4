#!/usr/bin/env python3
"""ONE finite whole-axis ZZ[a] Bareiss with exact patchcertificates."""
import argparse,gc,hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_BAREISS_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_BAREISS_SOURCE_20261105_v63';r=src(base+'/report.json');assert r['all_passed']and pins[base+'/report.json']=='f2331a0663ca70f068aa20a0df35a854cab795b848fb3b307a6ebc951a3505b7';g=src(base+'/whole-worker-source-plan.json');proof=src(base+'/integer-polynomial-source-proof.json');receipt=src('reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_BAREISS_WORKER_GATE_20261105_v63/report.json');assert receipt['all_passed'];axis=src('reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_SOURCE_V2_20261105_v63/whole-axis-literal-source.json')
 cpp=ROOT/'code/s61_nonsquare_K92_axis_integer_bareiss_20261105_v63.cpp';pins[str(cpp.relative_to(ROOT))]=sha(cpp);work=ROOT/'work/d-zero-v63/k92-axis-integer-bareiss';work.mkdir(parents=True,exist_ok=True);exe=work/'axis-bareiss';cp=subprocess.run(['g++','-std=c++17','-O2',str(cpp),'-lflint','-lgmp','-o',str(exe)],capture_output=True,text=True,timeout=25);(out/'compile.stdout').write_text(cp.stdout);(out/'compile.stderr').write_text(cp.stderr)
 if cp.returncode:dump(out/'report.json',dict(id=ID,all_passed=False,compile_only_failure=True,no_arithmetic_executed=True,input_sha256=pins));return
 hp=subprocess.run([str(exe),'--help'],capture_output=True,text=True,timeout=3);assert hp.returncode==0;(out/'native-help.txt').write_text(hp.stdout);deadline=min(45,56-(time.monotonic()-start));assert deadline>1
 infile=out/'native-input.txt'
 with infile.open('w')as f:
  f.write(f'K92AXIS001 11 88 864 1738 27812 {deadline}\n'+proof['every_minor_L1_upper']+'\n')
  for e in sorted(axis['axis_coefficients'],key=lambda e:(e['row'],e['column'])):
   pp=e['axis_integer_polynomial'];f.write(f"{e['row']} {e['column']} {len(pp)} "+' '.join(f'{k} {v}'for k,v in pp)+'\n')
 denominators=axis['axis_row_denominators'];del axis;gc.collect();dump(out/'phase.json',dict(phase='immutable_source_and_compile_precheck_complete; ONE_native_integer_axis_job',code_sha=pins[str(cpp.relative_to(ROOT))],input_sha256=sha(infile),guard_seconds=deadline))
 def limits():resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));resource.setrlimit(resource.RLIMIT_CPU,(55,55))
 rr=subprocess.run([str(exe),str(infile),str(out)],capture_output=True,text=True,timeout=min(deadline+5,59-(time.monotonic()-start)),preexec_fn=limits);(out/'native.stdout').write_text(rr.stdout);(out/'native.stderr').write_text(rr.stderr)
 if rr.returncode:
  dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,native_partial_or_failure=True,all_existing_state_and_patches_preserved=True,no_rank_or_determinant_result_claim=True,input_sha256=pins));return
 n=read(out/'native.json');assert n['all_passed']and 6<=n['rank_over_QQ_a']<=11 and n['whole_owned_with_input_and_metadata']<=g['all_owned_with_metadata15000']
 D=int(math.prod(map(int,denominators)));det={int(k):int(v)for k,v in n['determinant']};assert bool(det)==(n['rank_over_QQ_a']==11)
 original=[[k+5,str(Q(v,D))]for k,v in sorted(det.items())];assert all(5<=k<=72 for k,v in original)
 dump(out/'axis-rank-and-determinant-normalization.json',dict(rank_over_QQ_a=n['rank_over_QQ_a'],axis_detM_e0d0_QQ_coefficients=original,exact_identity='detM(a,0,0)=a^5*detCaxis(a)/prod(axis_rowLCD)',
  polynomial_determinant_sign_and_all_exact_divisions_bound=True,wholeaxis_fixed_minor_rows=n['row_permutation'][:n['rank_over_QQ_a']],wholeaxis_fixed_minor_columns=n['column_permutation'][:n['rank_over_QQ_a']],
  remaining_normal_dimension_if_detZERO=11-n['rank_over_QQ_a'],d0_OUTSIDE_FULLU=True,NO_all_a_pointwise_or_fullU_cubic57_claim=True,independent_peer_required=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new<g['total_possible_new_RAW_upper']and new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,ONE_finite_integer_polynomial_axis_job=True,rank_over_QQ_a=n['rank_over_QQ_a'],axis_determinant_NONZERO=bool(det),all_exact_integer_polynomial_divisions=True,
  integer_products=n['integer_products'],observed_raw_integer_bits=n['observed_raw_integer_bits'],whole_owned_slots=n['whole_owned_with_input_and_metadata'],no_point_prime_fieldRF_gaussian_or_coefwindow=True,
  no_fullU57_or_all_a_pointwise_rank_claim=True,independent_peer_pending=True,total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
