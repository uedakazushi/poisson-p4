#!/usr/bin/env python3
"""ONE guarded P02 continuation from the exact P01 cursor; no arithmetic replay."""
import argparse,gc,hashlib,importlib.util,json,resource,shlex,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_K92_70_TAIL_GUARDED_CONTINUATION_20261105_v63'
P01='reports/S61_NONSQUARE_PREPOWER_K92_70_LAPLACE_GUARDED_ATTEMPT_20261105_v63'
GATE='reports/S61_NONSQUARE_PREPOWER_K92_70_TAIL_ASSOCIATION_PREFLIGHT_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def source(rel):
  p=ROOT/rel;pins[rel]=sha(p);return read(p)
 dump(out/'phase.json',dict(phase='prepare ONEP02 continuation; no scalar work'))
 r=source(GATE+'/report.json');assert r['all_passed']and r['source_only']and pins[GATE+'/report.json']=='f42fcaf9f984c0f964b4be8aac5efa078f28007581cb7fd0f1fc8cc89917069d'
 resource_gate=source(GATE+'/resource-gate.json');resume=source(GATE+'/P01-exact-resume.json');recipe=source(GATE+'/P02-operator-recipe.json');profiles=source(GATE+'/groups70-tail-source.json')['groups'];assert len(profiles)==70
 old=source(P01+'/report.json');h=source(P01+'/restored-state-check.json')['header'];assert old['safety_and_restore_checks_passed']and not old['full_scalar_complete']and h['stage']=='global_accumulate'and h['group']==0 and h['term']==3
 state=ROOT/P01/'state.bin';statehash=sha(state);assert statehash=='b979ac727295fe6689bd8cd6805189899271ce7c124103e4479728268dce9aa7';pins[str(state.relative_to(ROOT))]=statehash
 norm=source('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/normalization.json')
 layout=out/'layout.txt'
 with layout.open('w')as f:
  f.write('11 8\n')
  for rr in h['matrix_ids']:f.write(' '.join(map(str,rr))+'\n')
  f.write(' '.join(map(str,h['weight_ids']))+'\n');f.write(norm['content_lower_divisor']+' 5\n')
  for i,p in enumerate(profiles):assert p['group']==i;f.write(p['tail_common_QQfactor']+'\n')
 dump(out/'source-binding.json',dict(origin_P01_state_sha256=statehash,origin_partial_is_NOT_undone=True,exact_resume=resume,all_tail_signs_and_common_contents_from=GATE+'/groups70-tail-source.json',new_layout_sha256=sha(layout),
  source_DAG_evaluations=0,old_completed_first3_terms_NOT_replayed=True,group0_tail='4..9 AFTER term3boundary',group1to69_tail='all10 BEFORE oneP4product',all700_original_terms_retained=True,all_57source_and_tail_effects_unchanged=True))
 cpp=Path(__file__).with_suffix('.cpp');inc=ROOT/'code/s61_nonsquare_prePOWER_K92_70_tail_p02_20261105_v63.inc';p01cpp=ROOT/'code/s61_nonsquare_prePOWER_K92_70_laplace_guarded_entry_20261105_v63.cpp';oldengine=ROOT/'code/s61_nonsquare_prePOWER_K92_det11_guarded_entry_20261104_v63.cpp';parser=ROOT/'code/s61_nonsquare_prePOWER_K92_restore_p02_sparse_state_20261105_v63.py';note=ROOT/'research/notes'/f'{ID}.md';exe=out/'native-worker'
 assert cpp.read_text()==p01cpp.read_text().split('int main(int argc,char**argv)')[0]+inc.read_text()
 for p in[Path(__file__),cpp,inc,p01cpp,oldengine,parser,note]:pins[str(p.relative_to(ROOT))]=sha(p)
 external_raw=sum((ROOT/rel).stat().st_size for rel in pins);assert external_raw<10*1024*1024
 maximum_binary=135000*144+512*1216+4096;assert external_raw+maximum_binary+2*1024*1024<32*1024*1024
 dump(out/'resource-gate.json',dict(resource_gate,actual_external_read_input_RAW=external_raw,max_new_binary=maximum_binary,all_read_input_plus_new_binary_plus2MiB=external_raw+maximum_binary+2*1024*1024,owned_numeric_hardcap=135000,metadata_and_scratch_reserve=15000,
  active_poly_analytic_upper=128,P01_input_restored_terms=97053,P01_pending_product_remainder=8094365,P02_group_products=70,clock_guard_seconds=45,total_RAW_including_readinput_cap=33554432))
 # Do not retain the 3.6MB source metadata tables in the parent during native work.
 del profiles,h,r,old,norm,recipe,resume,resource_gate;gc.collect()
 cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+shlex.quote(str(cpp))+' -o '+shlex.quote(str(exe))+' $(pkg-config --cflags --libs flint)';r=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,55-(time.monotonic()-began)));(out/'compile.log').write_text(r.stdout)
 if r.returncode:
  dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,scalar_STARTED=False,compile_failure_preserved=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()}));signal.alarm(0);return
 spec=importlib.util.spec_from_file_location('p02reader',parser);restore=importlib.util.module_from_spec(spec);spec.loader.exec_module(restore)
 fixture=out/'fixture';fixture.mkdir();cmd='source /workspace/cloud-setup/activate.sh && exec '+shlex.quote(str(exe))+' --fixture '+shlex.quote(str(fixture));r=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=4);(out/'fixture.log').write_text(r.stdout);assert r.returncode==0;fs=restore.read_state(fixture/'state.bin');restore.check_fixture(fs);restore.check_pending(fs);assert fs['header']['P01_resume_pairs_completed']==3
 dump(out/'fixture-restore.json',dict(all_passed=True,version=restore.VERSION,partial_product_resume_NO_oldpair_replay=True,exact_product_and_precommit_cell_bit_checks=True,state_sha256=sha(fixture/'state.bin')));del fs
 dump(out/'phase.json',dict(phase='compile/partialproductfixture PASS; ONEP02 scalar begins'))
 cmd='source /workspace/cloud-setup/activate.sh && exec '+shlex.quote(str(exe))+' '+shlex.quote(str(state))+' '+shlex.quote(str(out))+' '+statehash+' '+shlex.quote(str(layout))+' '+sha(layout)+' '+str(external_raw)
 r=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,58-(time.monotonic()-began)));(out/'native.log').write_text(r.stdout)
 nr=read(out/'native.json');s=restore.read_state(out/'state.bin');pendingcheck=restore.check_pending(s);full=bool(nr['full_scalar_complete']);assert full==s['header']['complete'];assert s['header']['origin_P01_sha256']==statehash and s['header']['layout_sha256']==sha(layout)
 newstatehash=sha(out/'state.bin');dump(out/'restored-state-check.json',dict(all_passed=True,version=restore.VERSION,state_sha256=newstatehash,complete=full,active_polynomials=len(s['polynomials']),saved_terms=s['total_terms'],pending=s['pending'],pending_cursor_check=pendingcheck,header=s['header'],arithmetic_NOT_resumed=True));del s
 assert sha(state)==statehash
 dump(out/'phase.json',dict(phase='FINAL ALL70GROUPScomplete'if full else'FINAL exactP02prefix only; no automatic continuation'))
 actual=sum(p.stat().st_size for p in out.rglob('*')if p.is_file());assert actual+external_raw<32*1024*1024
 report=dict(id=ID,all_passed=full,frozen=True,safety_and_roundtrip_checks_passed=True,ONE_changed_P02_continuation=True,full_scalar_complete=full,actual_det11_nonzero=nr['output_nonzero']if full else None,
  groups_completed=nr['groups_completed'],globally_completed_H3_terms=nr['globally_completed_H3_terms'],tail_terms_ready=nr['tail_terms_ready'],P01_term3_boundary_reached=nr['P01_boundary_reached'],P01_pending_pairs_completed=nr['P01_pending_pairs_completed'],guard_reason=nr['reason'],restored_P01_terms=nr['restored_P01_terms'],source_DAG_evaluations=0,coefficient_updates=nr['coefficient_updates'],peak_owned_cells=nr['peak_owned_cells'],peak_owned_plus_work=nr['peak_owned_plus_work'],state_sha256=newstatehash,
  restore_command=f'.venv/bin/python code/{parser.name} --check {out.relative_to(ROOT)}/state.bin',independent_pending_cursor_checked=True,old_state_unchanged=True,no_prefix_to_det_or_rank_claim=True,no_retry_otherentry_or_global57_claim=True,
  actual_new_RAW_before_report=actual,actual_read_input_RAW=external_raw,actual_new_plus_read_RAW=actual+external_raw,input_sha256=pins,evidence_sha256={str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},seconds=time.monotonic()-began)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
