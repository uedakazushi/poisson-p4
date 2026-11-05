#!/usr/bin/env python3
"""ONE changed-association guarded scalar attempt; no source-DAG reevaluation."""
import argparse,hashlib,importlib.util,json,resource,shlex,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_K92_70_LAPLACE_GUARDED_ATTEMPT_20261105_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def source(rel):
  p=ROOT/rel;pins[rel]=sha(p);return read(p)
 dump(out/'phase.json',dict(phase='compile/fixture preparation; no scalar execution'))
 old='reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63';gate='reports/S61_NONSQUARE_PREPOWER_K92_70_LAPLACE_SOURCE_PREFLIGHT_20261104_v63'
 gr=source(gate+'/report.json');assert gr['source_only']and gr['all_passed']and pins[gate+'/report.json']=='a5414501c48a0cb39f685010b226222a8229f7d57aa327ed6df8752a94ae3ca6'
 binding=source(old+'/source-binding.json');previous=source(old+'/report.json');assert previous['safety_and_restore_checks_passed']and not previous['full_scalar_complete']
 norm=source('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/normalization.json');oldstate=ROOT/old/'state.bin';oldsha=sha(oldstate);assert oldsha=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02';pins[str(oldstate.relative_to(ROOT))]=oldsha
 layout=out/'layout.txt'
 with layout.open('w')as f:
  f.write('11 8\n')
  for rr in binding['actual_pure_matrix_roots']:f.write(' '.join(map(str,rr))+'\n')
  f.write(norm['content_lower_divisor']+' '+str(norm['known_output_a_power'])+'\n')
 recipe=source(gate+'/operator-recipe.json');resources=source(gate+'/resource-gate.json');assert resources['fixed_input_actual_terms']==25266
 dump(out/'source-binding.json',dict(pinned_old_state_sha256=oldsha,source575_already_evaluated572=True,NO_source_DAG_replay=True,matrix11x8_roots=binding['actual_pure_matrix_roots'],normative_H3_signs=binding['normative_sign_schedule'],association='P4four-row Laplace70groups; each bottom4 then H3*bottom4 then P4*weightedbottom streamed to global',
  layout_sha256=sha(layout),old_metadata_provenance=old+'/source-binding.json',all_tail_effects_and_current15_labels_retained=True))
 dump(out/'resource-gate.json',dict(resources,active_poly_analytic_upper=128,metadata_scratch_budget=dict(QQdeterminant_factors=160,polynomial_factor_and_header_records=2048,source_ID_and_group_arrays=2000,pending_scalar_scratch=100,loop_record_and_serialization_scratch=1000,remaining_reserve=9692),
  no_runtime_completion_promise=True,no_other_candidate=True,coefficient_bit_cap=1024,QQfactor_bit_cap=4096,pending_bit_cap=16384,clock_seconds=45,total_RAW_hardcap=33554432))
 cpp=Path(__file__).with_suffix('.cpp');oldengine=ROOT/'code/s61_nonsquare_prePOWER_K92_det11_guarded_entry_20261104_v63.cpp';parser=ROOT/'code/s61_nonsquare_prePOWER_K92_restore70_sparse_state_20261105_v63.py';note=ROOT/'research/notes'/f'{ID}.md';exe=out/'native-worker'
 for p in[Path(__file__),cpp,oldengine,parser,note]:pins[str(p.relative_to(ROOT))]=sha(p)
 cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+shlex.quote(str(cpp))+' -o '+shlex.quote(str(exe))+' $(pkg-config --cflags --libs flint)';r=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,55-(time.monotonic()-began)));(out/'compile.log').write_text(r.stdout)
 if r.returncode:
  dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,arithmetic_STARTED=False,compile_failure_preserved=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()}));signal.alarm(0);return
 spec=importlib.util.spec_from_file_location('restore70',parser);restore=importlib.util.module_from_spec(spec);spec.loader.exec_module(restore)
 fixture=out/'fixture';fixture.mkdir();cmd='source /workspace/cloud-setup/activate.sh && exec '+shlex.quote(str(exe))+' --fixture '+shlex.quote(str(fixture));r=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=4);(out/'fixture.log').write_text(r.stdout);assert r.returncode==0;s=restore.read_state(fixture/'state.bin');restore.check_fixture(s);dump(out/'fixture-restore.json',dict(all_passed=True,exact_product_and_precommit_cell_bit_guards=True,state_sha256=sha(fixture/'state.bin'),parser_version=restore.VERSION));del s
 dump(out/'phase.json',dict(phase='compile/fixture PASS before ONE70group scalar attempt'))
 cmd='source /workspace/cloud-setup/activate.sh && exec '+shlex.quote(str(exe))+' '+shlex.quote(str(oldstate))+' '+shlex.quote(str(out))+' '+oldsha+' '+shlex.quote(str(layout))+' '+sha(layout)
 r=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,58-(time.monotonic()-began)));(out/'native.log').write_text(r.stdout)
 nr=read(out/'native.json');state=restore.read_state(out/'state.bin');h=state['header'];assert h['old_state_sha256']==oldsha and h['layout_sha256']==sha(layout)and h['source_DAG_evaluations']==0;full=bool(nr['full_scalar_complete']);assert full==h['complete']
 statehash=sha(out/'state.bin');dump(out/'restored-state-check.json',dict(all_passed=True,version=restore.VERSION,state_sha256=statehash,complete=full,active_polynomials=len(state['polynomials']),saved_terms=state['total_terms'],pending=state['pending'],header=h,arithmetic_NOT_resumed=True));del state
 assert sha(oldstate)==oldsha
 dump(out/'phase.json',dict(phase='FINAL ALL70GROUPScomplete'if full else'FINAL precise partial-prefix only; NO automatic retry'))
 total=sum(p.stat().st_size for p in out.rglob('*')if p.is_file());assert total<32*1024*1024
 report=dict(id=ID,all_passed=full,frozen=True,safety_and_restore_checks_passed=True,one_changed_association_attempt=True,full_scalar_complete=full,actual_det11_nonzero=nr['output_nonzero']if full else None,groups_completed=nr['groups_completed'],H3_terms_completed=nr['H3_terms_completed'],guard_reason=nr['reason'],reused_source_terms=nr['reused_source_terms'],source_DAG_evaluations=0,old_DP_terms_NOT_allocated=nr['old_terms_skipped_without_allocation'],coefficient_updates=nr['coefficient_updates'],peak_owned_cells=nr['peak_owned_cells'],peak_owned_plus_work=nr['peak_owned_plus_work'],state_sha256=statehash,
  restore_command=f'.venv/bin/python code/{parser.name} --check {out.relative_to(ROOT)}/state.bin',restore_NOT_arithmetic_continuation=True,no_partial_sum_to_nonzero_inference=True,all70x10_required=True,no_global57_or_otherentry_claim=True,total_RAW_bytes_before_report=total,
  input_sha256=pins,evidence_sha256={str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},seconds=time.monotonic()-began)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
