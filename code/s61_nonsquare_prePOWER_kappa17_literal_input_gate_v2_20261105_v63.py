#!/usr/bin/env python3
"""Copy only literal source atoms; repair phase ownership without new values."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_KAPPA17_LITERAL_INPUT_GATE_V2_20261105_v63'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def count(x):
    if isinstance(x,list):
        if len(x)==2 and isinstance(x[0],list) and len(x[0])==3 and all(type(t)==int for t in x[0]) and isinstance(x[1],str):
            try:Fraction(x[1]);return 1
            except ValueError:pass
        return sum(map(count,x))
    if isinstance(x,dict):return sum(map(count,x.values()))
    return 0
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False);pins={}
    def load(rel):
        p=ROOT/rel;pins[rel]=sha(p);return json.loads(p.read_text())
    base='reports/S61_NONSQUARE_PREPOWER_KAPPA17_SOURCE_GATE_20261105_v63/'
    br=load(base+'report.json');assert br['all_passed']
    plan=load(base+'literal-source-plan.json');gate=load(base+'source-support-and-worker-gate.json')
    paths=sorted({r['source']['path'] for r in plan['source_atom_references']})
    docs={p:load(p) for p in paths};counts={p:count(x) for p,x in docs.items()}
    assert sorted(counts.values())==[40,101604]
    atoms=[]
    for ref in plan['source_atom_references']:
        p=ref['source']['path'];x=docs[p]
        for k in ref['source']['json_pointer'].strip('/').split('/'):x=x[int(k)] if isinstance(x,list) else x[k]
        assert pins[p]==ref['source_reference_sha256'];assert len(x)==2
        # Exact original arrays, no coefficient rescaling or polynomial operation.
        atoms.append(dict(node=ref['node'],original_num_den=x,source_reference=ref['source'],source_sha256=pins[p]))
    copied=sum(count(x['original_num_den']) for x in atoms);assert copied==2858
    prep=sum(counts.values())+copied+35000;assert prep==139502 and prep<150000
    dump(out/'literal-atoms20.json',dict(atom_payloads=atoms,source_order=plan['source_root_order'],
      source49_definitions=[dict(node=n['source_node'],definition=n['definition']) for n in plan['source_nodes']],
      kappa_sources=plan['kappa_coefficients_in_exact_order'],
      no_new_polynomial_or_coefficient_values=True,all_original_num_den_records_preserved=True))
    supplement=dict(old_gate_reference=base+'report.json',old_files_unchanged=True,
      metadata_history_resolution='75390 bounds convolution/publication only; old111610 source-preparation estimate counted matrix25x32 only. Full literal source preparation now separated and counted.',
      source_prepare_all_original_coefficient_records=counts,copied_literal_coefficient_records=copied,
      source_prepare_ALL_metadata_reserve=35000,source_prepare_whole_owned_upper=prep,
      future_worker_input='literal-atoms20.json plus existing four kappa polynomials; do NOT load original checkpoint or Q2/pair JSON',
      future_worker_convolution_owned_upper=gate['whole_worker_owned_upper'],
      future_worker_owned_upper_including20_denominator_records=gate['whole_worker_owned_upper']+20,
      future_worker_integer_intermediate_bits_upper=2*(gate['tracked_integer_coefficient_bits_upper']+gate['tracked_denominator_bits_upper'])+16,
      output_support_upper=gate['cleared_kappa17_output']['support'],no_B_denominator=True,
      scalar_job_RUN=False,separate_parent_GO_required=True)
    dump(out/'phase-ownership-supplement.json',supplement)
    code=Path(__file__);pins[str(code.relative_to(ROOT))]=sha(code);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
    total=sum(p.stat().st_size for p in out.iterdir());assert total<2*1024*1024
    rep=dict(id=ID,all_passed=True,frozen=True,source_only_literal_copy=True,no_new_scalar_values=True,
      old_guard_phase_distinction_explicit=True,source_prepare_owned_upper=prep,
      future_convolution_owned_upper=gate['whole_worker_owned_upper']+20,temporary_integer_bits_upper=supplement['future_worker_integer_intermediate_bits_upper'],
      scalar_worker_not_run=True,requires_separate_parent_GO=True,total_new_RAW_before_report=total,
      input_sha256=pins,evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start)
    dump(out/'report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items() if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
