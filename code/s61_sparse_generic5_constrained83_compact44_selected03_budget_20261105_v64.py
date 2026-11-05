#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
ID="S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT44_SELECTED03_BUDGET_20261105_v64"
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main():
 d=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT44_UNIT_GATE_20261105_v64';r=d/'report.json';gp=d/'candidate-03-full-gate.json'
 assert sha(r)=='b2938534812cdaee4879f689bce730a5ec8469a34ed0a3b6e36e73d6939c5daf'
 g=json.loads(gp.read_text());assert g['FULL_safe']and g['source_indexed_pivot']['row']==32 and g['source_indexed_pivot']['column']==68
 old=json.loads(r.read_text());assert old['evidence_sha256'][gp.name]==sha(gp)
 out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
 files=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(d.iterdir())if p.is_file()]
 reserve=3500000;raw=g['actual_READ_input_RAW_bytes']+g['normalized_row_RAW_bound']+g['all_new_rows_RAW_bound']+reserve
 assert raw==33451152<32*1024*1024
 proof={'all_passed':True,'source_only':True,'selected_gate':str(gp.relative_to(ROOT)),'selected_gate_sha256':sha(gp),'selected':[32,68],
 'owned_algebraic_and_scratch_bound':14734,'integer_index_plan_receipt_slots':52919,'whole_owned_upper':67653,
 'all_four_forecasts_and_metadata_retained_and_charged':files,'old_gate_folder_bytes':sum(v['bytes']for v in files),
 'new_other_artifact_cap':reserve,'actual_READ_input_RAW_bytes':7432267,'row_RAW_cap':22473147,'normalization_RAW_cap':45738,'whole_RAW_bound':raw,
 'hard_RAW_cap_unchanged':32*1024*1024,'no_arithmetic_or_retry_of_unsafe_candidates':True,
 'reason':'The old 3MiB auxiliary reserve did not cover all four saved forecasts plus worker outputs. This new complete budget retains every file and uses3500000 auxiliary bytes within the same hard32MiB cap and replayable compact event records with a shared context. Compile/control/native receipts must pass actual per-file and total guards before the sole pivot.',
 'future_control_executable_events_manifests_receipts_budget':reserve-sum(v['bytes']for v in files)}
 (out/'whole-budget.json').write_text(json.dumps(proof,indent=2)+'\n')
 report={'id':ID,'all_passed':True,'source_only':True,'whole_RAW_bound':raw,'owned_bound':67653,
 'input_sha256':{**{f['path']:f['sha256']for f in files},str(Path(__file__).relative_to(ROOT)):sha(Path(__file__)),f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')},
 'evidence_sha256':{'whole-budget.json':sha(out/'whole-budget.json')}}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('report_sha256',sha(out/'report.json'),'RAW',raw,'other',reserve,'remaining',proof['future_control_executable_events_manifests_receipts_budget'])
if __name__=='__main__':main()
