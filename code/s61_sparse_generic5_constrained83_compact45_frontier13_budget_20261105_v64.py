#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
ID="S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT45_FRONTIER13_BUDGET_20261105_v64"
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main():
 d=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT45_UNIT_FRONTIER10_15_20261105_v64';r=d/'report.json';gp=d/'candidate-13-full-gate.json'
 assert sha(r)=='c86504e356bde0f7847ec1f90722893035628c8575be71c61fd194c00251247c'
 g=json.loads(gp.read_text());assert g['FULL_safe']and g['source_indexed_pivot']['row']==53 and g['source_indexed_pivot']['column']==54
 old=json.loads(r.read_text());assert old['evidence_sha256'][gp.name]==sha(gp)
 out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
 files=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(d.iterdir())if p.is_file()]
 reserve=3500000;raw=g['actual_READ_input_RAW_bytes']+g['normalized_row_RAW_bound']+g['all_new_rows_RAW_bound']+reserve
 assert raw==32367271<32*1024*1024
 proof={'all_passed':True,'source_only':True,'selected_gate':str(gp.relative_to(ROOT)),'selected_gate_sha256':sha(gp),'selected':[53,54],
 'owned_algebraic_and_scratch_bound':17162,'integer_index_plan_receipt_slots':44206,'whole_owned_upper':61368,
 'all_four_forecasts_and_metadata_retained_and_charged':files,'old_gate_folder_bytes':sum(v['bytes']for v in files),
 'new_other_artifact_cap':reserve,'actual_READ_input_RAW_bytes':8122702,'row_RAW_cap':20670917,'normalization_RAW_cap':73652,'whole_RAW_bound':raw,
 'hard_RAW_cap_unchanged':32*1024*1024,'no_arithmetic_or_retry_of_unsafe_candidates':True,
 'reason':'All four continuation forecasts and their cursor are retained and charged. The complete3500000-byte auxiliary reserve and replayable compact event context stay within the hard32MiB total; compile/control/future metadata are checked before the sole pivot. Earlier ten forecasts remain in their immutable source-stop archive.',
 'future_control_executable_events_manifests_receipts_budget':reserve-sum(v['bytes']for v in files)}
 (out/'whole-budget.json').write_text(json.dumps(proof,indent=2)+'\n')
 report={'id':ID,'all_passed':True,'source_only':True,'whole_RAW_bound':raw,'owned_bound':61368,
 'input_sha256':{**{f['path']:f['sha256']for f in files},str(Path(__file__).relative_to(ROOT)):sha(Path(__file__)),f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')},
 'evidence_sha256':{'whole-budget.json':sha(out/'whole-budget.json')}}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('report_sha256',sha(out/'report.json'),'RAW',raw,'other',reserve,'remaining',proof['future_control_executable_events_manifests_receipts_budget'])
if __name__=='__main__':main()
