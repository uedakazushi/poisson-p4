#!/usr/bin/env python3
"""Freeze completed source forecasts after the 45-second source-only stop."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT45_SOURCE_STOP_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    sd=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT45_UNIT_GATE_20261105_v64';short=sd/'metadata-shortlist.json';s=json.loads(short.read_text());assert len(s)==16
    fs=sorted(sd.glob('candidate-*-full-gate.json'));assert len(fs)==10 and not(sd/'report.json').exists()
    values=[]
    for k,p in enumerate(fs):
        g=json.loads(p.read_text());assert g['ordinal']==k and not g['FULL_safe']and g['RAW_TOTAL_bound']>32*1024*1024 and g['whole_owned_scalar_and_scratch_upper']<150000
        values.append({'ordinal':k,'pivot':[g['source_indexed_pivot']['row'],g['source_indexed_pivot']['column']],
                       'RAW_upper':g['RAW_TOTAL_bound'],'owned_upper':g['whole_owned_scalar_and_scratch_upper'],'file':str(p.relative_to(ROOT)),'sha256':sha(p)})
    cursor={'completed_exact_forecasts':values,'next_shortlist_ordinal':10,'next_candidate_partial_not_saved':{'row':s[10]['row'],'column':s[10]['column']},
            'remaining_shortlist_ordinals':list(range(10,16)),'no_field_or_coefficient_elimination':True,
            'reason':'explicit source45s guard raised RuntimeError at the beginning of a target row; prior complete candidates retained',
            'no_claim_all_units_unsafe':True,'current_generation45_unchanged':True}
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);p=out/'cursor.json';p.write_text(json.dumps(cursor,indent=2)+'\n')
    paths=[short,*fs,ROOT/'code/s61_sparse_generic5_constrained83_compact45_unit_gate_20261105_v64.py',Path(__file__),ROOT/f'research/notes/{ID}.md']
    report={'id':ID,'all_passed':True,'source_gate_completed':False,'bounded_source_stop_verified':True,'source_only':True,
            'completed_forecasts':10,'remaining_shortlist':6,'generation45_unchanged':True,'no_arithmetic_or_retry':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{'cursor.json':sha(p)}}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
