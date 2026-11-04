#!/usr/bin/env python3
"""Verify both immutable versions of one documented metadata phase artifact."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def verify(root=ROOT):
    history_rel=Path('reports/S61_NONSQUARE_PREPOWER_G115_PHASE_PROVENANCE_20261104_v63')
    primary_rel=Path('reports/S61_NONSQUARE_PREPOWER_G115_T19_LITERAL_PAIR_PREFLIGHT_20261104_v63')
    history_dir=root/history_rel; primary_dir=root/primary_rel
    history_path=history_dir/'report.json'; primary_path=primary_dir/'report.json'
    history=json.loads(history_path.read_text()); primary=json.loads(primary_path.read_text())
    old=history_dir/'phase-before-final.json'
    final=history_dir/'phase-final-unchanged.json'
    current=primary_dir/'phase.json'
    assert history['all_passed'] and history['metadata_only']
    assert digest(primary_path)==history['original_preflight_report_sha256']=='eb7495b6e5692bb84c9ca2b9df40922fc0a3132146a065e2651e0118776f98ca'
    assert primary['evidence_sha256']['phase.json']==digest(old)=='e88475844a64269e94ce3107061339deda98a9e3514527941cd8bea179bec3bd'
    assert digest(current)==digest(final)==history['current_phase_sha256']=='acfa1ab05f5cebab720c44abc02c0e303235f6e52bea7a11f11fa8a08b4bbee2'
    for name,want in history['preserved_files_sha256'].items():
        assert digest(history_dir/name)==want,name
    for name,want in primary['evidence_sha256'].items():
        if name!='phase.json': assert digest(primary_dir/name)==want,name
    for name,want in primary.get('input_sha256',{}).items():
        assert digest(root/name)==want,name
    return dict(original_report=str(primary_rel/'report.json'),
        original_report_sha256=digest(primary_path),metadata_path=str(primary_rel/'phase.json'),
        reported_historical_sha256=digest(old),current_sha256=digest(current),
        preserved_historical_path=str(history_rel/old.name),preserved_final_path=str(history_rel/final.name),
        provenance_receipt=str(history_rel/'report.json'),provenance_receipt_sha256=digest(history_path),
        all_other_primary_evidence_and_inputs_match_literally=True,
        primary_report_and_current_phase_unchanged=True,issue=history['issue'])

if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
