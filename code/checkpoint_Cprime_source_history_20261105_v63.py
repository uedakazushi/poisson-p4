#!/usr/bin/env python3
"""Verify the single immutable prepared-worker input recorded before a reader fix."""
import hashlib
import json
from pathlib import Path

HISTORY = Path('reports/AMBIENT_B18_CPRIME_LINEAR_INITIAL_SOURCE_HISTORY_20261105_v63/report.json')
HISTORY_SHA = '4538b10299661a210b86043ce1d702a21058ebe11f8b020d7a0a06a9909f44ef'
SOURCE_REPORT = 'reports/AMBIENT_B18_CPRIME_LINEAR_INITIAL_SOURCE_20261105_v63_REVIEW/report.json'
SOURCE_SHA = 'b205c1d56b3f2f3e1800d456ec4a7a0a2b1f3b9cd5c922ff5d17e86914e862f1'
WORKER = 'code/ambient_B18_Cprime_linear_initial_values_20261105_v63.py'
OLD_SHA = '1be2cda944e1dfc0a608e0ac03035168438f4240abc279ad34739146c474b028'
FINAL_SHA = '5f4df4c1ebbcd4c97ebfd7234961f71ed8300fa2dbd302046a5d88d4d99f0ccc'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def verify():
    assert sha(HISTORY) == HISTORY_SHA
    h = json.loads(HISTORY.read_text())
    assert h['all_passed'] and h['historical_input_path'] == WORKER
    assert h['prepared_source_sha256'] == OLD_SHA
    assert h['executed_final_source_sha256'] == FINAL_SHA
    assert sha(SOURCE_REPORT) == SOURCE_SHA
    s = json.loads(Path(SOURCE_REPORT).read_text())
    assert s['all_passed'] and s['input_sha256'][WORKER] == OLD_SHA
    old = Path(h['contentaddressed_prepared_source'])
    assert old == HISTORY.parent / 'source-history' / (OLD_SHA + '.py')
    assert sha(old) == OLD_SHA and sha(WORKER) == FINAL_SHA
    for p, digest in h['input_sha256'].items():
        assert sha(p) == digest, p
    for p, digest in h['evidence_sha256'].items():
        assert sha(HISTORY.parent / p) == digest, p
    return dict(
        original_report=SOURCE_REPORT, original_report_sha256=SOURCE_SHA,
        historical_input_path=WORKER, preserved_historical_path=str(old),
        historical_sha256=OLD_SHA, final_sha256=FINAL_SHA,
        provenance_receipt=str(HISTORY), provenance_receipt_sha256=HISTORY_SHA,
        resolution='Only the pinned SOURCE-review prepared-worker hash is historical; final code and value peer remain checked in full.',
    )

if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
