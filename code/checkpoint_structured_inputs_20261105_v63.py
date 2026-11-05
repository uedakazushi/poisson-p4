#!/usr/bin/env python3
"""Resolve literal structured report inputs, retaining the older alias manifest check."""
import hashlib,json,runpy
from pathlib import Path
LEGACY=Path('code/checkpoint_manifest_inputs_20261104_v63.py')
LEGACY_SHA256='30e1cbfa94c4a23f7afd9df50d9b22bd7a65bb4222f85614a923031796c88fa4'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(LEGACY)==LEGACY_SHA256
legacy_resolve=runpy.run_path(str(LEGACY))['resolve']
def resolve(report_path,alias,expected):
    if not isinstance(expected,dict):return legacy_resolve(report_path,alias,expected)
    source=expected.get('path');want=expected.get('sha256')
    if not isinstance(source,str)or not isinstance(want,str)or len(want)!=64:return None
    if any(c not in '0123456789abcdef'for c in want):return None
    report=Path(report_path);data=json.loads(report.read_text())
    assert data.get('input_sha256',{}).get(alias)==expected,(report,alias)
    p=Path(source)
    assert p.is_file() and sha(p)==want,(report,alias,source)
    return dict(path=source,sha256=want,report=str(report),alias=alias,
        manifest=str(report),manifest_sha256=sha(report),manifest_field='input_sha256',
        source_kind='literal_structured_input',legacy_checker=str(LEGACY),
        legacy_checker_sha256=LEGACY_SHA256)
