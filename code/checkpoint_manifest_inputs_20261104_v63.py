#!/usr/bin/env python3
"""Resolve report input aliases only through hash-bound evidence manifests."""
from functools import lru_cache
import hashlib
import json
from pathlib import Path

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

@lru_cache(maxsize=128)
def manifests(report_name, report_digest):
    report=Path(report_name)
    assert digest(report)==report_digest
    data=json.loads(report.read_text())
    result={}
    for source,want in data.get('evidence_sha256',{}).items():
        path=Path(source) if source.startswith(('reports/','code/','research/','papers/')) else report.parent/source
        if not path.is_file() or path.suffix!='.json' or digest(path)!=want:
            continue
        evidence=json.loads(path.read_text())
        if not isinstance(evidence,dict):continue
        for field in ('source_inputs','source_manifest'):
            items=evidence.get(field,{})
            if not isinstance(items,dict):continue
            for alias,item in items.items():
                if not isinstance(item,dict):continue
                if isinstance(item.get('path'),str) and isinstance(item.get('sha256'),str):
                    result.setdefault(alias,[]).append(dict(path=item['path'],sha256=item['sha256'],
                        manifest=str(path),manifest_sha256=want,manifest_field=field))
    return result

def resolve(report_path, alias, expected_digest):
    expected_path=None
    if isinstance(expected_digest,dict):
        expected_path=expected_digest.get('path')
        expected_digest=expected_digest.get('sha256')
        if not isinstance(expected_path,str) or not isinstance(expected_digest,str):
            return None
    candidates=manifests(str(report_path),digest(report_path)).get(alias,[])
    matches=[item for item in candidates if item['sha256']==expected_digest
             and (expected_path is None or item['path']==expected_path)]
    if not matches:return None
    assert len({item['path'] for item in matches})==1,(report_path,alias,matches)
    item=dict(matches[0])
    assert Path(item['path']).is_file(),item
    assert digest(Path(item['path']))==expected_digest,item
    item.update(report=str(report_path),alias=alias)
    return item
