#!/usr/bin/env python3
"""Restore published input caches needed by byte-identical historical workers."""
import argparse,gzip,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main(destination):
    destination.mkdir(parents=True,exist_ok=True)
    parent=ROOT/'reports/AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63'
    entries=[('point-setup.json.gz','CB-cubic-source-inversion-setup.json'),
             ('point-cubic.json.gz','CB-actualD-cubic-syzygy-probe.json')]
    restored=[]
    for source,target in entries:
        with gzip.open(parent/source,'rb')as f:data=f.read()
        assert json.loads(data)['all_passed']
        path=destination/target
        if path.exists():
            existing=json.loads(path.read_text());published=json.loads(data)
            assert {k:v for k,v in existing.items()if k!='seconds'}=={k:v for k,v in published.items()if k!='seconds'},str(path)
            restored.append(dict(path=str(path),action='existing equivalent cache retained'))
        else:
            path.write_bytes(data);restored.append(dict(path=str(path),action='restored from published gzip',sha256=hashlib.sha256(data).hexdigest()))
    print(json.dumps(dict(all_passed=True,restored=restored),indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--destination',type=Path,default=ROOT/'work/v62-review');a=p.parse_args();main(a.destination.resolve())
