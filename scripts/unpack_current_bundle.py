#!/usr/bin/env python3
"""Reconstruct the current research tree from its disjoint Library ZIP volumes."""
import argparse,hashlib,json
from pathlib import Path
import zipfile

def main(parts,destination):
    target=destination.resolve()
    if target.exists():raise FileExistsError(target)
    archives=[zipfile.ZipFile(p)for p in parts]
    try:
        manifests=[json.loads(z.read('BUNDLE_MANIFEST.json'))for z in archives]
        ids={m['bundle_id']for m in manifests};assert len(ids)==1
        manifest=manifests[0];assert len(archives)==manifest['part_count']
        chunks={}
        for z in archives:
            assert z.testzip()is None
            for name in z.namelist():
                if name.startswith('payload/'):
                    digest=name.split('/')[1]
                    if digest in chunks:raise ValueError('Duplicate chunk')
                    chunks[digest]=(z,name)
        assert set(chunks)==set(manifest['chunks'])
        target.mkdir(parents=True)
        for f in manifest['files']:
            path=target/f['path']
            if not path.resolve().is_relative_to(target):raise ValueError(f['path'])
            path.parent.mkdir(parents=True,exist_ok=True);h=hashlib.sha256();size=0
            with path.open('xb')as out:
                for digest in f['chunks']:
                    z,name=chunks[digest];data=z.read(name)
                    assert hashlib.sha256(data).hexdigest()==digest
                    out.write(data);h.update(data);size+=len(data)
            assert size==f['bytes']and h.hexdigest()==f['sha256']
        for name,link in manifest.get('symlinks',{}).items():
            path=target/name
            if not path.resolve().is_relative_to(target):raise ValueError(name)
            if not (path.parent/link).resolve().is_relative_to(target):raise ValueError(link)
            path.parent.mkdir(parents=True,exist_ok=True);path.symlink_to(link)
        for logical,physical in manifest.get('logical_paths',{}).items():
            assert (target/logical).resolve()==(target/physical).resolve()
        (target/'BUNDLE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
        print(json.dumps({'destination':str(target),'files':len(manifest['files']),'all_sha256':'PASS','source_commit':manifest['source_commit']},indent=2))
    finally:
        for z in archives:z.close()
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--parts',type=Path,nargs='+',required=True);p.add_argument('--destination',type=Path,required=True);a=p.parse_args();main(a.parts,a.destination)
