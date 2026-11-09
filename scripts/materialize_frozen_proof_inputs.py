#!/usr/bin/env python3
"""Prepare a new isolated input tree for one frozen certificate's producer.

This copies inputs only. It never runs code or overwrites the research tree.
Run the selected producer there with its original CLI and a new output path.
"""
from pathlib import Path
from collections import deque
import argparse,hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[1]
def main(report,destination):
    destination=destination.resolve()
    if destination.exists():raise FileExistsError(destination)
    if not report.is_relative_to(ROOT):raise ValueError('Report must be in this repository')
    manifest=json.loads((ROOT/'reports/required_proof_dependency_closure.json').read_text());assert manifest['all_required_inputs_available']
    name=str(report.relative_to(ROOT));todo=deque([name]);files=set()
    while todo:
        f=todo.popleft()
        if f in files:continue
        files.add(f);todo.extend(manifest['edges'].get(f,[]))
    aliases=manifest['root_input_aliases']
    frozen=json.loads((ROOT/'research/reproduction/frozen_producer_index.json').read_text())['records']
    index={(r['original_path'],r['sha256']):r['materialized_path']for r in frozen}
    destination.mkdir(parents=True)
    def copy(src,relative):
        target=destination/relative
        if not target.resolve().is_relative_to(destination):raise ValueError(relative)
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/src,target)
    links={}
    for f in sorted(files):
        source=ROOT/f;actual=source.resolve();assert actual.is_relative_to(ROOT)
        canonical=str(actual.relative_to(ROOT));copy(canonical,canonical)
        for i in range(1,len(Path(f).parts)+1):
            prefix=ROOT.joinpath(*Path(f).parts[:i])
            if prefix.is_symlink():
                key=str((prefix.parent.resolve()/prefix.name).relative_to(ROOT))
                links[key]=__import__('os').readlink(prefix)
    for name,link in links.items():
        p=destination/name;p.parent.mkdir(parents=True,exist_ok=True)
        assert (p.parent/link).resolve().is_relative_to(destination)
        p.symlink_to(link)
    for logical,actual in aliases.items():
        if actual in files:copy(actual,logical)
    bindings=[]
    def binding(path,expected):
        if not isinstance(path,str)or not isinstance(expected,str)or len(expected)!=64:return
        src=index.get((path,expected),aliases.get(path,path));p=ROOT/src
        if not p.is_file():return
        actual=hashlib.sha256(p.read_bytes()).hexdigest()
        if actual!=expected:
            if path.endswith(('.md','.tex'))or'/jobs/'in path:return
            raise ValueError(('Cannot bind exact requested input',path,expected,actual))
        copy(src,path);bindings.append({'path':path,'sha256':expected,'source':src})
    packet=json.loads(report.read_text())
    for path,expected in packet.get('input_sha256',{}).items():binding(path,expected)
    for row in packet.get('input_bindings',[]):binding(row.get('path'),row.get('sha256'))
    (destination/'REPLAY_INPUTS.json').write_text(json.dumps({'report':name,'bindings':bindings,'files_copied':len(files),'scope':__doc__},indent=2)+'\n')
    print(json.dumps({'destination':str(destination),'files_copied':len(files),'exact_top_level_bindings':len(bindings)},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report',type=Path,required=True);p.add_argument('--destination',type=Path,required=True);a=p.parse_args();main(a.report.resolve(),a.destination)
