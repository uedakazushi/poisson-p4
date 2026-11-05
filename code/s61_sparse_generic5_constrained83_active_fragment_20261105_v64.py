#!/usr/bin/env python3
"""Exact active subset of a retained binary snapshot; no polynomial change."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import MAGIC,decode,encode
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(numeric(x)for x in v)
    if isinstance(v,dict):return sum(numeric(x)for x in v.values())
    return 0
def tree(p):return sum(f.stat().st_size for f in p.rglob('*')if f.is_file())
def main(a):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    gate=f'S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM{a.generation}_ACTIVE_SOURCE_FRAGMENT_GATE_20261105_v64'
    result=f'S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM{a.generation}_ACTIVE_SOURCE_FRAGMENT_20261105_v64'
    if a.phase=='gate':
        pd=ROOT/'reports'/a.source;pr=pd/'report.json';mp=pd/'generation-1.json';assert sha(pr)==a.report_sha and sha(mp)==a.manifest_sha
        r=read(pr);m=read(mp);assert r['all_passed']and r['completed_total_on_reordered_path']==a.generation and m['complete_pivot']
        ip=ROOT/r['literal_input_index_path'];index=read(ip);assert sha(ip)==m['index_sha256']and index['encoding']=='primitive-int-exponent-v1'and'files'not in index
        original=ROOT/index['source_file'];assert sha(original)==index['source_file_sha256']
        ks=[k for row in m['rows']for k in row.get('indexed_source_cell_indices',[])];assert len(ks)==len(set(ks));ks.sort()
        untouched={row['original_row']for row in m['rows']if'indexed_source_cell_indices'in row};offset=len(MAGIC);cells={};terms=0;maxcell=0
        with original.open('rb')as f:
            assert f.read(len(MAGIC))==MAGIC
            for k in ks:
                c=index['all_literal_cells'][k];assert c['row']in untouched and c['column']in m['column_original_ids']
                f.seek(c['payload_offset']);payload=f.read(c['payload_bytes']);assert hashlib.sha256(payload).hexdigest()==c['payload_sha256']
                p=decode(payload);assert encode(p)==payload and len(p['terms'])==c['terms']
                cells[str(k)]={'source_offset':c['payload_offset'],'payload_offset':offset,'payload_bytes':c['payload_bytes'],'payload_sha256':c['payload_sha256'],'terms':c['terms'],'row':c['row'],'column':c['column']}
                offset+=c['payload_bytes'];terms+=c['terms'];maxcell=max(maxcell,c['terms']);del payload,p
        out=ROOT/'reports'/gate;out.mkdir(parents=True,exist_ok=False)
        plan={'source_generation':a.generation,'source_manifest':str(mp.relative_to(ROOT)),'source_manifest_sha256':a.manifest_sha,
              'source_index':str(ip.relative_to(ROOT)),'source_index_sha256':sha(ip),'original_snapshot':str(original.relative_to(ROOT)),'original_snapshot_sha256':sha(original),
              'cells':cells,'exact_fragment_bytes':offset,'selected_cells':len(cells),'selected_terms':terms,
              'scope':'Only original indexed source cells still referenced by the current manifest; all newly computed row files are unchanged and excluded from this extraction.'}
        pp=out/'fragment-plan.json';pp.write_bytes(enc(plan))
        paths=[pr,mp,ip,original,pp,Path(__file__),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/f'research/notes/{gate}.md']
        active=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in paths]
        own=2*numeric(index)+2*numeric(plan)+24*maxcell+8192;total=sum(x['bytes']for x in active)+offset+len(enc(plan))+16384+1024*1024
        budget={'all_passed':True,'source_only':True,'actual_READ_inputs':active,'actual_READ_input_RAW_bytes':sum(x['bytes']for x in active),
                'exact_fragment_bytes':offset,'replacement_index_RAW_upper':len(enc(plan))+16384,'extra_code_receipt_metadata_allowance':1024*1024,
                'fresh_RAW_upper':total,'owned_upper':own,'FULL_safe':own<=150000 and total<=32*1024*1024,'all_original_files_retained':True}
        (out/'whole-budget.json').write_bytes(enc(budget));report={'id':gate,'all_passed':True,'source_only':True,'source_generation':a.generation,'FULL_safe_fragment':budget['FULL_safe'],
            'selected_cells':len(cells),'selected_terms':terms,'fragment_bytes':offset,'fresh_RAW_upper':total,'owned_upper':own,
            'input_sha256':{x['path']:x['sha256']for x in active if x['path']!=str(pp.relative_to(ROOT))},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start}
        (out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
    else:
        gd=ROOT/'reports'/gate;gr=gd/'report.json';assert sha(gr)==a.gate_sha;r=read(gr);assert r['all_passed']and r['FULL_safe_fragment']
        bp=gd/'whole-budget.json';pp=gd/'fragment-plan.json'
        for p in(bp,pp):assert sha(p)==r['evidence_sha256'][p.name]
        b=read(bp);plan=read(pp);active=b['actual_READ_inputs']
        for x in active:assert sha(ROOT/x['path'])==x['sha256']
        out=ROOT/'reports'/result;out.mkdir(parents=True,exist_ok=False);source=ROOT/plan['original_snapshot'];fragment=out/'active-source.bin';tmp=out/'active-source.bin.partial'
        extra=[gr,bp,ROOT/f'research/notes/{result}.md'];ebytes=sum(p.stat().st_size for p in extra);base=b['actual_READ_input_RAW_bytes']
        with source.open('rb')as src,tmp.open('w+b')as f:
            f.write(MAGIC)
            for k,c in plan['cells'].items():
                assert time.monotonic()-start<55;src.seek(c['source_offset']);payload=src.read(c['payload_bytes']);assert hashlib.sha256(payload).hexdigest()==c['payload_sha256']
                assert f.tell()==c['payload_offset']and f.tell()+len(payload)<=b['exact_fragment_bytes'];f.write(payload);end=f.tell();f.flush();f.seek(c['payload_offset']);back=f.read(c['payload_bytes']);f.seek(end)
                assert back==payload and encode(decode(back))==payload and len(decode(back)['terms'])==c['terms']
                assert base+ebytes+tree(out)<b['fresh_RAW_upper']
            assert f.tell()==b['exact_fragment_bytes']
        tmp.rename(fragment);fh=sha(fragment)
        replacement={'source_generation':a.generation,'source_manifest_sha256':plan['source_manifest_sha256'],'source_index_sha256':plan['source_index_sha256'],
                     'original_snapshot':plan['original_snapshot'],'original_snapshot_sha256':plan['original_snapshot_sha256'],
                     'fragment_file':str(fragment.relative_to(ROOT)),'fragment_sha256':fh,'cells':plan['cells'],
                     'all_active_old_source_references_covered_exactly':True,'all_original_files_still_retained':True}
        ip=out/'replacement-index.json';ip.write_bytes(enc(replacement));assert ip.stat().st_size<=b['replacement_index_RAW_upper']
        assert sha(source)==plan['original_snapshot_sha256']
        report={'id':result,'all_passed':True,'source_generation':a.generation,'selected_cells':plan['selected_cells'],'selected_terms':plan['selected_terms'],
                'fragment_bytes':fragment.stat().st_size,'all_selected_values_and_bytes_roundtrip':True,
                'source_manifest_sha256':plan['source_manifest_sha256'],'source_index_sha256':plan['source_index_sha256'],
                'replacement_index':ip.name,'replacement_index_sha256':sha(ip),'old_snapshot_retained_unchanged':True,
                'not_a_whole_matrix_or_pivot_claim':True,'owned_upper':b['owned_upper'],'fresh_RAW_upper':b['fresh_RAW_upper'],
                'input_sha256':{**{x['path']:x['sha256']for x in active},**{str(p.relative_to(ROOT)):sha(p)for p in extra}},
                'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start}
        rp=out/'report.json';rp.write_bytes(enc(report));assert base+ebytes+tree(out)<=b['fresh_RAW_upper']<32*1024*1024
        print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp))
    signal.alarm(0)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--phase',choices=['gate','copy'],required=True);p.add_argument('--generation',type=int,required=True);p.add_argument('--source');p.add_argument('--report-sha');p.add_argument('--manifest-sha');p.add_argument('--gate-sha');main(p.parse_args())
