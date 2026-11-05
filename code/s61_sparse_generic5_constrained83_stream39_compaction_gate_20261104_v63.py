#!/usr/bin/env python3
"""Source-only current-generation compaction gate; no coefficient updates."""
from pathlib import Path
import hashlib,json,resource,signal,time
from s61_sparse_generic5_constrained83_stream38_source_gate_20261104_v63 import index_row
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM39_COMPACTION_GATE_20261104_v63'
PREV='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM38_PIVOT_20261104_v63'
PRIOR_INDEX='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM38_SOURCE_GATE_20261104_v63'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return (json.dumps(v,indent=2)+'\n').encode()
def numeric(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(numeric(z)for z in v)
    if isinstance(v,dict):return sum(numeric(z)for z in v.values())
    return 0
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False)
    pd=ROOT/'reports'/PREV;pr=pd/'report.json';mp=pd/'generation-1.json'
    assert sha(pr)=='8b9b6f54be58ca29f9a29ec8876b50b60dd02f0ac299a759aae199ca8ebab3bc'
    assert sha(mp)=='2e64cc4f0f79465a6aa9b04ed3b4874aa422c06a5733625869a6163a128dc232'
    report=read(pr);assert report['all_passed']and report['completed_total_on_reordered_path']==39
    m=read(mp);assert m['shape']==[197,44]and m['logical_total_terms']==234898
    oldip=ROOT/'reports'/PRIOR_INDEX/'mixed-literal-file-index.json';old=read(oldip);assert sha(oldip)==m['index_sha256']
    inputpaths={pr,mp,oldip,Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_sparse_generic5_constrained83_stream38_source_gate_20261104_v63.py'}
    filemap={};files=[];cells=[];rowcells={};colcells={};receipts=[]
    def fid(path,h):
        k=str(path.relative_to(ROOT))
        if k not in filemap:
            assert sha(path)==h;filemap[k]=len(files);files.append({'path':k,'sha256':h,'bytes':path.stat().st_size});inputpaths.add(path)
        else:assert files[filemap[k]]['sha256']==h
        return filemap[k]
    for row in m['rows']:
        i=row['original_row'];assert i!=15
        if'literal_row_file'in row:
            p=pd/row['literal_row_file'];h=row['sha256'];k=fid(p,h);new=index_row(p,i,k)
            assert sum(c['terms']for c in new)==row['terms']
        else:
            new=[]
            for a in row['indexed_source_cell_indices']:
                c=old['all_literal_cells'][a];f=old['files'][c['file_index']];k=fid(ROOT/f['path'],f['sha256'])
                new.append({**c,'file_index':k})
        for c in new:
            assert c['row']==i and c['column']in m['column_original_ids']and c['column']!=42
            a=len(cells);cells.append(c);rowcells.setdefault(i,[]).append(a);colcells.setdefault(c['column'],[]).append(a)
    # Every active generated row also has its immutable back-add receipt.
    for f in files:
        p=ROOT/f['path']
        if p.parent.name=='rows':
            receipt=p.parent.parent/'receipts'/p.name;q=read(receipt);inputpaths.add(receipt)
            assert q['literal_row_sha256']==f['sha256']and q['all_cell_back_add_identities_checked']and q['pivot_column_cleared']
            receipts.append({'path':str(receipt.relative_to(ROOT)),'sha256':sha(receipt)})
    assert sum(c['terms']for c in cells)==234898 and len({(c['row'],c['column'])for c in cells})==len(cells)
    maxcell=max(c['terms']for c in cells)
    index={'generation':39,'source_manifest':str(mp.relative_to(ROOT)),'source_manifest_sha256':sha(mp),
           'files':files,'all_literal_cells':cells,'row_cell_indices':rowcells,'column_cell_indices':colcells,
           'all_original_row_ids':[r['original_row']for r in m['rows']],'all_original_column_ids':m['column_original_ids'],
           'logical_shape':[197,44],'logical_coefficient_terms':234898,'max_cell_terms':maxcell,'active_row_receipts':receipts}
    ip=out/'current-generation-cell-index.json';ip.write_bytes(enc(index))
    active=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(inputpaths)]
    active.append({'path':str(ip.relative_to(ROOT)),'sha256':sha(ip),'bytes':ip.stat().st_size})
    # Future snapshot writes source payload bytes literally, with small headers.
    snapshot=sum(c['payload_bytes']+len(str(c['row']))+len(str(c['column']))+8 for c in cells)+32768
    # Each output-cell record removes file_index; even a 32-digit future byte
    # offset and extra top-level fields fit +64 bytes per current cell.
    newindex=ip.stat().st_size+64*len(cells)+32768
    metaupper=numeric(old)+numeric(index)+numeric(m)+numeric(active)+6*maxcell+2048
    # Compaction owns current+new metadata indices but no old prior index.
    compactowned=2*numeric(index)+numeric(active)+12*maxcell+4096
    activebytes=sum(f['bytes']for f in active)
    reserve=1024*1024
    whole=activebytes+snapshot+newindex+reserve
    bound={'fresh_job_actual_read_inputs':active,'actual_read_input_RAW_bytes':activebytes,
           'old_full_retained_history_unchanged_RAW_bytes':23923327,
           'retired_history_not_deleted_and_not_reclassified_as_current_input':True,
           'current_literal_files':len(files),'current_cells':len(cells),'current_logical_terms':234898,
           'max_owned_cell_terms':maxcell,'source_gate_owned_numeric_slots_upper':metaupper,
           'compaction_owned_numeric_and_cell_scratch_upper':compactowned,
           'new_full_snapshot_RAW_upper':snapshot,'new_full_index_RAW_upper':newindex,
           'all_source_gate_extra_metadata_future_code_IO_receipts_manifests_reserve':reserve,
           'fresh_compaction_job_RAW_TOTAL_upper':whole,
           'FULL_safe':max(metaupper,compactowned)<=150000 and whole<=32*1024*1024,
           'no_polynomial_math_or_pivot':True}
    (out/'source-and-whole-budget.json').write_bytes(enc(bound))
    rr={'id':ID,'all_passed':True,'source_only':True,'current_shape':[197,44],'logical_terms':234898,
        'current_cells':len(cells),'actual_read_input_files':len(active),'FULL_safe_compaction':bound['FULL_safe'],
        'fresh_RAW_upper':whole,'source_gate_owned_numeric_upper':metaupper,'compaction_owned_numeric_upper':compactowned,
        'all_current_source_cells_hash_bound':True,'all_active_row_receipts_checked':True,
        'retired_history_preserved_not_an_artificial_deletion':True,'no_actual_compaction_or_pivot':True,
        'input_sha256':{f['path']:f['sha256']for f in active if f['path']!=str(ip.relative_to(ROOT))},
        'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
        'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_bytes(enc(rr));print(json.dumps({k:v for k,v in rr.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
