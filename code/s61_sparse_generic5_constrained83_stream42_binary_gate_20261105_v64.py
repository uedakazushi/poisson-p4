#!/usr/bin/env python3
"""Changed representation SOURCE gate, no full snapshot or coefficient math."""
from pathlib import Path
import hashlib,json,resource,signal,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import MAGIC,encode,decode,canonical
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM42_BINARY_GATE_20261105_v64'
FUTURE='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM42_BINARY_COMPACTION_20261105_v64'
TEXT='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM42_COMPACTION_GATE_20261105_v64'
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
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);td=ROOT/'reports'/TEXT
    tr=td/'report.json';tb=td/'source-and-whole-budget.json';ip=td/'current-generation-cell-index.json'
    assert sha(tr)=='01d11d53dca151358ba68bf36986ff8368763938f43a5c4c653815ebc958778a'
    r=read(tr);assert r['all_passed']and not r['FULL_safe_compaction']
    for p in(tb,ip):assert sha(p)==r['evidence_sha256'][p.name]
    old=read(tb);index=read(ip);assert index['logical_shape']==[194,41]and index['logical_coefficient_terms']==316841
    active={f['path']:f for f in old['fresh_job_actual_read_inputs']}
    codec=ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py'
    for p in(tr,tb,Path(__file__),codec,ROOT/f'research/notes/{ID}.md'):
        active[str(p.relative_to(ROOT))]={'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}
    for f in active.values():assert sha(ROOT/f['path'])==f['sha256']
    cells=[];offset=len(MAGIC);termcount=0;maxcell=0;maxbits=0
    for k,c in enumerate(index['all_literal_cells']):
        assert time.monotonic()-start<52
        p=ROOT/index['files'][c['file_index']]['path']
        with p.open('rb')as f:f.seek(c['payload_offset']);payload=f.read(c['payload_bytes'])
        assert hashlib.sha256(payload).hexdigest()==c['payload_sha256']
        value=json.loads(payload);binary=encode(value);back=decode(binary)
        assert back==canonical(value)and len(back['terms'])==c['terms']
        maxbits=max(maxbits,max(abs(e).bit_length()for e in value['shift']),max(e.bit_length()for t in value['terms']for e in t[1:]))
        z={x:v for x,v in c.items()if x not in('file_index','payload_offset','payload_bytes','payload_sha256')}
        z.update(source_input_cell_index=k,payload_offset=offset,payload_bytes=len(binary),payload_sha256=hashlib.sha256(binary).hexdigest())
        cells.append(z);offset+=len(binary);termcount+=c['terms'];maxcell=max(maxcell,c['terms'])
        del payload,value,binary,back
    assert termcount==316841 and maxbits<64
    plan={'source_index':str(ip.relative_to(ROOT)),'source_index_sha256':sha(ip),'encoding':'primitive-int-exponent-v1',
          'shape':[194,41],'logical_terms':316841,'all_original_row_ids':index['all_original_row_ids'],
          'all_original_column_ids':index['all_original_column_ids'],'all_literal_cells':cells,
          'row_cell_indices':index['row_cell_indices'],'column_cell_indices':index['column_cell_indices'],
          'exact_binary_file_bytes':offset,'all_source_cells_roundtrip_exact':True,'max_shift_or_exponent_bits':maxbits}
    pp=out/'binary-cell-plan.json';pp.write_bytes(enc(plan));active[str(pp.relative_to(ROOT))]={'path':str(pp.relative_to(ROOT)),'sha256':sha(pp),'bytes':pp.stat().st_size}
    template={'encoding':'primitive-int-exponent-v1','source_file':f'reports/{FUTURE}/current-matrix.bin','source_file_sha256':'0'*64,
              'logical_shape':[194,41],'logical_coefficient_terms':316841,'all_literal_cells':cells,
              'row_cell_indices':index['row_cell_indices'],'column_cell_indices':index['column_cell_indices'],
              'all_original_row_ids':index['all_original_row_ids'],'all_original_column_ids':index['all_original_column_ids'],
              'source_plan':str(pp.relative_to(ROOT)),'source_plan_sha256':sha(pp)}
    newindex=len(enc(template))+4096
    owned=numeric(index)+numeric(plan)+numeric(list(active.values()))+24*maxcell+8192
    inputbytes=sum(f['bytes']for f in active.values());reserve=1024*1024;whole=inputbytes+offset+newindex+reserve
    proof={'all_passed':True,'source_only':True,'changed_method':'Uncompressed canonical signed magnitude integers, zigzag shifts, unsigned varint exponents; every source cell encode/decode exact. No full output snapshot yet.',
           'actual_READ_inputs':list(active.values()),'actual_READ_input_RAW_bytes':inputbytes,
           'binary_snapshot_exact_bytes':offset,'binary_output_index_RAW_upper':newindex,'all_extra_new_metadata_code_receipt_IO_reserve':reserve,
           'fresh_job_RAW_upper':whole,'owned_numeric_and_codec_scratch_upper':owned,'source_cells':len(cells),'logical_terms':316841,
           'FULL_safe':owned<=150000 and whole<=32*1024*1024,'text35MB_stop_preserved':True,'retired_history_retained_not_deleted':True}
    (out/'whole-budget.json').write_bytes(enc(proof))
    report={'id':ID,'all_passed':True,'source_only':True,'all_source_cell_binary_roundtrips':len(cells),'logical_terms':316841,
            'FULL_safe_binary_compaction':proof['FULL_safe'],'fresh_RAW_upper':whole,'owned_upper':owned,'binary_exact_bytes':offset,
            'input_sha256':{k:f['sha256']for k,f in active.items()if k!=str(pp.relative_to(ROOT))},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,
            'no_arithmetic_pivot_gzip_or_cap_increase':True}
    (out/'report.json').write_bytes(enc(report));signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
