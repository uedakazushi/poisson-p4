#!/usr/bin/env python3
"""ONE exact current-generation byte/value compaction; no polynomial math."""
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM40_COMPACTION_20261105_v64'
GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM40_COMPACTION_GATE_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def numeric(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(numeric(z)for z in v)
    if isinstance(v,dict):return sum(numeric(z)for z in v.values())
    return 0
def treebytes(p):return sum(q.stat().st_size for q in p.rglob('*')if q.is_file())
def jbytes(v):return json.dumps(v,separators=(',',':')).encode()
def streamjson(path,value,cap):
    tmp=path.with_name(path.name+'.partial');assert not path.exists();n=0
    with tmp.open('wb')as f:
        for s in json.JSONEncoder(separators=(',',':')).iterencode(value):
            b=s.encode();assert n+len(b)<=cap;f.write(b);n+=len(b)
        assert n+1<=cap;f.write(b'\n')
    tmp.rename(path)
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);gd=ROOT/'reports'/GATE
    gr=gd/'report.json';bp=gd/'source-and-whole-budget.json';ip=gd/'current-generation-cell-index.json'
    assert sha(gr)=='ccca04eef129b9a99a6ef14088f8d6016a770fdcd52d545d08660b4f06617d20'
    g=read(gr);assert g['all_passed']and g['FULL_safe_compaction']
    for p in(bp,ip):assert sha(p)==g['evidence_sha256'][p.name]
    b=read(bp);assert b['FULL_safe']and b['fresh_compaction_job_RAW_TOTAL_upper']==32394791
    active=b['fresh_job_actual_read_inputs'];inputbytes=sum(v['bytes']for v in active);assert inputbytes==19941661
    for f in active:
        p=ROOT/f['path'];assert p.stat().st_size==f['bytes']and sha(p)==f['sha256']
    index=read(ip);assert index['logical_coefficient_terms']==292652 and index['logical_shape']==[196,43]
    extra={gr,bp,Path(__file__),ROOT/f'research/notes/{ID}.md'}
    activepaths={f['path']for f in active};extra={p for p in extra if str(p.relative_to(ROOT))not in activepaths}
    extrabytes=sum(p.stat().st_size for p in extra)
    def whole():
        n=inputbytes+extrabytes+treebytes(out)
        assert n<=32394791 and n<=32*1024*1024
        return n
    pre={'id':ID,'source_gate_sha256':sha(gr),'fresh_actual_READ_input_bytes':inputbytes,
         'all_actual_READ_inputs':active,'all_additional_READ_inputs':{str(p.relative_to(ROOT)):sha(p)for p in extra},
         'additional_READ_bytes':extrabytes,'previous_fresh_job_RAW_receipt_bytes_unchanged':20889002,
         'hard_total_seconds':60,'address_space_limit_bytes':2000000000,'owned_slots_upper':94504,
         'fresh_RAW_upper':32394791,'no_old_pivot_replay_or_polynomial_arithmetic':True}
    streamjson(out/'preflight.json',pre,200000);whole()
    # Keep current metadata and one output metadata index; cell value buffers
    # are bounded by the largest 548-term source polynomial.
    newcells=[];rowcells={};colcells={};rowchecks=[];nterms=0;payloadreads=0;readbacks=0
    matrix=out/'current-matrix.json';tmp=out/'current-matrix.json.partial';limit=b['new_full_snapshot_RAW_upper']
    prefix={'generation':40,'shape':[196,43],'row_original_ids':index['all_original_row_ids'],'column_original_ids':index['all_original_column_ids']}
    initial=jbytes(prefix)[:-1]+b',"sparse_entries":['
    rowstats={i:[0,0,hashlib.sha256(),hashlib.sha256()]for i in index['all_original_row_ids']}
    with tmp.open('w+b')as f:
        def write(buf):
            assert f.tell()+len(buf)<=limit;f.write(buf)
        write(initial)
        for k,c in enumerate(index['all_literal_cells']):
            assert time.monotonic()-start<56
            p=ROOT/index['files'][c['file_index']]['path']
            with p.open('rb')as source:source.seek(c['payload_offset']);payload=source.read(c['payload_bytes'])
            assert len(payload)==c['payload_bytes']and hashlib.sha256(payload).hexdigest()==c['payload_sha256']
            payloadreads+=len(payload)
            value=json.loads(payload);assert len(value['terms'])==c['terms']and len(value['shift'])==5
            assert len({tuple(z[1:])for z in value['terms']})==c['terms']
            if k:write(b',')
            write(f"[{c['row']},{c['column']},".encode());offset=f.tell();write(payload);write(b']')
            end=f.tell();f.flush();f.seek(offset);copied=f.read(c['payload_bytes']);readbacks+=len(copied);f.seek(end)
            assert copied==payload and hashlib.sha256(copied).hexdigest()==c['payload_sha256']
            z={a:v for a,v in c.items()if a!='file_index'};z['payload_offset']=offset;z['source_input_cell_index']=k
            newcells.append(z);rowcells.setdefault(c['row'],[]).append(k);colcells.setdefault(c['column'],[]).append(k)
            rr=rowstats[c['row']];rr[0]+=1;rr[1]+=c['terms'];rr[2].update(payload);rr[3].update(copied)
            nterms+=c['terms'];del value,payload,copied
            if k%256==0:whole()
        suffix=b'],"total_terms":292652,"nonzero_entries":'+str(len(newcells)).encode()+b',"all_passed":true}\n';write(suffix);f.flush()
    assert nterms==292652 and len(newcells)==4949
    tmp.rename(matrix);matrixsha=sha(matrix)
    for i in index['all_original_row_ids']:
        n,t,a,z=rowstats[i];assert a.digest()==z.digest()
        rowchecks.append({'original_row':i,'cells':n,'terms':t,'source_payload_sequence_sha256':a.hexdigest(),'output_payload_sequence_sha256':z.hexdigest(),'all_cell_bytes_and_values_equal':True})
    compact={'source_file':str(matrix.relative_to(ROOT)),'source_file_sha256':matrixsha,
             'logical_shape':[196,43],'logical_coefficient_terms':292652,'all_literal_cells':newcells,
             'row_cell_indices':rowcells,'column_cell_indices':colcells,
             'all_original_row_ids':index['all_original_row_ids'],'all_original_column_ids':index['all_original_column_ids'],
             'source_input_index':str(ip.relative_to(ROOT)),'source_input_index_sha256':sha(ip),
             'all4949_payload_byte_equalities_checked':True,'all196_rows_preserved_including_zero':True}
    owned=numeric(index)+numeric(compact)+numeric(active)+12*548+4096
    assert owned<=94504<150000
    compactip=out/'literal-file-index.json';streamjson(compactip,compact,b['new_full_index_RAW_upper']);whole()
    # A second independent file read checks every output offset and literal hash.
    with matrix.open('rb')as f:
        for c in newcells:
            f.seek(c['payload_offset']);s=f.read(c['payload_bytes']);readbacks+=len(s)
            assert hashlib.sha256(s).hexdigest()==c['payload_sha256']
    receipt={'all4949_cells_literal_byte_and_value_equal':True,'all196_rows_equal':rowchecks,
             'current_generation_unchanged':40,'shape':[196,43],'terms':292652,'zero_rows_preserved':True,
             'source_payload_bytes_read':payloadreads,'output_payload_bytes_readback':readbacks,
             'owned_numeric_plus_cell_scratch_upper':owned,'no_coefficient_math_pivot_or_new_inverse':True,
             'all_prior_versions_remain_stored':True}
    streamjson(out/'cell-and-row-equality-receipt.json',receipt,100000)
    for f in active:assert sha(ROOT/f['path'])==f['sha256']
    for p in extra:assert sha(p)==pre['all_additional_READ_inputs'][str(p.relative_to(ROOT))]
    assert not list(out.glob('*.partial'));whole()
    report={'id':ID,'all_passed':True,'exact_current_generation_compacted':40,'shape':[196,43],'logical_terms':292652,
            'all4949_source_output_cell_equalities':True,'all196_row_equalities':True,
            'owned_numeric_and_scratch_upper':owned,'owned_cap':150000,'fresh_RAW_upper':32394791,
            'hard_address_space_limit_bytes':2000000000,'hard_seconds':60,'seconds':time.monotonic()-start,
            'retired_history_preserved_and_not_deleted':True,'previous_fresh_job_RAW_receipt_bytes':20889002,
            'current_snapshot':'current-matrix.json','current_snapshot_sha256':matrixsha,
            'current_index':'literal-file-index.json','current_index_sha256':sha(compactip),
            'input_sha256':{**{f['path']:f['sha256']for f in active},**pre['all_additional_READ_inputs']},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},
            'no_pivot_arithmetic_rank_or_exhaustion_claim':True,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    ledger={'actual_READ_input_bytes':inputbytes,'additional_READ_bytes':extrabytes,'new_artifact_bytes':0,'fresh_job_TOTAL_RAW_bytes':0,
            'upper_bytes':32394791,'all_old_versions_preserved':True,'no_gzip_or_artificial_deletion':True}
    lp=out/'final-RAW-ledger.json';report['evidence_sha256'][lp.name]='0'*64
    def size(v):return len((json.dumps(v,separators=(',',':'))+'\n').encode())
    existing=treebytes(out)
    for _ in range(10):
        artifacts=existing+size(ledger)+size(report);total=inputbytes+extrabytes+artifacts
        if ledger['fresh_job_TOTAL_RAW_bytes']==total:break
        ledger['new_artifact_bytes']=artifacts;ledger['fresh_job_TOTAL_RAW_bytes']=total
    else:raise AssertionError('RAW accounting length fixed point')
    assert total<=32394791<32*1024*1024
    streamjson(lp,ledger,10000);report['evidence_sha256'][lp.name]=sha(lp);streamjson(out/'report.json',report,200000)
    assert whole()==ledger['fresh_job_TOTAL_RAW_bytes'];signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('RAW',total,'report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
