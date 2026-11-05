#!/usr/bin/env python3
"""ONE exact binary current snapshot, every source cell roundtrip checked."""
from pathlib import Path
import hashlib,json,resource,signal,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import MAGIC,encode,decode,canonical
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM42_BINARY_COMPACTION_20261105_v64'
GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM42_BINARY_GATE_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def streamjson(p,v,cap):
    t=p.with_name(p.name+'.partial');assert not p.exists();n=0
    with t.open('wb')as f:
        for s in json.JSONEncoder(separators=(',',':')).iterencode(v):
            b=s.encode();assert n+len(b)<=cap;f.write(b);n+=len(b)
        assert n+1<=cap;f.write(b'\n')
    t.rename(p)
def tree(p):return sum(f.stat().st_size for f in p.rglob('*')if f.is_file())
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);gd=ROOT/'reports'/GATE;gr=gd/'report.json';bp=gd/'whole-budget.json';pp=gd/'binary-cell-plan.json'
    assert sha(gr)=='22a26cf1bdd8dcef4d31c4aa7b9dc93e1138ee3f34101768fea2c9e432558025'
    r=read(gr);assert r['all_passed']and r['FULL_safe_binary_compaction']
    for p in(bp,pp):assert sha(p)==r['evidence_sha256'][p.name]
    b=read(bp);assert b['FULL_safe']and b['fresh_job_RAW_upper']==29954386
    active=b['actual_READ_inputs'];base=sum(f['bytes']for f in active);assert base==b['actual_READ_input_RAW_bytes']
    for f in active:assert sha(ROOT/f['path'])==f['sha256']and(ROOT/f['path']).stat().st_size==f['bytes']
    plan=read(pp);ip=ROOT/plan['source_index'];index=read(ip);assert sha(ip)==plan['source_index_sha256']
    extra={gr,bp,Path(__file__),ROOT/f'research/notes/{ID}.md'};apaths={f['path']for f in active};extra={p for p in extra if str(p.relative_to(ROOT))not in apaths};ebytes=sum(p.stat().st_size for p in extra)
    def guard():
        n=base+ebytes+tree(out);assert n<=29954386 and n<32*1024*1024;return n
    pre={'id':ID,'all_actual_READ_inputs':active,'additional_READ_input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in extra},
         'actual_READ_input_RAW_bytes':base,'additional_READ_bytes':ebytes,'fresh_RAW_upper':29954386,
         'owned_upper':100852,'hard_seconds':60,'hard_address_space_bytes':2000000000,
         'method':'uncompressed exact signed magnitude integers and varint exponents, all current values retained'}
    streamjson(out/'preflight.json',pre,200000);guard()
    matrix=out/'current-matrix.bin';tmp=out/'current-matrix.bin.partial';rows={i:[0,0,hashlib.sha256(),hashlib.sha256()]for i in plan['all_original_row_ids']};count=0;total=0
    with tmp.open('w+b')as f:
        f.write(MAGIC)
        for c in plan['all_literal_cells']:
            assert time.monotonic()-start<55
            s=index['all_literal_cells'][c['source_input_cell_index']];p=ROOT/index['files'][s['file_index']]['path']
            with p.open('rb')as src:src.seek(s['payload_offset']);text=src.read(s['payload_bytes'])
            assert hashlib.sha256(text).hexdigest()==s['payload_sha256']
            source=json.loads(text);payload=encode(source);assert len(payload)==c['payload_bytes']and hashlib.sha256(payload).hexdigest()==c['payload_sha256']
            assert f.tell()==c['payload_offset']and f.tell()+len(payload)<=b['binary_snapshot_exact_bytes']
            f.write(payload);end=f.tell();f.flush();f.seek(c['payload_offset']);saved=f.read(c['payload_bytes']);f.seek(end)
            assert saved==payload;back=decode(saved);assert back==canonical(source)and len(back['terms'])==c['terms']
            row=rows[c['row']];row[0]+=1;row[1]+=c['terms'];row[2].update(enc(canonical(source)));row[3].update(enc(back));count+=1;total+=c['terms']
            del text,source,payload,saved,back
            if count%256==0:guard()
        assert f.tell()==b['binary_snapshot_exact_bytes']==4431853
    assert count==4734 and total==316841;tmp.rename(matrix);matrixsha=sha(matrix)
    checks=[]
    for i in plan['all_original_row_ids']:
        n,t,a,z=rows[i];assert a.digest()==z.digest();checks.append({'row':i,'cells':n,'terms':t,'source_value_hash':a.hexdigest(),'decoded_value_hash':z.hexdigest(),'equal':True})
    with matrix.open('rb')as f:
        assert f.read(len(MAGIC))==MAGIC
        for c in plan['all_literal_cells']:
            f.seek(c['payload_offset']);payload=f.read(c['payload_bytes']);assert hashlib.sha256(payload).hexdigest()==c['payload_sha256'];assert len(decode(payload)['terms'])==c['terms']
    compact={'encoding':'primitive-int-exponent-v1','source_file':str(matrix.relative_to(ROOT)),'source_file_sha256':matrixsha,
             'logical_shape':[194,41],'logical_coefficient_terms':316841,'all_literal_cells':plan['all_literal_cells'],
             'row_cell_indices':plan['row_cell_indices'],'column_cell_indices':plan['column_cell_indices'],
             'all_original_row_ids':plan['all_original_row_ids'],'all_original_column_ids':plan['all_original_column_ids'],
             'source_plan':str(pp.relative_to(ROOT)),'source_plan_sha256':sha(pp)}
    oi=out/'literal-file-index.json';streamjson(oi,compact,b['binary_output_index_RAW_upper']);guard()
    streamjson(out/'all-cell-row-roundtrip.json',{'all4734_cells_equal':True,'all194_rows':checks,'logical_terms':316841,'no_polynomial_math_or_pivot':True,'all_retired_text_bytes_preserved':True},100000)
    for f in active:assert sha(ROOT/f['path'])==f['sha256']
    assert not list(out.glob('*.partial'));guard()
    report={'id':ID,'all_passed':True,'exact_current_generation_compacted':42,'encoding':'primitive-int-exponent-v1','shape':[194,41],'logical_terms':316841,
            'all4734_cell_roundtrips':True,'all194_row_value_equalities':True,'owned_upper':100852,'fresh_RAW_upper':29954386,
            'current_snapshot':'current-matrix.bin','current_snapshot_sha256':matrixsha,'current_index':'literal-file-index.json','current_index_sha256':sha(oi),
            'all_prior_text_generations_and_stop_frontiers_preserved':True,'no_gzip_no_coefficient_loss_no_pivot':True,
            'input_sha256':{**{f['path']:f['sha256']for f in active},**pre['additional_READ_input_sha256']},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    lp=out/'final-RAW-ledger.json';report['evidence_sha256'][lp.name]='0'*64
    ledger={'actual_READ_input_bytes':base,'additional_READ_bytes':ebytes,'new_artifact_bytes':0,'fresh_job_TOTAL_RAW_bytes':0,'upper_bytes':29954386,'all_old_files_retained':True}
    existing=tree(out)
    for _ in range(10):
        new=existing+len(enc(ledger))+len(enc(report));t=base+ebytes+new
        if ledger['fresh_job_TOTAL_RAW_bytes']==t:break
        ledger['new_artifact_bytes']=new;ledger['fresh_job_TOTAL_RAW_bytes']=t
    else:raise AssertionError('RAW length fixed point')
    assert t<=29954386;streamjson(lp,ledger,10000);report['evidence_sha256'][lp.name]=sha(lp);streamjson(out/'report.json',report,200000)
    assert guard()==t;signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('RAW',t,'report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
