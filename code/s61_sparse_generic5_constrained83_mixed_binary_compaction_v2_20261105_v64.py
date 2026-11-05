#!/usr/bin/env python3
"""One source-gated exact mixed-source binary compaction; no pivot."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import MAGIC,encode,decode,canonical
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def stream(p,v,cap):
    t=p.with_name(p.name+'.partial');assert not p.exists();n=0
    with t.open('wb')as f:
        for s in json.JSONEncoder(separators=(',',':')).iterencode(v):
            b=s.encode();assert n+len(b)<=cap;f.write(b);n+=len(b)
        assert n+1<=cap;f.write(b'\n')
    t.rename(p)
def tree(p):return sum(x.stat().st_size for x in p.rglob('*')if x.is_file())
def numeric(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(numeric(x)for x in v)
    if isinstance(v,dict):return sum(numeric(x)for x in v.values())
    return 0
def main(a):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    gate=f'S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM{a.generation}_MIXED_BINARY_GATE_20261105_v64'
    id=f'S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM{a.generation}_MIXED_BINARY_COMPACTION_20261105_v64'
    gd=ROOT/'reports'/gate;gr=gd/'report.json';bp=gd/'whole-budget.json';pp=gd/'mixed-binary-cell-plan.json';assert sha(gr)==a.gate_sha
    r=read(gr);assert r['all_passed']and r['generation']==a.generation
    supplement_paths=[]
    if a.supplement:
        sd=ROOT/'reports'/a.supplement;sr=sd/'report.json';assert sha(sr)==a.supplement_sha
        sp=read(sr);assert sp['all_passed']and sp['source_only']and sp['generation']==a.generation
        bp=sd/'whole-budget.json';assert sha(bp)==sp['evidence_sha256'][bp.name]
        b=read(bp);assert b['original_budget_sha256']==r['evidence_sha256']['whole-budget.json']
        supplement_paths=[sr,ROOT/f'research/notes/{a.supplement}.md']
    else:
        assert r['FULL_safe_binary_compaction'];assert sha(bp)==r['evidence_sha256'][bp.name];b=read(bp)
    assert sha(pp)==r['evidence_sha256'][pp.name]
    plan=read(pp);assert b['FULL_safe']and plan['generation']==a.generation
    caps=b.get('writer_caps',{'preflight':200000,'row_receipt':100000,'report':200000,'ledger':10000})
    active=b['actual_READ_inputs'];base=sum(f['bytes']for f in active);assert base==b['actual_READ_input_RAW_bytes']
    for f in active:assert sha(ROOT/f['path'])==f['sha256']
    out=ROOT/'reports'/id;out.mkdir(parents=True,exist_ok=False)
    extra={gr,bp,Path(__file__),ROOT/f'research/notes/{gate}.md',ROOT/f'research/notes/{id}.md',*supplement_paths};paths={f['path']for f in active};extra={p for p in extra if str(p.relative_to(ROOT))not in paths};ebytes=sum(p.stat().st_size for p in extra)
    def guard():
        n=base+ebytes+tree(out);assert n<=b['fresh_job_RAW_upper']and n<=32*1024*1024;return n
    pre={'id':id,'all_actual_READ_inputs':active,'additional_READ_input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in extra},
         'actual_READ_input_RAW_bytes':base,'additional_READ_bytes':ebytes,'fresh_RAW_upper':b['fresh_job_RAW_upper'],
         'owned_upper':b['owned_numeric_codec_scratch_upper'],'hard_seconds':60,'hard_address_space_bytes':2000000000,'configuration':vars(a)}
    stream(out/'preflight.json',pre,caps['preflight']);guard()
    matrix=out/'current-matrix.bin';tmp=out/'current-matrix.bin.partial';rows={i:[0,0,hashlib.sha256(),hashlib.sha256()]for i in plan['all_original_row_ids']};terms=0
    with tmp.open('w+b')as f:
        f.write(MAGIC)
        for k,c in enumerate(plan['all_literal_cells']):
            assert time.monotonic()-start<55
            source=plan['files'][c['source_file_index']]
            with(ROOT/source['path']).open('rb')as s:s.seek(c['source_offset']);payload=s.read(c['source_bytes'])
            assert hashlib.sha256(payload).hexdigest()==c['source_sha256']
            value=decode(payload)if source['encoding']=='primitive-int-exponent-v1'else json.loads(payload)
            binary=encode(value);assert len(binary)==c['payload_bytes']and hashlib.sha256(binary).hexdigest()==c['payload_sha256']
            if source['encoding']=='primitive-int-exponent-v1':assert payload==binary
            assert f.tell()==c['payload_offset']and f.tell()+len(binary)<=b['binary_snapshot_exact_bytes']
            f.write(binary);end=f.tell();f.flush();f.seek(c['payload_offset']);backbytes=f.read(c['payload_bytes']);f.seek(end);assert backbytes==binary
            back=decode(backbytes);assert back==canonical(value)and len(back['terms'])==c['terms']
            row=rows[c['row']];row[0]+=1;row[1]+=c['terms'];row[2].update(enc(canonical(value)));row[3].update(enc(back));terms+=c['terms']
            del payload,value,binary,backbytes,back
            if k%256==0:guard()
        assert f.tell()==b['binary_snapshot_exact_bytes']
    assert terms==plan['logical_terms'];tmp.rename(matrix);matrixsha=sha(matrix);checks=[]
    for i in plan['all_original_row_ids']:
        n,t,s,d=rows[i];assert s.digest()==d.digest();checks.append({'row':i,'cells':n,'terms':t,'source_value_hash':s.hexdigest(),'decoded_value_hash':d.hexdigest(),'equal':True})
    cells=[{k:v for k,v in c.items()if not k.startswith('source_')or k=='source_input_cell_index'}for c in plan['all_literal_cells']]
    with matrix.open('rb')as f:
        assert f.read(len(MAGIC))==MAGIC
        for c in cells:
            f.seek(c['payload_offset']);z=f.read(c['payload_bytes']);assert hashlib.sha256(z).hexdigest()==c['payload_sha256'];assert len(decode(z)['terms'])==c['terms']
    compact={'encoding':'primitive-int-exponent-v1','source_file':str(matrix.relative_to(ROOT)),'source_file_sha256':matrixsha,
             'logical_shape':plan['shape'],'logical_coefficient_terms':terms,'all_literal_cells':cells,
             'row_cell_indices':plan['row_cell_indices'],'column_cell_indices':plan['column_cell_indices'],
             'all_original_row_ids':plan['all_original_row_ids'],'all_original_column_ids':plan['all_original_column_ids'],
             'source_plan':str(pp.relative_to(ROOT)),'source_plan_sha256':sha(pp)}
    own=numeric(plan)+numeric(compact)+numeric(active)+24*plan['max_cell_terms']+8192;assert own<=b['owned_numeric_codec_scratch_upper']<=150000
    oi=out/'literal-file-index.json';stream(oi,compact,b['binary_output_index_RAW_upper']);guard()
    stream(out/'all-cell-row-roundtrip.json',{'all_cells_equal':True,'all_rows':checks,'logical_terms':terms,'no_polynomial_math_or_pivot':True,'all_prior_versions_retained':True},caps['row_receipt'])
    for f in active:assert sha(ROOT/f['path'])==f['sha256']
    assert not list(out.glob('*.partial'));guard()
    report={'id':id,'all_passed':True,'exact_current_generation_compacted':a.generation,'encoding':'primitive-int-exponent-v1',
            'shape':plan['shape'],'logical_terms':terms,'all_cell_roundtrips':len(cells),'all_row_equalities':len(checks),'owned_upper':own,
            'fresh_RAW_upper':b['fresh_job_RAW_upper'],'current_snapshot':'current-matrix.bin','current_snapshot_sha256':matrixsha,
            'current_index':'literal-file-index.json','current_index_sha256':sha(oi),'all_prior_versions_preserved':True,
            'no_gzip_coefficient_loss_pivot_or_exhaustion':True,'input_sha256':{**{f['path']:f['sha256']for f in active},**pre['additional_READ_input_sha256']},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    lp=out/'final-RAW-ledger.json';report['evidence_sha256'][lp.name]='0'*64
    ledger={'actual_READ_input_bytes':base,'additional_READ_bytes':ebytes,'new_artifact_bytes':0,'fresh_job_TOTAL_RAW_bytes':0,'upper_bytes':b['fresh_job_RAW_upper'],'all_prior_files_retained':True};existing=tree(out)
    for _ in range(10):
        new=existing+len(enc(ledger))+len(enc(report));total=base+ebytes+new
        if ledger['fresh_job_TOTAL_RAW_bytes']==total:break
        ledger['new_artifact_bytes']=new;ledger['fresh_job_TOTAL_RAW_bytes']=total
    else:raise AssertionError('RAW length fixed point')
    assert total<=b['fresh_job_RAW_upper'];stream(lp,ledger,caps['ledger']);report['evidence_sha256'][lp.name]=sha(lp);stream(out/'report.json',report,caps['report'])
    assert guard()==total;signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('RAW',total,'report_sha256',sha(out/'report.json'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--generation',type=int,required=True);p.add_argument('--gate-sha',required=True);p.add_argument('--supplement');p.add_argument('--supplement-sha');main(p.parse_args())
