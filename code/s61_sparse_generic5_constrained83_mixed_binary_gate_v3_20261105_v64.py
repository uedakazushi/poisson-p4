#!/usr/bin/env python3
"""Current mixed JSON/binary exact source gate, one cell at a time."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
from s61_sparse_generic5_constrained83_stream38_source_gate_20261104_v63 import index_row
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import MAGIC,encode,decode,canonical
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
def main(a):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    mode='_ACTIVE_FRAGMENT'if a.active_fragment else''
    id=f'S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM{a.generation}{mode}_MIXED_BINARY_GATE_20261105_v64'
    future=f'S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM{a.generation}{mode}_MIXED_BINARY_COMPACTION_20261105_v64'
    out=ROOT/'reports'/id;out.mkdir(parents=True,exist_ok=False);pd=ROOT/'reports'/a.source
    pr=pd/'report.json';mp=pd/'generation-1.json';assert sha(pr)==a.report_sha and sha(mp)==a.manifest_sha
    r=read(pr);assert r['all_passed']and r['completed_total_on_reordered_path']==a.generation
    m=read(mp);assert m['complete_pivot']and m['logical_total_terms']==r['logical_total_coefficients']
    oldip=ROOT/r['literal_input_index_path'];old=read(oldip);assert sha(oldip)==m['index_sha256']
    oldfiles=old.get('files')or[{'path':old['source_file'],'sha256':old['source_file_sha256'],'encoding':old.get('encoding','primitive-json')}]
    paths={pr,mp,oldip,Path(__file__),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/'code/s61_sparse_generic5_constrained83_stream38_source_gate_20261104_v63.py'}
    overlay=None
    if a.active_fragment:
        od=ROOT/'reports'/a.active_fragment;orp=od/'report.json';assert sha(orp)==a.active_fragment_sha
        ore=read(orp);assert ore['all_passed']and ore['source_generation']==a.generation
        op=od/ore['replacement_index'];assert sha(op)==ore['replacement_index_sha256'];overlay=read(op)
        assert overlay['source_manifest_sha256']==a.manifest_sha and overlay['source_index_sha256']==sha(oldip)
        required={k for row in m['rows']for k in row.get('indexed_source_cell_indices',[])}
        assert set(map(int,overlay['cells']))==required
        paths.update([orp,op])
    files=[];filemap={};cells=[];rc={};cc={};receipts=[]
    def fid(p,h,encoding):
        name=str(p.relative_to(ROOT))
        if name not in filemap:
            assert sha(p)==h;filemap[name]=len(files);files.append({'path':name,'sha256':h,'bytes':p.stat().st_size,'encoding':encoding});paths.add(p)
        else:assert files[filemap[name]]['sha256']==h
        return filemap[name]
    for row in m['rows']:
        i=row['original_row']
        if'literal_row_file'in row:
            p=pd/row['literal_row_file'];k=fid(p,row['sha256'],'primitive-json');new=index_row(p,i,k);assert sum(c['terms']for c in new)==row['terms']
        elif 'literal_binary_row_file'in row:
            p=pd/row['literal_binary_row_file'];q=pd/row['paired_index_file'];assert sha(q)==row['paired_index_sha256'];paths.add(q);ix=read(q)
            assert ix['encoding']=='primitive-int-exponent-v1'and ix['source_binary_sha256']==row['sha256']and ix['original_row']==i and ix['source_binary_file']==p.name
            k=fid(p,row['sha256'],'primitive-int-exponent-v1');new=[];end=8
            for e in ix['entries']:
                col,off,nbytes,h,terms,mono,exps=e;assert off==end;end+=nbytes
                c={'file_index':k,'row':i,'column':col,'payload_offset':off,'payload_bytes':nbytes,'payload_sha256':h,'terms':terms}
                if terms==1:assert mono is not None and len(exps)==5;c.update(monomial_QQ=mono,monomial_exponent=exps)
                else:assert mono is None and exps is None
                new.append(c)
            assert end==p.stat().st_size and sum(c['terms']for c in new)==row['terms']==ix['total_terms']
        else:
            new=[]
            for k in row['indexed_source_cell_indices']:
                c=old['all_literal_cells'][k]
                if overlay:
                    e=overlay['cells'][str(k)];assert e['payload_sha256']==c['payload_sha256']and e['payload_bytes']==c['payload_bytes']and e['terms']==c['terms']
                    j=fid(ROOT/overlay['fragment_file'],overlay['fragment_sha256'],'primitive-int-exponent-v1')
                    new.append({**c,'file_index':j,'payload_offset':e['payload_offset']})
                else:
                    f=oldfiles[c.get('file_index',0)];j=fid(ROOT/f['path'],f['sha256'],f.get('encoding','primitive-json'));new.append({**c,'file_index':j})
        for c in new:
            assert c['row']==i and c['column']in m['column_original_ids'];k=len(cells);cells.append(c);rc.setdefault(i,[]).append(k);cc.setdefault(c['column'],[]).append(k)
    for f in files:
        p=ROOT/f['path']
        if p.parent.name=='rows':
            q=p.parent.parent/'receipts'/(p.stem+'.json'if p.suffix=='.bin'else p.name);x=read(q);assert x['literal_row_sha256']==f['sha256']and x['all_cell_back_add_identities_checked'];paths.add(q);receipts.append({'path':str(q.relative_to(ROOT)),'sha256':sha(q)})
    logical=m['logical_total_terms'];assert sum(c['terms']for c in cells)==logical
    overlay_slots=(numeric(overlay)+numeric(ore)+2*len(required)+512)if overlay else 0
    oldmeta=numeric(old)+numeric(m)+numeric(cells)+numeric(files)+numeric(rc)+numeric(cc)+overlay_slots
    del old,oldfiles,r
    offset=len(MAGIC);maxcell=0;maxbits=0
    for k,c in enumerate(cells):
        assert time.monotonic()-start<52
        f=files[c['file_index']];p=ROOT/f['path']
        with p.open('rb')as src:src.seek(c['payload_offset']);payload=src.read(c['payload_bytes'])
        assert hashlib.sha256(payload).hexdigest()==c['payload_sha256']
        value=decode(payload)if f['encoding']=='primitive-int-exponent-v1'else json.loads(payload)
        target=encode(value);back=decode(target);assert back==canonical(value)and len(back['terms'])==c['terms']
        if c['terms']==1:
            from fractions import Fraction
            assert c['monomial_QQ']==str(Fraction(value['content'])*int(value['terms'][0][0]))and c['monomial_exponent']==[value['shift'][v]+value['terms'][0][v+1]for v in range(5)]
        if f['encoding']=='primitive-int-exponent-v1':assert target==payload
        maxbits=max(maxbits,max(abs(e).bit_length()for e in value['shift']),max(e.bit_length()for t in value['terms']for e in t[1:]));maxcell=max(maxcell,c['terms'])
        c.update(source_file_index=c.pop('file_index'),source_offset=c.pop('payload_offset'),source_bytes=c.pop('payload_bytes'),source_sha256=c.pop('payload_sha256'),source_input_cell_index=k,
                 payload_offset=offset,payload_bytes=len(target),payload_sha256=hashlib.sha256(target).hexdigest())
        offset+=len(target);del payload,value,target,back
    assert maxbits<64
    plan={'generation':a.generation,'shape':m['shape'],'logical_terms':logical,'source_report':str(pr.relative_to(ROOT)),'source_report_sha256':a.report_sha,
          'source_manifest':str(mp.relative_to(ROOT)),'source_manifest_sha256':a.manifest_sha,'files':files,'all_literal_cells':cells,
          'row_cell_indices':rc,'column_cell_indices':cc,'all_original_row_ids':[x['original_row']for x in m['rows']],
          'all_original_column_ids':m['column_original_ids'],'exact_binary_file_bytes':offset,'max_cell_terms':maxcell,
          'active_row_receipts':receipts,'encoding':'primitive-int-exponent-v1'}
    pp=out/'mixed-binary-cell-plan.json';pp.write_bytes(enc(plan));paths.add(pp)
    active=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in sorted(paths)]
    futurecells=[{k:v for k,v in c.items()if not k.startswith('source_')or k=='source_input_cell_index'}for c in cells]
    template={'encoding':'primitive-int-exponent-v1','source_file':f'reports/{future}/current-matrix.bin','source_file_sha256':'0'*64,
              'logical_shape':m['shape'],'logical_coefficient_terms':logical,'all_literal_cells':futurecells,
              'row_cell_indices':rc,'column_cell_indices':cc,'all_original_row_ids':plan['all_original_row_ids'],'all_original_column_ids':m['column_original_ids'],
              'source_plan':str(pp.relative_to(ROOT)),'source_plan_sha256':sha(pp)}
    indexbytes=len(enc(template))+4096;own=max(oldmeta+12*maxcell+8192,numeric(plan)+numeric(template)+numeric(active)+24*maxcell+8192+overlay_slots)
    inputbytes=sum(f['bytes']for f in active);reserve=1024*1024;total=inputbytes+offset+indexbytes+reserve
    proof={'all_passed':True,'source_only':True,'generation':a.generation,'actual_READ_inputs':active,'actual_READ_input_RAW_bytes':inputbytes,
           'binary_snapshot_exact_bytes':offset,'binary_output_index_RAW_upper':indexbytes,'all_extra_worker_code_receipt_IO_reserve':reserve,
           'fresh_job_RAW_upper':total,'owned_numeric_codec_scratch_upper':own,'FULL_safe':own<=150000 and total<=32*1024*1024,
           'all_source_cells_roundtrip':len(cells),'logical_terms':logical,'all_old_text_and_binary_versions_retained':True,
           'no_pivot_coefficient_arithmetic_gzip_or_cap_increase':True}
    (out/'whole-budget.json').write_bytes(enc(proof))
    report={'id':id,'all_passed':True,'source_only':True,'generation':a.generation,'FULL_safe_binary_compaction':proof['FULL_safe'],
            'shape':m['shape'],'logical_terms':logical,'all_source_cell_roundtrips':len(cells),'owned_upper':own,'fresh_RAW_upper':total,'binary_exact_bytes':offset,
            'input_sha256':{f['path']:f['sha256']for f in active if f['path']!=str(pp.relative_to(ROOT))},
            'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,
            'configuration':vars(a)}
    (out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--generation',type=int,required=True);p.add_argument('--source',required=True);p.add_argument('--report-sha',required=True);p.add_argument('--manifest-sha',required=True);p.add_argument('--active-fragment');p.add_argument('--active-fragment-sha');main(p.parse_args())
