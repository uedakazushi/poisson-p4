#!/usr/bin/env python3
"""ONE authorized streamed pivot from compact39; no full matrix array."""
from pathlib import Path
import gc,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT42_PIVOT_20261105_v64'
GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_COMPACT42_UNIT_GATE_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def treebytes(p):return sum(x.stat().st_size for x in p.rglob('*')if x.is_file())
def prepare(out,work):
    gd=ROOT/'reports'/GATE;gr=gd/'report.json';gp=gd/'candidate-01-full-gate.json'
    assert sha(gr)=='6cf47c3b9f1ddf740e264e54224932b92afd6e6cbe59d1126866c9a2621eada2'
    g=read(gr);assert g['all_passed']and g['source_only']and g['FULL_safe_found']
    assert sha(gp)==g['evidence_sha256'][gp.name];gate=read(gp)
    assert gate['source_indexed_pivot']['row']==10 and gate['source_indexed_pivot']['column']==69
    assert gate['source_indexed_pivot']['monomial_QQ']=='-131072/3645'and gate['source_indexed_pivot']['monomial_exponent']==[2,0,0,0,-1]
    active=gate['actual_READ_inputs'];assert sum(f['bytes']for f in active)==gate['actual_READ_input_RAW_bytes']==5511096
    for f in active:
        p=ROOT/f['path'];assert p.stat().st_size==f['bytes']and sha(p)==f['sha256']
    ip=next(ROOT/f['path']for f in active if f['path'].endswith('/literal-file-index.json'))
    index=read(ip);assert index['encoding']=='primitive-int-exponent-v1';matrix=ROOT/index['source_file'];assert sha(matrix)==index['source_file_sha256']
    cells=index['all_literal_cells'];rows=gate['all_row_forecasts'];cpp=ROOT/f'code/{Path(__file__).stem}.cpp';note=ROOT/f'research/notes/{ID}.md'
    header=ROOT/'code/s61_sparse_generic5_row_stream_kernel_20261104_v63.hpp'
    assert 'using Matrix'not in cpp.read_text()
    external=treebytes(gd)+sum(p.stat().st_size for p in[Path(__file__),cpp,note,header])
    control=work/'control.txt'
    with control.open('w')as o:
        for p in(matrix,ip,gp):o.write(str(p)+' '+sha(p)+'\n')
        o.write(str(external)+'\n')
        for key in('all_original_row_ids','all_original_column_ids'):
            v=index[key];o.write(str(len(v))+' '+' '.join(map(str,v))+'\n')
        o.write(str(len(cells))+'\n')
        for c in cells:o.write(' '.join(map(str,[c['row'],c['column'],c['payload_offset'],c['payload_bytes'],c['terms'],c['payload_sha256']]))+'\n')
        o.write(str(len(rows))+'\n')
        for r in rows:
            o.write(' '.join(map(str,[r['row'],r['source_terms'],r['new_row_terms_upper'],r['new_immutable_row_file_bytes_upper'],r['factor_reference_cell_index'],len(r['updates'])]))+'\n')
            for u in r['updates']:o.write(' '.join(map(str,u))+'\n')
    meta=5*len(cells)+len(cells)+len(index['row_cell_indices'])+sum(len(index[k])for k in('all_original_row_ids','all_original_column_ids'))+5*len(rows)+8*sum(len(r['updates'])for r in rows)+3*len(rows)+512
    assert meta==57371 and gate['whole_owned_scalar_and_scratch_upper']==70032
    pre={'id':ID,'prior_completed_steps':42,'one_pivot':[10,69,'-131072*H^2/(3645*t0)'],'source_objects_released_before_native':True,
         'no_full_coefficient_matrix':True,'hard_seconds':60,'native_soft_seconds':42,'hard_address_space_bytes':2000000000,
         'owned_algebraic_and_scratch_upper':12661,'index_plan_receipt_integer_slots_upper':57371,'whole_owned_upper':70032,
         'RAW_actual_read_input_bytes':5511096,'RAW_other_external_charge':external,'RAW_rows_cap':20934023,
         'RAW_normalization_cap':17305,'RAW_other_cap':3145728,'RAW_TOTAL_bound':29608152,
         'expected_back_add_checks':sum(len(r['updates'])for r in rows)+len(index['row_cell_indices']['10']),
         'input_sha256':{**{f['path']:f['sha256']for f in active},**{str(p.relative_to(ROOT)):sha(p)for p in[gr,gp,Path(__file__),cpp,note,header]}},
         'control_sha256':sha(control),'retired_history_retained_not_deleted_or_read':True}
    dump(out/'preflight.json',pre);del index,cells,rows,gate,g,active
    return control
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;work=ROOT/'work'/ID;out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    control=prepare(out,work);gc.collect();cpp=ROOT/f'code/{Path(__file__).stem}.cpp';binary=work/'native';status={}
    try:
        with(out/'compile.log').open('wb')as log:
            subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint openssl) -lgmp','compile',str(cpp),str(binary)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,69-(time.monotonic()-start)))
        assert treebytes(out)+treebytes(work)+read(out/'preflight.json')['RAW_other_external_charge']<3145728
        with(out/'native.log').open('wb')as log:
            subprocess.run([str(binary),str(control),str(out),str(work)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,69-(time.monotonic()-start)))
        n=read(out/'native-summary.json');assert n['all_passed']
        status={'one_worker_finished_verified':True,'completed_pivot':n['completed_pivot'],'updated_rows':n['updated_rows'],'total_seconds':time.monotonic()-start}
    except Exception as e:
        status={'one_worker_finished_verified':False,'failure':repr(e),'no_automatic_retry':True,'total_seconds':time.monotonic()-start};raise
    finally:dump(out/'execution-status.json',status)
    signal.alarm(0);print(json.dumps(status))
if __name__=='__main__':main()
