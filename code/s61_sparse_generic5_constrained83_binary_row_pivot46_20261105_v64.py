#!/usr/bin/env python3
"""ONE resource-gated exact pivot with binary row outputs and paired indices."""
from pathlib import Path
import gc,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_BINARY_ROW_PIVOT46_20261105_v64'
PROOF='S61_SPARSE_GENERIC5_CONSTRAINED83_BINARY_ROW_OUTPUT_GATE46_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def tree(p):return sum(x.stat().st_size for x in p.rglob('*')if x.is_file())
def prepare(out,work):
    pd=ROOT/'reports'/PROOF;pr=pd/'report.json';bp=pd/'whole-budget.json';wp=pd/'binary-row-writer-bounds.json'
    assert sha(pr)=='72e1b346938cb6a8032c7716e0fec5fdcb547ea3cd572f82afc8efb9ad64d5e8'
    r=read(pr);assert r['all_passed']and r['FULL_safe_binary_row_pivot']
    for p in(bp,wp):assert sha(p)==r['evidence_sha256'][p.name]
    b=read(bp);assert b['FULL_safe']and b['fresh_RAW_upper']==31181871 and b['whole_owned_upper']==93447
    gp=ROOT/b['source_text_gate'];assert sha(gp)==b['source_text_gate_sha256'];g=read(gp)
    assert g['source_indexed_pivot']['row']==76 and g['source_indexed_pivot']['column']==66
    active=b['actual_READ_inputs'];assert sum(f['bytes']for f in active)==9593321
    for f in active:assert sha(ROOT/f['path'])==f['sha256']
    ip=next(ROOT/f['path']for f in active if f['path'].endswith('/literal-file-index.json'));index=read(ip);matrix=ROOT/index['source_file']
    assert index['encoding']=='primitive-int-exponent-v1'and sha(matrix)==index['source_file_sha256']
    bound=read(wp);assert bound['source_gate_sha256']==sha(gp);writer={x['row']:x for x in bound['rows']}
    cpp=ROOT/f'code/{Path(__file__).stem}.cpp';note=ROOT/f'research/notes/{ID}.md';header=ROOT/'code/s61_sparse_generic5_row_stream_kernel_20261104_v63.hpp'
    finalizer=ROOT/'code/s61_sparse_generic5_constrained83_streamed_one_finalize_v4_20261105_v64.py'
    activepaths={f['path']for f in active};extra={pr,bp,Path(__file__),cpp,note,header,finalizer};extra={p for p in extra if str(p.relative_to(ROOT))not in activepaths}
    external=sum(p.stat().st_size for p in extra);control=work/'control.txt';plans=g['all_row_forecasts'];cells=index['all_literal_cells']
    with control.open('w')as o:
        for p in(matrix,ip,gp,bp,wp):o.write(str(p)+' '+sha(p)+'\n')
        o.write(str(external)+'\n')
        for key in('all_original_row_ids','all_original_column_ids'):
            v=index[key];o.write(str(len(v))+' '+' '.join(map(str,v))+'\n')
        o.write(str(len(cells))+'\n')
        for c in cells:o.write(' '.join(map(str,[c['row'],c['column'],c['payload_offset'],c['payload_bytes'],c['terms'],c['payload_sha256']]))+'\n')
        o.write(str(len(plans))+'\n')
        for p in plans:
            w=writer[p['row']]
            o.write(' '.join(map(str,[p['row'],p['source_terms'],p['new_row_terms_upper'],w['binary_file_upper']+w['paired_index_upper'],p['factor_reference_cell_index'],len(p['updates'])]))+'\n')
            for u in p['updates']:o.write(' '.join(map(str,u))+'\n')
            o.write(' '.join(map(str,[w['binary_file_upper'],w['paired_index_upper'],w['writer_owned_upper'],len(w['cells'])]))+'\n')
            for c in w['cells']:o.write(' '.join(map(str,c))+'\n')
    meta=g['index_plan_receipt_integer_slots_upper']+4*sum(len(w['cells'])for w in writer.values())+4*len(plans)+512
    assert meta==71309
    eventupper=sum(len(f"[{p['row']},{u[0]},{u[3]},true]\n")for p in plans for u in p['updates'])
    future=eventupper+4096+2048*len(plans)+2*180000+12000+512*len(plans)+200000+10000+5000+20000
    pre={'id':ID,'prior_completed_steps':46,'one_pivot':[76,66,'-56*s0/38443359375'],'hard_seconds':60,'native_soft_seconds':42,
         'hard_address_space_bytes':2000000000,'no_full_coefficient_matrix':True,'source_objects_released_before_native':True,
         'new_binary_row_output':True,'compact_event_schema':True,'finalizer_already_charged_in_external':True,
         'owned_algebraic_and_scratch_upper':22138,'index_plan_receipt_integer_slots_upper':71309,'whole_owned_upper':93447,
         'RAW_actual_read_input_bytes':9593321,'RAW_other_external_charge':external,'RAW_rows_cap':18394424,'RAW_normalization_cap':48398,
         'RAW_other_cap':3145728,'RAW_TOTAL_bound':31181871,'future_auxiliary_write_upper':future,
         'expected_back_add_checks':sum(len(p['updates'])for p in plans)+len(index['row_cell_indices']['76']),
         'binary_output_proof_path':str(bp.relative_to(ROOT)),'binary_output_bounds_path':str(wp.relative_to(ROOT)),
         'input_sha256':{**{f['path']:f['sha256']for f in active},**{str(p.relative_to(ROOT)):sha(p)for p in extra}},
         'control_sha256':sha(control),'all_prior_versions_preserved':True}
    dump(out/'preflight.json',pre);del b,g,r,index,plans,cells,writer,bound,active
    return control
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;work=ROOT/'work'/ID;out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    control=prepare(out,work);gc.collect();cpp=ROOT/f'code/{Path(__file__).stem}.cpp';binary=work/'native';status={}
    try:
        with(out/'compile.log').open('wb')as log:
            subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint openssl) -lgmp','compile',str(cpp),str(binary)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,59-(time.monotonic()-start)))
        f=read(out/'preflight.json');assert tree(out)+tree(work)+f['RAW_other_external_charge']+f['future_auxiliary_write_upper']<3145728
        with(out/'native.log').open('wb')as log:
            subprocess.run([str(binary),str(control),str(out),str(work)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,59-(time.monotonic()-start)))
        n=read(out/'native-summary.json');assert n['all_passed']
        status={'one_worker_finished_verified':True,'completed_pivot':n['completed_pivot'],'updated_rows':n['updated_rows'],'total_seconds':time.monotonic()-start}
    except Exception as e:
        status={'one_worker_finished_verified':False,'failure':repr(e),'no_automatic_retry':True,'total_seconds':time.monotonic()-start};raise
    finally:dump(out/'execution-status.json',status)
    signal.alarm(0);print(json.dumps(status))
if __name__=='__main__':main()
