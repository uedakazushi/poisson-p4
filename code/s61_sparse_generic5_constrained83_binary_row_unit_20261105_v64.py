#!/usr/bin/env python3
"""ONE generic source-gated exact parameter-monomial pivot; no search."""
from pathlib import Path
import argparse,gc,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for z in iter(lambda:f.read(65536),b''):h.update(z)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def tree(p):return sum(q.stat().st_size for q in p.rglob('*')if q.is_file())
def prepare(a,out,work):
    pd=ROOT/'reports'/a.proof;pr=pd/'report.json';bp=pd/'whole-budget.json'
    assert sha(pr)==a.proof_sha
    r=read(pr);assert r['all_passed']and r['FULL_safe_binary_row_pivot']and r['source_generation']==a.prior
    assert sha(bp)==r['evidence_sha256'][bp.name]
    b=read(bp);assert b['all_passed']and b['FULL_safe']and b['source_generation']==a.prior
    gp=ROOT/b['source_text_gate'];wp=ROOT/b['binary_bounds_file']
    assert sha(gp)==b['source_text_gate_sha256']and sha(wp)==b['binary_bounds_sha256']
    g=read(gp);bound=read(wp);assert bound['source_gate_sha256']==sha(gp)
    active=b['actual_READ_inputs'];assert sum(x['bytes']for x in active)==b['actual_READ_input_RAW_bytes']
    for f in active:assert sha(ROOT/f['path'])==f['sha256']
    ip=next(ROOT/f['path']for f in active if f['path'].endswith('/literal-file-index.json'))
    index=read(ip);matrix=ROOT/index['source_file'];assert index['encoding']=='primitive-int-exponent-v1'and sha(matrix)==index['source_file_sha256']
    assert len(index['all_original_row_ids'])==236-a.prior and len(index['all_original_column_ids'])==83-a.prior
    pivot=g['source_indexed_pivot'];prrow,pc=pivot['row'],pivot['column'];assert b['candidate']==[prrow,pc]
    writer={w['row']:w for w in bound['rows']};plans=g['all_row_forecasts'];cells=index['all_literal_cells']
    cpp=Path(__file__).with_suffix('.cpp');note=ROOT/f'research/notes/{a.id}.md'
    finalizer=ROOT/'code/s61_sparse_generic5_constrained83_streamed_one_finalize_v4_20261105_v64.py'
    extra={pr,bp,Path(__file__),cpp,note,ROOT/'code/s61_sparse_generic5_row_stream_kernel_20261104_v63.hpp',finalizer}
    activepaths={f['path']for f in active};extra={p for p in extra if str(p.relative_to(ROOT))not in activepaths};external=sum(p.stat().st_size for p in extra)
    control=work/'control.txt'
    with control.open('w')as o:
        for p in(matrix,ip,gp,bp,wp):o.write(str(p)+' '+sha(p)+'\n')
        o.write(str(external)+'\n')
        settings=[prrow,pc,a.prior,index['logical_coefficient_terms'],g['pivot_row_terms'],g['logical_after_support_union_upper'],b['algebraic_and_writer_scratch_upper'],g['charged_candidate_scalar_metadata'],b['actual_READ_input_RAW_bytes'],b['new_binary_row_files_RAW_upper']+b['new_paired_indices_RAW_upper'],b['normalized_pivot_RAW_upper'],b['other_control_code_executable_events_receipts_manifests_upper'],b['fresh_RAW_upper'],b['unit_QQ'],*b['unit_exponent']]
        o.write(' '.join(map(str,settings))+'\n')
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
    nu=sum(len(p['updates'])for p in plans);nbc=sum(len(w['cells'])for w in writer.values())
    meta=6*len(cells)+len(index['row_cell_indices'])+len(index['all_original_row_ids'])+len(index['all_original_column_ids'])+5*len(plans)+8*nu+3*len(plans)+512+4*nbc+4*len(plans)+512
    assert meta==g['index_plan_receipt_integer_slots_upper']==b['integer_index_plan_writer_metadata_slots']
    assert meta+b['algebraic_and_writer_scratch_upper']==b['whole_owned_upper']<=150000
    eventupper=sum(len(f"[{p['row']},{u[0]},{u[3]},true]\n")for p in plans for u in p['updates'])
    future=eventupper+4096+2048*len(plans)+360000+12000+512*len(plans)+200000+10000+5000+20000
    pre={'id':a.id,'prior_completed_steps':a.prior,'one_pivot':[prrow,pc,b['unit_QQ'],b['unit_exponent']],
         'hard_seconds':60,'native_soft_seconds':42,'hard_address_space_bytes':2000000000,
         'no_full_coefficient_matrix':True,'source_objects_released_before_native':True,'new_binary_row_output':True,
         'compact_event_schema':True,'finalizer_already_charged_in_external':True,
         'owned_algebraic_and_scratch_upper':b['algebraic_and_writer_scratch_upper'],
         'index_plan_receipt_integer_slots_upper':meta,'whole_owned_upper':b['whole_owned_upper'],
         'RAW_actual_read_input_bytes':b['actual_READ_input_RAW_bytes'],'RAW_other_external_charge':external,
         'RAW_rows_cap':b['new_binary_row_files_RAW_upper']+b['new_paired_indices_RAW_upper'],
         'RAW_normalization_cap':b['normalized_pivot_RAW_upper'],'RAW_other_cap':b['other_control_code_executable_events_receipts_manifests_upper'],
         'RAW_TOTAL_bound':b['fresh_RAW_upper'],'future_auxiliary_write_upper':future,
         'expected_back_add_checks':nu+len(index['row_cell_indices'][str(prrow)]),
         'binary_output_proof_path':str(bp.relative_to(ROOT)),'binary_output_bounds_path':str(wp.relative_to(ROOT)),
         'input_sha256':{**{f['path']:f['sha256']for f in active},**{str(p.relative_to(ROOT)):sha(p)for p in extra}},
         'control_sha256':sha(control),'all_prior_versions_preserved':True}
    dump(out/'preflight.json',pre)
    return control
def main(a):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/a.id;work=ROOT/'work'/a.id;out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    control=prepare(a,out,work);gc.collect();cpp=Path(__file__).with_suffix('.cpp');binary=work/'native';status={}
    try:
        with(out/'compile.log').open('wb')as log:
            subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint openssl) -lgmp','compile',str(cpp),str(binary)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,59-(time.monotonic()-start)))
        f=read(out/'preflight.json');assert tree(out)+tree(work)+f['RAW_other_external_charge']+f['future_auxiliary_write_upper']<f['RAW_other_cap']
        with(out/'native.log').open('wb')as log:
            subprocess.run([str(binary),str(control),str(out),str(work)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,59-(time.monotonic()-start)))
        n=read(out/'native-summary.json');assert n['all_passed']
        status={'one_worker_finished_verified':True,'completed_pivot':n['completed_pivot'],'updated_rows':n['updated_rows'],'total_seconds':time.monotonic()-start}
    except Exception as e:
        status={'one_worker_finished_verified':False,'failure':repr(e),'no_automatic_retry':True,'total_seconds':time.monotonic()-start};raise
    finally:dump(out/'execution-status.json',status)
    signal.alarm(0);print(json.dumps(status))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('--proof',required=True);p.add_argument('--proof-sha',required=True);p.add_argument('--prior',type=int,required=True);main(p.parse_args())
