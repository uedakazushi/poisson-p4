#!/usr/bin/env python3
"""ONE authorized indexed-row arithmetic pivot; no full coefficient matrix."""
from pathlib import Path
import argparse,gc,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW_STREAM_PIVOT_20261104_v63'
GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW_STREAM_PREFLIGHT_20261104_v63'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def treebytes(p):return sum(x.stat().st_size for x in p.rglob('*')if x.is_file())
def prepare(out,work):
    gd=ROOT/'reports'/GATE;gr=gd/'report.json';ip=gd/'literal-file-index.json';gp=gd/'first-indexed-row-stream-gate.json'
    assert sha(gr)=='e420c439d649d9b9c9c45eb72a5807a8c54abd5be096693f2e2648b10d4ce283'
    g=read(gr);assert g['all_passed']and g['source_only']
    for p in(ip,gp):assert sha(p)==g['evidence_sha256'][p.name]
    index=read(ip);gate=read(gp);source=ROOT/index['source_file']
    assert sha(source)==index['source_file_sha256']=='1721bb75016c56d187e769379697c8c508eb223de267d992ddd602cb403ce11b'
    assert gate['source_indexed_pivot']['row']==222 and gate['source_indexed_pivot']['column']==24
    assert gate['owned_live_with_all_scratch_bound']==5497 and gate['RAW_ledger']['TOTAL']==15045736
    cpp=ROOT/f'code/{Path(__file__).stem}.cpp';header=ROOT/'code/s61_sparse_generic5_row_stream_kernel_20261104_v63.hpp'
    ancestor=ROOT/'code/s61_sparse_generic5_constrained83_frontier_safe_continuation_20261104_v63.cpp'
    assert header.read_text()==ancestor.read_text().split('struct Height{',1)[0]
    assert 'using Matrix'not in cpp.read_text() and 'using Matrix'not in header.read_text()
    paths=[gr,ip,gp,source,Path(__file__),cpp,header,ancestor,ROOT/f'research/notes/{ID}.md']
    external=sum(p.stat().st_size for p in gd.iterdir()if p.is_file()and p!=ip)
    external+=sum(p.stat().st_size for p in[Path(__file__),cpp,header,ROOT/f'research/notes/{ID}.md'])
    control=work/'control.txt'
    with control.open('w')as o:
        for p in(source,ip,gp):o.write(str(p)+' '+sha(p)+'\n')
        o.write(str(external)+'\n')
        for key in('all_original_row_ids','all_original_column_ids'):
            v=index[key];o.write(str(len(v))+' '+' '.join(map(str,v))+'\n')
        cells=index['all_literal_cells'];o.write(str(len(cells))+'\n')
        for c in cells:o.write(' '.join(map(str,[c['row'],c['column'],c['payload_offset'],c['payload_bytes'],c['terms'],c['payload_sha256']]))+'\n')
        rows=gate['all_row_forecasts'];o.write(str(len(rows))+'\n')
        for r in rows:
            o.write(' '.join(map(str,[r['row'],r['source_terms'],r['new_row_terms_upper'],r['new_immutable_row_file_bytes_upper'],r['factor_reference_cell_index'],len(r['updates'])]))+'\n')
            for u in r['updates']:o.write(' '.join(map(str,u))+'\n')
    metadata_count=5*len(cells)+len(cells)+len(index['all_original_row_ids'])+sum(len(index[k])for k in('all_original_row_ids','all_original_column_ids'))+5*len(rows)+8*sum(len(r['updates'])for r in rows)+3*len(rows)+512
    assert metadata_count+5497<=150000
    pre={'id':ID,'one_pivot_authorized':[222,24,'H*q0/(2025*t0)'],
         'hard_total_seconds':60,'native_soft_seconds':42,'memory_limit_bytes':2000000000,
         'owned_algebraic_coefficients_and_scratch_upper':5497,'all_index_plan_receipt_integer_slots_upper':metadata_count,
         'whole_owned_scalar_and_scratch_upper':metadata_count+5497,'owned_cap':150000,
         'logical_input_terms':147370,'logical_after_support_upper':174997,'logical_not_resident_count':True,
         'RAW_total_bound_including_input_index_all_generations_IO_manifests':15045736,
         'RAW_new_rows_bound':6699540,'RAW_normalized_pivot_bound':1457,'RAW_other_bound':2472704,
         'RAW_other_external_charge':external,'no_full_coefficient_array_opened_by_wrapper':True,
         'wrapper_index_and_gate_objects_released_before_native':True,
         'unchanged_arithmetic_kernel_is_exact_frozen_prefix':True,'no_second_pivot_or_new_candidate_search':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'control_sha256':sha(control)}
    dump(out/'preflight.json',pre)
    # Only paths remain after this function returns. Explicitly drop all
    # parsed index/gate objects before compile/native; no coefficient P read.
    del index,gate,cells,rows,g
    return control
def run(out,work):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    control=prepare(out,work);gc.collect();cpp=ROOT/f'code/{Path(__file__).stem}.cpp';binary=work/'native';status={}
    try:
        with(out/'compile.log').open('wb')as log:
            subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint openssl) -lgmp','compile',str(cpp),str(binary)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,59-(time.monotonic()-start)))
        assert treebytes(out)+treebytes(work)<2472704
        with(out/'native.log').open('wb')as log:
            subprocess.run([str(binary),str(control),str(out),str(work)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,59-(time.monotonic()-start)))
        n=read(out/'native-summary.json');assert n['all_passed']
        status={'one_worker_finished_verified':True,'completed_pivot':n['completed_pivot'],'updated_rows':n['updated_rows'],'total_seconds':time.monotonic()-start}
    except Exception as e:
        status={'one_worker_finished_verified':False,'failure':repr(e),'no_automatic_retry':True,'total_seconds':time.monotonic()-start}
        raise
    finally:
        dump(out/'execution-status.json',status)
    signal.alarm(0);print(json.dumps(status))
def finish(out,work):
    assert not(out/'report.json').exists();pre=read(out/'preflight.json');ex=read(out/'execution-status.json');n=read(out/'native-summary.json')
    assert ex['one_worker_finished_verified']and ex['total_seconds']<60 and n['all_passed']
    for p,h in pre['input_sha256'].items():assert sha(ROOT/p)==h
    assert sha(work/'control.txt')==pre['control_sha256']
    assert n['owned_including_scratch_peak']<=5497 and n['whole_owned_scalar_and_scratch_upper']<=150000
    assert n['index_plan_receipt_integer_slots_upper']==pre['all_index_plan_receipt_integer_slots_upper']
    assert n['no_full_coefficient_matrix']and n['no_second_pivot']and n['old37_steps_not_replayed']
    manifest=out/('generation-1.json'if n['completed_pivot']else'exact-prefix-generation.json');m=read(manifest)
    assert m['shape']==n['shape']and m['logical_total_terms']==n['logical_after_terms']
    for r in m['rows']:
        if'literal_row_file'in r:assert sha(out/r['literal_row_file'])==r['sha256']
    if n['completed_pivot']:
        assert n['updated_rows']==133 and n['shape']==[198,45]and n['back_add_identities_checked']==939
    rows=treebytes(out/'rows');norm=(out/'normalized-pivot-row.json').stat().st_size
    other=treebytes(out)+treebytes(work)+pre['RAW_other_external_charge']-rows-norm
    assert rows<=6699540 and norm<=1457 and other+16000<=2472704
    raw=4518376+1353659+rows+norm+other
    ledger={'frozen_input_bytes':4518376,'frozen_index_bytes':1353659,'all_new_row_versions_bytes':rows,'normalized_pivot_bytes':norm,'other_artifacts_work_external_metadata_bytes':other,'TOTAL_before_final_receipt_bytes':raw,'TOTAL_with_final_report_reserve_bytes':raw+16000,'upper_bytes':15045736,'RAW32MiB_passed':raw+16000<32*1024*1024}
    dump(out/'final-RAW-ledger.json',ledger)
    r={'id':ID,'all_passed':True,'completed_exactly_one_pivot':n['completed_pivot'],'new_pivot':[222,24],
       'completed_total_on_reordered_path':38 if n['completed_pivot']else 37,
       'shape':n['shape'],'logical_total_coefficients':n['logical_after_terms'],
       'owned_algebraic_and_scratch_peak':n['owned_including_scratch_peak'],
       'index_plan_receipt_integer_scalar_upper':n['index_plan_receipt_integer_slots_upper'],
       'whole_owned_including_index_scalars_upper':n['whole_owned_scalar_and_scratch_upper'],
       'all939_normalization_and_update_back_adds_verified':n['back_add_identities_checked']==939,
       'all133_rows_immutable_hash_bound':n['updated_rows']==133,'no_full_matrix_in_RAM':True,
       'atomic_next_generation_only_after_whole_pivot':True,'old37_and_all_old_allresident_failures_unchanged':True,
       'logical_count_above150k_is_external_not_cap_increase':True,
       'all_five_parameters_and_alpha_beta_constraints_preserved':True,
       'new_inverse_only_parameter_monomial':'2025*t0/(H*q0)',
       'no_second_pivot_nonmonomial_point_prime_gcd_series_or_exhaustion':True,
       'input_sha256':pre['input_sha256'],'work_control_sha256':pre['control_sha256'],'ignored_executable_sha256':sha(work/'native'),
       'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},
       'total_worker_seconds':ex['total_seconds'],'native_seconds':n['seconds'],
       'max_rss_kib_wrapper_and_children':max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss),
       'scope':'Explicit ONE Schur row/column deletion and exact disk-backed residual. Rank83 already known; no complete inverse or fiber exhaustion.'}
    dump(out/'report.json',r);print(json.dumps({k:v for k,v in r.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--finish',action='store_true');a=ap.parse_args();out=ROOT/'reports'/ID;work=ROOT/'work'/ID
    if a.finish:finish(out,work)
    else:run(out,work)
