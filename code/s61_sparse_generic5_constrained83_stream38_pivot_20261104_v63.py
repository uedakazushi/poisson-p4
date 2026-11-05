#!/usr/bin/env python3
"""ONE conditionally authorized mixed-file streamed pivot from exact38."""
from pathlib import Path
import gc,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM38_PIVOT_20261104_v63'
GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM38_SOURCE_GATE_20261104_v63'
PREV='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW_STREAM_PIVOT_20261104_v63'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def treebytes(p):return sum(x.stat().st_size for x in p.rglob('*')if x.is_file())
def prepare(out,work):
    gd=ROOT/'reports'/GATE;gr=gd/'report.json';gp=gd/'candidate-00-full-gate.json';ip=gd/'mixed-literal-file-index.json'
    assert sha(gr)=='1316ebae63acb4dc353819e7ba00d4dc0dab66233f9a4b390b37cabc11fd7069'
    g=read(gr);assert g['all_passed']and g['source_only']and g['FULL_safe_found']
    for p in(gp,ip):assert sha(p)==g['evidence_sha256'][p.name]
    gate=read(gp);index=read(ip);mp=ROOT/index['source_manifest']
    assert sha(mp)==index['source_manifest_sha256']=='e2a4db6346a784e150d77c50fa7f801baba476a7d35422c4112960f874d2716a'
    assert gate['source_indexed_pivot']['row']==15 and gate['source_indexed_pivot']['column']==42
    assert gate['source_indexed_pivot']['monomial_QQ']=='-16/405'
    assert gate['source_indexed_pivot']['monomial_exponent']==[0,0,0,0,1]
    assert gate['whole_owned_scalar_and_scratch_upper']==75719 and gate['RAW_TOTAL_bound']==27981532
    for f in index['files']:
        p=ROOT/f['path'];assert sha(p)==f['sha256']and p.stat().st_size==f['bytes']
    cpp=ROOT/f'code/{Path(__file__).stem}.cpp';note=ROOT/f'research/notes/{ID}.md'
    header=ROOT/'code/s61_sparse_generic5_row_stream_kernel_20261104_v63.hpp'
    assert 'using Matrix'not in cpp.read_text()and 'using Matrix'not in header.read_text()
    sourcegatecode=ROOT/'code/s61_sparse_generic5_constrained83_stream38_source_gate_20261104_v63.py'
    sourcegatenote=ROOT/f'research/notes/{GATE}.md'
    otherexternal=treebytes(gd)+sourcegatecode.stat().st_size+sourcegatenote.stat().st_size-gate['RAW_current_gate_index_source_bytes']
    otherexternal+=sum(p.stat().st_size for p in[Path(__file__),cpp,note])
    assert otherexternal>0
    control=work/'control.txt';rows=gate['all_row_forecasts'];cells=index['all_literal_cells'];files=index['files']
    with control.open('w')as o:
        for p in(mp,ip,gp):o.write(str(p)+' '+sha(p)+'\n')
        o.write(str(otherexternal)+'\n'+str(len(files))+'\n')
        for f in files:o.write(f"{ROOT/f['path']} {f['sha256']} {f['bytes']}\n")
        for key in('all_original_row_ids','all_original_column_ids'):
            v=index[key];o.write(str(len(v))+' '+' '.join(map(str,v))+'\n')
        o.write(str(len(cells))+'\n')
        for c in cells:o.write(' '.join(map(str,[c['file_index'],c['row'],c['column'],c['payload_offset'],c['payload_bytes'],c['terms'],c['payload_sha256']]))+'\n')
        o.write(str(len(rows))+'\n')
        for r in rows:
            o.write(' '.join(map(str,[r['row'],r['source_terms'],r['new_row_terms_upper'],r['new_immutable_row_file_bytes_upper'],r['factor_reference_cell_index'],len(r['updates'])]))+'\n')
            for u in r['updates']:o.write(' '.join(map(str,u))+'\n')
    idxslots=6*len(cells)+len(cells)+len(index['row_cell_indices'])+sum(len(index[k])for k in('all_original_row_ids','all_original_column_ids'))+5*len(rows)+8*sum(len(r['updates'])for r in rows)+3*len(rows)+4*len(files)+512
    assert idxslots==68162
    normalcells=len(index['row_cell_indices']['15']);expected=sum(len(r['updates'])for r in rows)+normalcells
    pre={'id':ID,'prior_completed_steps':38,'one_pivot':[15,42,'-16*t0/405'],'source_arrays_released_before_native':True,
         'no_full_coefficient_matrix_opened':True,'hard_total_seconds':60,'native_soft_seconds':42,'hard_address_space_bytes':2000000000,
         'owned_algebraic_and_scratch_upper':7557,'index_plan_receipt_integer_slots_upper':68162,
         'whole_owned_scalar_and_scratch_upper':75719,'owned_cap':150000,'RAW_TOTAL_bound':27981532,
         'RAW_history_base':12681245,'RAW_charged_source_gate_base':1484995,'RAW_other_external_charge':otherexternal,
         'RAW_rows_cap':10660591,'RAW_normalization_cap':8973,'RAW_other_cap':3145728,
         'logical_input_terms':174569,'logical_support_upper':235333,'expected_back_add_checks':expected,
         'all_current_input_files':len(files),'all_current_cells':len(cells),'single_pivot_target_rows':118,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in[gr,gp,ip,mp,Path(__file__),cpp,note,header,sourcegatecode,sourcegatenote]},
         'control_sha256':sha(control)}
    dump(out/'preflight.json',pre);del rows,cells,files,index,gate,g
    return control
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;work=ROOT/'work'/ID;out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    control=prepare(out,work);gc.collect();cpp=ROOT/f'code/{Path(__file__).stem}.cpp';binary=work/'native';status={}
    try:
        with(out/'compile.log').open('wb')as log:
            subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint openssl) -lgmp','compile',str(cpp),str(binary)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,59-(time.monotonic()-start)))
        assert treebytes(out)+treebytes(work)+read(out/'preflight.json')['RAW_other_external_charge']<3145728
        with(out/'native.log').open('wb')as log:
            subprocess.run([str(binary),str(control),str(out),str(work)],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,59-(time.monotonic()-start)))
        n=read(out/'native-summary.json');assert n['all_passed']
        status={'one_worker_finished_verified':True,'completed_pivot':n['completed_pivot'],'updated_rows':n['updated_rows'],'total_seconds':time.monotonic()-start}
    except Exception as e:
        status={'one_worker_finished_verified':False,'failure':repr(e),'no_automatic_retry':True,'total_seconds':time.monotonic()-start};raise
    finally:dump(out/'execution-status.json',status)
    signal.alarm(0);print(json.dumps(status))
if __name__=='__main__':main()
