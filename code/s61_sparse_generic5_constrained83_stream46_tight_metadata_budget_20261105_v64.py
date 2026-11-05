#!/usr/bin/env python3
"""Metadata-only re-costing of the immutable exact45 binary source plan."""
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM46_TIGHT_METADATA_BUDGET_20261105_v64'
GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM46_MIXED_BINARY_GATE_20261105_v64'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    gd=ROOT/'reports'/GATE;gr=gd/'report.json';ob=gd/'whole-budget.json';pp=gd/'mixed-binary-cell-plan.json'
    assert sha(gr)=='8022943091d57c9f4c9115bab71fdbf5631a65f4a45b7258f3f817a123b58960'
    r=read(gr);assert r['all_passed']and not r['FULL_safe_binary_compaction']and r['generation']==46
    for p in(ob,pp):assert sha(p)==r['evidence_sha256'][p.name]
    b=read(ob);plan=read(pp);assert b['fresh_job_RAW_upper']==33952204 and plan['exact_binary_file_bytes']==8028307
    active=b['actual_READ_inputs'];hashmap={f['path']:f['sha256']for f in active}
    # These are exact fixed-width hash placeholders and actual row counts.
    rows=[]
    for i in plan['all_original_row_ids']:
        cells=[plan['all_literal_cells'][k]for k in plan['row_cell_indices'].get(str(i),[])]
        rows.append({'row':i,'cells':len(cells),'terms':sum(c['terms']for c in cells),'source_value_hash':'0'*64,'decoded_value_hash':'0'*64,'equal':True})
    rowreceipt={'all_cells_equal':True,'all_rows':rows,'logical_terms':514213,'no_polynomial_math_or_pivot':True,'all_prior_versions_retained':True}
    caps={'preflight':len(enc(active))+8192,'row_receipt':len(enc(rowreceipt))+4096,'report':len(enc(hashmap))+16384,'ledger':4096}
    additional=160000;reserve=524288
    fixed_extra=additional+sum(caps.values())
    assert fixed_extra<reserve
    b={**b,'original_budget_sha256':sha(ob),'all_extra_worker_code_receipt_IO_reserve':reserve,
       'fresh_job_RAW_upper':b['actual_READ_input_RAW_bytes']+b['binary_snapshot_exact_bytes']+b['binary_output_index_RAW_upper']+reserve,
       'FULL_safe':True,'writer_caps':caps,'additional_READ_files_upper':additional,
       'complete_fixed_format_metadata_bound':fixed_extra,'old_one_MiB_stop_preserved':True,
       'source_only_no_coefficient_roundtrip_or_pivot_rerun':True}
    assert b['fresh_job_RAW_upper']==33427916<32*1024*1024 and b['owned_numeric_codec_scratch_upper']<150000
    out=ROOT/'reports'/ID;out.mkdir(parents=True,exist_ok=False);bp=out/'whole-budget.json';bp.write_bytes(enc(b))
    worker=ROOT/'code/s61_sparse_generic5_constrained83_mixed_binary_compaction_v2_20261105_v64.py'
    paths=[gr,ob,pp,worker,Path(__file__),ROOT/f'research/notes/{ID}.md']
    report={'id':ID,'all_passed':True,'source_only':True,'generation':46,'FULL_safe_compaction':True,
            'fresh_RAW_upper':33427916,'owned_upper':b['owned_numeric_codec_scratch_upper'],'metadata_reserve':reserve,
            'fixed_format_metadata_bound':fixed_extra,'writer_caps':caps,'no_raw_limit_increase_or_coefficient_recompute':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{'whole-budget.json':sha(bp)},'seconds':time.monotonic()-start}
    (out/'report.json').write_bytes(enc(report));signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
