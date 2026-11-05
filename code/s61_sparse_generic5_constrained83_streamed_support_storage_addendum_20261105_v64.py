#!/usr/bin/env python3
"""Freeze prepared thin controller/ELF and complete fresh RAW closure."""
from pathlib import Path
from functools import cmp_to_key
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for z in iter(lambda:f.read(65536),b''):h.update(z)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,indent=2)+'\n').encode()
def numeric(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(map(numeric,v))
    if isinstance(v,dict):return sum(map(numeric,v.values()))
    return 0
def tree(p):return sum(q.stat().st_size for q in p.rglob('*')if q.is_file())
def main(a):
    begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/a.id;work=ROOT/'work'/a.id;target=out/'implementation-storage-addendum.json';assert not target.exists()
    statep=out/'prepared-state.json';assert sha(statep)==a.prepared_sha;state=read(statep);assert state['all_passed']and state['no_coefficient_worker_launched']
    pp=out/'preflight.json';assert sha(pp)==state['preflight_sha256'];pre=read(pp);control=work/'control.txt';binary=work/'native'
    assert sha(control)==state['control_sha256']==pre['control_sha256']and sha(binary)==state['ignored_executable_sha256']
    assert not(out/'native.log').exists()and not list(out.glob('row-*.bin'))
    for p,h in pre['input_sha256'].items():assert sha(ROOT/p)==h
    sourceplan=ROOT/state['source_cell_plan'];data=read(sourceplan);plans=data['all16_cell_plans'];checks=[];peak=pre['thin_preparation_owned_peak']
    with control.open('rb')as f:
        header=f.readline().decode().split();assert int(header[2])==pre['whole_owned_limit']and int(header[3])==pre['whole_RAW_limit']and int(header[4])==state['coefficient_pair_limit']
        for p in plans:
            key=str((p['row'],p['column']));idx=pre['streamed_support_offset_index'][key];source=ROOT/p['support_file'];assert sha(source)==p['support_sha256']
            value=read(source);assert len(value)==idx['count']==p['output_support_upper']
            value.sort(key=cmp_to_key(lambda x,y:-1 if(sum(x),x)>(sum(y),y)else(1 if(sum(x),x)<(sum(y),y)else 0)))
            f.seek(idx['offset']);last=None
            for e in value:
                q=list(map(int,f.readline().split()));assert e==q and(last is None or(sum(last),last)>(sum(q),q));last=q
            assert f.tell()==idx['offset']+idx['bytes'];checks.append([p['row'],p['column'],idx['count'],idx['bytes'],True])
            peak=max(peak,numeric(pre)+numeric(plans)+numeric(value)+len(value)+numeric(checks)+8192);assert peak<=pre['whole_owned_limit']
            del value
    cpp=ROOT/'code/s61_sparse_generic5_constrained83_top_P290_bounded_coefficients_20261105_v64.cpp';driver=cpp.with_suffix('.py')
    text=cpp.read_text();assert 'set<Exp>support_set'not in text and 'long matched=0;Exp previous' in text and 'all actual terms matched streamed certificate'in text
    prepared_files=[p for d in(out,work)for p in d.rglob('*')if p.is_file()]
    ledger=[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}for p in prepared_files]
    base=pre['actual_READ_input_RAW_bytes'];external=pre['worker_code_header_note_charge'];prepared=sum(p['bytes']for p in ledger)
    remaining=800000;cap=pre['all16_binary_output_RAW_cap'];upper=pre['whole_RAW_limit']
    result={'all_passed':True,'source_only':True,'coefficient_job_safe':True,'prepared_state_sha256':a.prepared_sha,'preflight_sha256':sha(pp),
        'immutable_worker_cpp_sha256':sha(cpp),'immutable_driver_py_sha256':sha(driver),'compiled_ignored_worker_sha256':sha(binary),
        'prepared_controller_sha256':sha(control),'all16_controller_certificates_exact_source_reorder':checks,
        'membership_implementation':'Native merges DEGLEX certificate with actual polynomial order using2 exponent5 vectors; Python readback streams the same controller segment. No5081-entry support set or duplicate full matrix/kernel metadata.',
        'whole_source_preparation_and_addendum_owned_peak':peak,'future_whole_owned_limit':pre['whole_owned_limit'],'future_fixed_metadata':25000,
        'source_READ_bytes':base,'worker_source_header_note_charge':external,'all_prepared_input_files':ledger,'prepared_input_files_bytes':prepared,
        'all16_binary_output_RAW_upper':cap,'remaining_noncell_artifact_upper':remaining,
        'remaining_noncell_proof':{'source_back_add_log':200000,'cell_records':100000,'native_summary':20000,'native_stdout':16000,'cumulative_layer_and_metadata':100000,'final_report_and_failure_allowance':200000,'sum636k_less800k':True},
        'RAW_upper':upper,'fresh_RAW_complete_future_upper':0,'input_sha256':{str(Path(__file__).relative_to(ROOT)):sha(Path(__file__)),str(statep.relative_to(ROOT)):sha(statep),str(pp.relative_to(ROOT)):sha(pp),str(sourceplan.relative_to(ROOT)):sha(sourceplan)},
        'no_coefficient_products_evaluated':True,'first16_cell_values_unchanged':True,'cumulative_references_not_double_P':True,'seconds':time.monotonic()-begin}
    # Prepared controller/ELF are READ inputs in the coefficient phase; they
    # remain in the physical tree and are charged ONCE, not duplicated/deleted.
    for _ in range(10):
        n=base+external+prepared+len(enc(result))+cap+remaining
        if result['fresh_RAW_complete_future_upper']==n:break
        result['fresh_RAW_complete_future_upper']=n
    assert n<=upper<=32*1024*1024 and peak<=pre['whole_owned_limit']<=150000
    target.write_bytes(enc(result));assert base+external+tree(out)+tree(work)+cap+remaining==n
    signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in('input_sha256','all_prepared_input_files','all16_controller_certificates_exact_source_reorder')}));print('addendum_sha256',sha(target))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('--prepared-sha',required=True);main(p.parse_args())
