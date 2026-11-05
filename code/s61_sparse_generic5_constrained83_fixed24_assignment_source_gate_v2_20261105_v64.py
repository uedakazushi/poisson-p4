#!/usr/bin/env python3
"""Freeze finite combinatorial fixed24 worker before it runs."""
from pathlib import Path
import json,hashlib,time
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_ASSIGNMENT_SOURCE_GATE_V2_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
begin=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);sp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64/report.json';assert sha(sp)=='a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a';r=read(sp);assert r['all_passed'];ep=sp.parent/'actual-fixed24-entry-source-and-exponent-plan.json';assert sha(ep)==r['evidence_sha256'][ep.name];v=read(ep);assert len(v['actual_entry_sources576'])==576
worker=ROOT/'code/s61_sparse_generic5_constrained83_fixed24_assignment_v2_20261105_v64.py';inputs={sp:sha(sp),ep:sha(ep),worker:sha(worker),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_ASSIGNMENT_20261105_v64.md':sha(ROOT/'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_ASSIGNMENT_20261105_v64.md')}
for z in v['actual_entry_sources576']:
 if z['terms']:
  p=ROOT/z['source']['file'];h=r['input_sha256'][str(p.relative_to(ROOT))];assert sha(p)==h;inputs[p]=h
# Everything is fixed integer-cost Hungarian on24x24, then exact dual and
# graph certificates. No row/weight search. Leading scans stream one entry.
n=24;lo=v['global_entry_exponent_lower5'];hi=v['global_entry_exponent_upper5'];w=v['injective_Kronecker_weights5'];cost_bound=sum(max(abs(lo[j]),abs(hi[j]))*w[j]for j in range(5));arithmetic_bound=(32*n*n+100)*(cost_bound+1);maxterms=max(z['terms']for z in v['actual_entry_sources576']);gaussian_nodes_upper=2+576+276+2*4324+2*276+24+24
# Separate cases: matching phase metadata <=48k +two raw-entry coefficient
# buffers. If unique, discard raw entry buffers and retain onlythin refs;
# allRFnodes have<=7numericfields and no expanded scalarcoefficients.
matching_owned=25000+2*numeric(v)+2*maxterms+576*20+24*128
inverse_owned=25000+7*gaussian_nodes_upper+2*numeric(v)+576*20
owned=max(matching_owned,inverse_owned);assert owned<=150000
base=sum(p.stat().st_size for p in inputs);output_upper=6*1024*1024;raw=base+output_upper+131072;assert raw<32*1024*1024
b={'all_passed':True,'source_only':True,'ONE_fixed_assignment_job_safe':True,'prepared_gate_READ_reserve':131072,'original_gate_preserved_no_arithmetic':'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_ASSIGNMENT_SOURCE_GATE_20261105_v64/report.json','source_plan':str(ep.relative_to(ROOT)),'source_plan_report_sha256':sha(sp),'immutable_worker_sha256':sha(worker),'actual_READ_input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'actual_READ_input_RAW_bytes':base,'whole_owned_numeric_upper':owned,'matching_phase_owned_upper':matching_owned,'unique_inverse_phase_owned_upper':inverse_owned,'max_actual_entry_terms':maxterms,'all_actual_terms_scanned':v['logical_terms_not_resident'],'integer_cost_abs_upper':cost_bound,'all_Hungarian_integer_arithmetic_abs_upper':arithmetic_bound,'integer_bit_upper':arithmetic_bound.bit_length(),'worst_Hungarian_inner_comparisons':24**4,'RF_node_upper_if_unique':gaussian_nodes_upper,'raw_output_upper':output_upper,'fresh_RAW_upper':raw,'hard_seconds':60,'hard_address_space_bytes':2000000000,'fixed_row_set_no_other_rank_or_matching_search':True,'no_series_or_numeric_parameter_specialization':True,'streaming_24_selected_lead_recheck':True}
(out/'whole-worker-source-gate.json').write_bytes(enc(b));report={'id':ID,'all_passed':True,'source_only':True,'ONE_fixed_assignment_job_safe':True,'owned_upper':owned,'fresh_RAW_upper':raw,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{'whole-worker-source-gate.json':sha(out/'whole-worker-source-gate.json')},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
