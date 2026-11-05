#!/usr/bin/env python3
"""Close the source input-parse ownership bound, without rerunning any producer."""
from pathlib import Path
import json,hashlib,time
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_SOURCE_OWNERSHIP_AUDIT_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
begin=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);rp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_UNIT_COL76_SOURCE_20261105_v64/report.json';assert sha(rp)=='349d1606feb479950eb690241bfdb95c49586731465b9f456e6656d20ea467cc';r=read(rp);assert r['all_passed'];inputs={rp:sha(rp),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')};records=[]
paths=[ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64/actual-fixed24-entry-source-and-exponent-plan.json',ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_THIRD_UNIT_SOURCE_20261105_v64/actual-third-unit-and-full174x21-typed-operator.json',ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR3_COL48_COEFFICIENTS_20261105_v64/actual-residual-col48-and-structural-frontier.json',ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR_COL81_COEFFICIENTS_20261105_v64/actual-residual-col81-and-structural-frontier.json',ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROWS11_44_LITERAL_SOURCE_20261105_v64/whole-cached-rows11-44-exact-comparison.json']+list(rp.parent.glob('*.json'))
for p in paths:
 h=sha(p)
 if p.parent==rp.parent and p.name!='report.json':assert h==r['evidence_sha256'][p.name]
 elif p!=rp:assert h==r['input_sha256'][str(p.relative_to(ROOT))]
 inputs[p]=h;records.append({'path':str(p.relative_to(ROOT)),'numeric_slots_if_all_JSON_retained_together':numeric(read(p))})
bundle=read(rp.parent/'col76-fixed-source-operands.json')['source_operands'];nleaf=sum(len(z['original_literal']['terms'])for z in bundle.values());nlabels=len(bundle);plans=read(rp.parent/'all20-cleared-col76-support-height-plans.json')['all20_source_plans'];supportmax=max([p['exact_union_support_cardinality']for p in plans]+[z['exact_product_support_cardinality']for p in plans for z in p['products']]);parsed=sum(z['numeric_slots_if_all_JSON_retained_together']for z in records);extras={'E_key_and_row_index_slots':2*576+1024,'additional_source_profile_support_tuple_slots':6*nleaf+8*nlabels,'simultaneous_S_T_union_tuple_support_buffers':5*4*supportmax,'scalar_bigint_QQ_work_and_report_objects':4096};bound=25000+parsed+sum(extras.values());assert bound<150000
v={'all_passed':True,'source_only':True,'producer_report_unchanged_sha256':sha(rp),'producer_38348_value_scope':'Recorded post-parse/post-GC tracked phases only. Use this larger proven wholeSOURCE upper for the input-parse phase as well.','all_parsed_JSON_simultaneous_upper':records,'all_parsed_JSON_numeric_sum':parsed,'fixed_metadata':25000,'additional_duplicate_and_stream_scratch_bounds':extras,'source_whole_numeric_upper':bound,'source_hard_limit':150000,'proof':'Conservatively retain ALL input/outputJSON numericobjects together, including the39487-slot old rowdifference object that source actuallydeletes before leaf/supportwork. Add distinctE keys, copied profile supports, allS/T/union buffers, and scalar scratch. This dominates every source phase. No arithmeticproducer/source-value rerun.','future_coefficient_owned_bound_28968_UNCHANGED':True,'all_coefficient_height_support_and_RAW_evidence_UNCHANGED':True}
(out/'whole-source-parse-and-support-ownership-proof.json').write_bytes(enc(v));report={'id':ID,'all_passed':True,'source_only':True,'whole_source_numeric_upper':bound,'future_coefficient_numeric_upper':r['future_coeff_owned_upper'],'producer_bytes_unchanged':True,'no_family_coefficient_arithmetic':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
