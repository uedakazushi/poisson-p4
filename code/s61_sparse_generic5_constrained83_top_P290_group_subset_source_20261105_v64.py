#!/usr/bin/env python3
"""SOURCE-only subset of completed immutable support forecasts; no products."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def read(p):return json.loads(p.read_text())
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def main(a):
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 out=ROOT/'reports'/a.id;out.mkdir(exist_ok=False);(out/'supports').mkdir();(out/'cell-plans').mkdir()
 src=ROOT/'reports'/a.source_id;rp=src/'report.json';assert sha(rp)==a.source_sha;r=read(rp);assert r['all_passed'] and r['source_only']
 ep=src/'sixteen-source-support-height-cell-plans.json';bp=src/'whole-budget.json'
 for p in(ep,bp):assert sha(p)==r['evidence_sha256'][p.name]
 frame=read(ep);b=read(bp);rows=frame['selected_rows'][:a.rows];plans=[dict(p) for p in frame['all16_cell_plans'] if p['row'] in rows];assert len(plans)==4*a.rows
 inputs={ROOT/p:h for p,h in r['input_sha256'].items()};inputs.update({rp:sha(rp),ep:sha(ep),bp:sha(bp),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{a.id}.md':sha(ROOT/f'research/notes/{a.id}.md')})
 for p in plans:
  original=ROOT/p['support_file'];assert sha(original)==p['support_sha256'];inputs[original]=sha(original);blob=original.read_bytes();target=out/'supports'/original.name;target.write_bytes(blob);assert sha(target)==p['support_sha256'];p['support_file']=str(target.relative_to(ROOT));(out/'cell-plans'/original.name).write_bytes(enc(p))
 for p,h in inputs.items():assert sha(p)==h
 frame={**frame,'selected_rows':rows,'all16_cell_plans':plans,'subset_of_completed_source_gate':str(rp.relative_to(ROOT)),'subset_of_completed_source_gate_sha256':sha(rp),'no_support_product_or_coefficient_replay':True}
 (out/'sixteen-source-support-height-cell-plans.json').write_bytes(enc(frame))
 caps=sum(p['RAW_literal_polynomial_upper'] for p in plans);steps=[s for p in plans for s in p['steps']];owned=max([p['binary_roundtrip_phase_upper'] for p in plans]+[s['whole_owned_scratch_upper'] for s in steps]);prep=numeric(plans)+2000+6*max(p['output_support_upper']for p in plans)+8192;owned=max(owned,prep)
 ledger=[{'path':str(p.relative_to(ROOT)),'sha256':h,'bytes':p.stat().st_size}for p,h in sorted(inputs.items())];base=sum(z['bytes']for z in ledger);source=sum(p.stat().st_size for p in out.rglob('*')if p.is_file());reserve=3*1024*1024;raw=base+source+caps+reserve;safe=owned<=56473 and raw<=32*1024*1024
 budget={**b,'FULL_safe_sixteen_cell_coefficient_job':safe,'actual_READ_inputs':ledger,'actual_READ_input_RAW_bytes':base,'all_source_plan_and_support_files_RAW_bytes':source,'all16_output_RAW_upper':caps,'fresh_RAW_upper':raw,'whole_owned_with_fixed_metadata_upper':owned,'future_preparation_owned_upper':prep,'coefficient_term_pair_operations_upper':sum(s['n_a']*s['n_b']for s in steps),'thin_cell_descriptor_count':len({(p['row'],s['source_column'])for p in plans for s in p['steps']}),'completed_support_subset_no_recalculation':True}
 (out/'whole-budget.json').write_bytes(enc(budget))
 report={'id':a.id,'all_passed':True,'source_only':True,'selected_rows':rows,'cells':len(plans),'FULL_safe_coefficient_job':safe,'source_owned_peak':numeric(frame)+numeric(b)+25000,'future_owned_upper':owned,'future_RAW_upper':raw,'coefficient_pair_upper':budget['coefficient_term_pair_operations_upper'],'no_coefficient_products_evaluated':True,'all_prior_source_bytes_preserved':True,'input_sha256':{z['path']:z['sha256']for z in ledger},'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},'seconds':time.monotonic()-begin}
 assert report['source_owned_peak']<150000;(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('--source-id',required=True);p.add_argument('--source-sha',required=True);p.add_argument('--rows',type=int,choices=[1,2,3,4],required=True);main(p.parse_args())
