#!/usr/bin/env python3
"""Cached numerator/denominator degree and NF last-use metadata only."""
import collections,hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parents[1];S='S61_NONSQUARE_P100_LOW98_DEGREE_LASTUSE_PREFLIGHT_20261104_v63';O=R/'reports'/S
BASE=R/'reports/S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def get(p,hashes):hashes[str(p.relative_to(R))]=sha(p);return json.loads(p.read_text())
def profile(r):
 n,d=r['numerator'],r['denominator'];return {'zero':not n,'numerator_degree':max((m[0] for m,c in n),default=-1),'denominator_degree':max(m[0] for m,c in d),'stored_terms':len(n)+len(d)}
def main():
 O.mkdir(exist_ok=True);hashes={};binding=get(BASE/'source-binding.json',hashes);rows=[]
 for x in binding['nonpivot_sources10']:
  i=x['genuine_index'];d=get(BASE/f'pure-cubic84-{i}.json',hashes);rows.append({'role':'pure_kernel','genuine':i,'columns84':[profile(r) for r in d['pure84_values']]})
 p=get(R/'reports/S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63/cubic504.json',hashes)
 for b,row in zip(binding['six_pivot_bindings'],p['values6x84']):rows.append({'role':'pivot_tail','genuine':b['genuine_index'],'columns84':[profile(r) for r in row]})
 assert len(rows)==16 and all(len(r['columns84'])==84 for r in rows)
 cols=[{'column':j,'max_cached_numerator_degree':max(r['columns84'][j]['numerator_degree'] for r in rows),'max_cached_denominator_degree':max(r['columns84'][j]['denominator_degree'] for r in rows),'all16_cached_terms':sum(r['columns84'][j]['stored_terms'] for r in rows)} for j in range(84)]
 manifest=get(R/'reports/S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/manifest.json',hashes);fd=get(R/'reports/S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json',hashes);nodes=fd['linear_nodes'];need=manifest['NF_needed_linear_ids'];roots=manifest['NF_pure98_linear_roots'];uses=collections.Counter(roots)
 for i in need:
  if nodes[i][0]=='l':uses.update(j for c,j in nodes[i][1])
 live=set();peak=0;trace=[]
 for i in need:
  live.add(i);peak=max(peak,len(live))
  if nodes[i][0]=='l':
   for c,j in nodes[i][1]:
    uses[j]-=1
    if not uses[j]:live.discard(j)
  trace.append([i,len(live)])
 assert all(i in live for i in roots)
 e={'cached16x84_RF_profiles':rows,'per_column84_profiles':cols,'all_cached_RF_terms':sum(c['all16_cached_terms'] for c in cols),'max_numerator_degree':max(c['max_cached_numerator_degree'] for c in cols),'max_denominator_degree':max(c['max_cached_denominator_degree'] for c in cols),'NF_linear_nodes':len(need),'NF_last_use_peak_values':peak,'NF_schedule':trace,'NF_output_roots98':roots,'source_generation_may_stream_one_kernel_and_six_pivots':True,'NF_coefficient_degrees_or_global98_matrix_degrees_not_computed':True,'bareiss_polynomial_state_bound_still_missing':True,'no_new_coefficient_or_polynomial_operation':True}
 p=O/'metadata.json';p.write_text(json.dumps(e,separators=(',',':'))+'\n');report={'id':S,'source_only':True,'input_sha256':hashes,'evidence_sha256':{p.name:sha(p)},'max_cached_degrees':[e['max_numerator_degree'],e['max_denominator_degree']],'cached_RF_terms':e['all_cached_RF_terms'],'NF_last_use_peak_values':peak,'no_global98_rank_or_Bareiss_bound_claimed':True,'code_sha256':sha(Path(__file__))};(O/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
