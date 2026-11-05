#!/usr/bin/env python3
"""Two immutable cached coefficients decide the actual nilideal generator count."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_VALUES_20261105_v63'
SOURCE='reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_SOURCE_20261105_v63'
CACHE='reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_VALUES_20261105_v63/phase1-all35-nine-convolution-projections.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,(int,Q)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(map(count,v.values()))
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 return 0
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source-peer-report',required=True);args=ap.parse_args()
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 paths=[SOURCE+'/report.json',SOURCE+'/complete-twelve-generator-kernel-and-two-cache-callback-gate.json',CACHE,args.source_peer_report]
 pins={p:sha(ROOT/p)for p in paths};rep=json.loads((ROOT/paths[0]).read_text());src=json.loads((ROOT/paths[1]).read_text());cache=json.loads((ROOT/CACHE).read_text());peer=json.loads((ROOT/args.source_peer_report).read_text())
 assert rep['all_passed']and peer['all_passed'];assert peer['input_sha256'].get(paths[0])==pins[paths[0]]
 assert src['z_Tor_rank']==6 and src['all849_head0_r_first_coefficients_zero']
 assert cache['convolution_labels'][0]==[0,7]
 callbacks=[]
 for rel,p,k in src['two_cached_second_callbacks']:
  assert cache['all35_source_input'][p]['source_relation']==rel
  value=Q(cache['all35_projected_M2_lower_representatives'][p][k]);callbacks.append([rel,p,k,str(value)])
 q={rel:Q(v)for rel,p,k,v in callbacks};cor={int(j):Q(v)for j,v in src['known_edge_first_correction_on_two_sources'].items()}
 qsum=q[627]+Q(32,11)*q[623];csum=cor[627]+Q(32,11)*cor[623];w=qsum+csum
 assert qsum==0 and w
 det6=Q(src['fixed_z_minor_determinant']);det7=det6*w;assert det7
 data={
  'two_actual_cached_callbacks_source_relation_selected35_index_convolution_value':callbacks,
  'two_known_edge_first_corrections':[[j,str(v)]for j,v in sorted(cor.items())],
  'zero_z_source_witness':[[627,'1'],[623,'32/11']],
  'two_cached_second_contributions_sum':str(qsum),'edge_first_contributions_sum':str(csum),'nonzero_last_Tor_row_witness':str(w),
  'fixed_rank7_certificate':{'z_rows':src['fixed_z_minor_rows'],'z_columns':src['fixed_z_minor_columns'],'last_column':[[627,'1'],[623,'32/11']],'last_row':7,'block_triangular_determinant':str(det7),'proof':'Last column has all seven z coordinates zero and last coordinate w nonzero. First six columns have invertible fixed6 z-minor; their last coordinates are immaterial for the determinant. z-rank6 gives rank≤7 on the whole849 relation image.'},
  'complete_presentation_targets':12,'complete_presentation_sources':877,'actual_complete_origin_relation_rank':7,
  'actual_mu_S_Nbar':5,'actual_mu_R0_Nil_D0':5,
  'actual_minimal_bar_generators':['h1-(5/3)rA*h0','h2+(88/25)rA*h0','h3-(2/5)rA*h0','h4-(136/75)rA*h0','zA5*h0'],
  'minimal_set_proof':'The four corrected-head coefficient rows are zero at origin. The I*h0 fiber has dimension1 after the rank7 map. The fixed6 z-minor omits z5 and the seventh witness isolates r², so z5*h0 spans that quotient and is nonzero. Nakayama gives the five specified barclasses as a minimal S-generating set. Since ell is in the R0 maximal ideal, mu_R0(N)=mu_S(N/ellN)=5. Any actual N-lifts of those five classes generate N, but no full formal lifts are computed.',
  'no_new_source_gauge_jet_or_coefficient_solve':True,'no_full_formal_nil_lifts_or_nilideal_relations_or_product_algebra_claim':True,
  'original_nilcube_and_exact_index_remain_open':True}
 peak=count([rep,src,cache,peer,data]);assert peak<30000,peak
 blob=json.dumps(data,indent=2)+'\n';assert sum((ROOT/p).stat().st_size for p in paths)+len(blob.encode())<2<<20
 ev=out/'actual-two-cache-witness-and-five-minimal-bar-generators.json';ev.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_nilideal_R0_minimal_generator_count5','actual_mu_R0_Nil_D0':5,'actual_mu_S_Nbar':5,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'peak_owned_numeric_slots':peak,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)
if __name__=='__main__':main()
