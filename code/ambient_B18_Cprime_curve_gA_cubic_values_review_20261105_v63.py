#!/usr/bin/env python3
"""Independent eight rational products; no higher source values."""
import hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_VALUES_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def pairs(v):return {int(i):F(c)for i,c in v if F(c)}
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json',h=None):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p)
  if h:assert sha(p)==h
  return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_VALUES_20261105_v63';r=load(primary,h='4199daeceb3ea39f9277217684b9a980e5e1b93fa09ff094477cf482dd5e3179');assert r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 e=load(primary,'actual-curve-scalar-and-complete-positive-rank-certificate.json',h='f1ba364dbbcdc8f24bc93b7447a960e219d6e1cacd568d1b742a3db0b29beb83')
 assert load('AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_SOURCE_20261105_v63_REVIEW',h='a94f569de33cbf0ee236f759b08547d2223cb63f9e708aaf99d945ce8351f717')['single_scalar_worker_gate_ready']
 assert load('AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63_REVIEW',h='481862051b9e604ac92163b452f88e0cb46682d043393668b6648a89610e3559')['known_cubic_rank']==26
 source=load('AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_SOURCE_PREFLIGHT_20261105_v63','actual-curve-source-and-one-residual-cubic-functional.json');X=load('AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63','phase2-complete-streamed-cubic-upper-certificate.json');q=pairs(source['residual_compressed_dual27_sparse']);dual=pairs(source['residual_actual_cubic_dual495_sparse']);curve=pairs(source['actual_curve_relation495_sparse']);proj=[pairs(v)for v in X['all495_source_projection_columns27_sparse']]
 terms=[[i,c,dual.get(i,F(0)),c*dual.get(i,F(0))]for i,c in sorted(curve.items())];assert len(terms)==8 and terms==[[int(row[0])]+list(map(F,row[1:]))for row in e['all_eight_source_pairing_terms']]
 value=sum((row[3]for row in terms),F(0));compressed={}
 for i,c in curve.items():
  for j,a in proj[i].items():compressed[j]=compressed.get(j,F(0))+c*a
 compressed={j:c for j,c in compressed.items()if c};assert compressed==pairs(e['same_actual_curve_compressed_column27_sparse'])
 assert value==sum((c*compressed.get(i,F(0))for i,c in q.items()),F(0))==F(e['actual_residual_scalar'])==F(r['actual_residual_scalar'])==0
 assert all(sum((a*pairs(col).get(i,F(0))for i,a in q.items()),F(0))==0 for col in X['selected_original_compressed_columns'])
 assert F(e['old_fixed26_minor_determinant'])==F(-61,18000)and F(e['new_27minor_determinant_rows_old26_then25_columns_old26_then_curve'])==F(-61,18000)*value==0
 assert e['old2019_source_factorization_sha256']==X['all2019_fixed_factorization_sha256']and e['old2019_source_projection_sha256']==X['all2019_source_and_projection_stream_sha256']and not r['positive_Cprime_ell3_length9_and_long_exponent2']
 save('independent-eight-products-and-zero-scope.json',dict(all8_literal_products=terms,residual_scalar=value,compressed_column=compressed,full_and_compressed_values_equal=True,all26_known_columns_annihilated=True,old2019_source_and_factor_ancestry_pinned=True,actual_curve_source_initial_membership_accepted=True,candidate_adds_no_new_cubic_rank=True,partial_negative_only=True,no_actual_HF3_e_gt2_or_e2_conclusion=True,no_new_source_gauge_jet_or_rank_search=True))
 for p in[Path(__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir()if p.is_file());assert raw<32<<20 and time.monotonic()-st<60
 save('report.json',dict(id=ID,all_passed=True,frozen=True,actual_curve_residual_scalar='0',all8_products_full_compressed_and_source_ancestry_verified=True,only_this_candidate_negative=True,long_torsion_exponent_open=True,no_new_source_gauge_jet_or_rank_search=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},actual_RAW=raw,seconds=time.monotonic()-st))
 print(json.dumps(dict(all_passed=True,actual_scalar='0',partial_negative_only=True,seconds=time.monotonic()-st),indent=2))
if __name__=='__main__':main()
