#!/usr/bin/env python3
"""ONE fixed known-cubic closure rank; no source/gauge/jet or GB."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import pivots,invq,mv,enc
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63'
PATHS={
 'preflight':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63/report.json',
 'input':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63/literal-known-cubic-input-and-exact24-row-reduction.json',
 'quad_peer':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_REMAINING_20261105_v63_REVIEW/report.json',
 'length11':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json',
 'Tnil_length':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert all(D[k]['all_passed']for k in ['preflight','quad_peer','length11','Tnil_length'])
 for p,h in D['preflight']['input_sha256'].items():assert sha(ROOT/p)==h
 assert D['preflight']['aggregate_QQ_upper']==11610
 (out/'phase.json').write_text(json.dumps({'phase':'all_fixed_inputs_pinned_before_rank','input_sha256':{p:sha(ROOT/p)for p in PATHS.values()}},indent=2)+'\n')
 I=D['input'];cols=[]
 for sp in I['equivalent_integer24_by1308_sparse_columns']:
  v=[Q(0)]*24
  for row,c in sp:v[row]=Q(c)
  cols.append(v)
 assert len(cols)==1308
 cp,cr=pivots(cols);rank=len(cp);iv,det=invq([[cols[j][i]for j in cp]for i in cr])if rank else([],Q(1));M=[[cols[j][i]for j in cp]for i in cr]
 for A,B in [(M,iv),(iv,M)]:assert all(sum((A[i][k]*B[k][j]for k in range(rank)),Q(0))==Q(i==j)for i in range(rank)for j in range(rank))
 factorhash=hashlib.sha256()
 for col in cols:
  f=mv(iv,[col[i]for i in cr])if rank else[]
  assert all(sum((cols[cp[k]][i]*f[k]for k in range(rank)),Q(0))==col[i]for i in range(24))
  factorhash.update((json.dumps(enc(f),separators=(',',':'))+'\n').encode())
 upper=24-rank;assert 1<=upper<=3
 exact_HF=[5,5,3,1,1,1]if upper==1 else None
 data={'source_known_cubic_shape':[600,1284],'exact_reduced_known_cubic_shape':[24,1308],
  'integer_row_scaling_denominators24':I['row_denominator_units24'],
  'fixed_rank':rank,'fixed_pivot_columns':cp,'fixed_pivot_rows':cr,'fixed_minor':M,'fixed_minor_determinant':det,'fixed_minor_inverse':iv,
  'both_inverse_products_verified':True,'all1308_column_factorizations_verified':True,'all1308_factorization_stream_sha256':factorhash.hexdigest(),
  'source_full_reduction':'Each of1284 actual known source columns is explicitly q-projected; all24 Id-qE consistency columns remain. Id-Eq has its600 typed divisor witnesses in known quadratic*variable relations, so exact quotient equivalence—not a projection-only rank guess—supplies the upper bound.',
  'known_cubic_original_rank':600-upper,'actual_HF3_upper':upper,'no_new_actual_cubic_source_relations_assumed':True,
  'actual_total_length_Tnil':16,'actual_HF_through2':[5,5,3],
  'actual_full_HF_if_forced_by_length':exact_HF,
  'length_completion_proof':'If upper1, m^3Tnil nonzero because total16 exceeds13; HF3=1. A graded module with one standard monomial in degree3 has each later HF<=1. Exactly3 length remains, so degrees3,4,5 each1 and degree6zero. Nakayama gives m^6Tnil=0.',
  'all_remaining_tail_options_if_not_forced':[h for h in [[3],[2,1],[1,1,1]]if h[0]<=upper],
  'new_source_gauge_jet_or_GB':False,'nil_generators_or_nilcube_not_inferred':True}
 p=out/'fixed-known-cubic-rank-and-length-completion.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 r={'id':STEM,'all_passed':True,'status':'PROVED_known_cubic_closure_and_conditional_length_completion','known_cubic_original_shape':[600,1284],'exact_reduced_shape':[24,1308],'known_reduced_rank':rank,'actual_HF3_upper':upper,'actual_full_HF_if_forced':exact_HF,'all1308_factorizations_and_both_inverse_products':True,'new_source_gauge_jet_or_GB':False,'nilcube_or_generators_decided':False,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()
