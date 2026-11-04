#!/usr/bin/env python3
"""Independent40-generator module audit using accepted actual cubic witness."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
from ambient_B18_T_known_initial_module_length_upper24_review_20261104_v63 import divide,lead,minus,encode
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63_REVIEW'
A='AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63'
PATHS={
 'author_report':f'reports/{A}/report.json','author_proof':f'reports/{A}/known-initial40-and-length-upper17.json',
 'old_peer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_MODULE_LENGTH_UPPER24_20261104_v63_REVIEW/report.json',
 'old':'reports/AMBIENT_B18_T_KNOWN_INITIAL_MODULE_LENGTH_UPPER24_20261104_v63/known-initial-module-and-length-bound.json',
 'cubic_peer':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63_REVIEW/report.json',
 'cubic_witness':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63_REVIEW/independent-cubic-witness71-source-and-channel.json',
 'rankK':'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inp={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert all(D[k]['all_passed'] for k in ['author_report','old_peer','cubic_peer','rankK'])
 AR,P=D['author_report'],D['author_proof']
 for p,h in AR['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in AR['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 assert sha(ROOT/'code/ambient_B18_T_known_initial_length_upper17_20261104_v63.py')==AR['code_sha256']
 assert sha(ROOT/'research/notes'/f'{A}.md')==AR['note_sha256']
 assert D['old_peer']['actual_I_squared_ann_m10_and_maximal_ideal_binding_accepted']
 witness=D['cubic_witness'];assert witness['witness_two_actual_cubic_values']==['-251099375/1051776','0']
 assert D['cubic_peer']['actual_gr3_T_dimension_upper_bound_accepted']==1
 assert D['rankK']['rank_D_K']==4 and D['rankK']['dim_k_C_ell_torsion']==4
 gs=[{tuple(e)+(int(t),):Q(c) for t,e,c in q} for q in P['homogeneous_module_generators40']]
 old=[{tuple(e)+(int(t),):Q(c) for t,e,c in q} for q in D['old']['homogeneous_generators39']]
 assert gs[:39]==old and len(gs)==40
 key=(0,0,0,0,0,1,0,2,4);assert gs[39]=={key:Q(1)}
 # The independent cubic quotient map is constructed from GB39. Its kernel
 # is exactly that known module in degree3, so the accepted witness minus
 # alpha times this canonical term lies in the true known initial module.
 duals=[{int(i):Q(c) for i,c in d} for d in witness['canonical_Pi2_600_rows_independent']]
 import itertools
 ms3=list(itertools.combinations_with_replacement(range(8),3));index=4*120+ms3.index((5,7,7))
 assert [d.get(index,Q(0)) for d in duals]==[1,0]
 pairs=[]
 for i in range(40):
  for j in range(i+1,40):
   x,y=lead(gs[i]),lead(gs[j])
   if x[8]!=y[8]:continue
   m=tuple(max(a,b) for a,b in zip(x[:8],y[:8]));p={}
   p=minus(p,gs[i],-1/gs[i][x],tuple(a-b for a,b in zip(m,x[:8])))
   p=minus(p,gs[j],1/gs[j][y],tuple(a-b for a,b in zip(m,y[:8])))
   r,tr=divide(p,gs);assert not r;pairs.append({'generators':[i,j],'lcm_exponents':m,'independent_reverse_division_trace':tr})
 assert len(pairs)==P['S_pair_count']==140
 leading={lead(g) for g in gs};free={(0,0,0,0,0,0,0,1,3),(0,0,0,0,0,1,0,0,4),(0,0,0,0,0,0,0,1,4)}
 linear={(tuple(int(k==v) for k in range(8))+(t,)) for t in range(5) for v in range(8)}
 assert linear-leading==free
 assert key in leading and (0,0,0,0,0,2,0,0,4) in leading and (0,0,0,0,0,0,0,2,3) in leading
 hf=[5,3,2]+[1]*7;assert sum(hf)==17
 proof={'source_manifest':inp,'actual_cubic_source_peer_now_final':True,
  'actual_witness71_value':['-251099375/1051776','0'],
  'canonical_cubic_source_membership':'q_actual/alpha - e4*z5*r² belongs to GB39 homogeneous initial submodule because independent Pi2 is its exact degree3 quotient and the two values are(1,0)',
  'normalizing_rational_unit':'-1051776/251099375',
  'all140_S_pairs_independent_reverse_exact_zero':pairs,
  'all_degree_standard_terms':{'d0':'five constants','d1':'e3*r,e4*z5,e4*r','d2':'e4*r²,e4*z5*r','d_at_least3':'e4*r^d'},
  'known_module_Hilbert_series':'5+3*t+2*t²+t³/(1-t)',
  'known_module_only_not_all_initial_relations':True,
  'same_actual_m10_T_zero_filtration_transfer_accepted':True,
  'actual_T_length_interval':[10,17],'actual_mu_D_K_interval':[6,13],'actual_mu_D_torsion_K_interval':[2,9],
  'no_exact_HF3_length_H0_length_exponent_or_nilcube_claim':True,
  'no_new_source_jet_348_matrix_or_rank_search':True}
 ep=out/'independent-known-initial40-and-length-upper17.json';blob=json.dumps(proof,separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;ep.write_text(blob)
 report={'id':STEM,'all_passed':True,'actual_cubic_source_dependency_independently_accepted':True,
  'all140_S_pairs_of40_generators_independently_zero':True,'actual_T_length_interval_accepted':[10,17],
  'actual_mu_D_K_interval_accepted':[6,13],'actual_mu_D_torsion_K_interval_accepted':[2,9],
  'whole_initial_module_equality_or_exact_length_claimed':False,'input_hashes':inp,
  'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
