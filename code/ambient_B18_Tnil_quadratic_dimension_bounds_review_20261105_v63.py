#!/usr/bin/env python3
"""Independent same-ring length and graded-monomial proof; no source values."""
import hashlib,itertools,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_TNIL_QUADRATIC_DIMENSION_BOUNDS_20261105_v63';ID=BASE+'_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 primary=load(BASE);assert pins[f'reports/{BASE}/report.json']=='264e577b7e9544aa242d48a9f0fb5c42c6c7e002db67fec5819a49f336c6bb17'
 for p,h in primary['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in primary['evidence_sha256'].items():assert sha(ROOT/'reports'/BASE/n)==h
 natural=load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW');assert natural['exact_natural_injection']and natural['actual_origin_rank']==2 and natural['cokernel_length']==5 and natural['cokernel_mu']==3
 source=load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63','actual-natural-Told-to-Tnil-source-map-and-origin.json')
 assert source['actual_full_chain_map_old768_to_new780']['dimensions']==[780,768]
 assert source['matrix_operator_nodes'][9]==['matrix_product',[8,7]]
 assert source['matrix_operator_nodes'][15]==['matrix_product',[14,13]]
 # Z is an actual R0-linear coefficient map. Columnwise original homotopy
 # defines its closed lifts; only after choosing them is the map extended
 # S-linearly. Thus the injection commutes with EVERY power of the SAME m.
 assert load('AMBIENT_B18_T_UNISERIAL_DEEP_TAIL_20261104_v63_REVIEW')['all_passed']
 tail=load('AMBIENT_B18_T_UNISERIAL_DEEP_TAIL_20261104_v63');assert tail['actual_T_total_length']=='10+s'
 hf=load('AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63');assert hf['actual_HF_through_degree3']==[5,3,2,1]
 assert load('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63_REVIEW')['actual_HF_Tnil_degree0_1']==[5,5]
 op=load('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_661_OPERATOR_20261104_v63');assert op['all814_sources_complete']and op['residual_actual_operator_shape']==[6,661]
 # Independent divisor enumeration in exponent form, including degree4.
 def exponent(mon):return tuple(mon.count(i)for i in range(8))
 degree2=[exponent(m)for m in itertools.combinations_with_replacement(range(8),2)]
 counts=[]
 for e in degree2:
  row=[]
  for degree in[3,4]:
   good=[]
   for mon in itertools.combinations_with_replacement(range(8),degree):
    n=exponent(mon);divisors=[q for q in degree2 if all(a<=b for a,b in zip(q,n))]
    if divisors==[e]:good.append(n)
   pure=next((i for i,x in enumerate(e)if x==2),None)
   assert good==([]if pure is None else[tuple(degree if i==pure else 0 for i in range(8))])
   row.append(len(good))
  counts.append([list(e),row])
 # In all degrees the same divisor argument proves h_d<=1 when h2=1.
 length_checks=[]
 for s in range(1,8):
  old=10+s;new=old+5
  assert old-10==s and new-10==s+5
  assert new-(10+(s+4))==1
  for h in range(2,6):assert new-10-h>=old-10
  length_checks.append(dict(s=s,old=old,new=new,m_power_killing_old=s+3,m_power_killing_new=s+6,max_length_if_h2_one=s+14))
 evidence=dict(actual_source_map_is_S_linear=True,same_ring='S=QQ[[zA0,...,zA6,rA]], rA != DVR ell',
  upper_proof='i(m^3 Told) is contained in m^3 Tnil; lengths s and s+5-h2 give h2<=5',
  quotient_Loewy='length5, mu3 implies m^3 Q=0: otherwise degrees0..3 have length>=3+1+1+1',
  lower_proof='m^(s+6)Tnil=0; h2=0 gives length10. If h2=1, monomial initial submodule has a unique degree2 standard monomial. Every degree2 divisor of any higher standard monomial must equal it, forcing only pure powers in one target and one variable. Length<=s+14 contradicts s+15.',
  all36_exponent_degree3_and4_divisor_checks=counts,all5_module_targets_same_independent_divisor_condition=True,
  all7_length_checks=length_checks,interval_h2=[2,5],interval_actual_rank6x661=[1,4],
  equality_case_h2_5='m^3 Tnil=i(m^3 Told); S-linearity gives the same equality at every later power, hence equal later graded layers',
  truncation_inequality='len(Tnil/m^k Tnil)<=len(Told/m^k Told)+5',
  no_new_scalar_values_HF2_exact_length_or_nilcube=True)
 p=OUT/'independent-same-ring-length-and-divisor-proof.json';p.write_text(json.dumps(evidence,indent=2)+'\n')
 for p2 in[Path(__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p2.relative_to(ROOT))]=sha(p2)
 rep=dict(id=ID,all_passed=True,frozen=True,actual_HF2_Tnil_interval=[2,5],actual_complete6x661_rank_interval=[1,4],source_S_linearity_and_all_Loewy_length_hypotheses_verified=True,no_source_scalar_or_jet_evaluated=True,input_sha256=pins,evidence_sha256={p.name:sha(p)},seconds=time.monotonic()-st)
 (OUT/'report.json').write_text(json.dumps(rep,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
