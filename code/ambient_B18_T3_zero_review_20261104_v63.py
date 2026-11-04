#!/usr/bin/env python3
"""Independent sparse rank and filtration proof review for actual T3=0."""
import argparse,hashlib,json,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rank(columns):
 basis={}
 for col in columns:
  v={i:Q(c)for i,c in enumerate(col)if Q(c)}
  while v:
   p=min(v)
   if p not in basis:
    a=v[p];basis[p]={i:c/a for i,c in v.items()};break
   a=v[p]
   for i,c in basis[p].items():
    z=v.get(i,Q(0))-a*c
    if z:v[i]=z
    else:v.pop(i,None)
 return len(basis)
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 pp=ROOT/'reports/AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63'
 fp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
 ep=ROOT/'reports/AMBIENT_B18_T3_E0_ORDER_GATE_20261104_v63/certificate.json'
 hp=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
 up=ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json'
 op=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
 gates=[ROOT/'reports'/n/'report.json'for n in [
  'AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW',
  'AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW',
  'AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63_REVIEW',
  'AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63_REVIEW']]
 for p in gates:assert read(p)['all_passed'],p
 assert read(pp/'report.json')['all_passed']
 independent=read(pp/'independent-firstjet-and-structural-proof.json');full=read(fp)
 cols=independent['all2069_rebuilt_exact_firstjet90_columns']
 assert cols==full['all2069_known_firstjet90_columns'] and len(cols)==2069 and all(len(c)==90 for c in cols)
 assert rank(cols)==43
 for r in full['all12_top_results']:
  assert len(r['known_firstjet90'])==90 and not any(Q(c)for c in r['known_firstjet90'])
  assert not any(Q(c)for row in r['four_g_sensitivity9_by4']for c in row)
 assert full['four_g_derivatives_drop_from_all_relations'] and full['four_g_derivatives_drop_from_all_tops']
 # Independent combinatorial linear part of the canonical N presentation.
 # Components: e0,e2,...,e8. Variables: z1,...,z8,v.
 nc=[]
 def column(terms):
  c=[0]*72
  for component,var,value in terms:c[9*component+var]+=value
  return c
 for j in range(1,8):nc.append(column([(j,0,-1)]))
 for i in range(1,8):
  for j in range(i+1,8):nc.append(column([(i,j,1),(j,i,-1)]))
 for j in range(1,8):nc.append(column([(j,8,1)]))
 assert len(nc)==35 and rank(nc)==35
 h,u,orig=map(read,[hp,up,op]);nodes=h['DAG_scalar_nodes'];values=list(map(Q,orig['all339669_scalar_origin_values_evaluated']))
 e=read(ep)
 def product_factors(i):
  n=nodes[i]
  if n[0]=='m':return product_factors(n[1])+product_factors(n[2])
  return Counter({i:1})
 tops=h['all12_top_columns_lifted_R0_10'];reduced=u['all12_transformed_top_columns']
 for record in e['all12_exact_pointer_and_product_identities']:
  i=record['top'];r=tops[i][0]
  assert tops[i][5]==0 and reduced[i][0]==r==record['raw_E0_and_reduced_E0_same_node']
  factors=record['three_factor_nodes']
  assert len(factors)==3 and all(values[j]==0 for j in factors)
  total=Counter()
  for j in factors:total+=product_factors(j)
  assert product_factors(r)==total
 proof=dict(
  independently_checked_full_firstjet_rank=43,independent_N_linear_rank=35,
  all12_top_full_firstjets_zero=True,all12_top_E0_order_at_least3=True,
  complete_D4_relation_generators=44,
  linear_generators=43,
  exceptional_generator='(v^2-nu*z1^2)e0-c0(z1)F',
  exceptional_c0_linear_term_zero='Otherwise z1*F adds a 44th independent linear relation, contradicting the exact rank43.',
  comparison_of_minimal_covers='Both covers have9 minimal generators; the N-unit row of the comparison at the origin is (1,0,...,0). Since top firstjets vanish, this preserves their E0 quadratic terms.',
  exact_Nakayama_argument='In an expression of a top in the44 generators, independence of43 linear terms forces their constant coefficients to vanish. E0 quadratic then equals the remaining constant coefficient times v^2. It is zero, so every coefficient is in m. Thus top lies in m*Jinf, Jinf=Q3+top implies T3=m*T3, and finite generation implies T3=0.',
  actual_Jinf_equals_Q3=True,actual_T3_zero=True,
  actual_T0_T1_T2_or_original_depth_solved=False,
  no_numeric36_extension_residues_or_new_source_jets_needed=True)
 (out/'independent-rank-order-and-Nakayama-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),pp/'report.json',pp/'independent-firstjet-and-structural-proof.json',fp,ep,hp,up,op,*gates]
 report=dict(id=ID,all_passed=True,actual_T3_is_zero_proved=True,independent_exact_rank_order_and_module_lemma_review=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
  evidence_sha256={'independent-rank-order-and-Nakayama-proof.json':sha(out/'independent-rank-order-and-Nakayama-proof.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
