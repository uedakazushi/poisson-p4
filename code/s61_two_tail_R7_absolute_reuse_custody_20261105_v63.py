#!/usr/bin/env python3
"""Publish exact-decomposition peer or prepare ONLY14cached-root R7 candidates."""
from pathlib import Path
import argparse,ast,hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(kind,out):
 out.mkdir(parents=True,exist_ok=False);pins={};runtime={}
 def pin(s,h=None,run=True):
  p=ROOT/s;v=sha(p)
  if h:assert v==h
  pins[s]=v
  if run:runtime[s]=v
  return p
 b='reports/S61_TWO_TAIL_H_EXACT_QQ_ROOT_DENSE_REVIEW_20261105_v63/'
 peer=json.loads(pin(b+'report.json','b7a6b60c04a288e54128f0df3df8a601152e8d58f576f9fbc0d6b0eb94fc8625').read_text())
 assert peer['all_passed'] and peer['independent_exact_reverse_dense_product_all350'] and peer['S4_geometrically_irreducible_proved']
 for p,h in peer['input_sha256'].items():pin(p,h)
 for p,h in peer['evidence_sha256'].items():pin(b+p,h)
 pin('research/notes/S61_TWO_TAIL_H_EXACT_QQ_ROOT_DENSE_REVIEW_20261105_v63.md',run=False)
 if kind=='peer':
  pin(str(Path(__file__).relative_to(ROOT)),run=False)
  rep=dict(id=out.name,all_passed=True,frozen=True,publication_only=True,exact_H_R7_times_S4_and_S4_absolute_peer_accepted=True,
   R7_absolute_NOT_claimed=True,no_arithmetic_replay=True,no_Fitting_rank_or_exhaustion_claim=True,input_sha256=pins)
  p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');print(sha(p));return
 fac=json.loads((ROOT/'reports/S61_TWO_TAIL_H_EXACT_QQ_ENTRY_20261105_v63/canonical-factor-candidates.json').read_text())['factors']
 R=next(f for f in fac if f['total_degree']==7);assert len(R['terms'])==29 and R['multiplicity']==1
 leading={tuple(e):int(v) for e,v in R['terms'] if e[0]==4}
 assert leading=={(4,1,2):9281250,(4,0,2):-3515625,(4,0,1):11250000}
 R0={}
 for (a,e,d),v in R['terms']:R0[a,e]=R0.get((a,e),0)+int(v)
 R0={k:v for k,v in R0.items() if v};assert max(a for a,e in R0)==4 and max(e for a,e in R0)==6
 assert {e:v for (a,e),v in R0.items() if a==4}=={0:7734375,1:9281250}
 const=[0]*7
 for (a,e,d),v in R['terms']:
  if a==0:assert d==0;const[e]=int(v)
 assert const[0]==-6408
 rem=const[:];quo=[0]*6
 for k in range(6,0,-1):
  t=rem[k];rem=[6*v for v in rem];quo=[6*v for v in quo];quo[k-1]+=t;rem[k]-=6*t;rem[k-1]-=5*t
 assert rem[0] and rem[0]%32003 and not any(rem[1:])
 back=[0]*7
 for e,v in enumerate(quo):back[e]+=5*v;back[e+1]+=6*v
 back[0]+=rem[0];assert back==[6**6*v for v in const]
 old='reports/S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_ENTRY_20261105_v63/'
 report=json.loads(pin(old+'report.json','3028597f392c2c16d3a1aa51065b59cb5702401df3ed2f93f045742dc2ad4fa6').read_text())
 cache=json.loads(pin(old+'complete-six-root-certificate.json',report['evidence_sha256']['complete-six-root-certificate.json']).read_text())
 matched=cache['matched_polynomial_factor_pair'];assert matched['mask']==3 and matched['complement']==60 and matched['LC_allocation']==0
 Q=matched['Q'];assert len(Q)==50 and all(v[1]==0 for v in Q)
 binding=[]
 for a in range(5):
  for e in range(10):
   v=R0.get((a,e),0);oldv=Q[10*a+e][0];assert (v-140625*oldv)%32003==0
   if e>=7:assert oldv==0
   binding.append(dict(a=a,e=e,R7_integer=v,old_Q_real=oldv,imaginary=Q[10*a+e][1],cross_remainder_zero=True))
 roots=[[v for v in cache['roots'][i][:7]] for i in [2,3,4,5]];assert all(len(r)==7 for r in roots)
 q=cache['q'];assert q==[7499,25984] and cache['quadratic_seed_certificate']['Legendre']==32002
 pin('reports/S61_TWO_TAIL_H_D1_MODFACTOR_ROOT_REVIEW_20261105_v63/report.json','a41a96ded2bb3fe8d0b17dd98a38c98253b9510ee8618343f4bb5e07317a7f5c')
 pin('reports/S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_ROOT_CODE_REVIEW_20261105_v63/report.json','c68e2870124f99016edcda9d7a4b21b2e2b8946095a53a63480d9a929add117e')
 proof=dict(source_only=True,no_root_seed_Hensel_or_candidate_values=True,
  global_geometric_a_content_unit='LC=140625*d*(80-25d+66ed); a0 purelye nonzero; excludes d and absolutelylinearL',
  d1_a_degree=4,e_degree_upper=6,leading_a4='140625*(55+66e)',scalar140625_punit=bool(140625%32003),
  a0_integer=const,pseudo_divisor=[5,6],scale=6**6,quotient=quo,remainder=rem[0],remainder_modp=rem[0]%32003,
  primitive_modp_and_geometric=True,full_old_Q_binding=binding,old_root_indices=[2,3,4,5],precision=7,
  cache_is_INPUT_root_residuals_and_product_reverified=True,no_cached_Hensel_replay=True,
  root_identity_base_product_upper=4900,subset_upper=29400,candidates_upper=30870,other_upper=9000,total_upper=75000,
  hard_products=100000,hard_owned=12000,RAW_cap=2097152,field_temp_bits=30,existing_integer_input_bits=28,
  subset_max_cells=1120,all_live_with6000control_reserve_upper=11000,
  candidates=14,full_cells_each=65,not_only_witness_tested_by_producer=True,
  positive_PRIMARY_requires_independent_peer=True,negative_only_fixed_R7_modslice_factor=True,
  no_Fitting_or_rank_claim=True,soft_seconds=45,hard_seconds=60,RAM=2000000000)
 (out/'full-root-Q-binding-content-and-budget.json').write_text(json.dumps(proof,indent=2)+'\n')
 worker=pin('code/s61_two_tail_R7_absolute_reuse_worker_20261105_v63.py',run=False);ast.parse(worker.read_text())
 pin('research/notes/S61_TWO_TAIL_R7_ABSOLUTE_REUSE_SOURCE_20261105_v63.md',run=False);pin(str(Path(__file__).relative_to(ROOT)),run=False)
 h=subprocess.run([sys.executable,str(worker),'--help'],capture_output=True,text=True,timeout=10);assert h.returncode==0
 (out/'worker-help.txt').write_text(h.stdout);(out/'worker-help.stderr').write_text(h.stderr)
 rep=dict(id=out.name,all_passed=True,frozen=True,prepared_only=True,no_root_or_candidate_value_worker_run=True,
  worker_sha256=sha(worker),runtime_input_sha256=runtime,input_sha256=pins,
  R7_d1_integer_coefficients=[[[a,e],str(v)] for (a,e),v in sorted(R0.items())],
  old_Q_e_less7=[Q[10*a+e] for a in range(5) for e in range(7)],four_cached_root_jets=roots,quadratic_field=q,
  evidence_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()})
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n')
 raw=sum(x.stat().st_size for x in out.iterdir())+sum((ROOT/p).stat().st_size for p in pins);assert raw+1048576+131072<2097152
 print(json.dumps(dict(SHA256=sha(p),worker_SHA=sha(worker),RAW_prior=raw+1048576+131072,source_remainder=rem[0],remainder_modp=rem[0]%32003,values_run=False),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kind',choices=['peer','source'],required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.kind,a.output.resolve())
