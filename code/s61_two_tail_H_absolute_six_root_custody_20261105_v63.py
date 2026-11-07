#!/usr/bin/env python3
"""Prepare ONLY source pins, loop budgets and inert help for six-root absolute-H proof."""
from pathlib import Path
import argparse,ast,hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 out.mkdir(parents=True,exist_ok=False);pins={};runtime={}
 def pin(s,h=None,run=True):
  p=ROOT/s;v=sha(p)
  if h:assert v==h
  pins[s]=v
  if run:runtime[s]=v
  return p
 src='reports/S61_TWO_TAIL_H_ABSOLUTE_HENSEL_SOURCE_20261105_v63/'
 r=json.loads(pin(src+'report.json','6b1ebfdb365ee0b24e6ebcd14e07b4e4a4e7e95feeecf28c6e8ed5c3b1fa4c99').read_text())
 assert r['all_passed'] and r['no_factor_or_Hensel_values_executed']
 d=json.loads(pin(src+'fixed-d1-input-and-formal-content.json',r['evidence_sha256']['fixed-d1-input-and-formal-content.json']).read_text())
 assert d['primitive_in_a_over_Fp_e_and_closure'] and d['remainder_mod32003']==24204
 assert d['source_a_degree']==6 and d['source_e_degree']==9 and len(d['canonical_integer_coefficients'])==42
 p=pin('reports/S61_TWO_TAIL_H_ABSOLUTE_HENSEL_ROOT_SOURCE_REVIEW_20261105_v63/report.json','523ac8adfaec65208092de73dcb724693b2b084f4c07df8a4062a29f51865113')
 assert json.loads(p.read_text())['all_passed']
 pin('reports/S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63/report.json','49990187b3b0b8edaf8939e4d13ac297f4f2ecba7f8b9f5326639089fa043ce2')
 pin('reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/factor-sparse.txt','241488fafa9b45e5353f0b707f75902bc5bc32e85ee093b9f2e0a5f07c7d571c')
 pin('reports/S61_TWO_TAIL_H_FIXED_RABIN_ENTRY_20261105_v63/report.json','8628312b3b266d5dd12107db3f996492eb98f0607b781411c21225bc838fee58')
 worker=pin('code/s61_two_tail_H_absolute_six_root_worker_20261105_v63.py',run=False)
 ast.parse(worker.read_text())
 pin('research/notes/S61_TWO_TAIL_H_ABSOLUTE_HENSEL_SOURCE_20261105_v63.md',run=False)
 pin('research/notes/S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_CUSTODY_20261105_v63.md',run=False)
 pin(str(Path(__file__).relative_to(ROOT)),run=False)
 h=subprocess.run([sys.executable,str(worker),'--help'],text=True,capture_output=True,timeout=10);assert h.returncode==0
 (out/'worker-help.txt').write_text(h.stdout);(out/'worker-help.stderr').write_text(h.stderr)
 proof=dict(prepared_only=True,no_root_seed_factor_series_or_candidate_values_executed=True,
  same_d1_p32003_only=True,prime_certificate='all remainders P%div for2<=div<=178, sqrt bound178²<=P<179²',
  source_fixed_quartic=[2668,18535,17561,16979,1],seed_monic_f0=[5457,11750,555,10369,14038,10960,1],
  seed_complete_verification='exactf0/quartic; full4linearproduct; quadraticLegendre; six distinctKroots/f0residuals',
  quadratic_field='q monic irreducible degree2; Karatsuba cross sums reduced BEFORE multiply',
  root_Newton_precision=10,root_coefficients=54,derivative_inverses=6,root_identity_coefficients=60,
  subsets=64,unordered_proper_subsets=31,leading_allocations=2,full_candidate_products=62,full_product_cells_each=133,
  failure_witness_peer_recipe='reconstruct root/subset products; check one reported nonzero full-product coefficient per62cases',
  product_upper_components=dict(seed_quartic_Horner=128012,root_Newton=80640,root_final_checks=21000,
   subset_products=378000,full_candidates=496000,leading_allocations_and_other_checks=50000),
  source_base_product_upper=1300000,hard_base_products=2000000,hard_field_temporary_bits=30,
  existing_source_input_bits=40,source_formal_content_integer_bits=80,global384_unchanged=True,
  state_numeric_guard='6000+recursive numbers/tuple/dict/control cells <=24000 before periodic commits',
  subset_max_coefficient_capacity=8960,hard_owned=24000,
  all_live_bound='all64subsets8960 +root/input/seed/updates<2500 +candidate/doubleproduct<1500 +local/parser/serialization/control reserve<=10000 <24000',
  whole_RAW=4194304,stream_one_prefix=True,no_all62_candidates_resident=True,
  prefix_includes='pending scalar,field/basepoly/series/bivar/workspaces,allrootjets,allsubsets,currentcandidate,cursors,sourcepins',
  soft_seconds=45,hard_seconds=60,RAM_bytes=2000000000,
  positive_scope='PRIMARY modp absolute irreducibility and char0geometricH, requires independent finite identity and theorem peer',
  negative_scope='one actualfactor of fixed Fp2 bivariate image only; nochar0 reducibility claim/no newprime or sample',
  no_Fitting_fullrank_or_exhaustion_claim=True)
 (out/'whole-source-cursor-resource-and-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 rep=dict(id=out.name,all_passed=True,frozen=True,prepared_only=True,no_value_worker_run=True,
  worker_sha256=sha(worker),source_input=src+'fixed-d1-input-and-formal-content.json',runtime_input_sha256=runtime,input_sha256=pins,
  evidence_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()})
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n')
 raw=sum(x.stat().st_size for x in out.iterdir())+sum((ROOT/s).stat().st_size for s in pins)
 assert raw+2097152+262144<4194304
 print(json.dumps(dict(SHA256=sha(p),worker_SHA=sha(worker),prepared_input_code_RAW=raw,RAW_prior=raw+2097152+262144,values_executed=False),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
