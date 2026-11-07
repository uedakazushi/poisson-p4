#!/usr/bin/env python3
"""Prepare ONLY the one directQQ H94 factor backend/conditional identity pipeline."""
from pathlib import Path
import argparse,ast,hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 out.mkdir(parents=True,exist_ok=False);pins={};runtime={}
 def pin(s,h=None,run=False):
  p=ROOT/s;v=sha(p)
  if h:assert v==h
  pins[s]=v
  if run:runtime[s]=v
  return p
 H=pin('reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/factor-sparse.txt','241488fafa9b45e5353f0b707f75902bc5bc32e85ee093b9f2e0a5f07c7d571c',True)
 pin('reports/S61_TWO_TAIL_CRT_ROOT_DENSE_REVIEW_20261105_v63/report.json','6d7f852148000830f26d4f9ed709a7e7a5bc205e8c2678e4613648ea82b23bc7',True)
 pin('reports/S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63/report.json','49990187b3b0b8edaf8939e4d13ac297f4f2ecba7f8b9f5326639089fa043ce2',True)
 neg=pin('reports/S61_TWO_TAIL_H_D1_MODFACTOR_ROOT_REVIEW_20261105_v63/report.json','a41a96ded2bb3fe8d0b17dd98a38c98253b9510ee8618343f4bb5e07317a7f5c',True)
 nr=json.loads(neg.read_text());assert nr['all_passed'] and nr['all133_coefficients_match_original_H_d1_div_alpha']
 pin('research/notes/S61_TWO_TAIL_H_D1_MODFACTOR_ROOT_REVIEW_20261105_v63.md')
 lines=H.read_text().splitlines();assert int(lines[0])==94;rows=[tuple(map(int,s.split())) for s in lines[1:]];assert len(rows)==94
 assert max(abs(v).bit_length() for a,e,d,v in rows)<=40
 box=[max(r[j] for r in rows) for j in range(3)];assert box==[6,9,4]
 kron=max(50*a+5*e+d for a,e,d,v in rows);L1=sum(abs(v) for a,e,d,v in rows)
 assert kron==309 and L1.bit_length()<=47
 cpp=pin('code/s61_two_tail_H_exact_QQ_factor_20261105_v63.cpp');driver=pin('code/s61_two_tail_H_exact_QQ_entry_20261105_v63.py')
 ast.parse(driver.read_text());me=pin(str(Path(__file__).relative_to(ROOT)))
 note=pin('research/notes/S61_TWO_TAIL_H_EXACT_QQ_FACTOR_SOURCE_20261105_v63.md')
 command=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" $(pkg-config --cflags --libs flint) -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu','HQQ-build',str(cpp),str(out/'native-worker')]
 with (out/'compile.log').open('w') as f:r=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,timeout=30)
 if r.returncode:
  for p in [cpp,driver,me,note]:(out/('FAILED_PREPARATION_'+p.name)).write_bytes(p.read_bytes())
  (out/'failure.json').write_text(json.dumps(dict(stage='compile_only',no_factor_values_executed=True,returncode=r.returncode,input_sha256=pins),indent=2)+'\n')
  raise RuntimeError('compile-only failure preserved')
 for name,cmd in [('native-help.txt',[str(out/'native-worker'),'--help']),('driver-help.txt',[sys.executable,str(driver),'--help'])]:
  r=subprocess.run(cmd,capture_output=True,text=True,timeout=10);(out/name).write_text(r.stdout+r.stderr);assert r.returncode==0
 proof=dict(prepared_only=True,no_factor_call_or_values_executed=True,direct_three_variable_factor=True,
  API='FLINT3.1 fmpz_mpoly_factor, integer STATUS; primitive_part normalization also opaque library phase',
  source_H_terms=94,input_bits=40,box=box,total_degree=11,Kronecker_BOUND_ONLY=[50,5,1],Kronecker_degree309=kron,L1_bits=L1.bit_length(),
  Mignotte_primitive_output_bits_upper357=True,no_univariate_factor_decoding=True,canonical_emission_hard_bits=357,
  native_internal_bits_and_support_OPAQUE=True,native_AS_limit=1900000000,whole_RAM_limit=2000000000,
  internal_interrupt_frontier='source/code/phase/logs/any bounded canonical outputs preserved; NOT restartable internal API state',
  verification_STAGE_condition='actual partial_max_bits + next_factor_max_bits +18 <=384 BEFORE multiplication',
  verified_identity_hard_bits=384,conditional_frontier_on_failure_no_capraise=True,all_coefficients_and350box_zeros_checked=True,
  factor_count_cap=11,canonical_factor_terms_cap_each=350,aggregate_terms_cap=3850,hard_owned_native_and_Python_each=50000,
  source_pair_upper=11*350*350,hard_verification_pairs=1500000,
  RAW_cap=4194304,prefix_or_success_bundle_upper=2097152,report_log_reserve=262144,
  shared_soft_seconds=45,driver_hard_timeout58=True,hard_seconds=60,
  no_new_prime_d_value_beta_source_modular_factor_farm=True,
  exact_product_needed_before_factorization_claim=True,API_contract_not_multiplyback_alone_for_QQ_irreducibility=True,
  no_geometric_irreducibility_Fitting_fullrank_or_exhaustion_claim=True)
 (out/'source-height-lifetimes-opaque-boundary-and-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 pub=dict(accepted_independent_fixed_modular_factor_peer_only=True,report_sha256=sha(neg),Japanese_note_sha256=sha(ROOT/'research/notes/S61_TWO_TAIL_H_D1_MODFACTOR_ROOT_REVIEW_20261105_v63.md'),no_char0_factor_inferred=True,no_arithmetic_replay=True)
 (out/'negative-modular-peer-publication.json').write_text(json.dumps(pub,indent=2)+'\n')
 rep=dict(id=out.name,all_passed=True,frozen=True,prepared_only=True,no_factor_values_run=True,
  fixed_H_input=str(H.relative_to(ROOT)),native_executable=str((out/'native-worker').relative_to(ROOT)),native_executable_sha256=sha(out/'native-worker'),
  driver_sha256=sha(driver),runtime_input_sha256=runtime,input_sha256=pins,evidence_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()})
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n')
 raw=sum(x.stat().st_size for x in out.iterdir())+sum((ROOT/s).stat().st_size for s in pins)
 assert raw+2097152+262144<4194304
 print(json.dumps(dict(SHA256=sha(p),native_SHA=rep['native_executable_sha256'],driver_SHA=sha(driver),prepared_RAW=raw,RAW_prior=raw+2097152+262144,values_run=False),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
