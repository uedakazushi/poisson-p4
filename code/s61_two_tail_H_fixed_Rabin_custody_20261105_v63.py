#!/usr/bin/env python3
"""Prepare immutable one-specialization H degree6 Rabin custody; help only, no field arithmetic."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p,h=None):
  s=sha(p)
  if h:assert s==h
  pins[str(p.relative_to(ROOT))]=s;return p
 factor=pin(ROOT/'reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/factor-sparse.txt','241488fafa9b45e5353f0b707f75902bc5bc32e85ee093b9f2e0a5f07c7d571c')
 peer=pin(ROOT/'reports/S61_TWO_TAIL_CRT_ROOT_DENSE_REVIEW_20261105_v63/report.json','6d7f852148000830f26d4f9ed709a7e7a5bc205e8c2678e4613648ea82b23bc7')
 assert json.loads(peer.read_text())['exact_char0_gcd_proved_up_to_QQ_scalar']
 lines=factor.read_text().splitlines();assert int(lines[0])==94
 rows=[tuple(map(int,line.split())) for line in lines[1:]];assert len(rows)==94
 # Literal source coefficient selection only; NO e0,d1 evaluation or field operation.
 a6=sorted((e,d,v) for a,e,d,v in rows if a==6)
 assert a6==[(0,3,1406250000),(0,4,-439453125),(1,4,1160156250)]
 a0e0=sorted((d,v) for a,e,d,v in rows if a==0 and e==0)
 assert a0e0==[(0,-1025280),(1,-12816)]
 worker=pin(ROOT/'code/s61_two_tail_H_fixed_Rabin_20261105_v63.py')
 note=pin(ROOT/'research/notes/S61_TWO_TAIL_H_FIXED_RABIN_SOURCE_20261105_v63.md');pin(Path(__file__))
 r=subprocess.run([sys.executable,str(worker),'--help'],capture_output=True,text=True,timeout=10)
 (out/'worker-help.txt').write_text(r.stdout);(out/'worker-help.stderr').write_text(r.stderr);assert r.returncode==0
 bounds=dict(fixed_p=32003,fixed_e=0,fixed_d=1,degree=6,primitivity_literal_a6=a6,primitivity_literal_a0_e0=a0e0,
  primitivity_reason='only possible coefficient content primes d or 80-25d+66ed; d and L both fail on saved a0 e0 coefficient',
  fixed_leading_integer_expression='17578125*55; nonzero mod32003 required again at execution',
  source_modular_product_upper=13900,hard_modular_products=14000,
  product_derivation='6*30*66=11880; all normalization and inverse checks plus at most14 degree6 Euclid divisions add fewer than2020',
  max_new_field_temporary_bits=30,max_existing_input_bits=40,hard_owned=8192,hard_RAW=1048576,
  owned_proof='94*4 source numbers + <=7*6 Frobenius + two <=7-step Euclid chains with <=28 numbers each + <=4 transient length11 vectors + keys/control2048 + full failure reserve <8192',
  RAW_proof='all read sources/code/note/custody plus full state evidence/failure JSON <1MiB; phase includes current cursor and pending scalar',
  soft_seconds=45,hard_seconds=60,address_space_bytes=2000000000,
  no_actual_specialization_Rabin_gcd_or_product_executed=True,
  positive_logic='primitive in QQ[e,d][a], preserved degree6, irreducible fixed Fp image implies QQ[a,e,d] irreducible',
  no_extra_prime_or_specialization_on_negative=True)
 b=out/'source-and-whole-guard-proof.json';b.write_text(json.dumps(bounds,ensure_ascii=False,indent=2)+'\n')
 rep=dict(id=out.name,all_passed=True,frozen=True,preparation_only=True,no_target_values_or_Rabin_executed=True,
  worker_sha256=sha(worker),input_sha256=pins,evidence_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()})
 p=out/'report.json';p.write_text(json.dumps(rep,ensure_ascii=False,indent=2)+'\n');assert sum(x.stat().st_size for x in out.iterdir())+sum((ROOT/s).stat().st_size for s in pins)<1048576
 print(sha(p))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
