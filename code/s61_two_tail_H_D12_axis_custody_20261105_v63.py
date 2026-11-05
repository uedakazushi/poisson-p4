#!/usr/bin/env python3
"""SOURCE-only binding and help for fixed D12 determinant in a four-dimensional ring."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 out.mkdir(parents=True,exist_ok=False);pins={};runtime={}
 def pin(s,h,run=False):
  p=ROOT/s;assert sha(p)==h;pins[s]=h
  if run:runtime[s]=h
  return p
 factor=pin('reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/factor-sparse.txt','241488fafa9b45e5353f0b707f75902bc5bc32e85ee093b9f2e0a5f07c7d571c',True)
 peer=pin('reports/S61_NONSQUARE_PREPOWER_M12_AXIS_ONE_COEFFICIENT_20261105_v63_REVIEW/report.json','c5ebf88f59f93ab7a4a4158ed10f50c89618fa9d1ffb9f8d3e7fab5ef38d5709',True)
 r=json.loads(peer.read_text());assert r['all_passed'] and r['full_QQ_source_scope']
 mname='reports/S61_NONSQUARE_PREPOWER_M12_AXIS_ONE_COEFFICIENT_20261105_v63_REVIEW/own-integer-C12-normalization.json'
 mp=pin(mname,'1a834104707e4d7c85d964195b6cb768c7d249b497870888b2217932a9a4531a',True);m=json.loads(mp.read_text());M=m['matrix']
 assert m['all144_primary_integer_inputs_match'] and m['all23_actual_source_adapter_values_match']
 native=pin('reports/S61_NONSQUARE_PREPOWER_M12_AXIS_ENTRY_20261105_v63/native-input.txt','266b720d8cd7d2aa05b16fd6316e446a4d3a389bd16f2dabfefe66a314855328')
 with native.open() as f:
  assert f.readline().split()[0]=='M12AXIS001';f.readline() # old ALL-minor height bound, NOT the LCD
  for i in range(12):
   for j in range(12):
    x=list(map(int,f.readline().split()));assert x[:2]==[i,j] and len(x)==3+2*x[2]
    assert [(a,int(v)) for a,v in M[i][j]]==list(zip(x[3::2],x[4::2]))
  assert not f.read().strip()
 pin('reports/S61_NONSQUARE_PREPOWER_M12_AXIS_ENTRY_20261105_v63/source-normalization-and-guards.json','9b9c626b00f9422d8e31cc4e6da39745871b2c8f31e58fbbc1c367f1a820a978',True)
 pin('reports/S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63/report.json','49990187b3b0b8edaf8939e4d13ac297f4f2ecba7f8b9f5326639089fa043ce2',True)
 pin('reports/S61_TWO_TAIL_COMMON_FACTOR_SOURCE_20261105_v63/report.json','9a49a8f61d41faf0a2d2d68bfce6817d8c5d35b6d24aa570c33d54d703d844d9',True)
 rows=[tuple(map(int,s.split())) for s in factor.read_text().splitlines()[1:]]
 axis=sorted((a,v) for a,e,d,v in rows if e==d==0)
 h=[-1025280,35712000,-1148800000,26260000000,-96000000000];assert axis==list(enumerate(h))
 coeff=[int(v) for row in M for p in row for k,v in p];degree=max(k for row in M for p in row for k,v in p)
 assert len(coeff)==726 and degree==8 and max(abs(v).bit_length() for v in coeff)==110
 assert int(m['LCD_product'])%32003==30563 and h[-1]%32003==7157 and h[0]%32003==30819
 worker=ROOT/'code/s61_two_tail_H_D12_axis_quotient_20261105_v63.py';pins[str(worker.relative_to(ROOT))]=sha(worker)
 note=ROOT/'research/notes/S61_TWO_TAIL_H_D12_AXIS_QUOTIENT_SOURCE_20261105_v63.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 x=subprocess.run([sys.executable,str(worker),'--help'],capture_output=True,text=True,timeout=10)
 (out/'worker-help.txt').write_text(x.stdout);(out/'worker-help.stderr').write_text(x.stderr);assert x.returncode==0
 proof=dict(source_only=True,all144_original_integer_cells_literal_bound=True,no_67_or_214_or_572_replay=True,
  actual_source_terms=726,max_input_degree=8,input_coefficient_bits=110,LCD_bits=519,LCD_mod_p=30563,H_axis_lc_mod_p=7157,H_axis_const_mod_p=30819,
  exact_normalization='D12_axis=a5*detC12/LCD',axis_legal_denominators='a and B only, B_axis=-40',
  ring='Fp32003[a]/monic(H_axis)',finite_free_dimension=4,not_assumed_field_or_reduced=True,
  DP_ring_products=24576,ring_product_scalar_upper=28,total_scalar_product_upper=700000,hard_scalar_products=750000,
  DP_coefficients=16384,matrix_coefficients=576,hard_owned=28000,
  owned_bound='10000 fixed source/keys/control/failure reserve +16384DP+576matrix+36powers <=26996',
  field_temporary_bits=30,existing_ROM_coefficients=110,existing_ROM_LCD_bits=519,
  RAW_upper=1048576,RAW_proof='external union <200k; stream prefix/evidence at most28k numeric slots plus fixed64k textual reserve <400k; no DP copy',
  source_contradiction_if_positive='H coprime U and primitive; clear aB, Gauss at existing prime, legal axis retains nonzero Hdegree4; a and LCD are quotient units',
  zero_remainder_inconclusive=True,source_only_no_quotient_or_determinant_values=True,
  soft_seconds=45,hard_seconds=60,address_space=2000000000)
 (out/'source-frame-and-whole-budget.json').write_text(json.dumps(proof,indent=2)+'\n')
 rep=dict(id=out.name,all_passed=True,frozen=True,preparation_only=True,no_determinant_or_quotient_values_executed=True,
  worker_sha256=sha(worker),matrix_input=mname,H_axis_integer=h,runtime_input_sha256=runtime,input_sha256=pins,
  evidence_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()})
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');assert sum(x.stat().st_size for x in out.iterdir())+sum((ROOT/s).stat().st_size for s in pins)<200000
 print(sha(p))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
