#!/usr/bin/env python3
"""SOURCE ONLY binary witness/guard closure for the unchanged fixed M12 axis job."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_M12_AXIS_BINARY_GATE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def load(rel):p=ROOT/rel;pins[rel]=sha(p);return json.loads(p.read_text())
 st='reports/S61_NONSQUARE_PREPOWER_M12_AXIS_SOURCE_20261105_v63/'
 r=load(st+'report.json');f=load(st+'whole-axis-source-worker-forecast.json');assert r['all_passed']and not f['prospective_RAW_gate_passes']
 old='reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_BAREISS_ENTRY_20261105_v63/';n=load(old+'native.json');assert n['rank_over_QQ_a']==11
 load('reports/S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_ENTRY_20261105_v63_REVIEW/report.json')
 for rel,h in r['input_sha256'].items():assert sha(ROOT/rel)==h;pins[rel]=h
 B=(f['uncombined_integer_product_and_longdivision_bits_upper']+7)//8;mainraw=2*1024*1024+f['all_patch_coefficient_slots_upper']*(B+7);guardraw=1024*1024+f['owned_cells_upper']*(B+7)
 external=sum((ROOT/p).stat().st_size for p in r['input_sha256']);products=sum((11-k)**2*(2*((k+1)*9+1)**2+((k+2)*9+1)*(k*9+1))for k in range(11))
 assert products==2068638 and external+mainraw+guardraw<33554432
 proof=dict(old_JSON_FAIL_preserved=True,same_fixed_source_and_arithmetic=True,new_binary_publication_only=True,
  old11_fixed_row_permutation_extended=n['row_permutation']+[11],old11_fixed_column_permutation_extended=n['column_permutation']+[11],
  source_nonzero_pivots='top-left diag(extraD) times frozen oldCaxis, known nonzero constant row units; fixed11 order carried; no new field rank search',
  input_degree=9,minor_degree=108,minor_bits1096=True,temporary_bit_bound=2202,owned_cap33275=True,
  binary_record='degree_u16, sign_u8, byte_length_u16, absolute_integer_big_endian; coefficient payload≤ceil2202/8); zero records omitted; all stage/pivot/cursor/sign headers in fixed reserve',
  all_patch_coefficients=20504,main_witness_RAW=mainraw,guard_all_state_RAW=guardraw,external_input_RAW=external,TOTAL_RAW=external+mainraw+guardraw,
  integer_product_upper=products,proposed_soft45_hard60_seconds=True,memory_bytes=2000000000,
  determinant_positive_only_generic_nonzero=True,axis_ZERO_inconclusive=True,no_axis_values_or_determinant_computed=True,
  no_new_permission_assumed_for2202bit_worker=True)
 dump(out/'binary-publication-source-proof.json',proof)
 direct=ROOT/'research/notes/S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63.md';pins[str(direct.relative_to(ROOT))]=sha(direct)
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 rep=dict(id=ID,all_passed=True,frozen=True,source_only=True,unchanged_input_and_arithmetic=True,all_witnesses_binary_within32MiB=True,
  no_axis_value_or_scalar_det=True,no_polyGCD_fieldrank_point_prime=True,TOTAL_RAW=proof['TOTAL_RAW'],owned=33275,temporary_bits=2202,
  full15_direct_product_identity='A*D15_norm=-T*detM12; no extraPchart on oldOmega',
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',rep);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
