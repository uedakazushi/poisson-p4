#!/usr/bin/env python3
"""SOURCE ONLY order-specific polynomial L1 bounds for all exact division states."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_M12_AXIS_ORDERED_L1_GATE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def load(rel):p=ROOT/rel;pins[rel]=sha(p);return json.loads(p.read_text())
 f=load('reports/S61_NONSQUARE_PREPOWER_M12_AXIS_SOURCE_20261105_v63/whole-axis-source-worker-forecast.json')
 b=load('reports/S61_NONSQUARE_PREPOWER_M12_AXIS_BINARY_GATE_20261105_v63/binary-publication-source-proof.json')
 L=list(map(int,f['known_all_row_integer_L1_upper']));order=b['old11_fixed_row_permutation_extended'];assert len(L)==12 and order==list(range(12))
 pref=1;stages=[];maximum=0
 for k in range(11):
  w=max(pref*pref*L[order[k]]*L[order[i]]for i in range(k+1,12));bits=(3*w).bit_length();maximum=max(maximum,bits)
  stages.append(dict(stage=k,previous_row_L1_product=str(pref),pivot_row=order[k],product_poly_L1_upper=str(w),numerator_poly_L1_upper=str(2*w),ALL_partial_longdivision_L1_upper=str(3*w),integer_bits=bits));pref*=L[order[k]]
 assert maximum==1984 and f['all_minor_coefficient_bits_upper']==1096 and maximum<2048
 code=Path(__file__);pins[str(code.relative_to(ROOT))]=sha(code);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,completed_source_bound_bits=maximum,implementation_hard_guard_bits=2048,
  no_axis_value_or_determinant_computed=True,no_guard_cap_raise=True,old_uniform2202_retained_as_weaker=True,
  proof='whole polynomial L1 bounds quotient and denominator; partial convolution/longdivision use subsets and total subtraction L1≤W; remainder≤3W',
  owned_unchanged=33275,binary_RAW_unchanged=b['TOTAL_RAW'],stages=stages,input_sha256=pins,seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if k not in['input_sha256','stages']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
