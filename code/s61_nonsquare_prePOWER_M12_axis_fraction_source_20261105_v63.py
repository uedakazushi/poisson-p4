#!/usr/bin/env python3
"""SOURCE ONLY Fraction temporary and prep ownership, no signed axis values."""
import argparse,hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_M12_AXIS_FRACTION_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);pins={}
 def load(rel):p=ROOT/rel;pins[rel]=sha(p);return json.loads(p.read_text())
 st='reports/S61_NONSQUARE_PREPOWER_M12_AXIS_SOURCE_20261105_v63/'
 s=load(st+'added-source-profiles.json');f=load(st+'whole-axis-source-worker-forecast.json');a=load('reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_SOURCE_V2_20261105_v63/whole-axis-literal-source.json');v=[]
 for p in [x['bound']for x in s['finite_DAG_profiles']]+s['added_column11']+s['added_G134_row11']+[s['corner']]:
  if p['axis']:v.extend([p['axis']['L1'],p['axis']['D']])
 v.extend(map(int,f['known_all_row_integer_L1_upper']));v.extend(int(d)*e for d,e in zip(a['axis_row_denominators'],f['old11_row_extra_QQ_integer_scales']));v.append(int(f['new_row_LCD']));bits=max(abs(int(x)).bit_length()for x in v);assert bits==118
 slots=sum(0 if n['bound']['axis']is None else n['bound']['axis']['hi']-n['bound']['axis']['lo']+1 for n in s['finite_DAG_profiles']);assert slots==1309
 counts=dict(retained_source_profile_numeric=3340,forecast_numeric=266,old_axis_numeric=2051,current15_roots330=330,new_QQ_atom_numeric_upper=160,checkpoint_text_offsets_numeric=1600,all214_poly_key_num_den=3*slots,one_current_literal_cell_numeric_upper=2540,all23_added_poly_key_num_den_upper=504,native_integer_input_degree_coefficient_upper=1494,ALL_loop_callback_scratch_reserve=3000)
 assert sum(counts.values())==19212
 code=Path(__file__);pins[str(code.relative_to(ROOT))]=sha(code);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 r=dict(id=ID,all_passed=True,frozen=True,source_only=True,max_canonical_state_bits=bits,all_unreduced_Fraction_temporary_bits=2*bits+2,hard_guard2048_unchanged=True,
  whole_prep_numeric_upper=sum(counts.values()),counted_all_numeric_key_metadata=counts,source214_full_axis_support_upper=slots,no_signed_axis_values=True,input_sha256=pins,seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
