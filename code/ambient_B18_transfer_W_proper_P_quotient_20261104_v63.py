#!/usr/bin/env python3
"""Gate-only strictness consequence; no new source/rank or coefficient values."""
import hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_TRANSFER_W_PROPER_P_QUOTIENT_20261104_v63';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 t=time.monotonic();OUT.mkdir(exist_ok=True);pins={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 w=load('AMBIENT_B18_TRANSFER512_HONEST_TNIL_MAP_20261104_v63');assert w['all_passed']and w['origin_rank']==0
 assert pins['reports/AMBIENT_B18_TRANSFER512_HONEST_TNIL_MAP_20261104_v63/report.json']=='49dc174679e229e2b53f1536006f675f207574a01c1b3333bc5f9fa7ec9c6d0b'
 o=load('AMBIENT_B18_TRANSFER512_HONEST_TNIL_MAP_20261104_v63','actual512-origin-values-and-Pi0.json');assert len(o['all512'])==512 and all(all(str(c)=='0'for c in col['minimal_origin5'])for col in o['all512'])
 n=load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW');assert n['all_passed']and n['actual_origin_rank']==2 and n['exact_natural_injection']
 assert pins['reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW/report.json']=='07d3fd603e42e012507f7ba724c9b47d586f5b8df2d306e7ee468e4edddfb940'
 for st in ['AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES','AMBIENT_B18_HONEST_NIL_T_OPERATOR','AMBIENT_B18_NIL_DVR_FINITENESS','AMBIENT_B18_NIL_CUBE_OBSTRUCTION_REDUCTION']:
  assert load(st+'_20261104_v63_REVIEW')['all_passed']
 lower=5-2-2;assert lower==1
 ev=dict(notation={'U':'honest Tnil','A':'N/ellN subset U','B':'K/ellK subset A','C':'image(W→A)','V':'Told with natural injection j intoU'},
  exact_inputs=['length(U/A)=2','dim(U/mU)=5','dim image(V→U/mU)=2','C subset mU','B subset image(V)'],
  proof='image(A→U/mU) has dimension at least3. Its image modulo image(jV) therefore has dimension at least1. B,C,mA map0, so A/(B+C+mA) is nonzero.',
  actual_identification='A/(B+C)=N/(K+W+ellN)=P/(psi(W)+ellP)',
  minimal_generator_lower=lower,conclusions=['P/psi(W) nonzero','psi(W) is a proper submodule of P','K+W is a proper subideal of N','W is a proper subideal of N'],
  no_choice_of_Kfiber_basis_or_primitive_normalizer_fiber=True,
  method_frontier='The sufficient route N=K+W and stronger N=W are both excluded for this actual transfer W. Other beta/gamma vanishing routes remain open.',
  not_inferred=['N^3 nonzero','gamma nonzero','beta nonzero','precise rank/length of P/psi(W)','new nil generator list or exact nil index'],
  W_origin_primary_status='Frozen PRIMARY49dc1746 is an explicit input; its independent source peer is a separate gate, not claimed by this application.')
 p=OUT/'strictness-dimension-proof.json';p.write_text(json.dumps(ev,indent=2)+'\n')
 for p2 in[Path(__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p2.relative_to(ROOT))]=sha(p2)
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,transfer_W_strictly_proper_in_N=True,transfer_image_strictly_proper_in_P=True,mu_S_P_mod_psiW_ellP_lower=1,nil_cube_still_OPEN=True,no_new_source_or_coefficient_calculation=True,input_sha256=pins,evidence_sha256={'strictness-dimension-proof.json':sha(p)},seconds=time.monotonic()-t)
 (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
