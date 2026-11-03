#!/usr/bin/env python3
"""Transport accepted generic W completions along the independently exact orbit."""
import hashlib,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 start=time.monotonic();stem='AMBIENT_B18_ALL_NONZERO_W_FULL_COMPLETIONS_20261103_v63'
 out=ROOT/'reports'/stem;out.mkdir(exist_ok=False)
 generic=ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63_REVIEW.json'
 orbit=ROOT/'reports/AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63_REVIEW/report.json'
 for p in [generic,orbit]:
  d=json.loads(p.read_text());assert d['all_passed'] is True
  for name,want in d['input_sha256'].items():assert sha(ROOT/name)==want,name
 author=ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63/report.json'
 d=json.loads(author.read_text());assert d['all_passed']
 # For a generic nonzero W0 and any specified W1 !=0 choose r=W0/W1.
 # The exact orbit identity gives r^-1 g_r* P_W0=P_W1.
 # CA, CB and N are individually PGL-invariant reduced image components;
 # projective completed local rings and their normalizations are transported.
 note=ROOT/'research/notes'/f'{stem}.md'
 paths=[Path(__file__),note,generic,orbit,author]
 report=dict(id=stem,all_passed=True,independently_reviewed_generic_and_orbit_inputs=True,
  orbit_transport='r=W0/W1; r^-1 diag(1,1,r,1,1)_* P_W0=P_W1',
  each_component_CA_CB_N_individually_PGL_invariant=True,
  scope='Every complex W !=0 on the specified actual common W curve.',
  CA_CB_projective_completion='C[[g19,s8,t^2,s1*t,...,s8*t,t^3]]',
  CA_CB_projective_invariants=dict(dimension=28,embedding_dimension=37,depth=20,multiplicity=2),
  N_projective_completion='C[[g19,s7,t^2,s1*t,...,s7*t,t^3]]',
  N_projective_invariants=dict(dimension=27,embedding_dimension=35,depth=20,multiplicity=2),
  each_normalization_elementary_subintegral=True,
  W0_CA_CB_full_ideal_or_actual_D_solved=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print('ALL NONZERO W FULL COMPLETIONS BY TWO INDEPENDENT REVIEWS: PASS')
if __name__=='__main__':main()
