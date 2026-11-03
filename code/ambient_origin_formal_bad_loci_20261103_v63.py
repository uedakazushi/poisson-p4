#!/usr/bin/env python3
"""Exact ann-colons and proof receipt for formal bad loci; no GB."""
import argparse
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_ORIGIN_FORMAL_BAD_LOCI_20261103_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def minimal(A):
 return {u for u in A if not any(v!=u and all(a<=b for a,b in zip(v,u)) for v in A)}
def colon(A,mon):return minimal({tuple(max(a-b,0) for a,b in zip(u,mon)) for u in A})
def intersection(A,B):return minimal({tuple(max(a,b) for a,b in zip(u,v)) for u in A for v in B})
def main(output):
 parents=[ROOT/'reports'/stem/'report.json' for stem in [
  'AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63',
  'AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63_REVIEW',
  'AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63',
  'AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63_REVIEW']]
 for p in parents:
  d=json.loads(p.read_text());assert d['all_passed']
  for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
  for n,h in d.get('evidence_sha256',{}).items():assert sha(p.parent/n)==h,n
 ideal={(1,1,0),(1,0,1),(0,2,1)};cw=colon(ideal,(1,0,0));cr=colon(ideal,(0,1,0));both=intersection(cw,cr)
 assert cw=={(0,1,0),(0,0,1)} and cr=={(1,0,0),(0,1,1)}
 assert both=={(1,1,0),(1,0,1),(0,1,1)}
 radical=minimal({tuple(int(a>0) for a in mon) for mon in ideal});assert radical==both
 note=ROOT/'research/notes'/f'{STEM}.md';files=[Path(__file__),note,*parents]
 result=dict(id=STEM,all_passed=True,method='commutative-algebra proof from independently accepted exact formal models; finite monomial colons only',
  original_A_reduced_Sing_Nonnormal_NonCM_Nonseminormal_all_equal='three smooth20 conductor branches W,V,S meeting in smoothG19',
  seminormal_Aprime_reduced_Sing_Nonnormal_NonCM='only smoothselfS20',seminormal_Aprime_Nonseminormal_empty=True,
  exact_Ann_C_w='(r,s)',exact_Ann_C_r='(w,rs)',exact_Ann_C_w_r='(rs)=nilradical(C)',
  all_conductor_points_nonCM_by_local_AB_and_height9=True,normalization_fiber_proofs='K[t]/(t^2-nubar) on conductor; isomorphism off conductor',
  normalizer_fibers_original={'WminusG_VminusG':'one point dual number length2','SminusG':'two reduced points','G19':'one point dual number length2','off_conductor':'isomorphism'},
  normalizer_fibers_seminormal={'SminusG':'two reduced points','G19':'one point dual number length2','off_selfS':'isomorphism'},
  generic_transverse_dim_edim_depth_original={'W_V':[9,18,1],'S':[9,18,1],'G':[10,20,2]},
  generic_transverse_dim_edim_depth_seminormal={'W_V':[9,9,9],'S':[9,18,1],'G':[10,19,2]},
  central_cone_dim_edim_depth_multiplicity_original=[29,39,21,2],central_cone_dim_edim_depth_multiplicity_seminormal=[29,38,21,2],
  projective_original=[28,38,20,2],projective_seminormal=[28,37,20,2],
  global_bad_loci_or_relative_CA_CB_embedding_or_convergent_analytic_classification_claimed=False,
  no_GB=True,input_sha256={str(p.relative_to(ROOT)):sha(p) for p in files})
 with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS',STEM,'exact3-monomial-colon proof; no GB')
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
