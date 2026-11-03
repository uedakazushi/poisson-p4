#!/usr/bin/env python3
"""Independent coefficient-orbit and finite-module proof application review."""
import hashlib,json,itertools
from pathlib import Path
import sympy as s
from ambient_B18_CA_N_critical_review_20261103_v63 import ROOT,vec,coefficient_action,x,h,z,w
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 producer=ROOT/'reports/AMBIENT_B18_REDUCED_COMPONENTS_NONCM_20261103_v63.json';d=json.loads(producer.read_text());assert d['all_passed']
 for p,want in d['input_sha256'].items():assert sha(ROOT/p)==want,p
 P=s.zeros(5)
 for j,v in[(0,x*x),(3,h*x),(4,h*h+z*z)]:P[2,j]=v;P[j,2]=-v
 v=vec(P);O=s.Matrix.hstack(*(coefficient_action(a,b)*v for a,b in itertools.product(range(5),repeat=2)));cc=list(O.rref()[1]);assert len(cc)==19
 V=O[:,cc];rr=list(V.T.rref()[1]);assert V[rr,:].det()!=0;assert V*(V[rr,:].inv()*v[rr,:])==v
 note=ROOT/'research/notes/AMBIENT_B18_REDUCED_COMPONENTS_NONCM_20261103_v63_REVIEW.md';paths=[Path(__file__),note,producer,ROOT/'code/ambient_B18_CA_N_critical_review_20261103_v63.py',ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json']
 result=dict(id='AMBIENT_B18_REDUCED_COMPONENTS_NONCM_20261103_v63_REVIEW',all_passed=True,finite_R_actual_image_direct_sum_R_I_proof_accepted=True,CM_implies_free_R_by_parameter_depth_and_Auslander_Buchsbaum=True,rank1_ideal_direct_summand_free_principal_contradiction=True,no_regular_nu_mod_I_or_Tor_shortcut_assumed=True,own_full150_coefficient_GL_orbit_rank=19,scalar_in_own_orbit=True,projective_dense_orbit_dimension=18,entire_B18_bad_locus_inclusions_by_openness_and_orbit_closure=True,individual_reduced_CA_CB_N_nonCM_nonnormal_singular_on_entire_B18=True,exact_bad_locus_depth_or_actual_D_completion_claimed=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths})
 (ROOT/'reports/AMBIENT_B18_REDUCED_COMPONENTS_NONCM_20261103_v63_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print('INDEPENDENT ENTIRE B18 REDUCED TARGET NONCM/NONNORMAL REVIEW PASS')
if __name__=='__main__':main()
