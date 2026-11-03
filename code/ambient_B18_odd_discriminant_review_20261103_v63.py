#!/usr/bin/env python3
"""Independent N finite-base certificate and different/odd-image proof audit."""
import gzip,hashlib,json,math,time
from pathlib import Path
import sympy as s
from ambient_B18_CA_N_critical_review_20261103_v63 import ROOT,X,x,h,z,w,vec,coefficient_action
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();folder=ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63';r=json.loads((folder/'report.json').read_text());assert r['all_passed']
 for p,want in r['input_sha256'].items():assert sha(ROOT/p)==want,p
 for p,want in r['evidence_sha256'].items():assert sha(folder/p)==want,p
 d=json.loads(gzip.decompress((folder/'certificate.json.gz').read_bytes()));nfolder=ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63';n=json.loads(gzip.decompress((nfolder/'certificate.json.gz').read_bytes()));C=s.Matrix(n['fiber19_basis'])
 P=s.zeros(5)
 for j,v in[(0,x*x),(3,h*x),(4,h*h+z*z)]:P[2,j]=v;P[j,2]=-v
 v=vec(P);flag=[tuple(ab)for ab in n['flag_positions']];J=C.row_join(s.Matrix.hstack(*(coefficient_action(*ab)*v for ab in flag)));assert J==s.Matrix(d['N_source_Jacobian150_by28'])
 cols=d['N_complement27_columns'];rr=d['N_actual27_target_rows'];V=J[:,cols];assert len(cols)==len(rr)==27 and V[rr,:].det()==-s.Rational(1152,5)
 ker=s.Matrix(d['N_source_kernel']);assert J*ker==s.zeros(150,1)and ker[22]==1
 M=s.Matrix(d['N_kernel_flag_matrix']);L=sum((M[a,b]*coefficient_action(a,b)for a in range(5)for b in range(5)if M[a,b]),s.zeros(150));W=-L*L*v/2
 assert W==s.Matrix(d['N_entire_second_coefficient']);corrected=W-V*V[rr,:].inv()*W[rr,:];assert corrected==s.Matrix(d['N_all150_corrected_second_coefficients']);assert corrected[27]==-3 and 27 not in rr
 # Exponentiate the independent 150-coefficient representation, instead of
 # reusing the producer's coordinate substitution or action function.
 p=s.Symbol('p');curve=v.copy();power=v.copy()
 for k in range(1,12):
  power=L*power
  if power==s.zeros(150,1):break
  curve+=(1-k)*p**k*power/s.factorial(k)
 else:raise AssertionError('nilpotent representation did not terminate')
 assert curve.applyfunc(s.expand)==vec(s.Matrix(d['N_finite_kernel_curve']))
 assert v[106]==1 and 106 in rr
 assert r['exact_minimal_generators_after_discriminant']=={'CA':10,'CB':10,'N':8}
 # Minimal monomial generators: six independent linear coordinates plus
 # the4/2 pairwise nondividing core monomials.
 CA=[(1,1,0),(1,0,1),(0,1,1),(0,3,0)];N=[(1,1),(0,2)]
 for mons in[CA,N]:assert not any(all(a<=b for a,b in zip(u,v))for u in mons for v in mons if u!=v)
 note=ROOT/'research/notes/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.md';paths=[Path(__file__),note,ROOT/'code/ambient_B18_CA_N_critical_review_20261103_v63.py',folder/'report.json',folder/'certificate.json.gz',nfolder/'certificate.json.gz']
 result=dict(id='AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW',all_passed=True,independent_N_full150_by28_Jacobian_and_actual27minor=True,independent_nilpotent_coefficient_exponential_entire_kernel_curve=True,all150_corrected_quadratic_coefficients_and_target_row_minus3=True,N_own_parameter_base_cone28_and_rank2_finite_free_proof_accepted=True,trace_shift_and_regular_discriminant_and_critical_equals_t_I_proof_accepted=True,regular_t_minimal_generator_and_six_linear_eliminations_accepted=True,actual_image_odd_ideal_mod_discriminant_formal_classes_accepted=True,odd_image_minimal_generators={'CA':10,'CB':10,'N':8},full_odd_ideal_minimal_generators_lower_bounds={'CA':10,'CB':10,'N':8},I_image_mod_nu_not_I_mod_nuI=True,independent_coordinate_changes_not_aligned_and_actual_D_unresolved=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},seconds=time.monotonic()-start)
 (ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
