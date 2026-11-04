#!/usr/bin/env python3
"""Independent universal-polynomial proof and actual-source gate audit."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,resource,time
ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_QUADRATIC_NORMALIZER_NIL_UPPER12_20261104_v63'
STEM=AUTHOR+'_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def add(*ps):
 out={}
 for p in ps:
  for m,c in p.items():
   q=out.get(m,Q(0))+c
   if q:out[m]=q
   else:out.pop(m,None)
 return out
def scale(p,q):return {m:c*q for m,c in p.items() if c*q}
def mul(a,b):
 out={}
 for m,c in a.items():
  for n,d in b.items():out=add(out,{tuple(x+y for x,y in zip(m,n)):c*d})
 return out
def enc(p):return [[list(m),str(c)] for m,c in sorted(p.items())]
def main():
 t=time.monotonic();out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 hashes={}
 def get(p):hashes[p]=sha(p);return read(p)
 p='reports/'+AUTHOR+'/report.json';R=get(p)
 assert hashes[p]=='332701163619c4e3480e77e4240f604bb3200ef9c17a2131daa72d363bc0d77c' and R['all_passed']
 for p,h in R['input_sha256'].items():assert sha(p)==h;hashes[p]=h
 for n,h in R['evidence_sha256'].items():p='reports/'+AUTHOR+'/'+n;assert sha(p)==h;hashes[p]=h
 N=get('reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json')
 E=get('reports/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW/report.json')
 C=get('reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63/report.json')
 CR=get('reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63_REVIEW/report.json')
 T=get('reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json')
 QB=get('reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json')
 for d in (N,E,C,CR,T,QB):assert d['all_passed']
 assert N['actual_free4_48_and_product_conductor_comparison_accepted']
 assert E['core_and_normalized_E4_exact_nilradical_index3_accepted']
 assert C['all_positive_homogeneous_coefficient_maps_A_n_to_M_injective']
 assert C['tail_nilradical_square_zero'] and C['nil_square_nonzero_at_P0_by_Qb_coherent_closure']
 assert C['nilradical_index_at_P0_lower_bound']==3
 assert CR['mixed_positive_products_zero'] and CR['full_completed_nil_cube_equals_degree0_cube']
 assert T['Dalg_polynomial_Laurent_scope_correction_required_and_bound']
 assert T['tail_nil_principal_square_zero_and_cyclic_same_R0_regular_curve']
 assert QB['nilrank4_square_rank2_cubezero_Ass_N_and_embedded_SQ_depth20_e1_accepted']
 # Free commutative polynomial identity in (k,j0,j1,s,nu,tau).
 v=[{tuple(int(i==j) for i in range(6)):Q(1)} for j in range(6)]
 k,j0,j1,s,nu,tau=v;a=add(k,scale(j0,-1));j1tau=mul(j1,tau)
 relation=add(mul(tau,tau),scale(mul(s,tau),-1),scale(nu,-1))
 a_relation=add(mul(a,a),scale(mul(mul(s,j1),a),-1),scale(mul(nu,mul(j1,j1)),-1))
 witness=add(mul(add(a,scale(j1tau,-1)),add(a,j1tau,scale(mul(s,j1),-1))),mul(mul(j1,j1),relation))
 assert a_relation==witness
 # Rewriting k² modulo the source ideal J uses j0,j1 in J, without
 # dividing by nu, tau, the discriminant, or a source parameter.
 j_terms=add(scale(mul(k,j0),2),scale(mul(j0,j0),-1),mul(mul(s,j1),a),mul(nu,mul(j1,j1)))
 assert add(mul(k,k),scale(j_terms,-1))==a_relation
 u,v=v[0],v[1]
 assert add(mul(add(u,v),add(u,v)),scale(mul(u,u),-1),scale(mul(v,v),-1))==scale(mul(u,v),2)
 # No numerical ring instance is used: these polynomial identities prove
 # the kernel lemma over every injective A->B with the stated monic rule.
 proof={
  'variables':['k','j0','j1','s','nu','tau'],
  'universal_identity_left':enc(a_relation),'universal_identity_right':enc(witness),
  'k_square_minus_J_expression':enc(add(mul(k,k),scale(j_terms,-1))),
  'kernel_lemma':[
   'For an arbitrary A-algebra M choose a possibly infinitely generated polynomial presentation M=A[X]/J.',
   'If k maps to zero in B tensor_A M then k is in J B[X]. Since B=A+A tau, write k=j0+j1 tau with both ji in J.',
   'A[X] injects into B[X], although B need not be A-flat. The displayed monic identity therefore holds already in A[X].',
   'It gives k² in J. Polarization and 2 invertible give products of any two kernel classes zero.',
   'Neither 1,tau independence over A nor nu/discriminant invertibility is required.'],
  'actual_base_change':[
   'The accepted same-physical-R0 source inclusions are Ai=R0+Ii taui subset Bi=R0+R0 taui, taui²=nui.',
   'Both Ai and Bi are finite over the same complete Noetherian R0, since Ii is finitely generated; the source ambient maps include R0.',
   'Each tensor over the common zero ambient is a quotient of the tensor over R0, hence finite and already complete over R0.',
   'Consequently the actual completed tensors D0,D1,E4 agree with the finite ordinary tensors, and both comparison arrows are literal base changes.',
   'Apply the arbitrary-algebra lemma separately to D0->D1 and D1->E4; both kernel ideals have square zero.'],
  'power_bookkeeping':{'each_step_kernel_power':2,'composite_kernel_power':2*2,'whole_E4_nil_power':3,'original_D0_nil_power':3*2*2},
  'composition':[
   'For K=ker(D0->E4), its image in D1 is contained in the second kernel. Thus K² is contained in the first kernel; K⁴=0.',
   'Nil(D0)³ maps to zero in the whole E4, including its embedded localizations. Hence Nil(D0)^12 subset K⁴=0.'],
  'full_germ':[
   'The accepted finite-w-polynomial model is localized before w1-completion; no w=1 substitution in a w0-completed ring is used.',
   'All positive homogeneous pieces inject into the square-zero-nil stable tail. A product of two homogeneous nil elements with positive total weight is zero.',
   'The nil ideal is torus-stable and closed; weight projections modulo maximal-ideal powers and completeness extend this identity to arbitrary formal nil series.',
   'Thus for every q>=2, Nil(Dtrans)^q equals the degree-zero power. This transfers the power12 bound.',
   'Nil ideals extend unchanged under the actual G19 power-series product; the literal scalar chart removes one regular scalar factor without changing the index.'],
  'lower_bound':[
   'The actual Qb family through P0 has nonzero nil square on a dense unit open. The coherent nil-square support is closed.',
   'It therefore contains P0. Localization, faithful completion, and the actual group product preserve nonzero nil square.'],
  'scope':'Only actual B18 germ of CA intersect CB: 3<=nil index<=12. No exact index, cube nonzero/zero, upper8, full Poisson scheme, or global exhaustive nil inventory.'}
 ev=out/'independent-quadratic-kernel-and-power-proof.json';ev.write_text(json.dumps(proof,indent=2)+'\n')
 for p in ['code/ambient_B18_quadratic_normalizer_nil_upper12_review_20261104_v63.py','research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'arbitrary_algebra_base_change_quadratic_kernel_ideal_square_zero_proved':True,
  'trace_monic_variant_and_characteristic0_polarization_verified':True,
  'actual_same_R0_two_step_completed_tensor_application_accepted':True,
  'ordinary_finite_tensor_already_complete_no_flatness_assumption':True,
  'composite_comparison_kernel_power4_zero':True,'whole_E4_nil_power3_source_bound':True,
  'actual_D0_and_whole_B18_germ_nil_power12_zero':True,
  'actual_nil_index_interval':[3,12],'exact_index_upper8_or_nil_cube_claim':False,
  'no_new_source_jet_CAS_point_rank_or_ledger':True,'input_sha256':hashes,
  'evidence_sha256':{'independent-quadratic-kernel-and-power-proof.json':hashlib.sha256(ev.read_bytes()).hexdigest()},
  'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','evidence_sha256']}),flush=True)
if __name__=='__main__':main()
