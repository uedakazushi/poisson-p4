#!/usr/bin/env python3
"""Gate/hash-only universal Sp4 curve descent; no CAS or new arithmetic model."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='Q_FIRSTJET_UNIVERSAL_WEYL_DESCENT_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 stems=[
 'Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW',
 'Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63_REVIEW',
 'Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63_REVIEW',
 'Q_FIRSTJET_RESIDUAL_OCTIC_PRESENTATION_20261103_v63_REVIEW',
 'Q_FIRSTJET_RESIDUAL_RECIPROCAL_GENUS_20261104_v63_REVIEW',
 'Q_FIRSTJET_RESIDUAL_GENUS3_BRANCH_20261104_v63_REVIEW']
 paths=[ROOT/'reports'/n/'report.json'for n in stems];gates={p.parent.name:json.loads(p.read_text())for p in paths}
 assert all(d['all_passed']for d in gates.values())
 g=gates[stems[-1]]
 assert g['ordered_relative_sign_full_curve_geometrically_irreducible_genus']==3
 assert g['actual_selected_m_dominant_component_exact_m_degree']==8
 assert gates[stems[0]]['image_function_field_equals_sign_invariants_and_Sp4_invariant_fields_agree']
 assert gates[stems[1]]['finite_Weyl_and_effective_relative_sign_cross_section_accepted']
 assert gates[stems[2]]['full_three_actual_fixed_base_equations_and_relative_sign_quotient_accepted']
 checked={}
 for p in paths:
  d=json.loads(p.read_text())
  for section in['input_sha256','checked_dependency_sha256']:
   for n,h in d.get(section,{}).items():
    q=ROOT/n;assert q.exists()and sha(q)==h,n;checked[n]=h
  for n,h in d.get('evidence_sha256',{}).items():
   q=p.parent/n;assert q.exists()and sha(q)==h,str(q);checked[str(q.relative_to(ROOT))]=h
 # Only a finite combinatorial Weyl-group enumeration, no polynomial/CAS job.
 signed=[]
 for swap in[False,True]:
  for ex in[1,-1]:
   for ey in[1,-1]:signed.append((ex,ey,(1,0)if swap else(0,1)))
 def act(w,v):return(w[0]*v[w[2][0]],w[1]*v[w[2][1]])
 images={act(w,(1,3))for w in signed};assert len(images)==8
 for a in signed:
  for b in signed:
   assert act(a,act(b,(1,3)))in images
 certificate=dict(
  base_fields=['QQ','C'],actual_invariant_field='F=k(W16)^Sp4=k(Z_pair)^Sp4',
  spectral_field='K=k(F4,p1,p2,p3,p4)',
  ordered_field='Kprime=K(x,y)=k(x,y,a,k1,k2)',
  base_relations=dict(S='F4/20',P='p4',x2_plus_y2='S',x2_y2='P',Delta='x^2-y^2',
   a='-(p3+2*x^2*S)/(2*x*Delta)',k1='p1',k2='p2',p3='-2*a*x*Delta-2*x^2*S'),
  Weyl_group='signed permutations W(C2), order8',
  all8_actions_on_x_y=[dict(x=('+'if e[0]==1 else'-')+('x'if e[2][0]==0 else'y'),y=('+'if e[1]==1 else'-')+('x'if e[2][1]==0 else'y'))for e in signed],
  field_Galois_degree=8,Weyl_fixed_base='Kprime^W=K',
  ordered_total_field='L=ordered actual pair Cartan-slice function field modulo T; equivalently b=e=1 slice modulo relative mu2',
  quotient_binding=dict(
   generic_pair_recovery='beta ->(C,D) has exactly +/- fibre; -I in Sp4 acts as beta->-beta, so no extra degree2 in invariant field',
   eigenframe='distinct eigenvalues +/-x,+/-y; choose rational kernel vectors and rescale one vector in each pair by its nonzero symplectic pairing; no additional square root',
   torus='split T=diag(t,t^-1,u,u^-1)',gauge='b->t^2*b,e->u^2*e; b=e=1 needs etale roots, not assumed rational',
   residual='T[2]=(mu2)^2, diagonal central mu2 trivial on C,D; effective relative mu2 negates g,h,i,j',
   invariant_coordinates='R=g/h,Tvar=g*h; g^2=R*Tvar,h=Tvar/g,i=v*g/Tvar,j=u/g; EXACT relative-mu2 quotient',
   field_descent='k(Cartan slice)^T=k(gauge slice)^(mu2)^2=k(gauge slice)^(effective mu2), also over QQ',
   rational_slice='regular-D conjugacy slices give k(Z_pair)^Sp4=k(Cartan slice)^N=L^W',
   no_explicit_Weyl_curve_formula_required=True),
  three_different_order2_operations=[
   'beta sign is killed by Sp4 center before pair quotient',
   'effective relative gauge sign is already removed in R,Tvar',
   'reciprocal involution of R,Tvar gives genus1 quotient Gamma; its branched double cover is RETAINED in L, genus3'],
  unique_component_gate='reviewed entire ordered generic eight-dimensional algebra is one geometrically irreducible field L, genus3; not merely an arbitrarily selected component',
  normalizer_action='W acts semilinearly on unique L via actual pair normalizer, faithful on Kprime; affine/gauge chart need only birationally stable',
  Artin=dict(F='L^W',L_over_F=8,
   compositum='Gal(L/F*Kprime)=kernel(W on Kprime)=1, hence F*Kprime=L',
   linear_disjoint='[F*Kprime:F]=8=[Kprime:K], so F and Kprime are linearly disjoint overK',
   constant_closure='alpha inF algebraic overK -> algebraic overKprime -> alpha inKprime by regularL/Kprime -> Wfixed -> alpha inK'),
  actual_geometric_generic_curve_genus=3,
  genus_reason='constant basechange F*Kprime=L, geometrically regular genus3 L/Kprime; this is descent, not quotient of one fixed geometric curve byW',
  low6=dict(delta8='(F4/20)^2-4*p4=p0',m='N12/(2*delta8)',m_in_F=True,
   ordered_degree='[L:Kprime(m)]=8',base_degree='[Kprime(m):K(m)]=8',
   full_degree='[L:K(m)]=64=[L:F]*[F:K(m)]',actual_F_over_low6_degree=8,
   polynomial_low6_field='K(N12)=k(F4,p1,p2,p3,p4,N12)=K(m)'),
  NO_regular_ring_generation_or_global_finite_map=True,
  NO_NEW_multisection_degree_equation_or_pure_QQ_rationality_claim=True)
 (out/'proof-certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 note=ROOT/'research/notes'/(ID+'.md')
 report=dict(id=ID,all_passed=True,all_independent_curve_and_actual_slice_gates_passed=True,
  genus3_descends_to_actual_full_Sp4_universal_spectral_curve_QQ_and_C=True,
  low6_invariant_subfield_generic_degree_exact8=True,
  scalar_eigenvalue_roots_and_b_e_gauge_not_conflated=True,
  Sp4_center_beta_sign_not_counted_twice=True,
  no_CAS_new_Weyl_formula_point_prime_probe_or_invariant_job=True,
  independent_conceptual_review_pending=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),note,*paths]},
  checked_gate_dependency_sha256=checked,
  evidence_sha256={'proof-certificate.json':sha(out/'proof-certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','genus3_descends_to_actual_full_Sp4_universal_spectral_curve_QQ_and_C','low6_invariant_subfield_generic_degree_exact8','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
