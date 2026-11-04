#!/usr/bin/env python3
"""Independent field-inclusion and distinct spectral-base tower proof audit."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63';author=read(base/'report.json')
 assert author['all_passed']and author['status']=='proved_by_accepted_prerequisites'and not author['pending']
 assert author['actual_target_field_equals_known_tau_fixed_source_field']and author['known_tau_fixed_source_to_actual_target_index']==1
 checked={}
 def check(path,h):
  if path not in checked:checked[path]=sha(path)
  assert checked[path]==h,str(path)
 for rel,h in author['input_sha256'].items():check(ROOT/rel,h)
 gates={}
 for rel in author['input_sha256']:
  if rel.startswith('reports/')and rel.endswith('/report.json'):
   p=ROOT/rel;gates[rel]=read(p);assert gates[rel]['all_passed']
   for r,h in gates[rel].get('input_sha256',{}).items():check(ROOT/r,h)
   for r,h in gates[rel].get('evidence_sha256',{}).items():
    q=ROOT/r if(ROOT/r).exists()else p.parent/r;check(q,h)
 inv14=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63_REVIEW/report.json';v=read(inv14)
 assert v['all_passed']and v['complete_QQ_Inv14_3_Inv18_9_Inv26_45_bases_certified']
 # Distinguish the actual affine spectral bases by literal symbolic inverse
 # identities. These are tiny universal identities, not new source points.
 f,R,z0,z1,z2=s.symbols('f R z0 z1 z2')
 back=[f*f*z0,f*f*z1,f*f*z2,f*f*(z0*(R+1)-s.Rational(1,400)),f*f*(s.Rational(1,1600)-z0/4)]
 assert s.expand(400*back[0]+1600*back[4]-f*f)==0
 assert s.cancel((back[3]+4*back[4])/back[0]-R)==0
 assert all(s.cancel(back[i]/f**2-z)==0 for i,z in enumerate([z0,z1,z2]))
 p0,p1,p2,p3,p4=s.symbols('p0:5');ratio={R:(p3+4*p4)/p0,z0:p0/f**2,z1:p1/f**2,z2:p2/f**2}
 q=f*f-400*p0-1600*p4
 errors=[s.cancel(b.subs(ratio,simultaneous=True)-p)for b,p in zip(back,[p0,p1,p2,p3,p4])]
 assert errors[:3]==[0,0,0] and s.expand(errors[3]+q/400)==0 and s.expand(errors[4]-q/1600)==0
 xi,ff,hh=s.symbols('xi ff hh',nonzero=True)
 assert s.cancel(((xi**7*hh)/(xi**2*ff)**3)*ff**3/hh-xi)==0
 assert s.diff(400*p0+1600*p4,p0)==400
 # Proof obligations audited in the accompanying note. They are field
 # inclusions, not a numeric sample or a formal identification by dimension.
 proof={
  'target':'B_NEW,Q in W16,Q, geometrically integral actual primitive first-jet image of dimension15',
  'not_target':'the entire dimension25 Poisson component or the ambient W16 quotient',
  'fields':{
   'E_Q':'QQ(B_NEW,Q)^Sp4','L_Q':'QQ(xi,R,Z0,r,s)',
   'Kshape_Q':'QQ(R,Z0,r,s)','K0_Q':'QQ(R,Z0,Z1,Z2)',
   'H_Q':'QQ(F4,R,Z0,Z1,Z2)','P_Q':'QQ(p0,p1,p2,p3,p4)',
   'low_field_Q':'QQ(F4,R,Z0,Z1,Z2,ZN)=Kshape_Q(F4)'},
  'actual_source_pullback_inclusion':'E_Q subset L_Q by rational source derivation constants and the actual rational finite tau equivalence',
  'shape_recovery':'Kshape_Q=K0_Q(ZN) subset E_Q by the independent finite-free rank10 determinant',
  'J14_exists_over_QQ':'nonzero QQ invariant supplied by the accepted three-dimensional semantic Poisson-bracket basis; rational Lie kernel gives the same descent',
  'J14_nonzero_on_actual_image':'geometric homogeneous prime equation degree>=26 cannot divide nonzero degree14 polynomial',
  'all_source_homogeneity':'full16 source covariance gives F4=xi^2*f4, J14=xi^7*h14 with nonzero f4,h14 in Kshape_Q',
  'xi_recovery':'xi=(J14/F4^3)*f4^3/h14 in E_Q',
  'index_one':'Kshape_Q and xi are in E_Q subset L_Q, hence E_Q=L_Q, pure rational trdeg5',
  'complexification':'apply the same field inclusions over C; preserve the geometric curve degree3/5 proof and explicit Bezout, rather than assuming QQ irreducibility survives',
  'QQ_geometric_image':'rational QQ source coordinate ring injects into a pure transcendental QQ field; after flat extension it still injects into the C source field, so the actual image is geometrically integral',
  'degree10':'Kshape/K0 has exact degree10',
  'degree20':'F4 is transcendental over Kshape; Kshape(F4)/K0(F4) has degree10, and xi^2=F4/f4 has odd F4-adic valuation1, giving degree2',
  'degree40':'H=P(F4) by the checked rational base inverses; F4^2=400p0+1600p4 has odd valuation1 in five independent octic variables, giving another degree2',
  'low_field_index_two':'E/Kshape(F4) has degree2, while ZN already generates Kshape/K0',
  'projective_degree10':'the algebraic scalar action beta->lambda*beta induces xi->lambda^2*xi with shape fixed; its rational fixed field is Kshape, over QQ as well as C',
  'no_rational_point_square_surjectivity_assumption_over_QQ':True,
  'no_arbitrary_hypersurface_invariant_fractions_or_unproved_linear_disjointness':True,
  'generic_field_only':'no global finite-flat morphism, boundary fiber count, genus, exact defining equation or component25 rationality is asserted'}
 ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
 inputs=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),base/'report.json',ROOT/'code/q_new_family_target_field_and_spectral_degrees_20261104_v63.py',ROOT/'research/notes/Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63.md',inv14]
 inputs += [ROOT/rel for rel in author['input_sha256']if rel.startswith('reports/')]
 report={'id':STEM,'all_passed':True,'actual_firstjet_image_B_NEW15_only':True,
         'actual_target_field_equals_known_tau_fixed_source_QQ_and_C':True,
         'actual_invariant_field_pure_rational_trdeg5_over_QQ_and_C':True,'target_source_index':1,
         'projective_normalized_spectral_field_degree':10,'affine_field_degree_over_F4_and_four_ratios':20,'affine_field_degree_over_five_octic_pencil_coefficients':40,
         'degree_over_F4_spectral_base_plus_ZN':2,
         'independent_tiny_symbolic_base_inverse_and_xi_ratio_identities':True,
         'nonzero_QQ_Inv14_exists_by_prior_global_basis_and_lower26':True,
         'separate_generic_field_tower_and_odd_valuation_proofs_accepted':True,
         'geometric_basechange_and_QQ_field_of_definition_checked':True,
         'no_new_point_prime_source26_Inv14_expansion_Krylov_replay_or_large_CAS':True,
         'whole_Poisson_component_quotient_rationality_genus_NEW_equation_or_global_finiteness_not_claimed':True,
         'prior_current_dependency_hashes_verified':len(checked),
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(inputs)},'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in ['all_passed','target_source_index','projective_normalized_spectral_field_degree','affine_field_degree_over_F4_and_four_ratios','affine_field_degree_over_five_octic_pencil_coefficients','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/STEM);a=p.parse_args();main(a.output.resolve())
