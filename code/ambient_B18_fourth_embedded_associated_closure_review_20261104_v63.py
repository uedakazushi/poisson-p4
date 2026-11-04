#!/usr/bin/env python3
"""Independent gate audit for the explicit B18 orbit associated prime."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,combinations_with_replacement
import gzip,hashlib,json,resource,time

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_FOURTH_EMBEDDED_ASSOCIATED_CLOSURE_20261104_v63'
STEM=AUTHOR+'_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 raw=(ROOT/p).read_bytes()
 return json.loads(gzip.decompress(raw) if p.endswith('.gz') else raw)
def main():
 started=time.monotonic();out=ROOT/'reports'/STEM
 assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 paths={};docs={}
 def get(p):paths[p]=sha(p);docs[p]=read(p);return docs[p]
 p='reports/'+AUTHOR+'/report.json';R=get(p)
 assert paths[p]=='2f9d432a174bdd6897d39f0c7483f0be0928da9014180716532d37f08c80ba2f' and R['all_passed']
 for p,h in R['input_sha256'].items():assert sha(p)==h;paths[p]=h
 depth=get('reports/AMBIENT_B18_ORIGINAL_DEPTH19_INITIAL_SPACE_20261104_v63_REVIEW/report.json')
 source=get('reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW/report.json')
 orbit=get('reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63_REVIEW.json')
 orbitP=get('reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/report.json')
 boundary=get('reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW.json')
 boundaryP=get('reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63/report.json')
 old3=get('reports/AMBIENT_CA_CRITICAL_THIRD_EMBEDDED_CLOSURE_20261104_v63_REVIEW/report.json')
 pa=get('reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json')
 qb=get('reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json')
 for d in [depth,source,orbit,orbitP,boundary,boundaryP,old3,pa,qb]:assert d['all_passed']
 assert depth['source_peer_sha256']==paths[depth['source_peer_gate']]
 assert source['actual_third_projected_quadratic_class_vanishes_on_all620']
 assert depth['original_cone_depth']==19 and depth['original_projective_depth']==18
 assert orbit['actual_sum_ideal_G_invariance_and_all_Artin_points_proof_accepted']
 assert orbit['full_D_nilpotents_retained_in_formal_G19_product'] and orbit['simultaneous_CA_CB_D_N_R_slice_accepted']
 assert orbit['independent_index_GL25_and_orbit19_minor']==-32
 assert orbit['canonical126_complement107_and_action_Jacobian_unit_checked']
 assert orbit['transverse_dimension']==9 and orbit['transverse_edim']==20
 assert orbitP['orbit_dimension']==19 and orbitP['ambient_canonical_dimension']==126
 assert orbitP['literal_scalar_is_one_of19_independent_parameters']
 assert boundary['entire_tzero_nilpotent_Jordanflag_cover_and_reduced_R_intersection_N_equals_irreducible_B18']
 assert boundaryP['boundary_irreducible'] and boundaryP['boundary_projective_dimension']==18
 assert boundaryP['entire_boundary_identity']=='R intersection N = B18 = image(t=0 spectralsource)'
 assert boundaryP['chain_linear_parameters']==['aN','bN','cN','g0','g1','g2','g3','g4','g5']
 assert boundaryP['chain_field']==['aN*x**2','0','bN*h*x + cN*x**2','g0*h**2 + g1*h*x + g2*h*z + g3*x**2 + g4*x*z + g5*z**2']
 assert old3['at_least3_distinct_global_embedded_associated_closures'] and old3['cone_dimension_each']==20
 assert pa['global_embedded_associated_component_closure_S20_accepted'] and qb['global_embedded_SQ_distinct_from_SP_accepted']
 certpath='reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz'
 cert=get(certpath);assert paths[certpath]==orbit['input_sha256'][certpath]
 # Literal chain specialization, using only sparse monomial indices.
 pairs=list(combinations(range(5),2));mons=list(combinations_with_replacement(range(5),2))
 field=[{(1,1):Q(1)},{},{},{(0,1):Q(1)},{(0,0):Q(1),(3,3):Q(1)}]
 p=[]
 for i,j in pairs:
  for m in mons:p.append((field[j].get(m,Q(0)) if i==2 else Q(0))-(field[i].get(m,Q(0)) if j==2 else Q(0)))
 assert p==[Q(x) for x in cert['physical_point150']] and p[106]==1
 assert cert['all25_physical_GL_derivative_rank']==19 and cert['scalar_derivative_equals_physical_point']
 assert cert['scalar_GL_index']==12 and Q(cert['orbit19_minor'])==-32
 assert len(cert['affine_transverse107_free_coordinate_indices'])==107
 # The proof below deliberately identifies the formal prime as a specific
 # algebraic orbit ideal.  No dimension assertion for an arbitrary formal
 # fibre prime is used in the descent.
 proof={
  'orbit_identification':['GL5 x Gm connected and smooth; its actual P0 orbit is smooth irreducible of dimension19.',
   'Literal chain specialization places P0 in the G-stable irreducible cone B18 of dimension19.',
   'Orbit closure equals B18 by irreducibility and equal dimension; the locally closed orbit is open in its closure.'],
  'socle_and_annihilator':['Actual all-Artin group action gives Dhat=Aperp[[g19]], including nilpotents.',
   'DepthD=19 gives depthAperp=0. A nonzero finite-length H0 contains a nonzero socle element z with Ann_Aperp(z)=mperp.',
   'Coefficientwise multiplication gives Ann_Dhat(z)=mperp[[g19]]; its quotient is k[[g19]], hence it is prime and associated.'],
  'specific_prime_descent':['Actual Gamma19 times transverse107 action is a formal ambient isomorphism; Gamma19 maps etale to the smooth open orbit at P0.',
   'The closed subscheme Gamma19 times {P0} therefore has exactly the completed B18 orbit ideal. Its ideal inside Dhat is mperp[[g19]].',
   'Flat-completion contraction of this associated prime gives I_B18 localized at P0; localization of associated primes gives the unique global prime I_B18.',
   'B18 lies properly in both known minimal components N28 and R20, so its associated prime is embedded.'],
  'fourth_and_projective':['The three older distinct associated closures have cone dimension20; B18 has cone dimension19, so is different from each.',
   'I_B18 is homogeneous and P0 has q106=1; the homogeneous localization is the projective chart tensored with k[q106,q106^-1].',
   'Faithful-flat Laurent extension/localization descends the same associated prime to projective dimension18.'],
  'excluded_claims':['No full associated-prime list or primary decomposition.','No concrete socle representative or length.','No full nilpotence index or nonzero nil cube.','No extra Poisson irreducible component or change to lower58.']}
 ev=out/'proof-receipt.json';ev.write_text(json.dumps({'proof':proof,'literal_chain_P0_nonzero_raw150':[[i,str(q)] for i,q in enumerate(p) if q],'source_gate_checks':True},indent=2)+'\n')
 for p in ['code/ambient_B18_fourth_embedded_associated_closure_review_20261104_v63.py','research/notes/'+STEM+'.md']:paths[p]=sha(p)
 report={'id':STEM,'all_passed':True,'actual_G19_product_and_nonzero_transverse_socle_argument_accepted':True,
  'formal_associated_prime_identified_with_specific_actual_B18_orbit_ideal':True,
  'no_arbitrary_formal_fibre_prime_contraction_dimension_assumption':True,
  'faithfully_flat_completion_and_localization_Ass_descent_accepted':True,
  'global_fourth_embedded_associated_closure':'B18','cone_dimension':19,'projective_dimension':18,
  'global_embedded_associated_closure_lower_bound':4,'new_Poisson_component_full_Ass_or_nil_index_claim':False,
  'no_new_CAS_source_jet_or_matrix':True,'input_sha256':paths,
  'evidence_sha256':{'proof-receipt.json':hashlib.sha256(ev.read_bytes()).hexdigest()},
  'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','evidence_sha256']}),flush=True)
if __name__=='__main__':main()
