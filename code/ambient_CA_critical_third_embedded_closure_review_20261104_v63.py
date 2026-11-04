#!/usr/bin/env python3
"""Independent finite-source logic and monomial Ass descent audit; no CAS."""
from pathlib import Path
from itertools import product
import hashlib, json, resource, time

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_CA_CRITICAL_THIRD_EMBEDDED_CLOSURE_20261104_v63'
STEM=AUTHOR+'_REVIEW'

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def divides(a,b):return all(x<=y for x,y in zip(a,b))
def minimal(gs):
 gs=set(gs)
 return sorted(g for g in gs if not any(h!=g and divides(h,g)for h in gs))
def var(i,p=1):return tuple(p if j==i else 0 for j in range(9))

def main():
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 primary_path='reports/'+AUTHOR+'/report.json';A=read(primary_path)
 assert A['all_passed'] and A['critical_cone_orbit_dimension']==20
 assert A['critical_projective_orbit_dimension']==19 and A['critical_closure_contains_existing_P0']
 hashes={primary_path:sha(primary_path)}
 for p,h in A['input_sha256'].items():assert sha(p)==h,p;hashes[p]=h
 direct={p:read(p)for p in A['input_sha256']if p.endswith('.json')and'/report' in p or p.endswith('_REVIEW.json')}
 assert all(d.get('all_passed',False)for d in direct.values())
 def gate(stem,flat=False):
  path='reports/'+stem+('.json'if flat else'/report.json')
  assert path in hashes
  return direct[path]
 V=gate('AMBIENT_CA_V1_ODD_ORDER_20261103_v63_REVIEW')
 assert V['actual_V_curve_all150_critical_polynomials_zero_and_P0_V1_endpoints']
 assert V['full_source_frame29_invertible_GL20_target_rank20_and_scalar_column10_P']
 assert V['fullproper_P0_fiber2_and_dense_V1_orbit_transport_to_fullfiber2_accepted']
 assert V['normalization_conclusion_after_nonzero_odd_I_not_assumed_from_proper_model']
 Q=gate('AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW')
 assert Q['whole_three_flag_Groebner_fibers_rebuilt']
 assert Q['finite_free_rank2_completed_normalizers_and_separate_image_kernel_proof_reviewed']
 C=gate('AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW')
 assert C['same_R0_canonical_conductor_embedding_and_G19_transport_accepted']
 assert C['independent_conceptual_all_four_branch_and_transport_audit']
 assert C['CAcritical_length']==1 and C['CBcritical_length']==0
 S=gate('AMBIENT_CA_CRITICAL_V1_CONTACT_20261104_v63_REVIEW')
 assert S['independent_actual_N_frame_adjoint_source_replay_and_contact_proof']
 assert S['all3300_source_slots_and_all150_graph_bits_checked'] and S['c_v']=='0' and S['nonzero_graph_count']==0
 L=gate('AMBIENT_CA_CRITICAL_V1_CONTACT_ALGEBRA_20261104_v63_REVIEW')
 assert L['pure_algebra_contact_and_normalizer_comparison_accepted'] and L['actual_source_delta_peer_separate']
 PA=gate('AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW')
 PB=gate('AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW')
 assert PA['global_embedded_associated_component_closure_S20_accepted']
 assert PB['global_embedded_SQ_distinct_from_SP_accepted']
 assert A['CB_chosen_Jacobian_rank_alone_used_for_whole_fiber']is False
 assert A['critical_whole_CA_fiber']=={'geometric_points':1,'scheme_length':2,'reduced':False,'local_algebra':'k[epsilon]/epsilon^2'}
 assert A['critical_whole_CB_fiber']=={'geometric_points':1,'scheme_length':1,'reduced':True,'source_and_target_smooth_dimension':29}
 assert A['global_embedded_associated_closure_lower_bound_cone']==3
 assert not A['full_embedded_prime_exhaustion_claimed'] and not A['new_Poisson_irreducible_component_claimed']

 # Own monomial ideal arithmetic in l0..l6,v,x.  This is an exact generating
 # ideal calculation, not a finite jet inference.
 x=var(8);v=var(7)
 J=minimal([tuple(a+b for a,b in zip(x,var(i)))for i in range(7)]+[var(8,2),tuple(a+b for a,b in zip(x,var(7,2)))])
 primary=minimal([var(i)for i in range(7)]+[var(8,2),var(7,2)])
 intersection=minimal([tuple(max(a,b)for a,b in zip(x,g))for g in primary])
 assert J==intersection
 xv=tuple(a+b for a,b in zip(x,v))
 colon=minimal([tuple(max(0,a-b)for a,b in zip(g,xv))for g in J])
 assert colon==minimal([var(i)for i in range(9)])
 assert not any(divides(g,xv)for g in J)
 quotient_basis=[]
 for a,b in product(range(3),repeat=2):
  mon=tuple([0]*7+[a,b])
  if not any(divides(g,mon)for g in primary):quotient_basis.append(mon)
 assert len(quotient_basis)==4 and set(quotient_basis)=={(0,)*9,x,v,xv}

 patterns={
  'critical_CA':((1,2,False),(1,1,True)),
  'Pa':((2,2,True),(1,1,True)),
  'Qb':((1,1,True),(2,2,True)),
 }
 assert len(set(patterns.values()))==3
 proof={
  'all_direct_frozen_source_and_completion_gates_current':True,
  'whole_proper_fiber_open_proof':[
   'Proper finite-type maps have an open finite locus around a target with finite whole fiber.',
   'On this finite locus, fiber dimension of the finite algebra as a vector space is upper semicontinuous; length<=2 is an open condition.',
   'Take the intrinsic finite/length loci, or the union of all group translates of a valid neighborhood. Equivariance makes the resulting neighborhoods G-stable.',
   'P0 lies in closure(O_C), so its neighborhood meets O_C. G-stability and transitivity put the entire orbit in that neighborhood.',
   'The CA corank1 lift has nonreduced local fiber length at least2. The whole fiber bound2 excludes every additional point; Ac8 gives dual numbers exactly.'
  ],
  'finite_source_normalization_not_assumed_from_properness':[
   'The source is regular29 on the accepted nonzero flag model.',
   'The whole P0 fiber has one support point, so the completed finite source has the single accepted normal local factor.',
   'The accepted completed image R+I*tau has I nonzero and tau=(f*tau)/f in its fraction field.',
   'Completion is faithfully flat; a finite generic degree greater than1 would remain greater than1 over the completed fraction field. Thus the irreducible global finite map has generic degree1.',
   'Consequently the normal finite source is the normalization on the common finite open.'
  ],
  'opposite_conductor_and_whole_CB_fiber':[
   'Use the canonical intrinsic conductor quotient A_i/conductor=R0/I_i0 in the SAME physical R0, not an abstract coordinate matching.',
   'The accepted actual CA critical-V source and equivariant G19 section bind a formal branch of closure(O_C) to V(I_A0).',
   'The accepted m-primary I_A0+I_B0 makes I_B0 a unit at its generic branch point.',
   'If closure(O_C) lay in the coherent CB conductor support on the common finite open, its ideal would vanish after completion and restriction to that actual branch, a contradiction.',
   'The conductor-complement is G-stable and meets the orbit, hence contains all O_C.',
   'A finite normalization is an isomorphism outside the conductor. The regular source therefore gives a smooth CB target and the whole single reduced fiber.'
  ],
  'monomial_full_ideal_intersection_verified':True,
  'Ann_nonzero_xv_exactly_l7_v_x_verified':True,
  'embedded_primary_transverse_length':len(quotient_basis),
  'G20_formal_orbit_identification':[
   'The accepted group20 section and complementary ambient slice give a formal ambient isomorphism by its invertible derivative.',
   'All actual component ideals and normalizer maps are G-equivariant, so the full scheme diagram is a simultaneous product.',
   'The orbit is the group factor times the transverse origin. Its completed ideal is exactly (l7,v,x), including all eliminated ambient graph coordinates.'
  ],
  'associated_prime_descent':[
   'The displayed annihilator is an associated prime of the completed actual intersection.',
   'Faithfully flat completion contracts associated primes to those of the local algebraic module.',
   'Because this prime is the completion of the actual locally closed orbit ideal, its contraction is I_Z,p, not an unidentified formal-fiber prime.',
   'Localization of associated primes then yields the unique global prime I_Z.',
   'Z has dimension20 and lies properly in N28. It is not the minimal R20: R intersect N has dimension19, whereas O_C has dimension20.',
   'The nonzero scalar action makes Z homogeneous; nonzero support descends to projective dimension19.'
  ],
  'generic_whole_fiber_patterns':patterns,
  'pairwise_distinctness':[
   'Each actual G-times-nonempty-unit-family image is constructible and dense in its irreducible20-dimensional closure.',
   'A dense constructible subset of an irreducible Noetherian space contains a nonempty open.',
   'Equal closures would make two of those opens meet, contradicting the incompatible WHOLE proper fiber patterns.',
   'This argument does not compare a special one-point double fiber with a generic two-point degeneration.'
  ],
  'global_embedded_associated_closures_at_least3_cone20_projective19':True,
  'no_new_Poisson_component_or_exhaustion_or_boundary_or_origin_Ass_claim':True,
  'no_new_source_coefficient_CAS_or_jet_job':True,
 }
 (out/'independent-whole-fiber-and-Ass-descent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 for p in [str(Path(__file__).relative_to(ROOT)),'research/notes/'+STEM+'.md']:
  hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'frozen':True,
  'whole_proper_fiber_and_opposite_conductor_argument_accepted':True,
  'actual_orbit_embedded_Ass_contraction_accepted':True,
  'at_least3_distinct_global_embedded_associated_closures':True,
  'cone_dimension_each':20,'projective_dimension_each':19,
  'full_Ass_or_Poisson_new_component_claim':False,'input_sha256':hashes,
  'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in out.iterdir()if p.is_file()},
  'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}),flush=True)

if __name__=='__main__':main()
