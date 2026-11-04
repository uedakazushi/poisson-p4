#!/usr/bin/env python3
"""Static evidence binding for conditional full-Inv24 recognition logic; no coefficient job."""
import argparse
import gzip
import hashlib
import json
from math import comb
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_FULL_DEGREE24_LOGIC_20261104_v63_REVIEW'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    b=p.read_bytes()
    return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def accepted(p):
    r=read(p);assert r['all_passed']
    for key in ['input_sha256','checked_dependency_sha256']:
        for n,h in r.get(key,{}).items():assert sha(ROOT/n)==h,n
    for n,h in r.get('evidence_sha256',{}).items():
        q=ROOT/n if (ROOT/n).exists() else p.parent/n
        assert sha(q)==h,n
    return r

def main(out):
    start=time.monotonic();assert not out.exists();out.mkdir(parents=True)
    names=['Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW',
           'Q_NEW_FAMILY_SIX_SOURCE_20261103_v63_REVIEW',
           'Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW',
           'Q_NEW_FAMILY_DEGREE22_INVARIANT_EXCLUSION_20261103_v63_REVIEW',
           'Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW',
           'Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW']
    paths=[ROOT/f'reports/{n}/report.json' for n in names]
    image,six,constants,lower,char,low=[accepted(p) for p in paths]
    assert image['exact_primitive_image_dimension15_not_dominant']
    assert six['actual_six_parameter_raw_jet_full_coefficients_and_both_accepted_gauges']
    assert six['independent_original30x26_matrix_and_all20_QQ_pivot_steps']
    assert six['independent_explicit_nonzero_residual5minor_and_full26_polynomial_kernel']
    for k in ['complete_rational_source_field_inverse_and_algebraic_independence',
              'exact_source_constants_field_by_nonzero_derivative_on_rational_t',
              'symbolic_t2_section_and_all1032_nodes_composed_from_reviewed_source6']:
        assert constants[k]
    assert not constants['actual_NEW_beta_quotient_field_index_or_birationality_claimed']
    assert lower['NEW25_firstjet_hypersurface_equation_degree_lower_bound']==24
    certpath=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW/independent-certificate.json'
    primarypath=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
    own=read(certpath)['characters_through24'];primary=read(primarypath)['character_records_through24']
    for a,b in zip(own,primary):
        d=a['degree'];assert d==b['degree']
        assert a['total_weight_monomials']==comb(15+d,d)
        assert a['invariant_dimension']==b['invariant_dimension']==sum(v['sign']*v['multiplicity'] for v in b['Weyl_terms'])
        assert len(b['Weyl_terms'])==8
    dims={a['degree']:a['invariant_dimension'] for a in own}
    assert [dims[d] for d in [12,16,24]]==[9,20,67]
    assert all(dims[d]==0 for d in range(1,25,2))
    # Count formal weighted monomials only. No actual beta or source evaluation.
    def count(weights,target):
        n=[1]+[0]*target
        for w in weights:
            for d in range(w,target+1):n[d]+=n[d-w]
        return n[target]
    # Choose four Inv8 complements to F4**2; do not count that square twice.
    counts={'spectral':count([4]+[8]*4,24),
            'primitive12_times_spectral12':4*count([4]+[8]*4,12),
            'symmetric_square_primitive12':comb(4+2-1,2),
            'primitive16_times_Inv8':5}
    assert list(counts.values())==[35,20,10,5] and sum(counts.values())==70
    prepath=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_PREFLIGHT_20261104_v63/report.json';pre=accepted(prepath)
    assert pre['preflight_only_no_rank_or_invariant_expansion_job'] and pre['universal_rank67_required_before_completeness']
    extrapaths=[certpath,primarypath,prepath,ROOT/f'research/notes/{STEM}.md',Path(__file__),
                ROOT/'research/notes/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63.md',
                ROOT/'research/notes/Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW.md',
                ROOT/'research/notes/Q_NEW_FAMILY_SIX_SOURCE_20261103_v63_REVIEW.md',
                ROOT/'research/notes/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW.md']
    result=dict(id=STEM,all_passed=True,verification_kind='PURE conditional logic and frozen evidence/hash audit; no new coefficient job',
        accepted_hypersurface_irreducible_Sp4_stable_homogeneous_codimension1=True,
        accepted_full_characters={'Inv12':9,'Inv16':20,'Inv24':67},accepted_equation_degree_lower_bound=24,
        future_weight24_candidates=counts,future_candidate_count=70,
        ambient_p_integral_coefficient_minor_rank67_implies_QQ_basis_conditional=True,
        exact_raw30_kernel_plus_complete_constants_gives_invariant_section_detection=True,
        specialization_on_F5_injective_not_orbit_field_index1_required=True,
        any_nonzero_70_vector_is_not_sufficient_due_universal_relations=True,
        future_source_GF_rank67_implies_no_degree24_and_degree_lower26_conditional=True,
        source_GF_rank66_alone_does_not_imply_QQ_relation=True,
        exact_QQ_nonzero_basis_candidate_and_full_source_zero_implies_defining_degree24_conditional=True,
        actual_full_Inv24_basis_computed=False,actual_full_source_rank_computed=False,actual_G24_identity_computed=False,
        equation_found_or_degree_lower_bound_changed=False,
        no_new_beta_point_plane_expansion_native_solver_GB_or_source_coefficients=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths+extrapaths},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','actual_full_Inv24_basis_computed','actual_G24_identity_computed','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
