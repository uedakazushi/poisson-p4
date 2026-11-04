#!/usr/bin/env python3
"""Independent integer character convolution and CM/finite-map logic audit."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
BASE='Q_FIRSTJET_LOW6_NON_HSOP_20261104_v63'
STEM=BASE+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())


def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    primarypath=ROOT/'reports'/BASE/'report.json';primary=read(primarypath)
    assert primary['all_passed']
    assert primary['low6_denominator_times_Hilbert_t26']==-3
    assert primary['no_HSOP_of_degree_type_4_8_8_8_8_12_by_CM']
    assert primary['selected_generically_finite_low6_morphism_is_not_finite']
    assert primary['selected_low6_zero_fiber_in_affine_invariant_quotient_has_positive_dimension']
    assert not primary['NEW_equation_degree30_inferred_or_excluded']
    folder=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63'
    fullpath=folder/'report.json';full=read(fullpath)
    peerpath=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63_REVIEW/report.json'
    peer=read(peerpath)
    assert full['all_passed'] and full['complete_QQ_Inv26_basis']
    assert peer['all_passed'] and peer['independent_full_Weyl_DP_through26']
    assert peer['complete_QQ_Inv14_3_Inv18_9_Inv26_45_bases_certified']
    cp=folder/'character-records.json'
    assert sha(cp)==peer['input_sha256'][str(cp.relative_to(ROOT))]
    assert sha(cp)==full['evidence_sha256'][cp.name]
    # Independently reconstruct all eight C2 denominator shifts from four roots.
    roots=[(2,0),(0,2),(1,1),(1,-1)]
    shifts=Counter()
    for mask in range(16):
        weight=tuple(sum(roots[i][j] for i in range(4) if mask>>i&1) for j in range(2))
        shifts[weight]+=(-1)**mask.bit_count()
    shifts={w:c for w,c in shifts.items() if c}
    assert len(shifts)==8
    records=read(cp);assert [r['degree'] for r in records]==list(range(27))
    H=[]
    for r in records:
        terms=r['Weyl_terms']
        assert len(terms)==8 and {tuple(t['weight']):t['sign'] for t in terms}==shifts
        assert all(isinstance(t['multiplicity'],int) and t['multiplicity']>=0 for t in terms)
        value=sum(shifts[tuple(t['weight'])]*t['multiplicity'] for t in terms)
        assert value==r['invariant_dimension'] and value>=0
        H.append(value)
    dims={d:H[d] for d in [2,6,10,14,18,22,26]}
    assert list(dims.values())==[0,0,0,3,9,21,45]
    # Directly sum all 64 subsets, independently of sequential polynomial products.
    degrees=[4,8,8,8,8,12]
    D=Counter();coefficient=0
    for mask in range(64):
        degree=sum(degrees[i] for i in range(6) if mask>>i&1)
        sign=(-1)**mask.bit_count();D[degree]+=sign
        if degree<=26:coefficient+=sign*H[26-degree]
    D={d:c for d,c in D.items() if c}
    assert {d:c for d,c in D.items() if d<=26}=={0:1,4:-1,8:-4,12:3,16:7,20:-2,24:-8}
    assert coefficient==45-21-4*9+3*3==-3
    independence=ROOT/'reports/Q_FIRSTJET_SIX_LOW_INVARIANT_INDEPENDENCE_20261104_v63/report.json'
    ir=read(independence)
    assert ir['all_passed'] and ir['exact_selected_F4_p1_p2_p3_p4_N12_are_algebraically_independent']
    dimension=ROOT/'reports/Q_FIRSTJET_SHEPHERD_BARRON_RATIONALITY_20261104_v63/report.json'
    dr=read(dimension);assert dr['all_passed'] and dr['transcendence_degree']==6
    checked={}
    for p in [primarypath,fullpath,peerpath,independence,dimension]:
        data=read(p)
        for key in ['input_sha256','checked_dependency_sha256']:
            for name,want in data.get(key,{}).items():
                target=ROOT/name;assert sha(target)==want,str(target)
                checked[str(target.relative_to(ROOT))]=want
        for name,want in data.get('evidence_sha256',{}).items():
            target=p.parent/name;assert sha(target)==want,str(target)
            checked[str(target.relative_to(ROOT))]=want
    proof={
        'eight_Weyl_denominator_shifts_from_16_root_subsets':[[list(w),c] for w,c in sorted(shifts.items())],
        'full_Hilbert_dimensions_degree0_to26':H,
        'selected_dimensions':dims,'HSOP_degree_type':degrees,
        'denominator_coefficients_from64_subsets':D,
        'coefficient_t26':coefficient,
        'full_character_upper_bound_not_candidate_lower_bound':True,
        'CM_application':'C[W16]^Sp4 is CM; homogeneous SOP length6 is regular; quotient Hilbert coefficients cannot be negative',
        'specific_six_independence_chain':'accepted rank6 plus determinant 400*p0/F4 of exact differential change',
        'graded_finiteness_bridge':'finite quotient basis lifts; induction by positive degree proves module generation without assuming finiteness',
        'positive_dimensional_graded_zero_fiber':True,
        'ordered_boundary_scope':'geometric points F4=p4=0 force x=y=0; not an equality of scheme ideals',
        'new_nonnull_point_or_source26_or_equation_degree_claimed':False,
    }
    out.mkdir(parents=True,exist_ok=False)
    evidence=out/'independent-integer-proof.json'
    evidence.write_text(json.dumps(proof,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),
        ROOT/'research/notes'/(BASE+'.md'),primarypath,fullpath,peerpath,cp,independence,dimension]
    for p in inputs:checked[str(p.relative_to(ROOT))]=sha(p)
    result={'id':STEM,'all_passed':True,
        'independent_C2_Weyl_shifts_full_character_integer_sums_and64_subset_convolution':True,
        'full_Inv26_dimension45_not_merely_constructed_candidates':True,
        'denominator_times_Hilbert_t26':-3,
        'CM_forbids_any_HSOP_degree_type_4_8_8_8_8_12':True,
        'specific_six_generically_finite_but_not_module_finite':True,
        'graded_Nakayama_degree_induction_positive_dimensional_origin_fiber':True,
        'spectral_boundary_statement_is_geometric_not_scheme_ideal_equality':True,
        'NEW_degree30_proved_or_disproved':False,
        'nonnull_point_constructed_or_source26_test':False,
        'new_Sym_DP_CAS_GB_point_or_invariant_job':False,
        'input_sha256':checked,'evidence_sha256':{'independent-integer-proof.json':sha(evidence)},
        'seconds':time.monotonic()-started,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','denominator_times_Hilbert_t26','seconds']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
