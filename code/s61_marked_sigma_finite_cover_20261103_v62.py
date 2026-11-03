#!/usr/bin/env python3
"""Audit the pencil cover; keep the OLD double-junction bridge explicit."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
PRIME=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    h=s.Symbol('h');dependencies={}
    def load(name, passed_key='all_passed'):
        path=ROOT/'reports'/name;data=json.loads(path.read_text());assert data[passed_key]
        for filename,expected in data.get('input_sha256',{}).items():assert sha(ROOT/filename)==expected,filename
        for filename,expected in data.get('evidence_sha256',{}).items():assert sha(path.parent/filename)==expected,filename
        for filename,expected in data.get('output_sha256',{}).items():assert sha(path.parent/filename)==expected,filename
        dependencies[str(path.relative_to(ROOT))]=sha(path);return data
    def determinant(record):
        line=next(line for line in record['output'].splitlines()if line.startswith('DETERMINANT_COEFFICIENTS '))
        numbers=list(map(int,line.split()[1:]));assert numbers[0]==len(numbers)-1
        return s.Poly.from_list(list(reversed(numbers[1:])),h,modulus=PRIME).monic()
    full=load('S61_MARKED_SIGMA_PENCIL_FULL5_20261103/report.json')
    full5=next(r for r in full['results']if r['task']=='full:5');zero5=next(r for r in full['results']if r['task']=='zero:5')
    assert full5['shape']==[2603,1432]and zero5['shape']==[2022,932]
    expected=(h*h+9314)*(h*h+2)**2*(h*h-16)**2*(h*h-15105)**2*(h*h-12033)**4
    assert determinant(full5)==s.Poly(expected,h,modulus=PRIME).monic()
    assert determinant(zero5)==s.Poly(h*h+2,h,modulus=PRIME)
    heavy=load('S61_MARKED_SIGMA_PENCIL_RANK_20261103/report.json');heavy6=next(r for r in heavy['results']if r['task']=='heavy:6')
    assert heavy6['shape']==[140,23]
    assert determinant(heavy6)==s.Poly((h*h+2)**2*(h*h+9314),h,modulus=PRIME).monic()
    new=load('S61_MARKED_SIGMA_NEW_QUADRATIC_20261103_v62/report.json')
    assert new['orbit']=='T6'and new['transverse_length_upper']==2 and new['ambient_homogenizer_zero']
    assert new['zero_rho_reduced_support_exhausted_by_new_axis']
    extras=load('S61_MARKED_SIGMA_FINITE_EXTRAS_20261103_v62/report.json')
    assert {r['orbit']for r in extras['cases']}=={'T4','T7'}
    for r in extras['cases']:
        assert r['ambient_t_zero']and r['transverse_length_upper']==r['known_generic_same_incidence_length']
        num=52 if r['orbit']=='T4'else 16
        assert all((125*root*root+num)%PRIME==0 for root in r['residue_h_roots'])
    standard=load('S61_MARKED_SIGMA_STANDARD_FLAG_SPREAD_20261103_v62.json')
    assert standard['ambient_t_zero']and standard['all_h_included']and standard['mixed_characteristic_relative_free_rank']==5
    quadratics=load('S61_MARKED_SIGMA_RHO_ONE_QUADRATICS_20261103_v62/report.json')
    assert len(quadratics['cases'])==3
    assert all(r['length']==4 and r['origin_away_saturation_ideal_one']for r in quadratics['cases'])
    pure=load('S61_MARKED_SIGMA_PURE_PRIME_20261103_v62.json')
    assert pure['h']==4 and pure['three_points_remain_distinct_mod_prime']
    assert all(r['ambient_t_zero_over_V']for r in pure['cases'])
    infinity=load('S61_MARKED_SIGMA_PENCIL_INFINITY_PRIME_20261103_v62/report.json')
    assert infinity['core']==[0,0,1]and infinity['group_restoration_minor']=='-2'
    assert infinity['exact_reduced_special_support']=='standard rho axis'
    assert infinity['full_ambient_mixed_characteristic_homogenizer_nilpotent']
    assert infinity['actual_saturated_support_to_origin_weight_lemma']
    for root in [3083,28920]:
        old=load(f'S61_MARKED_SIGMA_OLD_FULL_SUPPORT_20261103_v62/root{root}/report.json','special_fiber_reduced_support_proved')
        assert old['root']==root and old['includes_rho']and old['returncode']==0
        assert old['reduced_support']=='rho-axis union weighted (2,3) plane'
    old_support_review=load('S61_MARKED_SIGMA_OLD_FULL_SUPPORT_20261103_v62_REVIEW.json')
    assert old_support_review['special_fiber_support_claim_accepted']
    old_curve=load('S61_MARKED_SIGMA_OLD_CURVE_COVER_20261103_v62/report.json')
    assert len(old_curve['results'])==2
    assert all(r['all_distinct_root_finite_k_upper2_proved']for r in old_curve['results'])
    old_curve_review=load('S61_MARKED_SIGMA_OLD_CURVE_COVER_20261103_v62_REVIEW.json')
    assert old_curve_review['mixed_characteristic_local_removal_off_special_discriminant_accepted']
    old_lower=load('S61_MARKED_SIGMA_OLD_DOUBLE_QUOTIENT_20261103_v62_FINAL/report.json')
    assert old_lower['exact_dual_number_thickening_with_t_zero']and old_lower['same_incidence_generic_length_lower']==2
    old_points=load('S61_MARKED_SIGMA_OLD_PLANE_AMBIENT_20261103_v62/report.json')
    old_distinct=[r for r in old_points['cases']if r['case']in ['distinct-u','distinct-v']and r['beta']=='actual']
    assert len(old_distinct)==2
    assert all(r['transverse_length_upper']==2 and r['quadratic_minor_norm_prime_unit']for r in old_distinct)
    old_note=ROOT/'research/notes/S61_MARKED_SIGMA_OLD_DISTINCT_ROOTS_20261103_v62.md'
    dependencies[str(old_note.relative_to(ROOT))]=sha(old_note)
    note=ROOT/'research/notes/S61_MARKED_SIGMA_PENCIL_FINITE_COVER_20261103_v62.md'
    return dict(id='S61_MARKED_SIGMA_FINITE_COVER_20261103_v62',all_passed=True,
                exceptional_radical_mod_prime=str(s.Poly((h*h+2)*(h*h+9314)*(h*h-16)*(h*h-15105)*(h*h-12033),h,modulus=PRIME).monic().as_expr()),
                all_unexceptional_finite_mod_prime_points_covered_by_standard_axis=True,
                new_quadratic_rho_zero_and_nonzero_covered=True,
                finite_T4_T7_extra_points_covered=True,h4_pure_extra_points_spread=True,
                old_quadratic_rho_nonzero_covered=True,
                old_quadratic_special_support_exhausted=True,
                old_plane_all_special_distinct_root_points_and_infinity_covered=True,
                pencil_projective_infinity_mixed_characteristic_covered=True,
                next_lemma='Apply the actual t-saturated-support finite-weight origin lemma, not saturation-specialization commutation.',
                remaining_conditions_not_claimed_here=[
                    'The OLD h^2+2=0 plane special-fiber discriminant-zero junction has a mixed-characteristic ambient nilpotent neighborhood.'],
                conditional_whole_pencil_implication_proved=True,whole_pencil_removal_claimed=False,
                input_sha256={**dependencies,str(Path(__file__).relative_to(ROOT)):sha(Path(__file__)),str(note.relative_to(ROOT)):sha(note)})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS',report['id'],'conditional; remaining bridges explicit')
