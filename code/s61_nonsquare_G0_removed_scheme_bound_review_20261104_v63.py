#!/usr/bin/env python3
"""Independent case coverage and ideal-theory consequence; no new solver."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from s61_nonsquare_raw_H_B_curve_components_review_20261104_v63 import poly,add,mul,scale,power
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_B_G0_REMOVED_SCHEME_BOUND_20261104_v63_REVIEW'
PREFIX='S61_NONSQUARE_RAW_H_B_';END='_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={};reports=ROOT/'reports'
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def audit(p,passed=False):
        r=read(p)
        if passed:assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():assert sha(p.parent/name)==h,name;hashes[str((p.parent/name).relative_to(ROOT))]=h
        return r
    primary=audit(reports/(PREFIX+'G0_REMOVED_SCHEME_BOUND'+END)/'report.json')
    assert primary['actual_B_rankbad_Fitting_scheme_length_upper']==468 and primary['all_G0_actual_Bbad_removed']
    gates={}
    for stem in ['MINORS_GCD','PAIR_NORM_FINITE112','ADDITIONAL_NORM_PENCIL','CURVE_COMPONENTS',
                 'Q_COMPONENT_FINITE','FINITE_EXCEPTION_REDUCTION','Q_LINEAR_COMPONENTS',
                 'Q_DEGREE6_FIELD_UNIT','Q_HIGHER_FIELD_6','Q_HIGHER_FIELD_7','Q_HIGHER_FIELD_8','OLD_RESIDUAL_LENGTH']:
        gates[stem]=audit(reports/(PREFIX+stem+END+'_REVIEW')/'report.json',True)
    assert gates['Q_LINEAR_COMPONENTS']['all_five_Q_linear_components_FULLU_empty_proved']
    assert gates['Q_DEGREE6_FIELD_UNIT']['whole_selected_degree6_closed_component_FULLU_empty']
    assert gates['Q_HIGHER_FIELD_6']['whole_selected_degree6_factor6_FULLU_empty']
    assert all(gates[f'Q_HIGHER_FIELD_{i}']['whole_selected_closed_component_FULLU_empty']for i in [7,8])
    assert gates['OLD_RESIDUAL_LENGTH']['exact_QQ_lengths148_117_468']
    # Explicit membership sources, not just set-theoretic zero-locus records.
    source=read(reports/(PREFIX+'35_MINORS'+END)/'certificate.json')
    pi=list(map(poly,source['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']));assert len(pi)==22
    assert sorted(set(i for i in source['minor35_to_distinct_generator_index']if i is not None))==list(range(22))
    Hcert=read(reports/(PREFIX+'PAIR_MONIC_COMPLETION'+END)/'certificate.json');H=poly(Hcert['monic_H_laurent'])
    assert add(mul(poly(Hcert['H_membership_f_cofactor']),pi[10]),mul(poly(Hcert['H_membership_g_cofactor']),pi[15]))==H
    assert {m:c for m,c in H.items()if m[2]==4}=={(0,0,4):F(1)}
    assert all(m[0]>=-5 and m[1]>=0 and m[2]>=0 for m in H)
    first=read(reports/(PREFIX+'PAIR_NORM'+END)/'native-certificate.json')
    second=read(reports/(PREFIX+'ADDITIONAL_MINOR_NORM'+END)/'native-certificate.json')
    common=read(reports/(PREFIX+'COMMON_CURVE_FACTOR'+END)/'certificate.json')
    G=poly(second['gcd_newK_oldh']);h=poly(first['remaining_h_coprime_U']);k=poly(second['remaining_h_coprime_U'])
    h1=poly(second['oldh_quotient_by_gcd']);k1=poly(second['newK_quotient_by_gcd'])
    assert mul(G,h1)==h and mul(G,k1)==k and G==poly(common['source_G'])
    assert F(common['factor_constant'])==1 and [f['exponent']for f in common['factors']]==[1,4]
    Q,L=map(lambda p:poly(p['factor']),common['factors']);assert mul(Q,power(L,4))==G
    old=read(reports/(PREFIX+'OLD_RESIDUAL_NORM16_INPUT'+END)/'certificate.json')
    hh=poly(old['h1_e16_monic']);kk=poly(old['k1_e25_monic']);units=list(map(F,old['source_original_QQleading_units']))
    assert units[0]and units[1]and scale(hh,units[0])==h1 and scale(kk,units[1])==k1
    assert {m:c for m,c in hh.items()if m[1]==16}=={(0,16,0):F(1)}
    assert {m:c for m,c in kk.items()if m[1]==25}=={(0,25,0):F(1)}
    assert all(min(m)>=0 and m[2]==0 for p in[hh,kk]for m in p)
    factor=read(reports/(PREFIX+'Q_EXCEPTION_FACTOR'+END)/'certificate.json')
    assert [f['degree_e']for f in factor['factor_inventory']]==[1,1,1,1,1,6,6,9,11]
    length=read(reports/(PREFIX+'OLD_RESIDUAL_LENGTH'+END)/'certificate.json')
    assert length['actual_det_degree']==148 and length['actual_det_a_valuation']==31
    assert length['exact_localized_projection_algebra_QQ_length']==148-31==117
    assert length['monic_H_d_degree']==4 and length['exact_old_containing_finite_cover_QQ_length']==117*4==468
    proof=dict(whole_G0_cases=dict(L_generic_boundary_and_all_finite_exceptions=True,
        Q_generic_root_norm_filter_preserves_all_denominator_exceptions=True,
        protected_Q_safe37_all9_closed_factors_eliminated=True,all_case_reviews_required=True),
        actual_ideal='I=I3(actual B7x3)=(22 original primitive source minors) in R=QQ[a,e,d,1/U]',
        source_H_membership_in_gen10_gen15_checked=True,
        source_norm_H_h_Gh1_K_Gk1_ideal_memberships_reused_from_independent_reviews=True,
        all_G_h1_k1_factor_and_monic_source_normalizations_rechecked=True,
        unit_upgrade='V_QQbar(I,G) empty; R finite-type and QQbar/QQ faithful imply (I,G)=R, hence G unit in R/I including nilpotents.',
        actual_scheme_surjection='QQ[a,a^-1,e,d]/(h1,k1,H)[1/U] -> R/I is surjective because H,h1,k1 belong to I.',
        exact_source_old_cover_QQ_length=468,actual_I3B_scheme_QQ_length_upper=468,
        geometric_point_upper=468,localization_and_quotient_cannot_increase_finite_QQ_dimension=True,
        no_new_membership_witness_for_G_unit_computed=True,
        no_Bbad_empty_coker_module_length_Poisson_length_full57_239_or_exhaustion_claim=True,
        no_solver_GCD_Gaussian_factorsearch_newpoint_prime_or_GB=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in[Path(__file__),ROOT/'code/s61_nonsquare_raw_H_B_curve_components_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    r=dict(id=ID,all_passed=True,all_G0_actual_Bbad_removed=True,actual_I3B_scheme_QQ_length_upper=468,
        all_nine_Q_closed_components_and_whole_L_independently_gated=True,
        input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'],report_sha256=sha(p))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
