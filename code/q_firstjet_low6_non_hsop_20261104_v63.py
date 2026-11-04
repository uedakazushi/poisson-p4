#!/usr/bin/env python3
"""Complete frozen characters plus tiny integer obstruction to low6 HSOP."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_LOW6_NON_HSOP_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    folder=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63'
    peerpath=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63_REVIEW/report.json'
    primary=read(folder/'report.json');peer=read(peerpath)
    assert primary['all_passed'] and primary['complete_QQ_Inv26_basis']
    assert peer['all_passed'] and peer['independent_full_Weyl_DP_through26']
    assert peer['complete_QQ_Inv14_3_Inv18_9_Inv26_45_bases_certified']
    cp=folder/'character-records.json'
    assert sha(cp)==primary['evidence_sha256'][cp.name]
    assert sha(cp)==peer['input_sha256'][str(cp.relative_to(ROOT))]
    records=read(cp);assert [r['degree'] for r in records]==list(range(27))
    H=[]
    for rec in records:
        assert len(rec['Weyl_terms'])==8
        value=sum(t['sign']*t['multiplicity'] for t in rec['Weyl_terms'])
        assert value==rec['invariant_dimension'] and value>=0
        H.append(value)
    assert [H[d] for d in [2,6,10,14,18,22,26]]==[0,0,0,3,9,21,45]
    independence=ROOT/'reports/Q_FIRSTJET_SIX_LOW_INVARIANT_INDEPENDENCE_20261104_v63/report.json'
    ir=read(independence)
    assert ir['all_passed'] and ir['exact_selected_F4_p1_p2_p3_p4_N12_are_algebraically_independent']
    dimension=ROOT/'reports/Q_FIRSTJET_SHEPHERD_BARRON_RATIONALITY_20261104_v63/report.json'
    dr=read(dimension);assert dr['all_passed'] and dr['transcendence_degree']==6
    degrees=[4,8,8,8,8,12];D={0:1}
    for d in degrees:
        new=dict(D)
        for j,c in D.items():new[j+d]=new.get(j+d,0)-c
        D={j:c for j,c in new.items() if c}
    numerator=[sum(c*H[n-j] for j,c in D.items() if j<=n) for n in range(27)]
    assert numerator[26]==-3
    terms=[{'denominator_degree':j,'denominator_coefficient':c,
            'Hilbert_degree':26-j,'complete_Hilbert_dimension':H[26-j],
            'contribution':c*H[26-j]}
           for j,c in sorted(D.items()) if j<=26]
    assert sum(t['contribution'] for t in terms)==45-21-4*9+3*3==-3
    out.mkdir(parents=True,exist_ok=False)
    cert={'complete_Hilbert_coefficients_degree0_through26':H,
          'degree26_is_complete_character_upper_not_only_candidate_rank':True,
          'degree_type':degrees,'all_denominator_polynomial_coefficients':D,
          'product_coefficients_through26':numerator,'degree26_exact_terms':terms,
          'standard_theorem_assumptions':{'group':'Sp4 over C, linearly reductive',
                                         'regular_ring':'C[W16]',
                                         'R':'C[W16]^Sp4, Cohen-Macaulay by Hochster-Roberts',
                                         'dimension':6,'selected_six_independent':True},
          'bibliography':{'authors':['M. Hochster','J. L. Roberts'],
                          'title':'Rings of invariants of reductive groups acting on regular rings are Cohen-Macaulay',
                          'journal':'Advances in Mathematics','volume':13,'issue':2,'year':1974,'pages':'115-175',
                          'doi':'10.1016/0001-8708(74)90067-X',
                          'author_metadata':'https://www-users.cse.umn.edu/~rober005/Publications.html'},
          'no_full_Hilbert_series_closed_form_claim':True}
    (out/'integer-certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),folder/'report.json',cp,
           peerpath,independence,dimension]
    result={'id':STEM,'all_passed':True,
            'complete_dimension_checks':{'14':H[14],'18':H[18],'22':H[22],'26':H[26]},
            'low6_denominator_times_Hilbert_t26':-3,
            'no_HSOP_of_degree_type_4_8_8_8_8_12_by_CM':True,
            'selected_generically_finite_low6_morphism_is_not_finite':True,
            'selected_low6_zero_fiber_in_affine_invariant_quotient_has_positive_dimension':True,
            'nonnull_low6_zero_point_constructed_or_evaluated':False,
            'NEW_equation_degree30_inferred_or_excluded':False,
            'full_Inv26_not_only_subset_and_no_source26_restriction':True,
            'no_R_rebuild_Sym_DP_CAS_GB_new_point_or_invariant_job':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
            'evidence_sha256':{'integer-certificate.json':sha(out/'integer-certificate.json')},
            'seconds':time.monotonic()-started,
            'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'numerator_t26':-3,'seconds':result['seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
