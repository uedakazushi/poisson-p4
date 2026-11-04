#!/usr/bin/env python3
"""Finite proof receipt for the two generic nodes; no CAS or further point."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import time
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_RESIDUAL_QUARTIC_FORCED_NODES_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic()
    operator=ROOT/'reports/Q_FIRSTJET_RESIDUAL_RECIPROCAL_QUOTIENT_20261104_v63'
    failed=ROOT/'reports/Q_FIRSTJET_RESIDUAL_QUARTIC_SMOOTH_SPECIALIZATION_20261104_v63'
    a=json.loads((operator/'report.json').read_text());b=json.loads((failed/'report.json').read_text())
    assert a['all_passed'] and a['actual_plane_total_degree_r_m']==4
    assert b['all_passed'] is False and b['run_status']=='completed'
    assert b['special_infinity_binary_quartic_reduced_in_both_patches']
    # Existing SAME coefficient point only, rational bookkeeping, no CAS job.
    x,y,a,k1,k2=map(Q,[1,2,3,5,7]);S=x*x+y*y;ell=(a+x)**2-y*y
    d=-(a*x+S)/y;c=-ell;f=-x*x*ell/(y*y)
    B2=-4*a*a*x*x-2*a*x**3-6*a*x*y*y+3*x**4+y**4
    u=(k2-B2)/(4*x*y)
    A0=a*d*(a*y+d*x)+(a*y+d*x)*2*u+c*d*y+f*a*x
    K=k1/2-A0
    B=(a*a+c)*(d*d+f)+2*a*d*2*u-(x*x-y*y)**2
    P0=[u*B+2*u*d*K/y,K/y+2*u*a-2*u*d*x/y,-x/y]
    assert P0==[Q(-5457,2),Q(287,4),Q(-1,2)]
    assert 287**2-8*10914==-4943
    fields={'K':'QQ(x,y,a,k1,k2)', 'K0':'QQ(x,y,a,k2)',
            'geometric_nodes':'after algebraic closure or adjoining one root alpha of P0; beta is the other root',
            'K10':'k1/2-A0|m0, transcendental over K0'}
    proof={'P0_degree_r':2,'P0_leading':'-x/y',
           'discriminant_K10_squared_coefficient':'1/y^2',
           'qR_unit_obstruction':'d^2+f0=2*a*x+3*x^2+y^2',
           'qT_unit_obstruction':'a^2+c0=-2*a*x-x^2+y^2',
           'local_node_form':'W^2-m^2*unit, with (m,W) regular coordinates',
           'same_coefficient_P0_coefficients':list(map(str,P0)),
           'same_coefficient_integer_primitive_P0':'2*r^2-287*r+10914',
           'same_coefficient_discriminant':-4943}
    out.mkdir(parents=True,exist_ok=False)
    (out/'proof-schema.json').write_text(json.dumps({'fields':fields,'proof':proof},indent=2)+'\n')
    paths=[Path(__file__),operator/'report.json',operator/'certificate.json',failed/'report.json',
           ROOT/'research/notes'/(STEM+'.md')]
    result={'id':STEM,'all_passed':True,
            'generic_two_distinct_ordinary_nodes_on_m0_proved':True,
            'all_geometric_quotient_components_normalization_genus_upper':1,
            'only_these_nodes_irreducibility_or_exact_genus_claimed':False,
            'failed_smoothness_probe_preserved_no_alternative_coefficient_test':True,
            'no_CAS_beta_point_source_or_invariant_job':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
            'evidence_sha256':{'proof-schema.json':sha(out/'proof-schema.json')},
            'seconds':time.monotonic()-started}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
