#!/usr/bin/env python3
"""ONE authorized coefficient specialization of the universal plane quartic.

Not a beta/source point.  Stop after this case, including if it is singular.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import resource
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_RESIDUAL_QUARTIC_SMOOTH_SPECIALIZATION_20261104_v63'
P=32003
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/Q_FIRSTJET_RESIDUAL_RECIPROCAL_QUOTIENT_20261104_v63'
    report=json.loads((source/'report.json').read_text())
    cert=json.loads((source/'certificate.json').read_text())
    assert report['all_passed'] and report['actual_plane_total_degree_r_m']==4
    point=[1,2,3,5,7]  # x,y,a,k1,k2: ONE fixed coefficient point, no alternatives.
    assert point[0]*point[1]*(point[0]**2-point[1]**2)%P
    coefficients={}
    for key,terms in cert['actual_primitive_plane_coefficients_r_m'].items():
        value=0
        for ex,co in terms:
            q=Fraction(co);assert q.denominator%P
            term=q.numerator*pow(q.denominator,-1,P)%P
            for v,e in zip(point,ex):term=term*pow(v,e,P)%P
            value=(value+term)%P
        coefficients[key]=value
    assert coefficients['4,0']==16*point[0]**4*point[1]**8%P
    expression=' + '.join(f'{co}*r^{i}*m^{j}*z^{4-i-j}'
                          for key,co in sorted(coefficients.items())
                          for i,j in [list(map(int,key.split(',')))] if co)
    script='''S=ZZ/32003[r,m,z,MonomialOrder=>GRevLex];
F=EXPRESSION;
J=ideal(diff(r,F),diff(m,F),diff(z,F));
Ar=(1_S % gb(J+ideal(r-1_S))==0_S);
Am=(1_S % gb(J+ideal(m-1_S))==0_S);
Az=(1_S % gb(J+ideal(z-1_S))==0_S);
print("PATCH_R "|toString Ar);
print("PATCH_M "|toString Am);
print("PATCH_Z "|toString Az);
H=substitute(F,{z=>0_S});
Ir=(1_S % gb(ideal(z,r-1_S,diff(r,H),diff(m,H)))==0_S);
Im=(1_S % gb(ideal(z,m-1_S,diff(r,H),diff(m,H)))==0_S);
print("INFINITY_R "|toString Ir);
print("INFINITY_M "|toString Im);
if Ar and Am and Az then print("PROJECTIVE_SMOOTH_PASS") else print("PROJECTIVE_SINGULAR_STOP");
if Ir and Im then print("INFINITY_TRANSVERSE_PASS") else print("INFINITY_REPEATED_STOP");
'''.replace('EXPRESSION',expression)
    m2=out/'smoothness.m2';m2.write_text(script)
    status='pending';elapsed=None;returncode=None
    try:
        t0=time.monotonic()
        proc=subprocess.run(['M2','--script',str(m2.resolve())],capture_output=True,text=True,timeout=55)
        elapsed=time.monotonic()-t0;returncode=proc.returncode
        log=proc.stdout+proc.stderr
        status='completed' if proc.returncode==0 else 'M2_failure_STOP_no_retry'
    except subprocess.TimeoutExpired as exc:
        log=(exc.stdout or b'')+(exc.stderr or b'')
        if isinstance(log,bytes):log=log.decode(errors='replace')
        status='timeout_STOP_no_retry'
    (out/'M2.log').write_text(log)
    smooth=returncode==0 and 'PROJECTIVE_SMOOTH_PASS' in log
    infinity=returncode==0 and 'INFINITY_TRANSVERSE_PASS' in log
    checks={key:(f'{key} true' in log) for key in ['PATCH_R','PATCH_M','PATCH_Z','INFINITY_R','INFINITY_M']}
    evidence={'coefficient_variable_order':['x','y','a','k1','k2'],
              'one_authorized_coefficient_specialization':point,'prime':P,
              'coefficient_point_is_NOT_beta_or_source_point':True,
              'all15_specialized_coefficients_r_m':coefficients,
              'homogeneous_special_F':expression,
              'good_model_r4_coefficient':coefficients['4,0'],
              'proper_flat_model':'Proj S[r,m,z]/F4, S=Z_(32003)[x,y,a,k1,k2,1/(x*y*(x^2-y^2))]; r4 unit gives graded free monic quotient',
              'patch_results':checks,'projectively_geometrically_smooth':smooth,
              'infinity_binary_quartic_separable':infinity,
              'no_other_coefficient_specialization_attempted':True}
    (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    result={'id':STEM,'all_passed':smooth and infinity,'run_status':status,
            'single_coefficient_specialization_only_not_beta_point':True,
            'special_projective_Jacobian_ideal_empty_in_all3_patches':smooth,
            'special_infinity_binary_quartic_reduced_in_both_patches':infinity,
            'generic_plane_quartic_geometrically_smooth_by_proper_flat_model':smooth,
            'generic_plane_quartic_geometrically_irreducible_and_genus3':smooth,
            'generic_four_distinct_infinity_points':infinity,
            'exact_full_cover_genus_claimed':False,
            'selected_actual_degree8_or_full_curve_irreducibility_claimed_here':False,
            'no_new_beta_source_invariant_job_or_alternative_probe':True,
            'input_sha256':{str(path.relative_to(ROOT)):sha(path) for path in [Path(__file__),source/'report.json',source/'certificate.json']},
            'evidence_sha256':{name:sha(out/name) for name in ['smoothness.m2','M2.log','certificate.json']},
            'M2_returncode':returncode,'M2_seconds':elapsed,
            'seconds':time.monotonic()-started,
            'max_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['all_passed','run_status','special_projective_Jacobian_ideal_empty_in_all3_patches','special_infinity_binary_quartic_reduced_in_both_patches','seconds']}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
