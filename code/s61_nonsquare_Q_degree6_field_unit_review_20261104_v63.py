#!/usr/bin/env python3
"""Independent fixed quotient-convolution verification of the entire degree6 component."""
import argparse
import hashlib
import json
import resource
import signal
import subprocess
import sys
import time
from fractions import Fraction as F
from pathlib import Path
from s61_nonsquare_finite_exception_review_20261104_v63 import ep, ua, um
from s61_nonsquare_raw_H_B_curve_components_review_20261104_v63 import poly

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_B_Q_DEGREE6_FIELD_UNIT_20261104_v63_REVIEW'
PREFIX='S61_NONSQUARE_RAW_H_B_'
END='_20261104_v63'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));sys.set_int_max_str_digits(0)
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    work=ROOT/'work/Q-degree6-independent-review';work.mkdir(parents=True,exist_ok=True)
    hashes={}
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def audit(p,passed=False):
        r=read(p)
        if passed:assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():assert sha(p.parent/name)==h,name;hashes[str((p.parent/name).relative_to(ROOT))]=h
        return r
    reports=ROOT/'reports'
    base=reports/(PREFIX+'Q_DEGREE6_FIELD_UNIT'+END)
    primary=audit(base/'report.json');assert primary['whole_selected_component_FULLU_empty_proved']
    assert primary['selected_factor_index']==5 and primary['field_degree']==6 and primary['field_ideal_matrix_shape']==[8,176]
    native=read(base/'native-certificate.json');assert native['membership_solved']
    audit(reports/(PREFIX+'Q_LINEAR_COMPONENTS'+END+'_REVIEW')/'report.json',True)
    factors=read(reports/(PREFIX+'Q_EXCEPTION_FACTOR'+END)/'certificate.json')
    product={0:F(factors['factor_constant'])}
    for f in factors['factor_inventory']:
        p=poly(f['factor']);assert all(a==d==0 for a,e,d in p)
        product=um(product,{e:c for (a,e,d),c in p.items()})
    assert product==ep(factors['source_safe_SF37'])
    selected=factors['factor_inventory'][5];factor=poly(selected['factor'])
    field={e:c for (a,e,d),c in factor.items()}
    assert selected['factor_index']==5 and max(field)==6 and field[6]==1
    qc=read(reports/(PREFIX+'COMMON_CURVE_FACTOR'+END)/'certificate.json')
    Q=poly(qc['factors'][0]['factor'])
    hc=read(reports/(PREFIX+'PAIR_MONIC_COMPLETION'+END)/'certificate.json')
    H={(a+5,e,d):c for (a,e,d),c in poly(hc['monic_H_laurent']).items()}
    source=read(reports/(PREFIX+'35_MINORS'+END)/'certificate.json')
    Pi=list(map(poly,source['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']))
    U=poly(source['FULLU_product']);assert len(Pi)==22
    paths=[Path(__file__),ROOT/'code/s61_nonsquare_Q_degree6_field_unit_review_20261104_v63.cpp',
        ROOT/'code/s61_nonsquare_finite_exception_review_20261104_v63.py',
        ROOT/'code/s61_nonsquare_raw_H_B_curve_components_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md']
    for p in paths:hashes[str(p.relative_to(ROOT))]=sha(p)
    inp=work/'fixed-quotient-input.txt'
    with inp.open('w') as f:
        def pwrite(p):
            p=ep(p) if isinstance(p,list) else p
            f.write(str(len(p))+''.join(f' {e} {c}'for e,c in sorted(p.items()))+'\n')
        def sourcewrite(p):f.write(str(len(p))+''.join(f' {a} {e} {d} {c}'for (a,e,d),c in sorted(p.items()))+'\n')
        def bwrite(p):
            for q in p:pwrite(q)
        def awrite(p):
            for q in p:bwrite(q)
        pwrite(field)
        for p in [Q,H,U]+Pi:sourcewrite(p)
        pwrite(native['q0_reduced']);pwrite(native['q1_reduced']);bwrite(native['a_inverse'])
        for p in native['Hclear_coefficients']:bwrite(p)
        for p in native['H_monic_coefficients']:bwrite(p)
        invs=native['all_field_inverse_certificates'];assert len(invs)==9
        f.write('9\n')
        for inv in invs:
            for key in ['element','inverse','Bezout_modulus','Bezout_element']:pwrite(inv[key])
        for p in native['all22_original_source_images']:awrite(p)
        awrite(native['U_image']);awrite(native['target_U_power8'])
        assert len(native['matrix_M'])==8 and all(len(row)==176 for row in native['matrix_M'])
        for row in native['matrix_M']:
            for p in row:pwrite(p)
        assert len(native['source22_cofactor_witness'])==22
        for p in native['source22_cofactor_witness']:awrite(p)
    assert inp.stat().st_size<16*1024*1024
    exe=work/'verify-native';cpp=paths[1]
    with (out/'compile.log').open('w') as log:
        subprocess.run(['c++','-O2','-std=c++17',str(cpp),'-o',str(exe),*subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
    with inp.open('r') as f:
        proc=subprocess.run([str(exe)],stdin=f,text=True,capture_output=True,check=True,timeout=max(1,55-(time.monotonic()-start)))
    (out/'native.log').write_text(proc.stdout);assert not proc.stderr
    check=json.loads(proc.stdout);assert check['all_passed']
    (out/'independent-verification.json').write_text(json.dumps(dict(
        selected_factor_index=5,field_factor_e=[[e,str(c)]for e,c in sorted(field.items())],
        exact_factor_product_and_source_pinning=True,all_used_inverses_full_QQe_Bezout1=True,
        source_ring='QQ[e,a,d]/(f6,Q,H) with a unit; each quotient is monic after the certified a-inverse',
        all_six_conjugates_and_all_a_d_roots_including_nilpotents=True,
        all176_source_ideal_columns_and_22cofactors=True,
        FULLU_empty_on_whole_selected_closed_component=True,
        no_Gaussian_rank_claim_needed_for_unit_identity=True,
        no_other_Q_components_old_residual_or_global_rank_solved=True,
        preflight_input_bytes=inp.stat().st_size,checker=check),indent=2)+'\n')
    report=dict(id=ID,all_passed=True,whole_selected_degree6_closed_component_FULLU_empty=True,
        input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},
        ignored_verifier_input_sha256=sha(inp),ignored_executable_sha256=sha(exe),
        no_author_worker_solver_inverse_GCD_Gaussian_points_primes_GB=True,
        no_scheme_length_globalB3_H7_cubic57_quartic239_claim=True,
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    p=out/'report.json';p.write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(p))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
