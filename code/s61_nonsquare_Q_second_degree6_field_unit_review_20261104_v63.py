#!/usr/bin/env python3
"""Fixed factor-index6 verification using the frozen independent degree6 checker."""
import argparse,hashlib,json,resource,signal,subprocess,sys,time
from fractions import Fraction as F
from pathlib import Path
from s61_nonsquare_finite_exception_review_20261104_v63 import ep,um
from s61_nonsquare_raw_H_B_curve_components_review_20261104_v63 import poly
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_B_Q_HIGHER_FIELD_6_20261104_v63_REVIEW'
PREFIX='S61_NONSQUARE_RAW_H_B_';END='_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));sys.set_int_max_str_digits(0)
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);work=ROOT/'work/Q-second-degree6-independent-review';work.mkdir(parents=True,exist_ok=True);hashes={}
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def audit(p,passed=False):
        r=read(p)
        if passed:assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():assert sha(p.parent/name)==h,name;hashes[str((p.parent/name).relative_to(ROOT))]=h
        return r
    reports=ROOT/'reports';base=reports/(PREFIX+'Q_HIGHER_FIELD_6'+END);r=audit(base/'report.json')
    assert r['selected_factor_index']==6 and r['field_degree']==6 and r['whole_selected_component_FULLU_empty_proved']
    audit(reports/(PREFIX+'Q_DEGREE6_FIELD_UNIT'+END+'_REVIEW')/'report.json',True)
    n=read(base/'native-certificate.json.model-checkpoint.json');main=read(base/'native-certificate.json');assert not(set(n)&set(main));n.update(main);assert n['membership_solved']
    factors=read(reports/(PREFIX+'Q_EXCEPTION_FACTOR'+END)/'certificate.json');product={0:F(factors['factor_constant'])}
    for factor in factors['factor_inventory']:
        p=poly(factor['factor']);assert all(a==d==0 for a,e,d in p);product=um(product,{e:c for (a,e,d),c in p.items()})
    assert product==ep(factors['source_safe_SF37']);factor=factors['factor_inventory'][6];p=poly(factor['factor']);field={e:c for (a,e,d),c in p.items()};assert max(field)==6 and field[6]==1
    qc=read(reports/(PREFIX+'COMMON_CURVE_FACTOR'+END)/'certificate.json');Q=poly(qc['factors'][0]['factor'])
    hc=read(reports/(PREFIX+'PAIR_MONIC_COMPLETION'+END)/'certificate.json');H={(a+5,e,d):c for (a,e,d),c in poly(hc['monic_H_laurent']).items()}
    sc=read(reports/(PREFIX+'35_MINORS'+END)/'certificate.json');Pi=list(map(poly,sc['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']));U=poly(sc['FULLU_product']);assert len(Pi)==22
    cpp=ROOT/'code/s61_nonsquare_Q_degree6_field_unit_review_20261104_v63.cpp'
    for p in [Path(__file__),cpp,ROOT/'code/s61_nonsquare_finite_exception_review_20261104_v63.py',ROOT/'code/s61_nonsquare_raw_H_B_curve_components_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    inp=work/'fixed-quotient-input.txt'
    with inp.open('w')as f:
        def pw(p):
            p=ep(p)if isinstance(p,list)else p;f.write(str(len(p))+''.join(f' {e} {c}'for e,c in sorted(p.items()))+'\n')
        def sw(p):f.write(str(len(p))+''.join(f' {a} {e} {d} {c}'for (a,e,d),c in sorted(p.items()))+'\n')
        def bw(p):
            for q in p:pw(q)
        def aw(p):
            for q in p:bw(q)
        pw(field)
        for p in[Q,H,U]+Pi:sw(p)
        pw(n['q0_reduced']);pw(n['q1_reduced']);bw(n['a_inverse'])
        for p in n['Hclear_coefficients']:bw(p)
        for p in n['H_monic_coefficients']:bw(p)
        assert len(n['all_field_inverse_certificates'])==9;f.write('9\n')
        for inv in n['all_field_inverse_certificates']:
            for key in ['element','inverse','Bezout_modulus','Bezout_element']:pw(inv[key])
        for p in n['all22_original_source_images']:aw(p)
        aw(n['U_image']);aw(n['target_U_power8']);assert len(n['matrix_M'])==8 and all(len(row)==176 for row in n['matrix_M'])
        for row in n['matrix_M']:
            for p in row:pw(p)
        assert len(n['source22_cofactor_witness'])==22
        for p in n['source22_cofactor_witness']:aw(p)
    assert inp.stat().st_size<16*1024*1024
    exe=work/'verify-native'
    with(out/'compile.log').open('w')as log:subprocess.run(['c++','-O2','-std=c++17',str(cpp),'-o',str(exe),*subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
    with inp.open()as f:proc=subprocess.run([str(exe)],stdin=f,text=True,capture_output=True,check=True,timeout=max(1,55-(time.monotonic()-start)))
    assert not proc.stderr;result=json.loads(proc.stdout);assert result['all_passed'];(out/'native.log').write_text(proc.stdout)
    proof=dict(selected_factor_index=6,field_degree6=True,factor_e=[[e,str(c)]for e,c in sorted(field.items())],
        source_Q_H_all22_U_bound=True,all9_full_QQe_Bezout_units=True,all176_columns_8448QQslots=True,
        original22_cofactor_dot_equals_U8=True,whole_closed_component_FULLU_empty=True,
        all_conjugates_roots_and_nilpotents_covered=True,no_other_factor_or_old_residual_claim=True,
        no_author_worker_solver_inverse_GCD_Gaussian_points_primes_or_GB=True,checker=result)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report=dict(id=ID,all_passed=True,whole_selected_degree6_factor6_FULLU_empty=True,input_sha256=hashes,
        evidence_sha256={p.name:sha(p)for p in out.iterdir()},ignored_verifier_input_sha256=sha(inp),ignored_executable_sha256=sha(exe),
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    p=out/'report.json';p.write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(p))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
