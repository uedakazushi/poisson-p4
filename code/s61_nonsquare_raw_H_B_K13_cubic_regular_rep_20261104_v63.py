#!/usr/bin/env python3
"""ONE positive-only p32003 regular-representation221x312 gate, no beta point."""
import argparse,hashlib,json,math,resource,signal,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_CUBIC_REGULAR_REP_20261104_v63';P=32003
sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    cb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_CUBIC17_COEFFICIENTS_20261104_v63';cp=cb/'report.json'
    assert sha(cp)=='b2f18c8ecbced0a2ad45280084849342a3760ae81d6fbe38d81fbcb07fdf87f3';c=read(cp);assert c['all_passed']
    mp=cb/'matrix17x24.jsonl';assert sha(mp)==c['evidence_sha256'][mp.name]
    dg=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63';g={int(i):int(v)for i,v in read(dg/'checkpoint.json')['g13']}
    rows=[json.loads(s)for s in mp.read_text().splitlines()];assert len(rows)==408
    lc=g[13]%P;denoms=[int(r['value']['denominator'])%P for r in rows]
    gate=dict(prime=P,g13_leading_coefficient_residue=lc,all408_QQ_common_denominator_residues=denoms,
        source_all4602_QQ_coefficient_denominators_covered=True,passed=bool(lc)and all(denoms),
        no_other_prime_point_or_numeric_beta_allowed=True)
    (out/'denominator-punit-gate.json').write_text(json.dumps(gate,indent=2)+'\n')
    status='punit_gate_FAILED_NO_retry';native=None
    if gate['passed']:
        gl=[g.get(i,0)%P*pow(lc,P-2,P)%P for i in range(13)]
        inp=out/'regular-input.txt'
        with inp.open('w')as f:
            f.write('32003 17 24 13\n'+' '.join(map(str,gl))+'\n')
            for index,r in enumerate(rows):
                assert r['index']==index and r['row17']==index//24 and r['column24']==index%24
                num={int(i):int(v)for i,v in r['value']['numerator']};assert all(0<=i<13 for i in num)
                invden=pow(denoms[index],P-2,P);coeff=[num.get(i,0)%P*invden%P for i in range(13)];f.write(' '.join(map(str,coeff))+'\n')
        pre=dict(full_regular_matrix_shape=[221,312],full_mod_slots=221*312,selected_combined_LU_slots=221*221,
            saved_QQ_input_terms=4602,peak_if_full_plus_ONE_LU_plus_QQ_inputs=221*312+221*221+4602,
            actual_native_full_freed_before_LU=True,no_two221square_copies=True,hard_seconds=60,hard_AS_bytes=2000000000,
            hard_JOB_TOTAL_artifact_bytes=32*1024*1024,hard_explicit_terms=150000,no_QQ_Gaussian_or_rank_failure_inference=True)
        assert pre['peak_if_full_plus_ONE_LU_plus_QQ_inputs']<150000
        (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
        cpp=Path(__file__).with_suffix('.cpp');exe=work/'regular-rep'
        with(out/'compile.log').open('w')as log:
            subprocess.run(['c++','-O2','-std=c++17',str(cpp),'-o',str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
        answer=out/'native-certificate.json'
        with(out/'native.log').open('w')as log:
            subprocess.run([str(exe),str(inp),str(answer)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,57-(time.monotonic()-start)))
        native=read(answer);assert native['prime']==P and native['regular_shape']==[221,312]
        if native['rank_modp']==221:
            assert native['actual_selected_minor_determinant']and native['all_selected_PLU_entry_checks']==221*221
            assert len(set(native['selected_columns']))==221
            piv=native['diagonal_U'];det=(-1 if native['row_swaps']%2 else 1)*math.prod(piv)%P;assert det==native['actual_selected_minor_determinant']
            status='positive_rank221_FULL_cubic57_ON_WHOLE_actual_Ka13'
        else:status='deficient_modular_rank_INCONCLUSIVE_no_retry'
    positive=status.startswith('positive_rank221')
    proof=dict(status=status,base_model='V=Z_(32003); T_V=V[a]/monic(g13), finiteFREE13 because lc(g13) is a V-unit.',
        actual_QQ_matrix='All408 canonicalcoefficients are V-integral by the complete denominator gate. C:T_V^24->T_V^17 is represented in basis1,a,...,a12 by the saved regular221x312 matrix.',
        regular_column='Column(j,k) is coefficients of C_ij*a^k modulo monicg13 for all17 outputblocks. This is multiplication, NOT evaluation ofa at a root.',
        positive_lift='Nonzero selected221minor modp means sameQQminor is nonzero. Thus C_QQ:Ka13^24->Ka13^17 is surjective as QQvectorspaces and as Ka13modules. Basechange gives rank17 atEVERY residuefield, without irreducibility or reducedness.',
        alternative_transpose='The source-row map free17->free24 is then splitinjective. No claim from a merely injective generic matrix is used.',
        source_full57_gate='C=cLCM*aB*M is a QQ-Ka13unitmultiple; accepted8unitprefix/product41/Cdelta criterion implies fullactualcubic57 at allpoints ofactualBbad. Prefixneednotbeclaimedunitmodp: the conclusion passes throughQQ.',
        full_cubic57_at_ALL_actual_Bbad_geometric_points_given_accepted_gates=positive,
        no_original_GF_cubic57_or_global_allFULLU57_or_quartic239_claim=True,
        deficient_case='Modular rankdeficiency is INCONCLUSIVE; no sourcecounterexample/QQGaussian/otherprime or automaticnextjob.',
        independent_review_of_new_coefficients_and_rank_certificate_pending=True)
    (out/'characteristic-zero-lift.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),Path(__file__).with_suffix('.cpp'),ROOT/f'research/notes/{ID}.md',cp,mp,dg/'report.json',dg/'checkpoint.json',ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_NONSQUARE_PRODUCT41_FULL_OPEN_20261103_v63/report.json']
    report=dict(id=ID,all_passed=bool(gate['passed']),status=status,prime=P,rank_modp=native['rank_modp']if native else None,
        selected221minor_determinant_modp=native.get('actual_selected_minor_determinant')if native else None,
        full_cubic57_on_ALL_actual_Bbad_given_accepted_source_gates=positive,
        no_beta_point_newprime_QQGaussian_global_allFULLU57_or_quartic_job=True,independent_review_pending=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
    try:main(a.output,a.work)
    except Exception as e:
        if a.output.exists():(a.output/'failure.json').write_text(json.dumps(dict(all_passed=False,error=str(e),no_auto_retry=True),indent=2)+'\n')
        raise
