#!/usr/bin/env python3
"""Independent whole-curve source identities; no source solve or NF search."""
import argparse,gzip,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_H1_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def native_matrix(p):
    lines=p.read_text().splitlines();n,m=map(int,lines[0].split());assert len(lines)==n*m+1;result=[]
    for line in lines[1:]:
        s=line.split();k=int(s[0]);assert len(s)==1+4*k
        poly={tuple(map(int,s[j+1:j+4])):F(s[j])for j in range(1,len(s),4)}
        assert len(poly)==k and all(min(e)>=0 and c for e,c in poly.items());result.append(poly)
    return n,m,result
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);work=ROOT/'work/P100-H1-independent-review';work.mkdir(parents=True,exist_ok=True);hashes={};base=ROOT/'reports'
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def audit(p,passed=False):
        r=read(p)
        if passed:assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():assert sha(p.parent/name)==h,name;hashes[str((p.parent/name).relative_to(ROOT))]=h
        for key,name in [('native_input_sha256','native-input.txt'),('native_output_sha256','native.json'),('certificate_sha256','certificate.json')]:
            if key in r and not passed:assert sha(p.parent/name)==r[key];hashes[str((p.parent/name).relative_to(ROOT))]=r[key]
        return r
    dual=base/'S61_NONSQUARE_P100_H1_DUAL_20261104_v63';right=base/'S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63'
    dr=audit(dual/'report.json');rr=audit(right/'report.json')
    assert dr['actual_rawH_dual_all16_verified']and rr['whole_P100_H1projection_surjective']
    for stem in ['S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63_REVIEW','S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63_REVIEW',
      'S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63_REVIEW','S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63_REVIEW',
      'S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63_REVIEW','S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63_REVIEW',
      'S61_NONSQUARE_H1_KNOWN_UNIT_ELIMINATION_20261104_v63_REVIEW']:audit(base/stem/'report.json',True)
    clarification=ROOT/'research/notes/S61_NONSQUARE_P100_H1_DUAL_ROW_INDEX_CLARIFICATION_20261104_v63.md';hashes[str(clarification.relative_to(ROOT))]=sha(clarification)
    for stem in ['S61_NONSQUARE_P100_H1_DUAL_20261104_v63','S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63']:
        p=ROOT/f'research/notes/{stem}.md';hashes[str(p.relative_to(ROOT))]=sha(p)
    D=read(dual/'native.json');R=read(right/'native.json');cert=read(right/'certificate.json')
    assert D['optional_branch']=='stopped_NF117_H5_descent_failure'and not D['H23_dual']and not D['optional_scalar_nonzero']
    raw=read(base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/preflight.json')
    h=read(base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json')
    cdir=base/'S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63';jdir=base/'S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63'
    for p in[cdir/'six-original-quadrics.txt',cdir/'source-metadata.txt',jdir/'coordinate-forward7x7.txt',jdir/'monomial-input.txt',jdir/'six-prime-quadrics-numerators.txt',jdir/'six-prime-quadrics-denominators.txt']:hashes[str(p.relative_to(ROOT))]=sha(p)
    n,m,qold=native_matrix(cdir/'six-original-quadrics.txt');assert(n,m)==(6,29)
    n,m,forward=native_matrix(jdir/'coordinate-forward7x7.txt');assert(n,m)==(7,7)
    n,m,prime=native_matrix(jdir/'six-prime-quadrics-numerators.txt');assert(n,m)==(6,29)
    n,m,den=native_matrix(jdir/'six-prime-quadrics-denominators.txt');assert(n,m)==(6,1)
    oldmeta=(cdir/'source-metadata.txt').read_text().splitlines();oldmons=[list(map(int,line.split()))for line in oldmeta[1:30]]
    newmeta=(jdir/'monomial-input.txt').read_text().splitlines();newmons=[list(map(int,line.split()))for line in newmeta[30:59]]
    assert len(oldmons)==len(newmons)==29
    hmons=[h['parent_joint_target239_monomials'][h['parent_target122_monomial_indices'][j]]for j in h['target122_Hdegree_classes']['1']];assert len(hmons)==23
    defs=h['all112_cubic_source319_definitions'];physical=[0,1,2,3,4,6,7]
    assert len(defs)==112 and [defs[7*i][2]for i in [1,6]]==[20,61]
    for i in range(16):assert[defs[7*i+j][3]for j in range(7)]==physical
    for i,rec in enumerate(cert['source_cubic_cofactor_rows']):
        assert rec['matrix_row']==i and rec['actual134_cache_row']==i+1 and rec['genuine_index']==i//7
        assert rec['physical_light']==physical[i%7]and rec['frozen_definition']==defs[i]
    pc=cert['actual_delta7_preconditioner'];assert pc['replace_rows']==list(range(43,50))and pc['formula']=='row(43+j)=delta*(oldrow(43+j)-oldrow(8+j))'
    expected=[({0:F(95,2)},{0:F(20),1:F(1)}),({0:F(-2865),1:F(-179,4)},{0:F(1200),1:F(80),2:F(1)}),
      ({0:F(-425,6)},{1:F(20),2:F(1)}),({0:F(1)},{0:F(1)}),({0:F(40)},{0:F(20),1:F(1)}),
      ({0:F(-6615),1:F(-429,4)},{0:F(1200),1:F(80),2:F(1)}),({0:F(2125,12)},{1:F(20),2:F(1)})]
    def uni(p):return{int(e[0]):F(c)for e,c in p if F(c)}
    assert[(uni(r['numerator']),uni(r['denominator']))for r in D['eta_physical7']]==expected
    assert len(R['right_inverse112x23'])==112 and all(len(row)==23 for row in R['right_inverse112x23'])
    assert len(R['linear_lift_coefficients7x16'])==7 and len(R['fixed_NF117_H5_normalized_rows'])==5
    for p in[Path(__file__),Path(__file__).with_suffix('.cpp'),ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    inp=work/'fixed-source-identities.txt'
    with inp.open('w')as f:
        def pw(p):f.write(str(len(p))+''.join(f' {e} {c}'for e,c in sorted(p.items()))+'\n')
        def sw(p):
            p={tuple(e):F(c)for e,c in p}if isinstance(p,list)else p
            f.write(str(len(p))+''.join(f' {a} {e} {d} {c}'for(a,e,d),c in sorted(p.items()))+'\n')
        def rw(r):pw(uni(r['numerator']));pw(uni(r['denominator']))
        def mw(ms):
            for m in ms:f.write(' '.join(map(str,m))+'\n')
        for c in [F(0),F(20),F(60),F(2700,203),F(12)]:pw({0:c,1:F(1)}if c else{1:F(1)})
        for row in raw['raw16x7']:
            for p in row:sw(p)
        mw(oldmons)
        for p in qold:sw(p)
        for p in forward:sw(p)
        mw(newmons)
        for p in prime:sw(p)
        for p in den:sw(p)
        mw(raw['constant_coefficient_column_frame7x7']);mw(raw['constant_inverse7x7']);mw(hmons)
        for r in D['eta_physical7']:rw(r)
        for r in D['quadratic_evaluations']:rw(r)
        for row in R['linear_lift_coefficients7x16']:
            for r in row:rw(r)
        for row in R['fixed_NF117_H5_normalized_rows']:
            for r in row:rw(r)
        for row in R['right_inverse112x23']:
            for r in row:rw(r)
        for r in R['q1_cofactors23']:rw(r)
        sw(raw['source_delta_change']['delta'])
    assert inp.stat().st_size<8*1024*1024
    exe=work/'verify-native'
    with(out/'compile.log').open('w')as log:subprocess.run(['c++','-O2','-std=c++17',str(Path(__file__).with_suffix('.cpp')),'-o',str(exe),*subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()],check=True,stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
    with inp.open()as f:proc=subprocess.run([str(exe)],stdin=f,text=True,capture_output=True,check=True,timeout=max(1,55-(time.monotonic()-start)))
    (out/'native.log').write_text(proc.stdout);assert not proc.stderr;result=json.loads(proc.stdout);assert result['all_passed']
    proof=dict(actual_raw16_eta_annihilations_and_W1=True,all6_original_quadratic_evaluations_and_q1_unit=True,
      square_eta_does_not_descend_NF23=True,all7_explicit_source_lifts_of_v_minus_etaW=True,
      raw_source_module_equals_kernel_eta_over_whole_curve=True,
      constant_frame_inverse_and_joint_lambda_orientation=True,all5_corrected_qNF_rows_and_original_source_transforms=True,
      all23_preNF_full_quadratic_identities_in_old_and_actual_delta7_frames=True,
      actual_cache_rows1_to112_and_preconditioner43_to49_offset_bound=True,
      old_unexecuted_optional_row112_error_explicitly_not_used=True,
      H1_projection23_globally_split_surjective_on_P100=True,
      accepted117_and_H2_unit_gates_give_pointwise_original_quartic_rank_atleast141=True,
      no_fixed_original141_minor_or_free_source_complement_claim=True,
      pure98_and_all21_quartics_retained_without_query=True,
      no_H_times_linear_form_or_H_W2_zero_in_full_ideal_inferred=True,
      no_full134Gaussian_POWER_newpoint_prime_GB_full57_full239_or_exhaustion=True,checker=result)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    r=dict(id=ID,all_passed=True,H1_projection23_surjective_whole_P100=True,pointwise_original_quartic_rank_lower141=True,
      dual_fails_NF23_descent=True,input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},
      ignored_verifier_input_sha256=sha(inp),ignored_executable_sha256=sha(exe),seconds=time.monotonic()-start,
      max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'],report_sha256=sha(p))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
