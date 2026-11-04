#!/usr/bin/env python3
"""Fixed independent P80 data substitution and all H1 right-inverse identities."""
import argparse,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
from s61_nonsquare_P100_H1_review_20261104_v63 import native_matrix
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P80_H1_RESCUE_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);work=ROOT/'work/P80-H1-independent-review';work.mkdir(parents=True,exist_ok=True);hashes={};base=ROOT/'reports'
    def read(p):hashes[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
    def audit(p,passed=False):
        r=read(p)
        if passed:assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for field in ['evidence_sha256','native_evidence_sha256']:
            for name,h in r.get(field,{}).items():assert sha(p.parent/name)==h,name;hashes[str((p.parent/name).relative_to(ROOT))]=h
        if not passed and 'certificate_sha256'in r:assert sha(p.parent/'certificate.json')==r['certificate_sha256'];hashes[str((p.parent/'certificate.json').relative_to(ROOT))]=r['certificate_sha256']
        return r
    b=base/'S61_NONSQUARE_P80_H1_RESCUE_20261104_v63';primary=audit(b/'report.json')
    assert primary['wholeP80_H1_rightinverse112x23_verified']and primary['chosen_eliminated_unit_quadric']==1
    assert primary['native_stored_terms']==17982
    audit(base/'S61_NONSQUARE_P100_H1_20261104_v63_REVIEW/report.json',True)
    audit(base/'S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63_REVIEW/report.json',True)
    D=read(b/'dual.json');R=read(b/'right-inverse.json');cert=read(b/'certificate.json')
    assert D['curve']=='P80'and D['no_cache_row_query']and not D['all5_eliminated_Hquad_directions_annihilated']
    raw=read(base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/preflight.json')
    h=read(base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json')
    cdir=base/'S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63';jdir=base/'S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63'
    for p in[cdir/'six-original-quadrics.txt',cdir/'source-metadata.txt',jdir/'coordinate-forward7x7.txt',jdir/'monomial-input.txt',jdir/'six-prime-quadrics-numerators.txt',jdir/'six-prime-quadrics-denominators.txt']:hashes[str(p.relative_to(ROOT))]=sha(p)
    n,m,qold=native_matrix(cdir/'six-original-quadrics.txt');assert(n,m)==(6,29)
    n,m,forward=native_matrix(jdir/'coordinate-forward7x7.txt');assert(n,m)==(7,7)
    n,m,prime=native_matrix(jdir/'six-prime-quadrics-numerators.txt');assert(n,m)==(6,29)
    n,m,den=native_matrix(jdir/'six-prime-quadrics-denominators.txt');assert(n,m)==(6,1)
    oldmons=[list(map(int,line.split()))for line in(cdir/'source-metadata.txt').read_text().splitlines()[1:30]]
    newmons=[list(map(int,line.split()))for line in(jdir/'monomial-input.txt').read_text().splitlines()[30:59]]
    hmons=[h['parent_joint_target239_monomials'][h['parent_target122_monomial_indices'][j]]for j in h['target122_Hdegree_classes']['1']];assert len(hmons)==23
    defs=h['all112_cubic_source319_definitions'];assert defs==cert['source_definitions112']and len(defs)==112
    assert [defs[7*i][2]for i in[1,6]]==[20,61]
    assert all([defs[7*i+j][3]for j in range(7)]==[0,1,2,3,4,6,7]for i in range(16))
    pc=cert['source_preconditioner'];assert pc['replace_rows']==list(range(43,50))and pc['formula']=='row(43+j)=delta*(oldrow(43+j)-oldrow(8+j))'
    assert pc['delta'][1]==[[[0,0,0],'1']]and {tuple(e):F(c)for e,c in pc['delta'][0]}=={tuple(e):F(c)for e,c in raw['source_delta_change']['delta']}
    expected=[({0:F(220,7)},{0:F(80,7),1:F(1)}),({0:F(-1500),1:F(-3265,98)},{0:F(24000,49),1:F(380,7),2:F(1)}),
      ({0:F(-1000,21)},{1:F(80,7),2:F(1)}),({0:F(1)},{0:F(1)}),({0:F(200,7)},{0:F(80,7),1:F(1)}),
      ({0:F(-166500,49),1:F(-7605,98)},{0:F(24000,49),1:F(380,7),2:F(1)}),({0:F(2500,21)},{1:F(80,7),2:F(1)})]
    def uni(p):return{int(e[0]):F(c)for e,c in p if F(c)}
    assert[(uni(r['numerator']),uni(r['denominator']))for r in D['eta_physical7']]==expected
    assert len(R['right_inverse112x23'])==112 and all(len(row)==23 for row in R['right_inverse112x23'])
    nonzero=sum(bool(r['numerator'])for row in R['right_inverse112x23']for r in row);assert nonzero==904
    positions=[i for i,row in enumerate(R['right_inverse112x23'])if any(r['numerator']for r in row)]
    assert positions==[7*s+j for s in[1,6,8,9,11,14]for j in range(7)]
    cpp=Path(__file__).with_suffix('.cpp')
    for p in[Path(__file__),cpp,ROOT/'code/s61_nonsquare_P100_H1_review_20261104_v63.py',ROOT/'code/s61_nonsquare_P100_H1_review_20261104_v63.cpp',ROOT/'research/notes/S61_NONSQUARE_P80_H1_RESCUE_20261104_v63.md',ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    # Audit exactly which parts of the accepted arithmetic were changed.
    old=(ROOT/'code/s61_nonsquare_P100_H1_review_20261104_v63.cpp').read_text();new=cpp.read_text()
    old=old.replace('// named linear units.','// six named P80 linear units.').replace('std::array<P,5> factors','std::array<P,6> factors').replace('std::array<std::vector<P>,5> powers','std::array<std::vector<P>,6> powers').replace('std::array<int,5>','std::array<int,6>').replace('for(int k=0;k<5;k++)','for(int k=0;k<6;k++)')
    old=old.replace('R a(pm(constant("-10/9"),zp20),{2,0,0,0,0}),e(pm(constant("1/9"),pa(z,constant("100"))),{1,0,0,0,0});','R a(pm(constant("-8/45"),pa(pm(constant("7"),z),constant("80"))),{2,0,0,0,0,0}),e(pm(constant("2/9"),pa(z,constant("40"))),{1,0,0,0,0,0});').replace('constant("-3/100")','constant("-3/80")').replace('constant("137500/243"),{2,0,0,0,0}','constant("41600/81"),{2,0,0,0,0,0}')
    assert old==new,'unexpected arithmetic change beyond P80 parameters and sixth unit factor'
    inp=work/'fixed-source-identities.txt'
    with inp.open('w')as f:
        def pw(p):f.write(str(len(p))+''.join(f' {e} {c}'for e,c in sorted(p.items()))+'\n')
        def sw(p):
            p={tuple(e):F(c)for e,c in p}if isinstance(p,list)else p;f.write(str(len(p))+''.join(f' {a} {e} {d} {c}'for(a,e,d),c in sorted(p.items()))+'\n')
        def rw(r):pw(uni(r['numerator']));pw(uni(r['denominator']))
        def mw(ms):
            for m in ms:f.write(' '.join(map(str,m))+'\n')
        for c in[F(0),F(20),F(300,7),F(80,7),F(216,17),F(60,7)]:pw({0:c,1:F(1)}if c else{1:F(1)})
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
        for key in['linear_lift_coefficients7x16','fixed_NF117_H5_normalized_rows','right_inverse112x23']:
            for row in R[key]:
                for r in row:rw(r)
        for r in R['q1_cofactors23']:rw(r)
        sw(raw['source_delta_change']['delta'])
    assert inp.stat().st_size<8*1024*1024;exe=work/'verify-native'
    with(out/'compile.log').open('w')as log:subprocess.run(['c++','-O2','-std=c++17',str(cpp),'-o',str(exe),*subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
    with inp.open()as f:proc=subprocess.run([str(exe)],stdin=f,text=True,capture_output=True,check=True,timeout=max(1,55-(time.monotonic()-start)))
    (out/'native.log').write_text(proc.stdout);assert not proc.stderr;result=json.loads(proc.stdout);assert result['all_passed']
    proof=dict(curve='P80',actual_P80_parameters_and_six_unit_factors_checked=True,q1_eta='41600/(81*zeta^2)',
      independent_arithmetic_reused_with_exact_change_audit=True,all16_raw_dual_annihilations_and_W1=True,
      all6_original_q_evaluations_and_q1_whole_curve_unit=True,all7_fixed_source_lifts_49slots=True,
      all5_corrected_prime_NF_rows_and_original_source_transforms=True,all1288_source_coefficients_both_frames=True,
      rightinverse_shape=[112,23],rightinverse_nonzero_coefficients=904,actual_source_positions=positions,
      actual_cache_rows1_to112_q5Hrow0_and_delta7_rows43_to49_bound=True,
      H1_projection23_split_surjective_whole_P80=True,original_pointwise_quartic_rank_atleast141=True,
      pure98_all21quartics_retained_no_full134_POWER_Gaussian_sourcequery_points_primes_or_GB=True,
      no_fixed141_minor_free_source_complement_full57_full239_or_exhaustion_claim=True,checker=result)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    r=dict(id=ID,all_passed=True,curve='P80',H1_projection23_surjective_whole_P80=True,pointwise_original_quartic_rank_lower141=True,
      input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},ignored_verifier_input_sha256=sha(inp),ignored_executable_sha256=sha(exe),
      seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'],report_sha256=sha(p))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
