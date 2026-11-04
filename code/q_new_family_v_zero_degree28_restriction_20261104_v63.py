#!/usr/bin/env python3
"""One capped exactQQ full108 source restriction; no source points."""
import argparse,gzip,json,resource,shlex,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_V_ZERO_DEGREE28_RESTRICTION_20261104_v63'
def selected_minor(data):
    M=[r[:]for r in data['selected_coefficient_matrix']];n=len(M);cols=[];det=1;p=32003
    for j in range(108):
        k=next((i for i in range(len(cols),n)if M[i][j]),None)
        if k is None:continue
        row=len(cols)
        if k!=row:M[k],M[row]=M[row],M[k];det=-det%p
        c=M[row][j];det=det*c%p;iv=pow(c,-1,p);M[row]=[a*iv%p for a in M[row]]
        for i in range(row+1,n):
            c=M[i][j];M[i]=[(a-c*b)%p for a,b in zip(M[i],M[row])]
        cols.append(j)
    assert len(cols)==data['rank']and det
    return dict(prime=p,rank=len(cols),pivot_columns=cols,free_columns=sorted(set(range(108))-set(cols)),minor=det,
        coefficient_rows=data['coefficient_monomials_r0_B0_s0_q0'],matrix=[[r[j]for j in cols]for r in data['selected_coefficient_matrix']])
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    old=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63';univ=ROOT/'reports/Q_NEW_FAMILY_DEGREE28_PRODUCT_BASIS_20261104_v63'
    gatepaths=[ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_DEGREE28_PRODUCT_BASIS_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed']for p in gatepaths)
    for p in[old/'report.json',univ/'report.json']:
        for n,h in read(p)['input_sha256'].items():assert sha(ROOT/n)==h,n
    full24=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63/native-certificate.json.gz';s24=read(full24)['degree24_selected_columns'];ub=read(univ/'native-certificate.json.gz');s28=ub['selected_columns'];assert ub['rank']==108
    inputfile=work/'input.txt'
    with inputfile.open('w')as f:
        with gzip.open(old/'native-input.txt.gz','rt')as g:f.write(g.read())
        for name,count in[('low-values19.json',19),('ambient-gradients128.json',128),('trace-values8.json',8)]:
            ps=read(old/name);assert len(ps)==count;f.write(str(count)+'\n')
            for p in ps:
                f.write(str(len(p['terms']))+'\n');c=Q(p['content'])
                for t in p['terms']:f.write(str(c*int(t[0]))+' '+' '.join(map(str,t[1:]))+'\n')
        f.write(' '.join(map(str,s24))+'\n'+' '.join(map(str,s28))+'\n')
    with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inputfile.read_bytes())
    defs=read(univ/'definitions.json');recipe=dict(universal_basis=defs,source_polynomial_variables=['H','L','s0','q0'],source_original_variables=['xi0','r0','s0','q0','B0'],
        source_generators14_order=['F4','p0','p1','p2','p3','p4','N12','CD000101','CD000111','theta5H','theta7H','Inv14_0','Inv14_1','Inv14_2'],
        original_scalar_xi0_power=14,finite_descent=['r0_exponent=(2Hexp+3Lexp)/4','B0_exponent=Lexp/2'],
        frozen128_gradients_are_ambient16_before_substitution=True,new_low_value_only='theta7H16 via two extra16x16 transports',
        all108_columns_exactQQ=True,GF32003_coefficient_rank_only=True,no_new_point_prime_higher_degree_or_candidate_growth=True)
    (out/'recipe.json').write_text(json.dumps(recipe,indent=2)+'\n')
    native=ROOT/'code/q_new_family_v_zero_degree28_restriction_20261104_v63.cpp';helper=ROOT/'code/q_new_family_v_zero_degree26_restriction_20261104_v63.cpp';binary=work/'source28-native'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(1,60-(time.monotonic()-start)))
    data=None;error=None
    try:
        with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inputfile),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,60-(time.monotonic()-start)))
        data=read(out/'native-certificate.json');assert data['all_passed']and data['exactQQ_columns_saved']==108 and data['selected_universal_candidates']==s28
        assert data['ambient128_gradient_Euler_and14_brackets_verified']and data['physical_theta5H_matches_saved']and data['physical_theta7H_new']and data['all108_parity_and_support']
        (out/'modular-selected-minor.json').write_text(json.dumps(selected_minor(data),indent=2)+'\n')
        assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    except Exception as e:error=repr(e)
    paths=[Path(__file__),native,helper,ROOT/f'research/notes/{ID}.md',old/'report.json',old/'native-input.txt.gz',old/'low-values19.json',old/'ambient-gradients128.json',old/'trace-values8.json',univ/'report.json',univ/'native-certificate.json.gz',univ/'definitions.json',full24,*gatepaths]
    report=dict(id=ID,all_passed=error is None,error=error,status='all108_exactQQ_and_GF_coefficient_rank'if data else'bounded_frontier',rank_mod32003=data['rank']if data else None,minor_mod32003=data['minor_mod32003']if data else None,
        exactQQ_columns_saved=108 if data else None,one_dense_v0_source_restriction_only=True,
        no_new_point_prime_higher_degree_or_QQ_kernel_inference=True,
        inference='rank108 gives QQ injectivity of completeInv28 after peer review; rankdefect alone is not QQ identity',
        prior_lower28_unchanged_until_review=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],rank=report['rank_mod32003'],minor=report['minor_mod32003'],seconds=report['seconds'],error=error)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
