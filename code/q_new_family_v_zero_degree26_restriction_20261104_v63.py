#!/usr/bin/env python3
"""One capped exactQQ, point-free full45 restriction on the dense v0 source."""
import argparse,gzip,hashlib,json,resource,shlex,subprocess,time
from collections import Counter
from pathlib import Path
import sympy as s
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63'
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    ip=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    sp=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    pp=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_SUPPORT_PREFLIGHT_20261104_v63/certificate.json'
    bp=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/native-certificate.json.gz'
    dp=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/definitions.json'
    gates=[ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_SUPPORT_PREFLIGHT_20261104_v63/report.json',ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_SOURCE_RELATION_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed']for p in gates)
    for p in gates:
        z=read(p)
        for n,h in z.get('input_sha256',{}).items():assert sha(ROOT/n)==h,(p,n)
    info,spin,pre,basis,defs=map(read,[ip,sp,pp,bp,dp]);assert pre['Inv26_support_bound_per_column']==2874 and len(defs['Inv26_basis'])==45
    Omega=s.Matrix(info['symplectic_form16']);Oi=Omega.inv();rho=[s.Matrix(x)for x in info['all_ten_representation_matrices']]
    J=s.diag(s.Matrix([[0,1],[-1,0]]),s.Matrix([[0,1],[-1,0]]));gs=[]
    for i in range(4):
        for j in range(i,4):
            A=s.zeros(4);A[i,j]=A[j,i]=1;gs.append(J*A)
    Gi=s.Matrix([[s.trace(A*B)for B in gs]for A in gs]).inv();mon=list(map(tuple,spin['quadratic_monomials']));mi={m:i for i,m in enumerate(mon)}
    mu=[]
    for A in rho:
        S=Omega*A;assert S==S.T;mu.append([S[i,j]/(2 if i==j else 1)for i,j in mon])
    rflat=s.Matrix([[A[i,j]for A in rho]for i in range(16)for j in range(16)]);Theta=2*rflat*Gi*s.Matrix(mu)
    hess=s.zeros(256,136)
    for indices,c in spin['quartic_trace_moment_square']:
        counts=Counter(indices);c=s.Rational(c)/2
        for i,ei in counts.items():
            for j,ej in counts.items():
                a=ei*(ej-(i==j))
                if a:
                    rest=list(indices);rest.remove(i);rest.remove(j);hess[16*i+j,mi[tuple(sorted(rest))]]+=c*a
    Hess=s.zeros(256,136)
    for i in range(16):
        for k in range(16):
            if Oi[i,k]:
                for j in range(16):Hess[16*i+j,:]+=2*Oi[i,k]*hess[16*k+j,:]
    assert all(c.is_Integer for A in[Theta,Hess]for c in A)
    selected=basis['degree26_selected_columns'];assert len(selected)==45
    inputfile=work/'input.txt'
    with inputfile.open('w')as f:
        f.write('16 136\n')
        for p in pre['beta16_sparse']:
            f.write(str(len(p))+'\n')
            for m,c in p:f.write(c+' '+' '.join(map(str,m))+'\n')
        for m in mon:f.write(' '.join(map(str,m))+'\n')
        for key in['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
            for row in spin[key]:
                cc=[4*s.Rational(c)for c in row];assert all(c.is_Integer for c in cc);f.write(' '.join(map(str,cc))+'\n')
        for A in[Theta,Hess]:
            for row in A.tolist():f.write(' '.join(map(str,row))+'\n')
        f.write(' '.join(map(str,Oi))+'\n');f.write(' '.join(map(str,selected))+'\n')
    with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inputfile.read_bytes())
    recipe=dict(source_preflight_sha256=sha(pp),variables=['H','L','s0','q0'],original_source_variables=['xi0','r0','s0','q0','B0'],
        original_scalar_xi0_power=13,descent_Ha_Lb_to_r0_B0=['(2a+3b)/4','b/2'],
        definitions=defs,selected45_universal_candidate_indices=selected,
        integer_input_normalizations=dict(CD='4*C,4*D',theta_Hessian='2*theta,2*H16',Omega_inverse='literalQQ'),
        all_ambient16_derivatives_before_substitution=True,Poisson_convention='grad(f)^T Omega^-1 grad(g)',
        low_values19_order=['Inv14[0..2]','primitive18[0..5]','F4','Inv12[0..8]'],
        coefficient_serialization='content QQ times primitive integer polynomial; exponent order H,L,s0,q0; no variable denominators',
        all45_source_columns_exactQQ_not_samples=True,GF_prime=32003,
        caps=dict(seconds_total=60,bytes_total=32*1024*1024,memory_bytes=2_000_000_000,intermediate_terms_per_entry=150000,output_terms_per_column=100000),
        no_new_point_prime_basis_growth_or_higher_degree=True)
    (out/'recipe.json').write_text(json.dumps(recipe,indent=2)+'\n')
    native=ROOT/'code/q_new_family_v_zero_degree26_restriction_20261104_v63.cpp';binary=work/'restriction-native'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
    subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(1,60-(time.monotonic()-start)))
    error=None;data=None
    try:
        with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inputfile),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,60-(time.monotonic()-start)))
        data=read(out/'native-certificate.json');assert data['all_passed']and data['exactQQ_columns_saved']==45
        assert data['all32_quadratic_and_all8_trace_ambient_Euler']and data['all45_finite_parity_and_support']
        assert data['selected_universal_candidates']==selected
        assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<=32*1024*1024
    except Exception as e:error=repr(e)
    paths=[Path(__file__),native,ROOT/f'research/notes/{ID}.md',ip,sp,pp,bp,dp,*gates]
    evidence={p.name:sha(p)for p in out.iterdir()if p.is_file()}
    report=dict(id=ID,all_passed=error is None,status='all45_exactQQ_and_GFcoefficient_rank'if error is None else'bounded_frontier',error=error,
        rank_mod32003=data['rank']if data else None,minor_mod32003=data['coefficient_minor_mod32003']if data else None,
        one_actual_v0_section_restriction_only=True,exactQQ_columns_saved=45 if data else None,
        no_new_point_prime_degree_growth_or_QQ_kernel_inference=True,
        inference='rank45 certifies QQ injectivity on complete Inv26 after independent coefficient review; rank<45 is modular evidence only',
        prior_lower26_unchanged_until_review=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256=evidence,seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],rank=report['rank_mod32003'],minor=report['minor_mod32003'],seconds=report['seconds'],error=error)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
