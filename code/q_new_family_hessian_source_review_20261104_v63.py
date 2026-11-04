#!/usr/bin/env python3
"""One bounded independent exact source4 audit via the full quartic Hessian."""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
import gzip,hashlib,json,resource,shlex,signal,subprocess,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW'
NAMES=['CD000101','CD000111','theta5H','theta7H']
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def main(out,work):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63';ar=read(folder/'report.json');assert ar['all_passed']
    source_path=ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/source-five-circuit.json.gz'
    spin_path=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    info_path=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    source,spin,info=map(read,[source_path,spin_path,info_path]);assert len(source['nodes'])==1032
    gates=[ROOT/'reports/Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW/report.json',
        ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW/report.json',
        ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed']for p in gates)
    R,a,X,W,B=ring('a,X,W,B',s.QQ);K=R.to_field();H=147456*X*X-(3072*W-1)**2
    replacement=[K(a*X),K.one,K(X),K(W),K(B)];cache={}
    def value(i):
        if i in cache:return cache[i]
        op,arg=source['nodes'][i]
        if op=='q':z=K(s.Rational(arg))
        elif op=='v':z=replacement[arg]
        elif op=='inv':z=1/value(arg)
        elif op=='+':z=sum((value(j)for j in arg),K.zero)
        else:
            assert op=='*';z=K.one
            for j in arg:z*=value(j)
        cache[i]=z;return z
    Omega=s.Matrix(info['symplectic_form16']);Oi=Omega.inv();rho=[s.Matrix(r)for r in info['all_ten_representation_matrices']]
    J=s.diag(s.Matrix([[0,1],[-1,0]]),s.Matrix([[0,1],[-1,0]]));gs=[]
    for i in range(4):
        for j in range(i,4):
            S=s.zeros(4);S[i,j]=S[j,i]=1;gs.append(J*S)
    Gi=s.Matrix([[s.trace(A*B)for B in gs]for A in gs]).inv();mon=[tuple(x)for x in spin['quadratic_monomials']];mi={v:i for i,v in enumerate(mon)}
    mu=[]
    for A in rho:
        S=Omega*A;assert S==S.T
        mu.append([S[i,j]/(2 if i==j else 1)for i,j in mon])
    rflat=s.Matrix([[A[i,j]for A in rho]for i in range(16)for j in range(16)])
    Theta=2*rflat*Gi*s.Matrix(mu)
    assert all(x.is_Integer for x in Theta)
    hess=s.zeros(256,136)
    for indices,coef in spin['quartic_trace_moment_square']:
        counts=Counter(indices);coef=s.Rational(coef)/2
        for i,ei in counts.items():
            for j,ej in counts.items():
                mult=ei*(ej-(i==j))
                if mult:
                    rest=list(indices);rest.remove(i);rest.remove(j);hess[16*i+j,mi[tuple(sorted(rest))]]+=coef*mult
    Hess=s.zeros(256,136)
    for i in range(16):
        for k in range(16):
            if Oi[i,k]:
                for j in range(16):Hess[16*i+j,:]+=2*Oi[i,k]*hess[16*k+j,:]
    assert all(x.is_Integer for x in Hess)
    # Original cocharacter, not a source-specific scaling guess.
    h2=s.diag(0,0,1,-1);sym=J.inv()*h2;co=[sym[i,j]for i in range(4)for j in range(i,4)]
    hrep=sum((c*A for c,A in zip(co,rho)),s.zeros(16));weights=[w[1]for w in info['Cartan_weights']]
    assert hrep==s.diag(*weights)
    for tensor in [Theta,Hess]:
        for row in range(256):
            i,j=divmod(row,16)
            assert all(not tensor[row,k]or weights[u]+weights[v]==weights[i]-weights[j]for k,(u,v)in enumerate(mon))
    # The character identities used by the author's alternative shortcut are
    # checked independently, although our worker does not use that shortcut.
    x,y=s.symbols('x y');s1=x*x+y*y;s2=x*x*y*y
    assert s.expand(sum((w[0]*x+w[1]*y)**6 for w in info['Cartan_weights'])-264*s1**3-408*s1*s2)==0
    assert s.expand(sum((w[0]*x+w[1]*y)**8 for w in info['Cartan_weights'])-1032*s1**4-3488*s1*s1*s2+4208*s2*s2)==0
    beta=[]
    for i,node in enumerate(source['outputs']['beta16']):
        p=value(node)*K(1080*B*H);assert p.denom==1
        terms=[]
        for (ea,ex,ew,eb),c in p.numer.items():
            assert s.denom(s.Rational(c))==1
            et=2*eb+weights[i]-1;assert et in [0,1]
            terms.append([int(c),ea,ex,ew,et])
        beta.append(terms)
    # Bind every author's transformed beta and integral tensor constant while
    # preserving the independent quartic-Hessian entrance to the source words.
    with gzip.open(folder/'native-input.txt.gz','rt')as f:tokens=iter(f.read().split())
    assert [int(next(tokens))for _ in range(4)]==[4,16,136,1080]
    for want in beta:
        n=int(next(tokens));got={}
        for _ in range(n):
            c,eu,ex,ew,et=[int(next(tokens))for _ in range(5)];got[eu,eu+ex,ew,et]=c
        assert got=={tuple(t[1:]):t[0]for t in want}
    assert [[int(next(tokens)),int(next(tokens))]for _ in range(136)]==[list(x)for x in mon]
    CD=[s.Matrix(spin[k])*4 for k in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']]
    assert all(x.is_Integer for A in CD for x in A)
    assert [[int(next(tokens))for _ in range(136)]for _ in range(32)]==CD[0].tolist()+CD[1].tolist()
    assert [int(next(tokens))for _ in range(256)]==[int(x)for x in Omega]
    assert [int(next(tokens))for _ in range(100)]==[int(2*x)for x in Gi]
    assert [[int(next(tokens))for _ in range(256)]for _ in range(10)]==[[int(x)for x in A]for A in rho]
    assert next(tokens,None)is None
    inp=out/'own-input.txt'
    with inp.open('w')as f:
        f.write('16 136\n')
        for terms in beta:
            f.write(str(len(terms))+'\n')
            for row in terms:f.write(' '.join(map(str,row))+'\n')
        for ij in mon:f.write(' '.join(map(str,ij))+'\n')
        for A in [*CD,Theta,Hess]:
            for row in A.tolist():f.write(' '.join(str(int(c))for c in row)+'\n')
        counts=[]
        for name,degree,factor in zip(NAMES,[12,12,12,16],[4096,4096,64,256]):
            p=read(folder/(name+'.json.gz'));assert p['source_c_power']==degree//2 and p['denominator_D_power']==degree
            assert p['T_squared_equals_B_preconditioned'];counts.append(len(p['terms']))
            f.write(f"{name} {degree} {factor} {p['denominator_B_power']} {p['denominator_integer']}\n")
            f.write(str(len(p['terms']))+'\n')
            for c,ea,ey,ew,eb in p['terms']:f.write(f'{c} {ea} {2*ey} {ew} {2*eb}\n')
    assert counts==[3449,3618,3948,10819]
    native=ROOT/'code/q_new_family_hessian_source_review_20261104_v63.cpp';binary=work/'direct-hessian-coefficients'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
    with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=20)
    try:
        with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(out/'own-result.json')],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=50)
    except Exception as exc:
        (out/'blocked.json').write_text(json.dumps(dict(error=repr(exc),input_sha256=sha(inp),one_direct_worker_only=True,no_budget_growth=True),indent=2)+'\n');raise
    result=read(out/'own-result.json');assert result['all_passed']and result['all4_exact_coefficient_maps']
    paths=[Path(__file__),native,ROOT/f'research/notes/{STEM}.md',source_path,spin_path,info_path,*gates,
        ROOT/'code/q_new_family_hessian_low_invariant_expand_20261104_v63.py',
        ROOT/'code/q_new_family_hessian_low_invariant_expand_20261104_v63.cpp',
        ROOT/'research/notes/Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63.md',
        ROOT/'reports/Q_NEW_FAMILY_HESSIAN_SOURCE_VECTOR_PREFLIGHT_20261104_v63/report.json',
        *[folder/p for p in ['report.json','recipe.json','native-input.txt.gz',*[n+'.json.gz'for n in NAMES]]]]
    report=dict(id=STEM,all_passed=True,actual1032_source_to_all16_preprocessed_beta_polynomials=True,
        exact_original_global_torus_cocharacter_and_all_theta_H_quadratic_weights=True,
        independent_integral_2theta_and_2H16_from_full_quartic_Hessian=True,
        independent_Cartan_character_pure_trace6_and8_identities=True,
        direct16_Hessian_column_recurrence_not_author_rank10_or_character_shortcut=True,
        all4_complete_QQ_coefficient_maps_and_exact_word_normalizations_verified=True,
        output_terms=counts,total_verified_coefficients=sum(counts),one_direct_native_worker=True,
        no_NEW_point_source67_relation_test_or_high_degree_growth=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k]for k in ['all_passed','total_verified_coefficients','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output.resolve(),a.work.resolve())
