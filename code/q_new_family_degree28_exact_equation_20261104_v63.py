#!/usr/bin/env python3
"""ONE exact kernel of frozen source28 columns and exhaustive QQ identity."""
import argparse,gzip,json,math,resource,shlex,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_DEGREE28_EXACT_EQUATION_20261104_v63'
def monomial24(z):
    e=[0]*13;k=z['kind']
    if k=='spectral':e[:5]=z['F4_p1_p2_p3_p4_exponents']
    elif k=='primitive12_times_spectral12':
        e[5+z['primitive']]+=1
        if z['spectral']==0:e[0]+=3
        else:e[0]+=1;e[z['spectral']]+=1
    elif k=='primitive12_product':
        for i in z['indices']:e[5+i]+=1
    else:
        assert k=='primitive16_times_Inv8';e[9]+=1
        if z['index']==0:e[0]+=2
        else:e[z['index']]+=1
    return e
def monomial28(z,defs24):
    e=[0]*13;k=z['kind']
    if k=='F4_times_Inv24_basis':e=monomial24(defs24[z['Inv24_basis_index']]);e[0]+=1
    elif k=='primitive12_times_p_p':
        e[5+z['primitive12_index']]+=1
        for j in z['p_indices']:e[j]+=1
    elif k=='primitive16_times_primitive12':e[9]+=1;e[5+z['primitive12_index']]+=1
    else:
        assert k=='Inv14_product'
        for i in z['indices']:e[10+i]+=1
    assert sum(a*b for a,b in zip(e,[4,8,8,8,8,12,12,12,12,16,14,14,14]))==28
    return e
def main(out,work):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE28_RESTRICTION_20261104_v63';univ=ROOT/'reports/Q_NEW_FAMILY_DEGREE28_PRODUCT_BASIS_20261104_v63';sr=read(source/'report.json');assert sr['all_passed']and sr['rank_mod32003']==107
    for n,h in sr['input_sha256'].items():assert sha(ROOT/n)==h,n
    piv=read(source/'modular-selected-minor.json');assert piv['rank']==107 and len(piv['free_columns'])==1;free=piv['free_columns'][0]
    defs=read(univ/'definitions.json');cols=[read(source/f'column{i}.json')['HLsq_polynomial']for i in range(108)]
    inp=work/'input.txt'
    with inp.open('w')as f:
        f.write('108 107\n')
        for p in cols:
            f.write(p['content']+' '+str(len(p['terms']))+'\n')
            for t in p['terms']:f.write(' '.join(map(str,t))+'\n')
        for rr,bb,ss,qq in piv['coefficient_rows']:
            aa=2*rr-3*bb;assert aa>=0;f.write(f'{aa} {2*bb} {ss} {qq}\n')
        f.write(' '.join(map(str,piv['pivot_columns']))+'\n'+str(free)+'\n')
    with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inp.read_bytes())
    native=ROOT/'code/q_new_family_degree28_exact_equation_20261104_v63.cpp';helper=ROOT/'code/q_new_family_v_zero_degree26_restriction_20261104_v63.cpp';binary=work/'exact-kernel-native';raw=work/'native-certificate.json'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(1,60-(time.monotonic()-start)))
    data=None;error=None;identity=False;equation=None
    try:
        with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(raw)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,60-(time.monotonic()-start)))
        assert raw.stat().st_size<32*1024*1024;data=read(raw);assert data['all_passed']and data['all_source_coefficients_zero']and data['exact_QQ_solve107']
        coeff=list(map(int,data['primitive_integer_kernel108']));assert len(coeff)==108 and any(coeff)and math.gcd(*coeff)==1
        # Separate exhaustive Fraction check, using only the frozen QQ columns.
        total={}
        for a,p in zip(coeff,cols):
            c=a*Q(p['content'])
            if c:
                for t in p['terms']:
                    e=tuple(t[1:]);total[e]=total.get(e,Q(0))+c*int(t[0])
        assert not any(total.values());identity=True
        terms=[]
        for i,a in enumerate(coeff):
            if a:terms.append(dict(coefficient=str(a),exponents=monomial28(defs['basis'][i],defs['degree24_basis_definitions']),selected_basis_index=i,universal_candidate_index=defs['selected_basis_indices'][i]))
        names=['F','p1','p2','p3','p4','Q0','Q1','Q2','Q3','T','A0','A1','A2']
        meanings=dict(F='trace(mu^2), mu^t=-C+2D',pencil='det(tC+zD)=sum p_i*t^(4-i)*z^i',Q0='trace(C D^2 C D^2)-trace(D^2)*trace(C^2 D^2)/2+det(D)*trace(C^2)',Q1='trace(C^3 D C D)',Q2='trace(C^3 D^3)',Q3='trace(theta^5 H16)',T='trace(theta^7 H16)',A0='{trace(C^4),trace(C^3 D)}_Omega',A1='{trace(C^4),trace(C^2 D^2)}_Omega',A2='{trace(C^3 D),trace(C^2 D^2)}_Omega',theta='rho16(mu)',H16='Omega^-1 Hess_beta(F/2)',Poisson='ambient16 gradient(f)^t Omega^-1 gradient(g)')
        equation=dict(name='G28',variables=names,weights=[4,8,8,8,8,12,12,12,12,16,14,14,14],definitions=meanings,terms=terms,nonzero_selected_QQ_basis=True,all_terms_weight28=True,exact_QQ_source_identity=True,python_Fraction_all_coefficients_zero=True,coefficient_union_size=len(total),source_terms_checked=sum(len(p['terms'])for p in cols),primitive_integer_kernel108=list(map(str,coeff)))
        (out/'equation.json').write_text(json.dumps(equation,indent=2)+'\n')
        assert all(not any(z['exponents'][10:])for z in terms)
        tt=[z for z in terms if z['exponents'][9]];rr=[z for z in terms if not z['exponents'][9]]
        assert [(z['coefficient'],z['exponents'][:5])for z in tt]==[('-80000',[3,0,0,0,0]),('9600000',[1,1,0,0,0]),('6400000',[1,0,1,0,0])]
        (out/'equation-linear-T.json').write_text(json.dumps(dict(formula='G28=-80000*F*(F^2-120*p1-80*p2)*T+R75',R75_terms=rr,active_variables=names[:10],unused_variables=names[10:]),indent=2)+'\n')
        text='G28 =\n'+'\n'.join(f" {z['coefficient']} * "+'*'.join(n+(f'^{e}'if e!=1 else'')for n,e in zip(names,z['exponents'])if e)for z in terms)+'\n'
        (out/'equation.txt').write_text(text)
        with gzip.open(out/'native-certificate.json.gz','wb')as f:f.write(raw.read_bytes())
        assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    except Exception as e:error=repr(e)
    gates=[ROOT/'reports/Q_NEW_FAMILY_DEGREE28_PRODUCT_BASIS_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE28_RESTRICTION_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed']for p in gates)
    paths=[Path(__file__),native,helper,ROOT/f'research/notes/{ID}.md',source/'report.json',source/'modular-selected-minor.json',source/'native-certificate.json',univ/'report.json',univ/'definitions.json',*gates,*[source/f'column{i}.json'for i in range(108)]]
    report=dict(id=ID,all_passed=error is None,exact_QQ_identity=identity,error=error,nonzero_universal_QQ_invariant=identity,homogeneous_degree=28 if identity else None,
        equation_term_count=len(equation['terms'])if equation else None,
        candidate_actual_irreducible_defining_equation28_after_source_peer=identity,prior_lower28=True,
        no_new_source_point_prime_degree_or_fit=True,one_exact107solve_and_allQQcoefficient_identity_only=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=report['all_passed'],identity=identity,terms=report['equation_term_count'],seconds=report['seconds'],error=error)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
