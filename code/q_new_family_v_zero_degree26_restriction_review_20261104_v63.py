#!/usr/bin/env python3
"""Independent integer-polynomial restriction, ambient16 gradients and fixed minor."""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
import gzip
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
BASE='Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63';STEM=BASE+'_REVIEW'
P=32003;RADIX=64;POW=[RADIX**i for i in range(4)]
STAGE='initial';PRODUCTS=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def mod(q):
    q=F(str(q));assert q.denominator%P
    return q.numerator*pow(q.denominator,-1,P)%P
def pack(e):return sum(a*b for a,b in zip(e,POW))
def unpack(m):return tuple((m//b)%RADIX for b in POW)
def add(a,b,c=1):
    if not b or c%P==0:return a.copy()
    z=dict(a)
    for m,v in b.items():
        q=(z.get(m,0)+c*v)%P
        if q:z[m]=q
        else:z.pop(m,None)
    return z
def scale(a,c):return {m:v*c%P for m,v in a.items() if v*c%P}
def mul(a,b):
    global PRODUCTS
    PRODUCTS+=1
    if not a or not b:return {}
    z=defaultdict(int)
    for m,x in a.items():
        for n,y in b.items():z[m+n]+=x*y
    return {m:v%P for m,v in z.items() if v%P}
def total(values):
    out={}
    for p in values:out=add(out,p)
    return out
def mm(a,b,n):
    c=[{} for _ in range(n*n)]
    for i in range(n):
        for k in range(n):
            if not a[n*i+k]:continue
            for j in range(n):
                if b[n*k+j]:c[n*i+j]=add(c[n*i+j],mul(a[n*i+k],b[n*k+j]))
    return c
def tracepair(a,b,n):return total(mul(a[n*i+j],b[n*j+i]) for i in range(n) for j in range(n))
def inverse(M):
    n=len(M);a=[row.copy()+[int(i==j) for j in range(n)] for i,row in enumerate(M)]
    for j in range(n):
        k=next(i for i in range(j,n) if a[i][j]);a[j],a[k]=a[k],a[j]
        q=pow(a[j][j],-1,P);a[j]=[x*q%P for x in a[j]]
        for i in range(n):
            if i!=j and a[i][j]:
                q=a[i][j];a[i]=[(x-q*y)%P for x,y in zip(a[i],a[j])]
    return [row[n:] for row in a]
def det(M):
    a=[row.copy() for row in M];answer=1
    for j in range(len(a)):
        k=next(i for i in range(j,len(a)) if a[i][j])
        if k!=j:a[k],a[j]=a[j],a[k];answer=-answer%P
        q=a[j][j];answer=answer*q%P;iv=pow(q,-1,P)
        for i in range(j+1,len(a)):
            c=a[i][j]*iv%P
            if c:a[i]=[(x-c*y)%P for x,y in zip(a[i],a[j])]
    return answer
def saved_poly(record):
    c=mod(record['content']);out={}
    for row in record['terms']:
        q=mod(row[0])*c%P
        if q:out[pack(row[1:])]=q
    return out
def sparse_records(poly):return [[list(unpack(e)),v] for e,v in sorted(poly.items())]


def main(out):
    global STAGE
    started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    def timeout(*_):raise TimeoutError('independent review60s cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    folder=ROOT/'reports'/BASE;primary=read(folder/'report.json');native=read(folder/'native-certificate.json');recipe=read(folder/'recipe.json')
    assert primary['all_passed'] and native['rank']==45 and native['coefficient_minor_mod32003']==7069
    assert recipe['all_ambient16_derivatives_before_substitution'] and native['selected_universal_candidates']==recipe['selected45_universal_candidate_indices']
    info=read(ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json')
    spin=read(ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz')
    prepath=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_SUPPORT_PREFLIGHT_20261104_v63/certificate.json';pre=read(prepath)
    v0path=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63/certificate.json';v0=read(v0path)
    gates=[ROOT/'reports'/n/'report.json' for n in [
        'Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_FULL_DEGREE24_SOURCE_RELATION_20261104_v63_REVIEW']]
    gd=[read(p) for p in gates];assert all(d['all_passed'] for d in gd)
    assert gd[1]['complete_QQ_Inv14_3_Inv18_9_Inv26_45_bases_certified']
    STAGE='actual v0 finite torus and scalar binding'
    xi,r0,ss,qq,B0,R,T,H,L=s.symbols('xi0 r0 s0 q0 B0 R T H L')
    weights=info['Cartan_weights'];beta=[]
    for text,(w1,w2),stored in zip(v0['beta_bar0_16'],weights,pre['beta16_sparse']):
        z=s.sympify(text)
        assert s.expand(2*xi*s.diff(z,xi)-(1-w1-w2)*z)==0
        p=s.Poly(s.cancel(z.subs({xi:1,r0:R**4,B0:T**2})*R**(-w1+2*w2)*T**(w1+2*w2)),R,T,ss,qq)
        exact={}
        for (er,et,es,eq),co in p.terms():
            if not co:continue
            assert er>=3*et and (er-3*et)%2==0
            exact[((er-3*et)//2,et,es,eq)]=F(str(co))
        assert exact=={tuple(e):F(c) for e,c in stored}
        beta.append({pack(e):mod(c) for e,c in exact.items()})
    assert sum(bool(p) for p in beta)==9 and sum(map(len,beta))==12
    assert all((w1+w2)%2 for w1,w2 in weights)
    mon=list(map(tuple,spin['quadratic_monomials']));mi={m:i for i,m in enumerate(mon)}
    quad=[mul(beta[i],beta[j]) for i,j in mon]
    cofs=[[mod(c) for c in row] for key in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136'] for row in spin[key]]
    STAGE='all32 global quadratics and512 ambient derivatives'
    matrices=[];derivatives=[]
    for a in range(2):
        values=[];der=[[{} for _ in range(16)] for _ in range(16)]
        for pos,row in enumerate(cofs[16*a:16*a+16]):
            value={}
            for k,((i,j),c) in enumerate(zip(mon,row)):
                if not c:continue
                value=add(value,quad[k],c)
                # Both branches retained even if beta_i vanishes on this section.
                der[i][pos]=add(der[i][pos],beta[j],c)
                der[j][pos]=add(der[j][pos],beta[i],c)
            assert total(mul(beta[k],der[k][pos]) for k in range(16))==scale(value,2)
            values.append(value)
            assert value=={pack(e):mod(c) for e,c in pre['CD_quadratic32'][16*a+pos] if mod(c)}
            for k in range(16):
                expected={pack(e):mod(c) for e,c in pre['CD_ambient_derivative512'][(16*a+pos)*16+k] if mod(c)}
                assert der[k][pos]==expected
        matrices.append(values);derivatives.append(der)
    Omega=[[mod(c) for c in row] for row in info['symplectic_form16']];Oi=inverse(Omega)
    rho=[[[mod(c) for c in row] for row in A] for A in info['all_ten_representation_matrices']]
    STAGE='independent global theta/H16 coefficient reconstruction'
    J4=[[0]*4 for _ in range(4)]
    J4[0][1]=J4[2][3]=1;J4[1][0]=J4[3][2]=P-1
    gs=[]
    for i in range(4):
        for j in range(i,4):
            S=[[0]*4 for _ in range(4)];S[i][j]=S[j][i]=1
            gs.append([[sum(J4[r][k]*S[k][c] for k in range(4))%P for c in range(4)] for r in range(4)])
    Gram=[[sum(A[i][j]*B[j][i] for i in range(4) for j in range(4))%P for B in gs] for A in gs];Gi=inverse(Gram)
    mu=[]
    for A in rho:
        S=[[sum(Omega[i][k]*A[k][j] for k in range(16))%P for j in range(16)] for i in range(16)]
        assert all(S[i][j]==S[j][i] for i in range(16) for j in range(16))
        mu.append([S[i][j]*(pow(2,-1,P) if i==j else 1)%P for i,j in mon])
    gmu=[[sum(Gi[a][b]*mu[b][k] for b in range(10))%P for k in range(136)] for a in range(10)]
    tco=[[sum(rho[a][i][j]*gmu[a][k] for a in range(10))%P for k in range(136)] for i in range(16) for j in range(16)]
    hc=[[0]*136 for _ in range(256)]
    for indices,c in spin['quartic_trace_moment_square']:
        counts=Counter(indices);c=mod(c)*pow(2,-1,P)%P
        for i,ei in counts.items():
            for j,ej in counts.items():
                multiplicity=ei*(ej-int(i==j))
                if multiplicity:
                    rest=list(indices);rest.remove(i);rest.remove(j)
                    k=mi[tuple(sorted(rest))];hc[16*i+j][k]=(hc[16*i+j][k]+c*multiplicity)%P
    hco=[[sum(Oi[i][k]*hc[16*k+j][m] for k in range(16) if Oi[i][k])%P for m in range(136)] for i in range(16) for j in range(16)]
    # Bind every normalized input coefficient, while never invoking the native solver.
    tokens=iter(gzip.decompress((folder/'native-input.txt.gz').read_bytes()).decode().split())
    assert (int(next(tokens)),int(next(tokens)))==(16,136)
    for expected in beta:
        found={}
        for _ in range(int(next(tokens))):
            c=mod(next(tokens));e=[int(next(tokens)) for _ in range(4)]
            if c:found[pack(e)]=c
        assert found==expected
    assert [tuple(int(next(tokens)) for _ in range(2)) for _ in range(136)]==mon
    assert [[int(next(tokens))%P for _ in range(136)] for _ in range(32)]==[[4*c%P for c in row] for row in cofs]
    assert [[int(next(tokens))%P for _ in range(136)] for _ in range(512)]==[[2*c%P for c in row] for row in tco+hco]
    assert [mod(next(tokens)) for _ in range(256)]==[x for row in Oi for x in row]
    selected=recipe['selected45_universal_candidate_indices'];assert [int(next(tokens)) for _ in range(45)]==selected
    assert next(tokens,None) is None
    theta=[total(scale(quad[k],c) for k,c in enumerate(row) if c) for row in tco]
    Hess=[total(scale(quad[k],c) for k,c in enumerate(row) if c) for row in hco]
    C,D=matrices
    cache={'':[{0:1} if i==j else {} for i in range(4) for j in range(4)],'0':C,'1':D}
    def word(w):
        if w not in cache:cache[w]=mm(word(w[:-1]),matrices[int(w[-1])],4)
        return cache[w]
    words=['0000','0001','0011','0101','0111','000101','000111','011011']
    traces=[];grads=[]
    STAGE='all8 traces and128 ambient gradients'
    for w in words:
        Wm=word(w);value=total(Wm[4*i+i] for i in range(4));parts=[[{} for _ in range(16)] for _ in range(2)]
        for i,a in enumerate(w):
            rest=word(w[i+1:]+w[:i])
            parts[int(a)]=[add(x,y) for x,y in zip(parts[int(a)],rest)]
        gradient=[total(tracepair(parts[a],derivatives[a][k],4) for a in range(2)) for k in range(16)]
        assert total(mul(beta[k],gradient[k]) for k in range(16))==scale(value,2*len(w))
        traces.append(value);grads.append(gradient)
    assert traces==[saved_poly(p) for p in read(folder/'trace-values8.json')]
    assert [p for row in grads for p in row]==[saved_poly(p) for p in read(folder/'ambient-gradients128.json')]
    def bracket(f,g):return total(scale(mul(grads[f][i],grads[g][j]),Oi[i][j]) for i in range(16) for j in range(16) if Oi[i][j])
    Q14=[bracket(0,1),bracket(0,2),bracket(1,2)]
    Q18=[bracket(i,j) for i in range(2) for j in range(5,8)]
    moment=[add(scale(c,-1),scale(d,2)) for c,d in zip(C,D)]
    F4=tracepair(moment,moment,4);c2=tracepair(C,C,4);d2=tracepair(D,D,4)
    assert scale(c2,5)==scale(d2,10)==scale(tracepair(C,D,4),-10)==F4
    D2=mm(D,D,4);C2=mm(C,C,4);CD2=mm(C,D2,4)
    determinantD=add(scale(mul(d2,d2),pow(8,-1,P)),scale(tracepair(D2,D2,4),-pow(4,-1,P)))
    N12=add(add(tracepair(CD2,CD2,4),scale(mul(d2,tracepair(C2,D2,4)),-pow(2,-1,P))),mul(determinantD,c2))
    STAGE='full16x16 theta fifth power Hessian trace'
    level=Hess
    for _ in range(5):level=mm(theta,level,16)
    H12=total(level[16*i+i] for i in range(16))
    I12=[mul(F4,p) for p in traces[:5]]+[N12,traces[5],traces[6],H12]
    low=Q14+Q18+[F4]+I12
    assert low==[saved_poly(p) for p in read(folder/'low-values19.json')]
    semantics=[dict(kind='Inv14_times_Inv12',Inv14_basis_index=i,Inv12_basis_index=j) for i in range(3) for j in range(9)]
    semantics +=[dict(kind='primitive18_times_Inv8',primitive18_basis_index=i,Inv8_index=j) for i in range(6) for j in range(5)]
    assert semantics==recipe['definitions']['degree26']
    assert [semantics[i] for i in selected]==recipe['definitions']['Inv26_basis']
    STAGE='all45 columns and2025 fixed coefficient entries'
    columns=[];descended=[];QQterms=0
    for j,index in enumerate(selected):
        p=mul(Q14[index//9],I12[index%9]) if index<27 else mul(Q18[(index-27)//5],traces[(index-27)%5])
        saved=read(folder/f'column{j}.json')
        assert saved['universal_candidate_index']==index and saved['xi0_power']==13
        assert p==saved_poly(saved['HLsq_polynomial']);QQterms+=len(saved['HLsq_polynomial']['terms'])
        assert len(saved['HLsq_polynomial']['terms'])<=2874
        q={}
        for row in saved['HLsq_polynomial']['terms']:
            aa,b,cc,dd=row[1:]
            assert b%2==0 and (2*aa+3*b)%4==0 and aa+b<=26 and b<=8 and cc+dd<=(26-3*b)//2
        for m,c in p.items():
            aa,b,cc,dd=unpack(m);ex=((2*aa+3*b)//4,b//2,cc,dd)
            assert ex not in q;q[ex]=c
        columns.append(p);descended.append(q)
    assert QQterms==22674
    coefficient_rows=[tuple(e) for e in native['coefficient_monomials_r0_B0_s0_q0']]
    ownminor=[[p.get(e,0) for p in descended] for e in coefficient_rows]
    assert len(coefficient_rows)==len(set(coefficient_rows))==45
    assert ownminor==native['selected_coefficient_matrix'] and det(ownminor)==7069
    checked={}
    for p in [folder/'report.json',*gates]:
        d=read(p)
        for key in ['input_sha256','checked_dependency_sha256']:
            for name,want in d.get(key,{}).items():
                target=ROOT/name;assert sha(target)==want,str(target);checked[name]=want
        for name,want in d.get('evidence_sha256',{}).items():
            target=p.parent/name;assert sha(target)==want,str(target);checked[str(target.relative_to(ROOT))]=want
    out.mkdir(parents=True,exist_ok=False)
    ep=out/'independent-polynomials.json.gz'
    evidence={'prime':P,'variables':['H','L','s0','q0'],'packed_base64_no_carry_degree_bound':26,
        'beta16':[sparse_records(p) for p in beta],
        'all128_ambient_trace_gradients':[sparse_records(p) for row in grads for p in row],
        'trace_values8':[sparse_records(p) for p in traces],
        'low_values19':[sparse_records(p) for p in low],
        'all45_source_columns':[sparse_records(p) for p in columns],
        'coefficient_monomials':coefficient_rows,'independent45_by45_minor':ownminor,'minor_mod32003':7069}
    with gzip.open(ep,'wt') as f:json.dump(evidence,f);f.write('\n')
    for p in [Path(__file__),ROOT/'research/notes'/(STEM+'.md'),ROOT/'research/notes'/(BASE+'.md'),folder/'report.json',folder/'recipe.json',prepath,v0path,*gates]:
        checked[str(p.relative_to(ROOT))]=sha(p)
    report={'id':STEM,'all_passed':True,
        'independent_actual_dense_v0_torus_pullback12_terms_and_xi13_binding':True,
        'all32_global_quadratics512_ambient_derivatives128_trace_gradients_verified':True,
        'global_theta_H16_Omega_and_all_native_input_coefficients_reconstructed':True,
        'full16_Hessian_not_confused_with4_by4_D':True,
        'all45_complete_universal_basis_source_coefficients_reconstructed_mod32003':True,
        'all2025_fixed_minor_entries_verified':True,'coefficient_minor_mod32003':7069,
        'QQ_restriction_injective_and_NEW_defining_degree_lower_bound':28,
        'degree28_equation_or_source28_claim':False,
        'author_native_replay_newpoint_newprime_or_high_basis_growth':False,
        'input_sha256':checked,'evidence_sha256':{'independent-polynomials.json.gz':sha(ep)},
        'seconds':time.monotonic()-started,'polynomial_products':PRODUCTS,
        'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','coefficient_minor_mod32003','seconds','max_rss_kib']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:main(a.output)
    except (TimeoutError,MemoryError) as e:
        a.output.mkdir(parents=True,exist_ok=True)
        (a.output/'blocked-review.json').write_text(json.dumps({'all_passed':False,'stage':STAGE,'reason':str(e),'no_automatic_larger_retry':True})+'\n')
        raise
