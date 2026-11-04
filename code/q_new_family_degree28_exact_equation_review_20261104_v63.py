#!/usr/bin/env python3
"""Exact independent QQ polynomial source binding and degree28 identity."""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
import gzip
import itertools
import json
import math
from pathlib import Path
import resource
import signal
import time
import q_new_family_v_zero_degree26_restriction_review_20261104_v63 as old

ROOT=Path(__file__).resolve().parents[1]
BASE='Q_NEW_FAMILY_DEGREE28_EXACT_EQUATION_20261104_v63'
STEM=BASE+'_REVIEW'
STAGE='initial'
PRODUCTS=0

class QP:
    """Canonical integer coefficient dictionary / one positive denominator."""
    __slots__=('a','d')
    def __init__(self,a=None,d=1):
        a={m:int(c) for m,c in (a or {}).items() if c}
        d=int(d);assert d
        if d<0:a={m:-c for m,c in a.items()};d=-d
        g=math.gcd(d,math.gcd(*a.values())) if a else d
        self.a={m:c//g for m,c in a.items()};self.d=d//g
    def __eq__(self,b):return self.d==b.d and self.a==b.a
    def __bool__(self):return bool(self.a)

ZERO=QP();ONE=QP({0:1})
def add(a,b):
    if not a:return b
    if not b:return a
    g=math.gcd(a.d,b.d);x=b.d//g;y=a.d//g
    c={m:v*x for m,v in a.a.items()}
    for m,v in b.a.items():c[m]=c.get(m,0)+v*y
    return QP(c,a.d*x)
def scale(a,c):
    c=F(c)
    return QP({m:v*c.numerator for m,v in a.a.items()},a.d*c.denominator)
def mul(a,b):
    global PRODUCTS
    PRODUCTS+=1
    if not a or not b:return ZERO
    c=defaultdict(int)
    for m,x in a.a.items():
        for n,y in b.a.items():c[m+n]+=x*y
    return QP(c,a.d*b.d)
def total(values):
    z=ZERO
    for p in values:z=add(z,p)
    return z
def product(values):
    z=ONE
    for p in values:z=mul(z,p)
    return z
def power(a,n):return product(a for _ in range(n))
def mm(a,b,n):
    c=[ZERO for _ in range(n*n)]
    for i in range(n):
        for k in range(n):
            if not a[n*i+k]:continue
            for j in range(n):
                if b[n*k+j]:c[n*i+j]=add(c[n*i+j],mul(a[n*i+k],b[n*k+j]))
    return c
def tracepair(a,b,n):return total(mul(a[n*i+j],b[n*j+i]) for i in range(n) for j in range(n))
def invert_numeric(M):
    n=len(M);a=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(M)]
    for j in range(n):
        k=next(i for i in range(j,n) if a[i][j]);a[j],a[k]=a[k],a[j]
        p=a[j][j];a[j]=[v/p for v in a[j]]
        for i in range(n):
            if i!=j and a[i][j]:
                c=a[i][j];a[i]=[v-c*w for v,w in zip(a[i],a[j])]
    return [row[n:] for row in a]
def saved_poly(z):
    c=F(z['content']);return QP({old.pack(t[1:]):c.numerator*int(t[0]) for t in z['terms'] if t[0]},c.denominator)
def from_records(records):
    cs=[F(c) for _,c in records];den=math.lcm(*(c.denominator for c in cs)) if cs else 1
    return QP({old.pack(e):c.numerator*(den//c.denominator) for (e,_),c in zip(records,cs)},den)
def record(p):return {'denominator':str(p.d),'integer_terms':[[str(c),*old.unpack(m)] for m,c in sorted(p.a.items())]}
def fraction_entry(p,e):return F(p.a.get(old.pack(e),0),p.d)

def main(out):
    global STAGE
    started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    def timeout(*_):raise TimeoutError('exact independent QQ review60s cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    folder=ROOT/'reports'/BASE
    source=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE28_RESTRICTION_20261104_v63'
    previous=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63'
    univ=ROOT/'reports/Q_NEW_FAMILY_DEGREE28_PRODUCT_BASIS_20261104_v63'
    gates=[ROOT/'reports'/n/'report.json' for n in [
        'Q_NEW_FAMILY_DEGREE28_PRODUCT_BASIS_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_V_ZERO_DEGREE28_RESTRICTION_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63_REVIEW',
        'Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63_REVIEW']]
    primary=old.read(folder/'report.json');gd=[old.read(p) for p in gates]
    assert primary['all_passed'] and primary['exact_QQ_identity'] and all(d['all_passed'] for d in gd)
    assert gd[0]['complete_QQ_Inv28_basis'] and gd[2]['QQ_restriction_injective_and_NEW_defining_degree_lower_bound']==28
    checked={}
    for reportpath in [folder/'report.json',*gates]:
        d=old.read(reportpath)
        for key in ['input_sha256','checked_dependency_sha256']:
            for name,want in d.get(key,{}).items():assert old.sha(ROOT/name)==want,name;checked[name]=want
        for name,want in d.get('evidence_sha256',{}).items():
            p=reportpath.parent/name;assert old.sha(p)==want,str(p);checked[str(p.relative_to(ROOT))]=want
    spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    infopath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    prepath=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_SUPPORT_PREFLIGHT_20261104_v63/certificate.json'
    spin,info,pre=map(old.read,[spinpath,infopath,prepath])
    beta=[from_records(z) for z in pre['beta16_sparse']]
    assert len(beta)==16 and sum(len(z.a) for z in beta)==12
    mon=list(map(tuple,spin['quadratic_monomials']));mi={m:i for i,m in enumerate(mon)}
    quad=[mul(beta[i],beta[j]) for i,j in mon]
    STAGE='QQ global C,D and all512 ambient partial derivatives'
    cofs=[list(map(F,row)) for key in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136'] for row in spin[key]]
    matrices=[];derivatives=[]
    for a in range(2):
        values=[];der=[[ZERO for _ in range(16)] for _ in range(16)]
        for pos,row in enumerate(cofs[16*a:16*a+16]):
            value=ZERO
            for k,((i,j),c) in enumerate(zip(mon,row)):
                if not c:continue
                value=add(value,scale(quad[k],c))
                der[i][pos]=add(der[i][pos],scale(beta[j],c))
                der[j][pos]=add(der[j][pos],scale(beta[i],c))
            assert total(mul(beta[k],der[k][pos]) for k in range(16))==scale(value,2)
            assert value==from_records(pre['CD_quadratic32'][16*a+pos])
            for k in range(16):assert der[k][pos]==from_records(pre['CD_ambient_derivative512'][(16*a+pos)*16+k])
            values.append(value)
        matrices.append(values);derivatives.append(der)
    C,D=matrices
    STAGE='QQ global moment theta and true16 Hessian reconstruction'
    Omega=[list(map(F,row)) for row in info['symplectic_form16']];Oi=invert_numeric(Omega)
    rho=[[list(map(F,row)) for row in A] for A in info['all_ten_representation_matrices']]
    J4=[[F(0) for _ in range(4)] for _ in range(4)];J4[0][1]=J4[2][3]=F(1);J4[1][0]=J4[3][2]=F(-1)
    gs=[]
    for i in range(4):
        for j in range(i,4):
            S=[[0 for _ in range(4)] for _ in range(4)];S[i][j]=S[j][i]=1
            gs.append([[sum(J4[r][k]*S[k][c] for k in range(4)) for c in range(4)] for r in range(4)])
    Gram=[[sum(A[i][j]*B[j][i] for i in range(4) for j in range(4)) for B in gs] for A in gs];Gi=invert_numeric(Gram)
    mu=[]
    for A in rho:
        S=[[sum(Omega[i][k]*A[k][j] for k in range(16) if Omega[i][k] and A[k][j]) for j in range(16)] for i in range(16)]
        assert all(S[i][j]==S[j][i] for i in range(16) for j in range(16))
        mu.append([S[i][j]/(2 if i==j else 1) for i,j in mon])
    gmu=[[sum(Gi[a][b]*mu[b][k] for b in range(10) if Gi[a][b] and mu[b][k]) for k in range(136)] for a in range(10)]
    tco=[[sum(rho[a][i][j]*gmu[a][k] for a in range(10) if rho[a][i][j] and gmu[a][k]) for k in range(136)] for i in range(16) for j in range(16)]
    hc=[[F(0) for _ in range(136)] for _ in range(256)]
    for indices,c in spin['quartic_trace_moment_square']:
        counts=Counter(indices);c=F(c)/2
        for i,ei in counts.items():
            for j,ej in counts.items():
                multiplicity=ei*(ej-int(i==j))
                if multiplicity:
                    rest=list(indices);rest.remove(i);rest.remove(j)
                    hc[16*i+j][mi[tuple(sorted(rest))]]+=c*multiplicity
    hco=[[sum(Oi[i][k]*hc[16*k+j][m] for k in range(16) if Oi[i][k] and hc[16*k+j][m]) for m in range(136)] for i in range(16) for j in range(16)]
    # Exact coefficient-level binding to the existing normalized integer entrance.
    tok=iter(gzip.decompress((previous/'native-input.txt.gz').read_bytes()).decode().split())
    assert (int(next(tok)),int(next(tok)))==(16,136)
    for expected in beta:
        records=[]
        for _ in range(int(next(tok))):
            c=next(tok);records.append(([int(next(tok)) for _ in range(4)],c))
        assert from_records(records)==expected
    assert [tuple(int(next(tok)) for _ in range(2)) for _ in range(136)]==mon
    assert [[F(next(tok)) for _ in range(136)] for _ in range(32)]==[[4*c for c in row] for row in cofs]
    assert [[F(next(tok)) for _ in range(136)] for _ in range(512)]==[[2*c for c in row] for row in tco+hco]
    assert [F(next(tok)) for _ in range(256)]==[x for row in Oi for x in row]
    theta=[total(scale(quad[k],c) for k,c in enumerate(row) if c) for row in tco]
    Hess=[total(scale(quad[k],c) for k,c in enumerate(row) if c) for row in hco]
    STAGE='QQ all8 trace values; selected108 basis needs no Poisson gradients'
    cache={'':[ONE if i==j else ZERO for i in range(4) for j in range(4)],'0':C,'1':D}
    def word(w):
        if w not in cache:cache[w]=mm(word(w[:-1]),matrices[int(w[-1])],4)
        return cache[w]
    words=['0000','0001','0011','0101','0111','000101','000111','011011']
    traces=[total(word(w)[4*i+i] for i in range(4)) for w in words]
    assert traces==[saved_poly(p) for p in old.read(previous/'trace-values8.json')]
    moment=[add(scale(c,-1),scale(d,2)) for c,d in zip(C,D)]
    F4=tracepair(moment,moment,4)
    D2=mm(D,D,4);C2=mm(C,C,4);CD2=mm(C,D2,4)
    d2=tracepair(D,D,4);c2=tracepair(C,C,4)
    determinantD=add(scale(mul(d2,d2),F(1,8)),scale(tracepair(D2,D2,4),F(-1,4)))
    N12=add(add(tracepair(CD2,CD2,4),scale(mul(d2,tracepair(C2,D2,4)),F(-1,2))),mul(determinantD,c2))
    STAGE='QQ sixteen-dimensional Hessian traces through7 and pencil'
    level=Hess;H12=None
    for degree in range(1,8):
        level=mm(theta,level,16)
        if degree==5:H12=total(level[16*i+i] for i in range(16))
    H16=total(level[16*i+i] for i in range(16))
    pencil=[ZERO for _ in range(5)]
    for perm in itertools.permutations(range(4)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        for letters in itertools.product(range(2),repeat=4):
            value=product((C,D)[letters[i]][4*i+perm[i]] for i in range(4));k=sum(letters)
            pencil[k]=add(pencil[k],scale(value,sign))
    assert mul(F4,F4)==add(scale(pencil[0],400),scale(pencil[4],1600))
    primitive12=[N12,traces[5],traces[6],H12]
    gens11=[F4,*pencil,*primitive12,H16]
    assert gens11==[saved_poly(p) for p in old.read(source/'source-generators14.json')[:11]]
    STAGE='QQ semantic108 basis columns'
    defs=old.read(univ/'definitions.json');spectral=[F4,*pencil[1:]]
    assert all(z['kind']!='Inv14_product' for z in defs['basis'])
    def inv24(z):
        k=z['kind']
        if k=='spectral':return product(power(p,e) for p,e in zip(spectral,z['F4_p1_p2_p3_p4_exponents']))
        if k=='primitive12_times_spectral12':
            j=z['spectral'];v=power(F4,3) if j==0 else mul(F4,pencil[j])
            return mul(primitive12[z['primitive']],v)
        if k=='primitive12_product':return product(primitive12[i] for i in z['indices'])
        assert k=='primitive16_times_Inv8'
        return mul(H16,power(F4,2) if z['index']==0 else pencil[z['index']])
    basis24=[inv24(z) for z in defs['degree24_basis_definitions']]
    columns=[];termcount=0
    for i,z in enumerate(defs['basis']):
        k=z['kind']
        if k=='F4_times_Inv24_basis':p=mul(F4,basis24[z['Inv24_basis_index']])
        elif k=='primitive12_times_p_p':p=product([primitive12[z['primitive12_index']],*[pencil[j] for j in z['p_indices']]])
        elif k=='primitive16_times_primitive12':p=mul(H16,primitive12[z['primitive12_index']])
        else:raise AssertionError('selected complete basis must contain no Inv14 product')
        saved=old.read(source/f'column{i}.json')
        assert p==saved_poly(saved['HLsq_polynomial'])
        assert saved['universal_candidate_index']==defs['selected_basis_indices'][i] and saved['xi0_power']==14
        termcount+=len(p.a);columns.append(p)
    assert termcount==96453
    STAGE='QQ primitive kernel, selected matrix and all source coefficients'
    equation=old.read(folder/'equation.json');native=old.read(folder/'native-certificate.json.gz')
    coeff=list(map(int,equation['primitive_integer_kernel108']))
    assert coeff==list(map(int,native['primitive_integer_kernel108'])) and len(coeff)==108 and math.gcd(*coeff)==1
    assert any(coeff) and total(scale(p,c) for p,c in zip(columns,coeff))==ZERO
    fixed=old.read(source/'modular-selected-minor.json');rows=fixed['coefficient_rows'];pivots=fixed['pivot_columns'];free=fixed['free_columns'][0]
    hrows=[(2*r-3*b,2*b,s,q) for r,b,s,q in rows]
    matrix=[[fraction_entry(columns[j],e) for j in pivots] for e in hrows]
    rhs=[-fraction_entry(columns[free],e) for e in hrows]
    assert matrix==[[F(c) for c in row] for row in native['selected_matrix_QQ']]
    assert rhs==list(map(F,native['rhs_QQ']))
    solution=list(map(F,native['solution_QQ']))
    assert len(solution)==107 and all(sum(a*b for a,b in zip(row,solution))==c for row,c in zip(matrix,rhs))
    assert coeff[free] and [F(coeff[j],coeff[free]) for j in pivots]==solution
    # Independently derive the low-generator monomials, with p0 eliminated as in basis.
    names=['F','p1','p2','p3','p4','Q0','Q1','Q2','Q3','T','A0','A1','A2']
    weights=[4,8,8,8,8,12,12,12,12,16,14,14,14]
    # Zero placeholders are never used: the actual selected108 and G28 omit A0..2.
    low13=[F4,*pencil[1:],*primitive12,H16,ZERO,ZERO,ZERO]
    def exp24(z):
        e=[0]*13;k=z['kind']
        if k=='spectral':e[:5]=z['F4_p1_p2_p3_p4_exponents']
        elif k=='primitive12_times_spectral12':
            e[5+z['primitive']]+=1;j=z['spectral'];e[0]+=3 if j==0 else 1
            if j:e[j]+=1
        elif k=='primitive12_product':
            for j in z['indices']:e[5+j]+=1
        else:
            assert k=='primitive16_times_Inv8';e[9]+=1;j=z['index'];e[0]+=2 if j==0 else 0
            if j:e[j]+=1
        return e
    def exp28(z):
        e=[0]*13;k=z['kind']
        if k=='F4_times_Inv24_basis':e=exp24(defs['degree24_basis_definitions'][z['Inv24_basis_index']]);e[0]+=1
        elif k=='primitive12_times_p_p':
            e[5+z['primitive12_index']]+=1
            for j in z['p_indices']:e[j]+=1
        elif k=='primitive16_times_primitive12':e[9]+=1;e[5+z['primitive12_index']]+=1
        else:
            for j in z['indices']:e[10+j]+=1
        return e
    expected=[]
    for i,(c,z) in enumerate(zip(coeff,defs['basis'])):
        e=exp28(z);assert sum(a*b for a,b in zip(e,weights))==28 and e[10:]==[0,0,0]
        assert product(power(p,j) for p,j in zip(low13,e))==columns[i]
        if c:expected.append({'coefficient':str(c),'exponents':e,'selected_basis_index':i,'universal_candidate_index':defs['selected_basis_indices'][i]})
    assert names==equation['variables'] and weights==equation['weights'] and expected==equation['terms'] and len(expected)==78
    STAGE='freeze exact QQ review'
    out.mkdir(parents=True,exist_ok=False)
    evidence={'low11_exactQQ':[record(p) for p in gens11],
        'all108_source_columns_exactQQ':[record(p) for p in columns],
        'primitive_integer_kernel108':list(map(str,coeff)),'all_QQ_coefficients_sum_zero':True,
        'all107_selected_QQ_equations_verified':True,'weighted28_terms78':expected}
    ep=out/'exact-polynomial-proof.json.gz'
    with gzip.open(ep,'wt') as f:json.dump(evidence,f);f.write('\n')
    assert ep.stat().st_size<32*1024*1024
    for p in [Path(__file__),ROOT/'code/q_new_family_v_zero_degree26_restriction_review_20261104_v63.py',ROOT/'research/notes'/(STEM+'.md'),ROOT/'research/notes'/(BASE+'.md'),folder/'report.json',folder/'equation.json',folder/'native-certificate.json.gz',spinpath,infopath,prepath,*gates]:checked[str(p.relative_to(ROOT))]=old.sha(p)
    report={'id':STEM,'all_passed':True,
        'exact_QQ_original_global_CD_Theta_Hessian_coefficient_binding':True,
        'all512_CD_partial_derivatives_exact_QQ':True,
        'all11_needed_low_invariant_values_exactQQ_independently_reconstructed':True,
        'selected108_complete_basis_and_G28_omit_all_Inv14_products_and_gradients':True,
        'all108_columns96453_trueQQ_coefficients_independently_equal':True,
        'all107_QQ_selected_equations_and_primitive_kernel_verified_without_solve':True,
        'all_source_QQ_coefficients_zero':True,'semantic_13low78term_degree28_equation_verified':True,
        'actual_reduced_primitive_firstjet_B15_irreducible_equation_degree':28,
        'QQ_and_complex_defining_equation_up_to_nonzero_scalar':True,
        'full_S61_ideal_or_nilstructure_or_invariant_ring_presentation_claim':False,
        'no_native_replay_newpoint_newprime_higherdegree_or_enlarged_retry':True,
        'input_sha256':checked,'evidence_sha256':{ep.name:old.sha(ep)},
        'seconds':time.monotonic()-started,'integer_polynomial_products':PRODUCTS,
        'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','seconds','integer_polynomial_products','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:main(a.output)
    except (TimeoutError,MemoryError) as e:
        a.output.mkdir(parents=True,exist_ok=True)
        (a.output/'blocked-review.json').write_text(json.dumps({'all_passed':False,'stage':STAGE,'reason':str(e),'no_automatic_larger_retry':True})+'\n')
        raise
