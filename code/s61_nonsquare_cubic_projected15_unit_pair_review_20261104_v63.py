#!/usr/bin/env python3
"""Independent actual57 H-coefficients, Fraction Bezout and chart rank scope."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
BASE='S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63'
STEM=BASE+'_REVIEW'
Z=(0,0,0)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def add(a,b,c=F(1)):
    r=a.copy()
    for m,v in b.items():
        w=r.get(m,F(0))+c*v
        if w:r[m]=w
        else:r.pop(m,None)
    return r
def scale(a,c):return {m:v*c for m,v in a.items() if v*c}
def mul(a,b):
    r={}
    for m,v in a.items():
        for n,w in b.items():
            e=tuple(x+y for x,y in zip(m,n));r[e]=r.get(e,F(0))+v*w
    return {m:v for m,v in r.items() if v}
def total(values):
    p={}
    for v in values:p=add(p,v)
    return p
def dec(records):return {tuple(m):F(v) for m,v in records if F(v)}
def enc(p):return [[list(m),str(v)] for m,v in sorted(p.items())]
ONE={Z:F(1)}
def const(v):return {Z:F(v)} if v else {}
def inverse_constant(A):
    n=len(A);r=[[F(c) for c in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for i in range(n):
        k=next(k for k in range(i,n) if r[k][i]);r[i],r[k]=r[k],r[i]
        p=r[i][i];r[i]=[c/p for c in r[i]]
        for k in range(n):
            if k!=i and r[k][i]:p=r[k][i];r[k]=[x-p*y for x,y in zip(r[k],r[i])]
    return [row[n:] for row in r]
def radd(a,b):
    if not a[0]:return b
    if not b[0]:return a
    return add(mul(a[0],b[1]),mul(b[0],a[1])),mul(a[1],b[1])
def rmul(a,b):return (mul(a[0],b[0]),mul(a[1],b[1])) if a[0] and b[0] else ({},ONE)
def rtotal(vs):
    r=({},ONE)
    for v in vs:r=radd(r,v)
    return r

def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    folder=ROOT/'reports'/BASE;primary=read(folder/'report.json');certificate=read(folder/'certificate.json')
    assert primary['all_exact_QQ_heavy_factor_and_Bezout_identities_passed'] and primary['frozen']
    checked={}
    pre=ROOT/'reports/S61_NONSQUARE_CUBIC_PROJECTED15_HEAVY_PREFLIGHT_20261104_v63'
    for p in [folder/'report.json',pre/'report.json']:
        d=read(p)
        for name,h in d['input_sha256'].items():assert sha(ROOT/name)==h,name;checked[name]=h
        assert sha(p.parent/'certificate.json')==d['certificate_sha256']
    srpath=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json'
    brpath=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json'
    source,basis=map(read,[srpath,brpath]);mons=source['normal_monomials']
    raw=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63/frozen57-source.txt'
    assert sha(raw)==basis['remaining_matrix_sha256']
    retained=basis['remaining_original_row_indices'];assert retained[42]==43 and retained[47]==49
    assert source['labels'][43]==['cubic',20] and source['labels'][49]==['cubic',61]
    assert all(source['labels'][i][0]=='product' for i in retained[:41])
    assert all(source['labels'][i][0]=='cubic' for i in retained[41:])
    npath=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/normalization.txt'
    normal=npath.read_text().splitlines();assert normal[0]=='63'
    assert normal[44].split()==normal[50].split()==['1','1/625000000','0','0','0']
    light=[0,1,2,3,4,6,7];columns=[]
    for j in light:
        e=[0]*8;e[5]=e[j]=1;columns.append(mons.index(e))
    assert len(set(columns))==7 and len(mons)==91 and sum(bool(m[5]) for m in mons)==7
    entries={}
    with raw.open() as f:
        assert f.readline().split()==['57','91','57']
        for i in range(57):
            for j in range(91):
                line=f.readline();assert line
                if j not in columns:continue
                ts=line.split();n=int(ts[0]);assert len(ts)==1+4*n
                p={}
                for k in range(n):
                    v=F(ts[1+4*k]);e=tuple(map(int,ts[2+4*k:5+4*k]));p[e]=p.get(e,F(0))+v
                p={e:c for e,c in p.items() if c}
                if i<41:assert not p
                if i in [42,47]:entries[i,j]=p
        assert not f.read().strip()
    difference=[add(entries[47,j],entries[42,j],F(-1)) for j in columns]
    expected=read(pre/'certificate.json')['unscaled_G61_minus_G20_heavy_coefficients']
    assert difference==[dec(expected[str(j)]) for j in light]
    # Independently derive the coefficient change from the actual variable substitution.
    # Physical(v0,v1,v2,v3,v4,v6,v7)=(x0,x1,x2,W,x3,Z-3x0+x1,x4).
    physical=[(1,0,0,0,0,0,0),(0,1,0,0,0,0,0),(0,0,1,0,0,0,0),
        (0,0,0,0,0,0,1),(0,0,0,1,0,0,0),(-3,1,0,0,0,1,0),(0,0,0,0,1,0,0)]
    H=[list(c) for c in zip(*physical)];Hi=inverse_constant(H)
    assert H==certificate['constant_heavy_coefficient_change7x7']
    assert Hi==certificate['constant_heavy_inverse7x7']
    transformed=[total(scale(p,F(c)) for p,c in zip(difference,row)) for row in H]
    assert not any(transformed[:5]);kZ,kW=transformed[5:]
    a={(1,0,0):F(1)};e={(0,1,0):F(1)};d={(0,0,1):F(1)}
    C=add(add(scale(a,100),scale(mul(e,e),18)),add(scale(e,-9),const(-2)))
    FF=add(add(scale(a,500),scale(mul(e,e),90)),add(scale(e,3),const(-24)))
    B=add(mul(add(C,scale(a,-50)),d),add(scale(e,60),const(-40)))
    g1=add(mul(add(scale(e,3),const(-2)),d),const(20))
    g2=add(mul(add(scale(e,6),const(1)),d),const(40))
    h=add(mul(C,g1),B,F(-1));ell=add(scale(e,6),const(1))
    M=add(add(scale(a,600),scale(mul(e,e),126)),add(scale(e,-63),const(-14)))
    K=add(M,mul(ell,C),F(-1))
    assert scale(kZ,10000000)==scale(mul(mul(mul(mul(a,FF),B),g1),g2),-1)
    assert scale(kW,10000000)==scale(mul(mul(FF,B),h),-1)
    assert add(scale(mul(ell,h),2),mul(K,g2))==scale(C,240)
    assert add(mul(mul(mul(K,C),g1),g2),mul(add(scale(mul(ell,B),2),mul(K,g2),F(-1)),h))==scale(mul(B,C),240)
    delta=scale(mul(C,C),-64000)
    cZ=mul(K,C);cW=mul(a,add(scale(mul(ell,B),2),mul(K,g2),F(-1)))
    eta=scale(mul(mul(mul(mul(a,FF),mul(B,B)),mul(C,C)),C),F(192,125))
    z,w=mul(delta,kZ),mul(delta,kW)
    assert add(mul(cZ,z),mul(cW,w))==eta
    for key,p in [('g1',g1),('g2',g2),('h',h),('Lscalar',ell),('Kscalar',K),('Bezout_cZ',cZ),('Bezout_cW',cW),('known_unit_eta',eta)]:assert p==dec(certificate[key])
    # This uses no inverse of g1,g2,h,Delta34 or an unknown minor.
    T=[[(cZ,eta),(cW,eta)],[(scale(w,-1),ONE),(z,ONE)]]
    Ti=[[(z,ONE),(scale(cW,-1),eta)],[(w,ONE),(cZ,eta)]]
    for k,name in [(T,'pair_target_coefficient_change2x2'),(Ti,'pair_target_inverse2x2')]:
        for i in range(2):
            for j in range(2):
                p0,p1=map(dec,certificate[name][i][j]);assert mul(k[i][j][0],p1)==mul(p0,k[i][j][1])
    for left,right in [(T,Ti),(Ti,T)]:
        for i in range(2):
            for j in range(2):
                num,den=rtotal(rmul(left[i][k],right[k][j]) for k in range(2))
                assert num==(den if i==j else {})
    for i in range(2):
        num,den=rtotal(rmul(T[i][j],((z,w)[j],ONE)) for j in range(2));assert num==(den if i==0 else {})
    assert add(mul(T[0][0][0],T[1][1][0]),mul(T[0][1][0],T[1][0][0]),F(-1))==eta
    qpath=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json'
    quadratics=read(qpath)['quadratic_rows_sparse'];assert len(quadratics)==6
    for row in quadratics:
        assert all(not m[5] or all(not F(v) for _,v in cs) for m,cs in row)
    gate_names=['S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW',
        'S61_NONSQUARE_PRODUCT41_FULL_OPEN_20261103_v63','S61_NONSQUARE_RAW24_PRODUCT41_BRIDGE_20261103_v63_REVIEW',
        'S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW','S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63_REVIEW',
        'S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63_REVIEW','S61_NONSQUARE_ELL_COLON_RANK10_20261104_v63_REVIEW',
        'S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63_REVIEW']
    gatepaths=[ROOT/'reports'/n/'report.json' for n in gate_names]
    assert all(read(p)['all_passed'] for p in gatepaths)
    mpath=ROOT/'reports/S61_NONSQUARE_ELL_CUBIC_RANK10_20261104_v63/certificate.json'
    minor=read(mpath)['parametric_nonzero34minor']
    assert certificate['projected15_filtration']['Delta34_definition']==minor and F(minor['oldpoint_value'])!=0
    assert len(set(minor['rows']))==len(set(minor['columns']))==34
    assert certificate['coefficient_domain']=='FULLU=QQ[a,e,d,1/(adFBCJ)]'
    assert certificate['global_cubic_unit_pivot_recipe']['same_cubic_cokernel_shape']==[56,90]
    out.mkdir(parents=True,exist_ok=False)
    proof={'actual_rows':[42,47],'source_labels':[[43,20],[49,61]],'common_actual_row_normalization':'1/625000000',
        'actual_seven_heavy_columns':columns,'independent_heavy_substitution_matrix':H,
        'actual_raw_heavy_difference':[enc(p) for p in difference],
        'actual_oldX_heavy_difference':[enc(p) for p in transformed],
        'kZ':enc(kZ),'kW':enc(kW),'cZ':enc(cZ),'cW':enc(cW),'delta':enc(delta),'eta':enc(eta),
        'unit_factorization_eta':'192/125*a*F*B^2*C^3','all8_two_sided_inverse_entries_and_unit_pivot_verified':True,
        'all41_product_source_rows_have_zero_all7_H_coefficients':True,
        'Delta34_is_additional_open_not_global_unit':True,
        'global_actual_same_cokernel_shape':[56,90],
        'full_actual_cubic_rank57_on_FULLU_intersect_Delta34':True,
        'no_wholeFULLU_rank57_or_quartic239_or_U78_CI_claim':True}
    ep=out/'independent-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),ROOT/'research/notes'/(BASE+'.md'),folder/'report.json',folder/'certificate.json',pre/'report.json',pre/'certificate.json',srpath,brpath,raw,npath,qpath,mpath,*gatepaths]
    checked.update({str(p.relative_to(ROOT)):sha(p) for p in paths})
    result={'id':STEM,'all_passed':True,'actual_frozen57_source7H_and_true_frame_independently_verified':True,
        'exact_Fraction_polynomial_Bezout_eta_known_FULLU_unit_verified':True,
        'target_change_det1_inverse_and_global_actual56_to90_samecokernel_accepted':True,
        'full_actual_cubic_rank57_on_FULLU_Delta34_only':True,
        'whole_FULLU_rank57_quartic239_U78_CI_or_new_kernel_claim':False,
        'input_sha256':checked,'evidence_sha256':{ep.name:sha(ep)},
        'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
