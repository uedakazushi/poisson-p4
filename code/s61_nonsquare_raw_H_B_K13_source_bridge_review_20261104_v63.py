#!/usr/bin/env python3
"""Independent integer-polynomial replay of four fixed Ka13 source certificates.
No author import, inverse/GCD search, Gaussian solve, new norm, point or prime.
"""
import gzip, hashlib, json, math, resource, signal, sys, time
from fractions import Fraction as F
from pathlib import Path
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_B_K13_SOURCE_BRIDGE_20261104_v63_REVIEW'
COUNTS={'integer_products':0,'exact_divisions':0,'identity_checks':0}
def load(p):
    data=gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes()
    return json.loads(data)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def up(p):return {int(i):F(c)for i,c in p if F(c)}
def ip(p):return {int(i):int(c)for i,c in p if int(c)}
def clean(p):return {i:c for i,c in p.items()if c}
def plus(p,q,k=1):
    r=p.copy()
    for i,c in q.items():
        v=r.get(i,0)+k*c
        if v:r[i]=v
        else:r.pop(i,None)
    return r
def scalar(p,c):return {i:c*v for i,v in p.items()if c*v}
def shift(p,s):return {i+s:c for i,c in p.items()}
def times(p,q):
    COUNTS['integer_products']+=1;r={}
    for i,c in p.items():
        for j,v in q.items():r[i+j]=r.get(i+j,0)+c*v
    return clean(r)
def equal(p,q):
    COUNTS['identity_checks']+=1
    assert clean(p)==clean(q)
def powers(p,n):
    r=[{0:1}]
    for _ in range(n):r.append(times(r[-1],p))
    return r
def integral(p):
    d=math.lcm(*(F(v).denominator for v in p.values()))
    return {i:int(F(v)*d)for i,v in p.items()},d
def exact_division(p,g):
    """Fixed divisibility check, integer long division; no gcd/inverse search."""
    COUNTS['exact_divisions']+=1;r=p.copy();q={};dg=max(g);lc=g[dg]
    while r and max(r)>=dg:
        k=max(r)-dg;c=r[max(r)]
        assert c%lc==0
        c//=lc;q[k]=c;r=plus(r,shift(g,k),-c)
    assert not r
    equal(times(g,q),p)
    return q
def congruent(p,q,g):
    a,da=integral(p);b,db=integral(q);d=math.lcm(da,db)
    return exact_division(plus(scalar(a,d//da),scalar(b,d//db),-1),g)
def ordinary(p):
    """Source sparse(a,e,d) -> grouped ordinary univariates in a."""
    ans={}
    for m,c in p:
        a,e,d=map(int,m);assert d==0
        ans[e]=plus(ans.get(e,{}),{a:F(c)})
    return ans
def common(ps):return math.lcm(*(F(c).denominator for p in ps for m,c in p))
def evaluate_cleared(source,tp,td,max_e,D,d_degree=0):
    out={}
    for m,c in source:
        a,e,d=map(int,m)
        if d!=d_degree:continue
        assert a>=0 and e<=max_e
        n=F(c)*D;assert n.denominator==1
        out=plus(out,shift(tp[e],a),int(n)*td**(max_e-e))
    return out
def check_parent(stem,paths):
    p=ROOT/'reports'/stem;r=load(p/'report.json');assert r['all_passed']
    paths.add(p/'report.json')
    for rel,h in r.get('input_sha256',{}).items():assert sha(ROOT/rel)==h,rel
    for rel,h in r.get('evidence_sha256',{}).items():assert sha(p/rel)==h,rel;paths.add(p/rel)
    return p,r

def main(out):
    started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(exist_ok=False,parents=True);paths=set()
    base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_LENGTH_K13_PREFLIGHT_20261104_v63'
    source=load(base/'K13-source-input.json');pref=load(base/'mixed-minor-W-source-preflight.json')
    paths.update([base/'report.json',base/'K13-source-input.json',base/'mixed-minor-W-source-preflight.json'])
    g=ip(source['g13']);assert max(g)==13 and g[0]!=0 and math.gcd(*map(abs,g.values()))==1
    ownmixed=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_NORM16_20261104_v63_REVIEW/report.json'
    assert load(ownmixed)['all_passed'];paths.add(ownmixed)
    mb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_NORM16_20261104_v63'
    ng=load(mb/'native-certificate.json')['gcd']
    equal(g,ip(ng));paths.add(mb/'native-certificate.json')
    mi=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_NORM16_INPUT_20261104_v63'
    oi=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_INPUT_20261104_v63'
    new,old=load(mi/'certificate.json'),load(oi/'certificate.json');cost=load(mi/'cost-preflight.json')
    paths.update([mi/'certificate.json',mi/'cost-preflight.json',oi/'certificate.json'])
    # Bind all source coefficient reductions, not just the W/E outputs.
    origH=ordinary(new['H_h1_monic16']);assert ordinary(old['h1_e16_monic'])==origH
    coeffs=source['coefficient_arrays_by_e_low_to_high']
    for j in range(17):congruent(origH.get(j,{}),up(coeffs[0][j]),g)
    for k,inp in enumerate([old,new],1):
        for j in range(16):
            p=inp['matrix_rows_source_order'][j][0]
            a={int(m[0]):F(c)for m,c in p};assert all(m[1:]==[0,0]for m,c in p)
            congruent(a,up(coeffs[k][j]),g)
    ainv=up(source['a_inverse']);congruent(shift(ainv,1),{0:F(1)},g)
    # Original accepted adjugate column: retain ONE common QQ*a^24 normalization.
    oldcanon=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_CANONICAL_20261104_v63/native-certificate.json.gz'
    Vsource=load(oldcanon)['only16_solution_cofactors'];paths.add(oldcanon)
    content=int(pref['V_common_QQ_content']);ashift=pref['V_common_known_a_power'];assert content!=0 and ashift==24
    Vbefore=list(map(up,pref['normalized_V_before_reduction']));V=list(map(up,pref['normalized_V_modg13']))
    for i,p in enumerate(Vsource):
        assert all(m[1:]==[0,0]for m,c in p)
        equal({int(m[0]):int(c)for m,c in p},scalar(shift(Vbefore[i],ashift),content))
        congruent(Vbefore[i],V[i],g)
    del Vsource
    D=int(cost['matrix_uniform_QQ_denominator']);assert D==int(pref['source_Nm_constant_clear'])
    A=[[up(p)for p in row]for row in pref['Nm1_modg13']]
    for i in range(16):
        for j in range(16):
            p=new['matrix_rows_source_order'][i][j];assert all(m[1:]==[0,0]for m,c in p)
            congruent({int(m[0]):F(c)*D for m,c in p},A[i][j],g)
    DA=math.lcm(*(v.denominator for row in A for p in row for v in p.values()))
    DV=math.lcm(*(v.denominator for p in V for v in p.values()))
    Ni=[[{i:int(c*DA)for i,c in p.items()}for p in row]for row in A]
    Vi=[{i:int(c*DV)for i,c in p.items()}for p in V]
    assert [DA.bit_length(),DV.bit_length()]==pref['common_numeric_denominator_bits']
    wp,wr=check_parent('S61_NONSQUARE_RAW_H_B_MIXED_MINOR_W_20261104_v63',paths)
    w=load(wp/'native-certificate.json.W-checkpoint.json');wn=load(wp/'native-certificate.json')
    assert len(w['rows'])==16
    for i,r in enumerate(w['rows']):
        raw={}
        for j in range(16):raw=plus(raw,times(Ni[i][j],Vi[j]))
        equal(raw,ip(r['raw_integer_product']))
        equal(scalar(raw,g[13]**r['leading_coefficient_power']),plus(times(g,ip(r['g13_quotient'])),ip(r['remainder'])))
        assert not ip(r['remainder'])
    assert wn['all_W_components_zero'] and not wn['gcd_chain'];equal(ip(wn['final_gcd']),g)
    timings={'W_and_original_source_seconds':time.monotonic()-started}
    # Certified V15 unit and all companion EV=tV equations.
    ep,er=check_parent('S61_NONSQUARE_RAW_H_B_K13_E_RECONSTRUCTION_20261104_v63',paths)
    e=load(ep/'native-certificate.json');j=e['selected'];assert j==15 and e['pivot_unit']
    equal(Vi[j],scalar(ip(e['pivot_primitive']),int(e['pivot_content'])))
    equal(plus(times(ip(e['unit_S']),ip(e['pivot_primitive'])),times(ip(e['unit_T']),g)),{0:int(e['unit_constant'])})
    assert int(e['unit_constant'])!=0
    tn=ip(e['e_numerator']);td=int(e['e_denominator']);assert td!=0 and max(tn)<13
    H=[up(p)for p in coeffs[0]];DH=math.lcm(*(c.denominator for p in H[:16]for c in p.values()))
    Hclear=[{i:int(c*DH)for i,c in p.items()}for p in H[:16]]
    for i,r in enumerate(e['EV_rows']):
        raw=plus(scalar(Vi[i-1],DH)if i else{},times(Hclear[i],Vi[15]),-1)
        equal(raw,ip(r['raw']))
        rem=ip(r['remainder']);equal(scalar(raw,g[13]**r['L_power']),plus(times(g,ip(r['quotient'])),rem))
        equal(plus(scalar(times(tn,Vi[i]),DH*g[13]**r['L_power']),scalar(rem,td),-1),times(g,ip(r['EV_equals_tV_quotient'])))
    equal(ip(e['e_raw_numerator']),times(ip(e['EV_rows'][j]['remainder']),ip(e['unit_S'])))
    lp=e['e_reduction_L_power']
    equal(scalar(ip(e['e_raw_numerator']),g[13]**lp),plus(times(g,ip(e['e_reduction_quotient'])),ip(e['e_reduction_remainder'])))
    assert int(e['e_raw_denominator'])==DH*g[13]**(e['EV_rows'][j]['L_power']+lp)*int(e['pivot_content'])*int(e['unit_constant'])
    equal(scalar(tn,int(e['e_raw_denominator'])),scalar(ip(e['e_reduction_remainder']),td))
    timings['E_total_seconds']=time.monotonic()-started
    # Full source Laurent H4 via closed a inverse; no unverified denominator.
    hp,hr=check_parent('S61_NONSQUARE_RAW_H_B_K13_H4_INPUT_20261104_v63',paths)
    h=load(hp/'native-certificate.json');hfields=load(hp/'H4-monic-field-coefficients.json')
    hsourcep=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63/certificate.json'
    Hsource=load(hsourcep)['monic_H_laurent'];paths.add(hsourcep);assert len(Hsource)==74
    anum={i-1:-v for i,v in g.items()if i};equal(ip(h['a_inverse_numerator']),anum);assert int(h['a_inverse_denominator'])==g[0]
    equal(plus(shift(anum,1),{0:-g[0]}),scalar(g,-1))
    DHsrc=common([Hsource]);hc=DHsrc*g[0]**5*td**8;assert hc==int(h['common_source_denominator'])
    ap=powers(anum,5);tp=powers(tn,23);rawH=[{}for _ in range(5)]
    for m,c in Hsource:
        a,b,d=map(int,m);cc=F(c)*DHsrc;assert cc.denominator==1 and -5<=a<=0 and b<=8
        rawH[d]=plus(rawH[d],times(ap[-a],tp[b]),int(cc)*g[0]**(5+a)*td**(8-b))
    for i,r in enumerate(h['H4_coefficients']):
        equal(rawH[i],ip(r['raw_common_numerator']))
        rem=ip(r['remainder']);equal(scalar(rawH[i],g[13]**r['L_power']),plus(times(g,ip(r['g13_quotient'])),rem))
        equal(scalar(ip(r['numerator']),hc*g[13]**r['L_power']),scalar(rem,int(r['denominator'])))
        equal(up(hfields['H4_coefficients_d0_to4'][i]),{a:F(v,int(r['denominator']))for a,v in ip(r['numerator']).items()})
    equal(up(hfields['H4_coefficients_d0_to4'][4]),{0:F(1)})
    timings['H4_total_seconds']=time.monotonic()-started
    # m2,m3 literal original source substitution; every one of the 84 B coefficients.
    mp,mr=check_parent('S61_NONSQUARE_RAW_H_B_K13_MIXED23_SCALARS_20261104_v63',paths)
    m=load(mp/'native-certificate.json');bc=load(mp/'native-certificate.json.B21.json')
    pp=ROOT/'reports/S61_NONSQUARE_RAW_H_B_NORM_PENCIL_20261104_v63/native-certificate.json'
    bp=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MINORS_COST_20261104_v63/certificate.json';paths.update([pp,bp])
    pencil=load(pp)['mixed3_quotients_by_final'];B=load(bp)['prepared_integer7x3'];assert len(B)==7 and sum(len(p)for row in B for p in row)==345
    for r in m['scalars']:
        sourcep=pencil[r['mixed_index']-1];Dsrc=common([sourcep]);be=max(ex[1]for ex,c in sourcep)
        raw=evaluate_cleared(sourcep,tp,td,be,Dsrc);cd=Dsrc*td**be
        equal(raw,ip(r['raw_numerator']));assert cd==int(r['common_source_denominator'])
        rem=ip(r['remainder']);equal(scalar(raw,g[13]**r['L_power']),plus(times(g,ip(r['g13_quotient'])),rem))
        equal(scalar(ip(r['numerator']),cd*g[13]**r['L_power']),scalar(rem,int(r['denominator'])))
        assert not ip(r['numerator']) and r['scalar_zero']
    DB=common([p for row in B for p in row]);bterms=bnz=0;assert len(bc['entries'])==21
    for entry,p in enumerate([p for row in B for p in row]):
        assert len(bc['entries'][entry])==4
        for d,r in enumerate(bc['entries'][entry]):
            raw=evaluate_cleared(p,tp,td,4,DB,d);num=ip(r['numerator']);den=int(r['denominator']);assert den!=0 and max(num,default=-1)<13
            exact_division(plus(scalar(raw,den),scalar(num,DB*td**4),-1),g)
            bterms+=len(num);bnz+=bool(num)
    assert bterms==789 and bnz==69
    proof={'all_passed':True,'arithmetic':'Own Python integer convolution, exact integer long division and integer-denominator cross-products; no author functions imported.',
      'original_source_V_single_common_normalization_and_all256_Nm_reductions':True,'W_all16_raw_products_and_zero_remainders':True,
      'V15_unit_Bezout_and_all16_companion_eigen_identities':True,'H4_all74_original_Laurent_terms_five_coefficients_monic':True,
      'm2_m3_original_source_substitutions_zero':True,'B21_all345_original_terms_all84_images':True,'B21_nonzero_coefficients':bnz,'B21_QQ_terms':bterms,
      'exact_projection_argument':'Let T=QQ[a,a^-1,e]/(h1,k1,m1). Accepted Smith/Fitting proof gives QQlength(T)<=13. The accepted old adjugate column and new W=0 show k1(E)V=m1(E)V=0 in Ka13^16. Companion h1(E)=0 and EV=tV, together with the unit V15, imply h1(t)=k1(t)=m1(t)=0. Thus T surjects onto Ka13=QQ[a]/g13, dimension13. The length bound forces an isomorphism. No irreducibility or reducedness is used.',
      'triple_projection_exact_QQ_length':13,'monic_H4_cover_QQ_length':52,
      'scope':'Only the finite source model. Neither actual Bbad nonemptiness nor emptiness, full cubic57, full quartic239, nor cokerB-module length follows.',
      'new_solve_inverse_GCD_norm_point_prime_GB_jobs':0,'counts':COUNTS,'timings':timings}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    paths.update([Path(__file__),ROOT/'research/notes'/(ID+'.md')])
    r={'id':ID,'all_passed':True,'four_frozen_packets_independently_replayed':True,'exact_projection_length13_H4cover52_accepted':True,
      'actual_Bbad_existence_or_empty_not_claimed':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},
      'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},'seconds':time.monotonic()-started,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    signal.alarm(0);(out/'report.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:v for k,v in r.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(p.parse_args().output.resolve())
