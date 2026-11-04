#!/usr/bin/env python3
"""Degree28 product-method preflight; no degree28 column multiplication/rank."""
import argparse,gzip,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha,add,mul,scale,l1,support,ceilq
from q_new_family_degree16_invariant_probe_20261103_v63 import character
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_DEGREE28_PRODUCT_PREFLIGHT_20261104_v63'
def decode(p):return {tuple(m):Q(c)for m,c in p}
def decode_native(p):
    content=Q(p['content']);return {tuple(t[1:]):content*int(t[0])for t in p['terms']}
def normden(p):return l1(p),math.lcm(*(c.denominator for c in p.values()))if p else 1
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    hp=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63/native-certificate.json.gz'
    dp=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63/definitions.json'
    cp=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_UNIVERSAL_BASIS_20261104_v63/native-certificate.json.gz'
    qp=ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/native-certificate.json.gz'
    pp=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_SUPPORT_PREFLIGHT_20261104_v63/certificate.json'
    lp=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63/low-values19.json'
    inp=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63/native-input.txt.gz'
    ip=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    gates=[ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63/report.json']
    assert all(read(p)['all_passed']for p in gates)
    h,d,c,q,pre,low,info=map(read,[hp,dp,cp,qp,pp,lp,ip]);assert len(d['degree24_basis'])==67
    chars=character(list(map(tuple,info['Cartan_weights'])),28);assert chars[28]['invariant_dimension']==108
    (out/'character-records.json').write_text(json.dumps(chars,indent=2)+'\n')
    definitions=[dict(kind='F4_times_Inv24_basis',Inv24_basis_index=i)for i in range(67)]
    definitions +=[dict(kind='primitive12_times_p_p',primitive12_index=i,p_indices=[j,k])for i in range(4)for j in range(1,5)for k in range(j,5)]
    definitions +=[dict(kind='primitive16_times_primitive12',primitive12_index=i)for i in range(4)]
    definitions +=[dict(kind='Inv14_product',indices=[i,j])for i in range(3)for j in range(i,3)]
    assert len(definitions)==117
    formal=dict(coefficient_order='i=0..degree,j=0..degree-i: S^i T^j',
        F4=c['F4'],p1_p2_p3_p4=c['E5'][1:5],
        primitive12=[h['degree12_columns'][i]for i in h['degree12_selected_columns'][5:]],
        primitive16=h['degree16_columns'][h['degree16_selected_columns'][-1]],
        Inv14=[q['degree14_columns'][i]for i in q['degree14_selected_columns']],
        Inv24=[h['degree24_columns'][i]for i in h['degree24_selected_columns']])
    assert[len(formal[k])for k in['F4','p1_p2_p3_p4','primitive12','primitive16','Inv14','Inv24']]==[15,4,4,153,3,67]
    with gzip.open(out/'saved-formal-low-columns.json.gz','wt')as f:json.dump(formal,f,separators=(',',':'))
    # Exact low-degree<=8 entrance only, to sharpen coefficient-height bounds.
    CD=[decode(z)for z in pre['CD_quadratic32']];C,D=CD[:16],CD[16:];pencil=[{}for _ in range(5)]
    for perm in itertools.permutations(range(4)):
        sign=(-1)**sum(perm[i]>perm[j]for i in range(4)for j in range(i+1,4))
        for letters in itertools.product(range(2),repeat=4):
            p={(0,0,0,0):Q(sign)}
            for i in range(4):p=mul(p,(C,D)[letters[i]][4*i+perm[i]])
            pencil[sum(letters)]=add(pencil[sum(letters)],p)
    lows=list(map(decode_native,low));F=lows[9];assert add(scale(pencil[0],400),scale(pencil[4],1600))==mul(F,F)
    # Bound the sole new low value theta^7 H16 by its already-bound TWO-degree entrance.
    with gzip.open(inp,'rt')as f:tok=iter(f.read().split())
    assert[int(next(tok)),int(next(tok))]==[16,136]
    for _ in range(16):
        for _ in range(int(next(tok))):
            for _ in range(5):next(tok)
    mon=[(int(next(tok)),int(next(tok)))for _ in range(136)]
    for _ in range(32*136):next(tok)
    hc=[[int(next(tok))for _ in range(136)]for _ in range(512)]
    beta=[decode(z)for z in pre['beta16_sparse']];quad=[mul(beta[i],beta[j])for i,j in mon];hm=[]
    for row in hc:
        p={}
        for z,t in zip(row,quad):
            if z:p=add(p,scale(t,Q(z,2)))
        hm.append(l1(p))
    theta_norm=max(hm[:256]);H_norm=max(hm[256:]);H8_norm=16*(16*theta_norm)**7*(16*H_norm);H8_den=2**8*54000**16
    P=[normden(p)for p in pencil];ff=normden(F);pr=[normden(lows[i])for i in[15,16,17,18]];hh=(H8_norm,H8_den);qq=[normden(lows[i])for i in range(3)]
    def product(xs):
        n=Q(1);de=1
        for a,b in xs:n*=a;de*=b
        return n,de
    def spectral(e):return product([ff]*e[0]+[P[i]for i in range(1,5)for _ in range(e[i])])
    def inv24(z):
        k=z['kind']
        if k=='spectral':return spectral(z['F4_p1_p2_p3_p4_exponents'])
        if k=='primitive12_times_spectral12':
            x=z['spectral'];return product([pr[z['primitive']],product([ff]*3)if x==0 else product([ff,P[x]])])
        if k=='primitive12_product':return product([pr[i]for i in z['indices']])
        assert k=='primitive16_times_Inv8';return product([hh,product([ff,ff])if z['index']==0 else P[z['index']]])
    bound=[]
    for z in definitions:
        k=z['kind']
        if k=='F4_times_Inv24_basis':a=product([ff,inv24(d['degree24_basis'][z['Inv24_basis_index']])])
        elif k=='primitive12_times_p_p':a=product([pr[z['primitive12_index']],*[P[j]for j in z['p_indices']]])
        elif k=='primitive16_times_primitive12':a=product([hh,pr[z['primitive12_index']]])
        else:a=product([qq[i]for i in z['indices']])
        n,de=a;bound.append(dict(denominator=str(de),coefficient_bits=ceilq(n*de).bit_length()))
    support28=len(support(28,(0,0)));assert support28==3696
    topbits=sorted([z['coefficient_bits']for z in bound],reverse=True)[:108]
    binary=sum(support28*((b+7)//8+16)for b in topbits);decimal=sum(support28*(math.ceil(b*math.log10(2))+32)for b in topbits)
    cert=dict(character_dimension28=108,candidate_count=117,candidates=definitions,degree24_basis_definitions=d['degree24_basis'],
        source_coordinate_binding='new xi0,r0,s0,q0,B0 from v0; polynomial H,L,s0,q0 after Sp4 finite gauge; restorexi0^14',
        source_support_bound3696=True,source_support_by_L={str(b):sum(m[1]==b for m in support(28,(0,0)))for b in range(10)},
        new_low_value_only='tr(theta^7 H16), degree16; ambient gradients do not need extending',new_low_value_support_bound=594,
        actual_pencil_low8_term_counts=list(map(len,pencil)),actual_low_values28_only_not_expanded=True,
        all117_height_denominator_bounds=bound,worst_selected108_binary_bytes_bound=binary,worst_selected108_decimal_bytes_bound=decimal,
        formal_support435=True,formal_117_coefficients_upper_bound=117*435,formal_rank_gate108_required=True,
        no_degree28_polynomial_column_or_rank_or_source_point=True,
        proposed_universal_job='ONE same saved formal2plane 117 coefficient products; rank108 gate and character108; stop if deficient, no candidate growth',
        proposed_source_job_after_universal_gate='Separate budget only: extend direct Hessian columns from5 to7, build selected108 exactQQ values, stream gzip with hard32MiB; ONE GF32003 coefficient rank; no source point',
        warning='No assertion that products span Inv28 until coefficient rank108. Worstcase uncompressed height bound may exceed32MiB; runtime/output caps remain real stopping conditions.')
    (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',hp,dp,cp,qp,pp,lp,inp,ip,*gates]
    report=dict(id=ID,all_passed=True,preflight_only=True,degree28_invariant_dimension=108,proposed_candidate_count=117,
        source_support_bound=3696,universal_or_source28_rank_not_run=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],source_terms=3696,maxbits=max(topbits),binary_bytes=binary,decimal_bytes=decimal)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
