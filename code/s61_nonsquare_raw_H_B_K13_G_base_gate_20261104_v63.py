#!/usr/bin/env python3
"""ONE exact Ka13 base gate: source G=Q*L^4, no d-variable computation."""
import argparse, hashlib, json, math, resource, signal, subprocess, sys, time
from fractions import Fraction as QQ
from pathlib import Path
from s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63 import Ring, pa, pm

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_B_K13_G_BASE_GATE_20261104_v63'
sys.set_int_max_str_digits(0)
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def ser(p): return [[i,str(v)] for i,v in sorted(p.items(),reverse=True)]
def exactdiv(p,g):
    p=dict(p); q={}; degree=max(g); lc=g[degree]
    while p and max(p)>=degree:
        k=max(p)-degree; v=p[max(p)]; assert v%lc==0
        c=v//lc; q[k]=c; p=pa(p,{i+k:c*x for i,x in g.items()},-1)
    assert not p
    return q
def sp(p): return {tuple(m):QQ(v) for m,v in p}
def sm(p,q):
    z={}
    for m,c in p.items():
        for n,d in q.items():
            k=tuple(a+b for a,b in zip(m,n)); z[k]=z.get(k,QQ(0))+c*d
    z={m:c for m,c in z.items() if c}; assert len(z)<=150000
    return z

def main(out,work):
    start=time.monotonic(); resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    signal.alarm(60); out=out.resolve(); work=work.resolve()
    out.mkdir(parents=True,exist_ok=False); work.mkdir(parents=True,exist_ok=False)
    fb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_COMMON_CURVE_FACTOR_20261104_v63'
    hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_H4_INPUT_20261104_v63'
    gp=ROOT/'reports/S61_NONSQUARE_RAW_H_B_G0_REMOVED_SCHEME_BOUND_20261104_v63'
    gr=ROOT/'reports/S61_NONSQUARE_RAW_H_B_G0_REMOVED_SCHEME_BOUND_20261104_v63_REVIEW'
    br=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_SOURCE_BRIDGE_20261104_v63_REVIEW'
    er=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_E_RECONSTRUCTION_20261104_v63'
    pre=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_LENGTH_K13_PREFLIGHT_20261104_v63'
    factor=read(fb/'certificate.json'); h=read(hb/'H4-monic-field-coefficients.json')
    assert read(gp/'report.json')['all_G0_actual_Bbad_removed']
    assert read(br/'report.json')['all_passed']
    # The accepted source factorization is reused, with no new factor search.
    Qsrc=sp(factor['factors'][0]['factor']); Lsrc=sp(factor['factors'][1]['factor'])
    Gsrc=sp(factor['source_G']); assert factor['factor_constant']=='1'
    assert [x['exponent'] for x in factor['factors']]==[1,4]
    assert Lsrc=={(1,0,0):QQ(1),(0,2,0):QQ(9,50),(0,1,0):QQ(3,50),(0,0,0):QQ(-3,25)}
    assert len(Qsrc)==9 and len(Gsrc)==49 and all(m[2]==0 for p in [Qsrc,Lsrc,Gsrc] for m in p)
    prod=Qsrc
    for _ in range(4): prod=sm(prod,Lsrc)
    assert prod==Gsrc
    g={int(i):int(v) for i,v in h['g13']}; ring=Ring(g)
    tn={int(i):int(v) for i,v in h['e_numerator']}; td=int(h['e_denominator'])
    a=ring.make({1:1}); t=ring.make(tn,td); one=ring.c(1)
    A=[one]; T=[one]; tp=[{0:1}]
    for _ in range(6): A.append(ring.mul(A[-1],a))
    for _ in range(12): T.append(ring.mul(T[-1],t)); tp.append(pm(tp[-1],tn))
    assert sum(map(len,tp))<150000
    records=[]; values=[]
    for name,p in [('L',Lsrc),('Q',Qsrc),('G',Gsrc)]:
        maxe=max(m[1] for m in p); D=math.lcm(*(c.denominator for c in p.values()))
        raw={}; val=ring.c(0)
        for m,c in p.items():
            val=ring.add(val,ring.scale(ring.mul(A[m[0]],T[m[1]]),c.numerator,c.denominator))
            raw=pa(raw,{i+m[0]:v for i,v in tp[m[1]].items()},int(c*D)*td**(maxe-m[1]))
        common=D*td**maxe
        diff=pa({i:c*val[1] for i,c in raw.items()},{i:c*common for i,c in val[0].items()},-1)
        quot=exactdiv(diff,g)
        assert len(raw)+len(diff)+len(quot)<150000
        records.append(dict(name=name,source=[[list(m),str(c)] for m,c in sorted(p.items(),reverse=True)],
            common_source_denominator=str(common),raw_numerator=ser(raw),value=ring.enc(val),
            exact_g13_divisibility_quotient=ser(quot),
            identity='value_denominator*raw-common_source_denominator*value_numerator=g13*quotient',
            raw_degree=max(raw,default=-1),raw_terms=len(raw),value_terms=len(val[0])))
        values.append(val)
    L,Q,G=values
    product=Q
    for _ in range(4): product=ring.mul(product,L)
    assert product==G
    status='G_zero_on_whole_Ka13' if not G[0] else 'G_nonzero_unit_gate_pending'
    unit=None
    if G[0]:
        # Only the authorized scalar G unit gate, never a d-Euclid or field guess.
        cpp=ROOT/'code/s61_nonsquare_raw_H_B_K13_named_unit_20261104_v63.cpp'
        exe=work/'unit-check'; inp=work/'unit-input.txt'; result=out/'G-unit-certificate.json'
        content=math.gcd(*(abs(v) for v in G[0].values())); primitive={i:v//content for i,v in G[0].items()}
        with inp.open('w') as f:
            for p in [g,primitive]:
                f.write(str(len(p))+'\n')
                for i,v in sorted(p.items(),reverse=True): f.write(str(i)+' '+str(v)+'\n')
        with (out/'compile.log').open('w') as f:
            subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
        subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','G-unit',str(exe),str(inp),str(result)],check=True,timeout=max(1,57-(time.monotonic()-start)))
        unit=read(result)
        if unit['unit']:
            S={int(i):int(v) for i,v in unit['S']}; TT={int(i):int(v) for i,v in unit['T']}
            assert pa(pm(S,primitive),pm(TT,g))=={0:int(unit['R'])}
            status='G_unit_whole_Ka13_cover_retained'
        else: status='G_proper_nonunit_base_components_retained'
    cert=dict(variable_order=['a','e','d'],coefficient_algebra='Ka13=QQ[a]/g13',g13=h['g13'],
        e_numerator=h['e_numerator'],e_denominator=h['e_denominator'],records=records,
        source_G_equals_Q_L4_all49_coefficients=True,whole_Ka13_G_equals_Q_L4=True,
        L_zero=not bool(L[0]),Q_zero=not bool(Q[0]),G_zero=not bool(G[0]),status=status,
        no_d_Euclid_Bminor_Norm4_point_prime_or_GB=True,field_integer_operations=ring.operations)
    (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    empty=not bool(G[0])
    proof=dict(status=status,actual_ring='R=QQ[a,e,d,1/(adFBCJ)]',I='I_3(actual B7x3)',
        accepted_Gunit='The accepted entire Q/L exclusion gives (I,G)=R, hence G is a unit in R/I.',
        accepted_projection='H,h1,k1,m1 belong to I after G inversion; accepted source bridge identifies QQ[a,a^-1,e]/(h1,k1,m1)=Ka13 via e=t(a). Thus Ka13[d]/H4 localized at FULLU surjects onto R/I.',
        exact_whole_base_gate=('The saved scalar source division certifies G(a,t(a))=0 in Ka13.' if empty else 'The saved source division and integer Bezout certify that G(a,t(a)) is a unit in the whole Ka13 algebra.' if unit and unit['unit'] else 'The saved source division and gcd retain the proper nonunit base components.'),
        consequence_if_G_zero='The image of G in R/I is both zero and a unit. Therefore R/I=0 as a scheme, and actual B has rank3 at every geometric FULLU point.',
        whole_actual_Bbad_empty_given_accepted_gates=empty,all_FULLU_B_rank3_given_accepted_gates=empty,
        not_global_full_cubic57_quartic239_or_generic_exhaustion=True,
        if_G_nonzero='A unit retains the cover; a proper nonunit retains/splits its proper base components. No unseen component is inverted away.',
        independent_review_of_this_new_gate_pending=True)
    (out/'scheme-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    if G[0]:
        cert['G_unit_primitive_numerator']=ser(primitive)
        cert['G_unit_numerator_content']=str(content)
        cert['G_unit_certificate_identity']='S*primitive_G_numerator+T*g13=R; nonzeroQQ content and denominators are units.'
        (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'code/s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63.py',ROOT/'code/s61_nonsquare_raw_H_B_K13_named_unit_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md',
        fb/'report.json',fb/'certificate.json',hb/'report.json',hb/'H4-monic-field-coefficients.json',
        gp/'report.json',gp/'certificate.json',gr/'report.json',br/'report.json',er/'report.json',
        er/'projection-isomorphism.json',pre/'report.json',pre/'fitting-length-lemma.json']
    report=dict(id=ID,all_passed=True,status=status,L_zero=not bool(L[0]),Q_zero=not bool(Q[0]),G_zero=empty,
        all_source_and_scalar_product_division_identities_checked=True,
        actual_Bbad_empty_given_accepted_gates=empty,all_FULLU_B_rank3_given_accepted_gates=empty,
        independent_review_pending=True,no_d_Euclid_Bminor_Norm4_point_prime_or_GB=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},
        evidence_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()},seconds=time.monotonic()-start)
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file())<31*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n'); signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
    try: main(a.output,a.work)
    except Exception as e:
        if a.output.exists(): (a.output/'failure.json').write_text(json.dumps(dict(all_passed=False,error=str(e),no_auto_retry=True),indent=2)+'\n')
        raise
