#!/usr/bin/env python3
"""Fixed QQ coefficient identities and saved L Euclid traces; no new division/search."""
import argparse, hashlib, json, resource, signal, time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_CURVE_COMPONENTS_20261104_v63_REVIEW'
Z=(0,0,0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(p):return {tuple(m):F(c) for m,c in p if F(c)}
def add(a,b):
    c=dict(a)
    for m,v in b.items():
        c[m]=c.get(m,F(0))+v
        if not c[m]:del c[m]
    return c
def scale(a,c):return {m:c*v for m,v in a.items() if c*v}
def mul(a,b):
    c={}
    for m,x in a.items():
        for n,y in b.items():
            k=tuple(i+j for i,j in zip(m,n));c[k]=c.get(k,F(0))+x*y
    return {m:v for m,v in c.items() if v}
def power(a,n):
    o={Z:F(1)}
    for _ in range(n):o=mul(o,a)
    return o
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);hashes={};base=ROOT/'reports'
    def bind(p):hashes[str(p.relative_to(ROOT))]=sha(p)
    def audit(p):
        d=read(p);bind(p)
        for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n;hashes[n]=h
        for n,h in d.get('evidence_sha256',{}).items():assert sha(p.parent/n)==h,n;bind(p.parent/n)
        return d
    prefix='S61_NONSQUARE_RAW_H_B_';suffix='_20261104_v63'
    names=['COMMON_CURVE_FACTOR','CURVE_COMPONENTS_INPUT','COMPONENT_COMMON_ROOT_L']
    reports=[audit(base/(prefix+n+suffix)/'report.json') for n in names]
    assert reports[0]['independent_factorization_and_reconstruction'] and reports[2]['all_whole_field_identities_checked']
    for name in ['MINORS_GCD_20261104_v63_REVIEW','PAIR_NORM_FINITE112_20261104_v63_REVIEW']:
        r=audit(base/(prefix+name)/'report.json');assert r['all_passed']
    factor,components,Lcert=[read(base/(prefix+n+suffix)/'certificate.json') for n in names]
    G=poly(read(base/(prefix+'NORM_PENCIL'+suffix)/'native-certificate.json')['final_common_curve_gcd'])
    assert G==poly(factor['source_G']) and len(G)==49 and F(factor['factor_constant'])==1
    factors=[(poly(r['factor']),r['exponent']) for r in factor['factors']]
    assert [n for p,n in factors]==[1,4]
    Q,L=factors[0][0],factors[1][0];assert mul(Q,power(L,4))==G
    assert L=={(1,0,0):F(1),(0,2,0):F(9,50),(0,1,0):F(3,50),Z:F(-3,25)}
    # The quadratic discriminant is squarefree/nonconstant in e, so Q is
    # irreducible over QQ(e), even after extending the constant field to C.
    qa={m:c for m,c in Q.items() if m[0]==2};assert qa=={(2,0,0):F(1)}
    q1={(0,m[1],0):c for m,c in Q.items() if m[0]==1};q0={m:c for m,c in Q.items() if m[0]==0}
    disc=add(mul(q1,q1),scale(q0,-4))
    assert disc=={(0,2,0):F(32481,810000),(0,1,0):F(30546,810000),Z:F(8881,810000)}
    assert disc[(0,1,0)]**2-4*disc[(0,2,0)]*disc[Z]!=0
    hcert=read(base/(prefix+'PAIR_MONIC_COMPLETION'+suffix)/'certificate.json')
    pair=read(base/(prefix+'COPRIME_PAIR_PREFLIGHT'+suffix)/'certificate.json')
    kcert=read(base/(prefix+'ADDITIONAL_MINOR_NORM_INPUT'+suffix)/'certificate.json')
    minors=read(base/(prefix+'35_MINORS'+suffix)/'certificate.json')
    pp=pair['selected_pair_source'];assert pp['indices']==[10,15]
    f,g=list(map(poly,pp['polynomials']));gens=list(map(poly,minors['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']))
    assert f==gens[10] and g==gens[15]
    k=poly(kcert['source_k']);assert kcert['selected_generator_index']==20 and k==gens[20]
    H=poly(hcert['monic_H_laurent']);Hn={(m[0]+5,m[1],m[2]):c for m,c in H.items()};assert all(m[0]>=0 for m in Hn)
    assert {m:c for m,c in Hn.items() if m[2]==4}=={(5,0,4):F(1)}
    # Membership and monicity use only a!=0; this also controls finite
    # exceptional fibers without imposing any generic field LC as a unit.
    a={(1,0,0):F(1)};e={(0,1,0):F(1)};d={(0,0,1):F(1)};one={Z:F(1)}
    hq=poly(hcert['q'])
    assert scale(Hn,-11)==scale(add(mul(a,f),scale(mul(hq,g),-1)),9)
    sources={'a5H':Hn,'g':g,'k':k};records={r['factor_name']:r for r in components['component_records']}
    for name,Fc in [('Q',Q),('L',L)]:
        rec=records[name];assert poly(rec['factor'])==Fc
        n=2 if name=='Q' else 1
        for p in rec['polynomials']:
            q,rem=poly(p['quotient']),poly(p['remainder'])
            assert add(mul(q,Fc),rem)==sources[p['source_label']]
            assert all(m[0]<n for m in rem)
    assert Lcert['factor_source']==records['L']['factor']
    assert Lcert['source_reduced_polynomials']==[p['remainder'] for p in records['L']['polynomials']]
    # Rational functions are pairs of UNREDUCED QQ[e] polynomial dictionaries.
    # Equality uses cross multiplication; no rational-field GCD/Euclid is run.
    ez=lambda p:{(int(i),):F(c) for i,c in p if F(c)}
    oz={(0,):F(1)};zero=({},oz);unit=(oz,oz);denominators=set()
    def ra(x,y):return (add(mul(x[0],y[1]),mul(y[0],x[1])),mul(x[1],y[1]))
    def rn(x):return (scale(x[0],-1),x[1])
    def rm(x,y):return (mul(x[0],y[0]),mul(x[1],y[1]))
    def req(x,y):return mul(x[0],y[1])==mul(y[0],x[1])
    def coef(data):
        val=zero
        for r in data:
            assert r['a_power']==0
            p=r['fraction'];nu,de=ez(p['numerator_e']),ez(p['denominator_e']);assert de
            deg=max(m[0] for m in de)
            if deg>0:
                lc=de[(deg,)];denominators.add(tuple(sorted((i[0],str(v/lc)) for i,v in de.items())))
            val=ra(val,(nu,de))
        return val
    def dp(data):return {r['d_power']:coef(r['coefficient']) for r in data}
    def pa(x,y):
        z=dict(x)
        for i,v in y.items():z[i]=ra(z.get(i,zero),v)
        return {i:v for i,v in z.items() if v[0]}
    def pn(x):return {i:rn(v) for i,v in x.items()}
    def pm(x,y):
        z={}
        for i,v in x.items():
            for j,w in y.items():z[i+j]=ra(z.get(i+j,zero),rm(v,w))
        return {i:v for i,v in z.items() if v[0]}
    def peq(x,y):return all(req(x.get(i,zero),y.get(i,zero)) for i in x.keys()|y.keys())
    def poly_to_dp(p):
        out={}
        for m,c in p.items():
            assert m[0]==0
            out[m[2]]=ra(out.get(m[2],zero),({(m[1],):c},oz))
        return out
    P=[poly_to_dp(poly(p)) for p in Lcert['source_reduced_polynomials']]
    final=Lcert['final'];GG=dp(final['common_gcd']);qq3=list(map(dp,final['all3_source_quotients_by_gcd']));bez=list(map(dp,final['Bezout3']))
    assert all(peq(pm(q,GG),p) for p,q in zip(P,qq3))
    assert peq(pa(pa(pm(bez[0],P[0]),pm(bez[1],P[1])),pm(bez[2],P[2])),GG)
    traces=Lcert['exact_Euclid_traces'];last=[];steps=0
    for t,(r0,r1) in zip(traces,[(P[0],P[1]),(dp(traces[0]['gcd_monic']),P[2])]):
        original=(r0,r1)
        for s in t['steps']:
            lc=coef(s['divisor_LC']);lci=coef(s['divisor_LC_inverse']);assert req(rm(lc,lci),unit)
            assert max(r1)==s['divisor_degree'] and req(lc,r1[max(r1)])
            quotient,remainder=dp(s['quotient']),dp(s['remainder'])
            assert peq(r0,pa(pm(quotient,r1),remainder)) and (not remainder or max(remainder)<max(r1))
            r0,r1=r1,remainder;steps+=1
        assert not r1
        gcd=dp(t['gcd_monic']);B=list(map(dp,t['source2_Bezout']))
        assert peq(pa(pm(B[0],original[0]),pm(B[1],original[1])),gcd)
        assert req(gcd[max(gcd)],unit)
        # The last nonzero remainder is a scalar times this monic gcd.
        assert peq(r0,pm({0:r0[max(r0)]},gcd));last.append(gcd)
    assert peq(last[-1],GG)
    root={1:unit,0:({(0,):F(20)},{(1,):F(3),(0,):F(-2)})}
    assert peq(GG,pm(root,root))
    # Original U and B reduce along a=-(9e²+3e-6)/50 by direct substitution.
    aval=scale({m:c for m,c in L.items() if m[0]==0},-1)
    def onL(p):
        out={}
        for m,c in p.items():out=add(out,mul({(0,m[1],m[2]):c},power(aval,m[0])))
        return out
    U=dp(final['U_reduced_on_component']);assert peq(U,poly_to_dp(onL(poly(minors['FULLU_product']))))
    B=add(mul(add(add(scale(a,50),scale(power(e,2),18)),add(scale(e,-9),scale(one,-2))),d),add(scale(e,60),scale(one,-40)))
    em=add(scale(e,3),scale(one,-2));assert onL(B)==mul(em,add(mul(em,d),scale(one,20)))
    rest=GG
    for s in Lcert['U_support_removal_steps']:
        removed,nxt,uq=map(dp,[s['removed_generic_gcd_with_U'],s['remaining_after'],s['U_quotient']])
        assert peq(rest,pm(removed,nxt)) and peq(U,pm(removed,uq));rest=nxt
    assert peq(rest,{0:unit}) and peq(rest,dp(final['remaining_common_factor_after_U_support']))
    listed={tuple(sorted((int(i),str(F(v))) for i,v in p)) for p in Lcert['finite_e_exception_polynomials']}
    assert denominators==listed and len(listed)==13 and sum(max(i for i,c in p) for p in listed)==65
    proof=dict(exact_G_equals_Q_times_L4=True,source_G49_bound_to_frozen_norm_pencil=True,Q_irreducible_by_squarefree_discriminant=True,
        all6_monic_a_division_quotient_remainder_identities=True,source_g_gen15_source_k_gen20_and_H_membership_verified=True,
        L_all_saved_Euclid_steps=steps,L_final_Bezout3_and_source_divisibilities=True,L_common_gcd='(d+20/(3e-2))^2',
        L_root_lies_on_actual_B_zero=True,all_U_support_removal_identities=True,
        all13_nonunit_denominator_exceptions_retained=True,exception_degree_sum_upper=65,
        L_FULLU_common_root_locus_is_finite_or_empty=True,Q_common_root_result_not_audited=True,
        old_finite_h_over_G_k_over_G_residual_preserved=True,
        no_all_FULLU_rank3_H7_cubic57_quartic239_or_exhaustion=True,no_new_Euclid_GCD_factor_search_GB_point_prime_or_author_worker=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/f'research/notes/{ID}.md'):bind(p)
    report=dict(id=ID,all_passed=True,reviewed_packets=[prefix+n+suffix for n in names],input_sha256=hashes,
        evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Exact G factor, all six component reductions, and generic L common roots only on B=0 outside finite denominator exceptions. Q common-root analysis and old finite residual remain.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
