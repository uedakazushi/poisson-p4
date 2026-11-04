#!/usr/bin/env python3
"""Strip only the six certified linear P80 units from the saved exact minor."""
import argparse,collections,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from s61_nonsquare_p_branch_chain_review_20261104_v63 import ROOT,read,sha
STEM='S61_NONSQUARE_THETA_P80_CUBIC_MINOR_UNITS_20261104_v63'
def enc(p):return [[i,str(c)] for i,c in sorted(p.items()) if c]
def mul(p,q):
    r=collections.defaultdict(F)
    for i,c in p.items():
        for j,d in q.items():r[i+j]+=c*d
    return {i:c for i,c in r.items() if c}
def linpower(a,b,n):return {i:F(math.comb(n,i)*a**i*b**(n-i)) for i in range(n+1) if a**i*b**(n-i)}
def divlin(p,a,b):
    r=dict(p);q={}
    for i in range(max(p),0,-1):
        c=r.pop(i,F(0))/a
        if c:q[i-1]=c;r[i-1]=r.get(i-1,F(0))-b*c
    return q,r.get(0,F(0))
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports/S61_NONSQUARE_THETA_P80_CUBIC_GENERIC_20261104_v63';report=read(base/'report.json');assert report['all_passed'] and report['full_cubic57_generic_P80_proved']
    for n,h in report['evidence_sha256'].items():assert sha(base/n)==h
    result=read(base/'native-rank.json');source={int(i):F(c) for i,c in result['nonzero_selected_minor']};assert max(source)==447 and len(source)==152 and all(c.denominator==1 for c in source.values())
    chartdir=ROOT/'reports/S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63';chart=read(chartdir/'certificate.json');cr=read(chartdir/'report.json');assert sha(chartdir/'certificate.json')==cr['certificate_sha256']
    curve=next(p for p in chart['curve_profiles'] if p['m']==80)
    factors=[('zeta',1,0),('zeta_plus20',1,20),('7zeta_plus300',7,300),('7zeta_plus80',7,80),('17zeta_plus216',17,216),('7zeta_plus60',7,60)]
    saved=[{int(m[0]):F(c) for m,c in p if F(c)} for p in curve['exact_univariate_open_factors']]
    assert saved==[{i:c for i,c in [(1,F(a)),(0,F(b))] if c} for name,a,b in factors]
    content=math.gcd(*(abs(c.numerator) for c in source.values()));current={i:c/content for i,c in source.items()};records=[];reconstruct={0:F(content)}
    for name,a,b in factors:
        before=current;n=0
        while max(current)>0:
            q,remainder=divlin(current,a,b)
            if remainder:break
            current=q;n+=1
        power=linpower(a,b,n);assert mul(power,current)==before;reconstruct=mul(reconstruct,power)
        q,remainder=divlin(current,a,b) if max(current)>0 else ({},current.get(0,F(0)))
        assert remainder
        records.append({'certified_unit':name,'linear_a_b':[a,b],'power':n,'whole_power_polynomial':enc(power),'one_final_quotient':enc(current),'next_linear_division_nonzero_remainder':str(remainder)})
    assert mul(reconstruct,current)==source
    constant=set(current)=={0} and bool(current[0]);evidence={'source_exact_minor_degree':447,'source_exact_minor_terms':152,'positive_QQ_content':str(content),'all6_known_unit_power_divisions':records,'remaining_polynomial':enc(current),'remaining_nonzero_QQ_constant':constant,'complete_source_reconstruction_exact':True,'no_unknown_factorization_or_inverse':True}
    ep=out/'exact-known-unit-strip.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',base/'report.json',base/'native-rank.json',chartdir/'report.json',chartdir/'certificate.json']
    r={'id':STEM,'all_passed':True,'frozen':True,'independent_review_pending':True,'only_saved_minor_and_certified_P80_units_used':True,'all_saved_power_divisions_and_source_reconstruction_exact':True,'remaining_degree':max(current),'remaining_terms':len(current),'remaining_nonzero_QQ_constant':constant,'ALL_P80_full_cubic57_proved_primary':constant,'no_P100_newrank_point_prime_or_unknown_factor':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'all_passed':True,'remaining_degree':max(current),'remaining_terms':len(current),'powers':{p['certified_unit']:p['power'] for p in records},'ALL_P80_full57':constant,'seconds':r['seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
