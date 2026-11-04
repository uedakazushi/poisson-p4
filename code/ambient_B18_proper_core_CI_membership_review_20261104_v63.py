#!/usr/bin/env python3
"""Independent positive membership/source proof and geometric CI gates; no GB."""
import hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def integral(p):
    d=math.lcm(*(c.denominator for c in p.values()));return {m:int(c*d)for m,c in p.items()},d
def check(products,expected):
    terms=[]
    for p,q in products:
        a,da=integral(p);b,db=integral(q);terms.append((a,b,da*db))
    ep,de=integral(expected);D=math.lcm(de,*(v[2]for v in terms));out={m:-c*(D//de)for m,c in ep.items()}
    for p,q,den in terms:
        c0=D//den
        for m,c in p.items():
            for n,v in q.items():
                k=tuple(x+y for x,y in zip(m,n));t=out.get(k,0)+c*v*c0
                if t:out[k]=t
                else:out.pop(k,None)
    assert not out
def enc(p):return [[list(m),str(c)]for m,c in sorted(p.items(),reverse=True)]
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    helper=ROOT/'code/ambient_B18_proper_rational_elimination_review_20261104_v63.py';spec=importlib.util.spec_from_file_location('own_sparse_QQ_review_helpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
    add,mul,scale,power,same=ar.add,ar.mul,ar.scale,ar.power,ar.same
    pb=ROOT/'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63';r=read(pb/'report.json');assert r['all_passed']
    paths={Path(__file__),helper,ROOT/'research/notes'/(ID+'.md'),pb/'report.json'}
    for rel,h in r['input_sha256'].items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in r['evidence_sha256'].items():assert sha(pb/rel)==h,rel;paths.add(pb/rel)
    source=ROOT/'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json'
    prior=ROOT/'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63_REVIEW/report.json';assert read(prior)['all_passed'];paths.update([source,prior])
    K=[ar.decode(p['actual_core_polynomial'])for p in read(source)['remaining7_actual_core_polynomials']]
    prep=read(pb/'preparation.json');e=read(pb/'membership-and-source-cofactors.json');assert prep['variable_order']==['a','d','v2','v6']
    a,d,v,t=[ar.var(4,i)for i in range(4)];D=add(ar.one(4),mul(d,d),-1);P=add(mul(d,d),scale(ar.one(4),5))
    units=[]
    for rec in prep['unit_stripping']:
        u=scale(mul(power(D,rec['D']),power(P,rec['P'])),Q(rec.get('QQ_sign',1)));assert u.get(ar.zero(4),Q())!=0;units.append(u)
    assert len(units)==7
    q0,q1,cubic=list(map(ar.decode,prep['generators3']));rr=add(add(v,a,-1),d,-1);ss=add(a,v)
    same(cubic,mul(ss,mul(rr,rr)));q2=scale(add(q0,mul(P,mul(ss,rr))),Q(1,2))
    same(add(q0,scale(q2,2),-1),scale(mul(P,mul(ss,rr)),-1))
    same(K[0],mul(units[0],q0));same(K[1],mul(units[1],q1));same(K[2],mul(units[2],mul(rr,q2)))
    targets=list(map(ar.decode,prep['stripped_targets4']))
    for j,target in enumerate(targets):same(K[j+3],mul(units[j+3],target))
    U0,U1,U2=units[:3];M=mul(P,mul(U0,mul(U1,U2)))
    same(M,scale(mul(P,power(D,3)),-1));assert M[ar.zero(4)]==-5
    same(mul(mul(P,mul(U0,U2)),cubic),add(scale(mul(U0,K[2]),2),mul(mul(rr,U2),K[0]),-1))
    C=[[ar.decode(p)for p in row]for row in e['cofactor3_by4']];assert len(C)==3 and all(len(row)==4 for row in C)
    assert all(not ar.decode(p)for p in e['remainder4'])
    for j in range(4):check([(q0,C[0][j]),(q1,C[1][j]),(cubic,C[2][j])],targets[j])
    source_records=e['all4_source_unit_cleared_identities'];assert len(source_records)==4
    for j,record in enumerate(source_records):
        assert record['source_index']==j+3;same(ar.decode(record['multiplier']),M)
        A,B,Cj=[C[i][j]for i in range(3)];uj=units[j+3]
        expected=[mul(mul(uj,mul(U1,U2)),add(mul(P,A),mul(rr,Cj),-1)),
          mul(mul(uj,mul(P,mul(U0,U2))),B),scale(mul(mul(uj,mul(U0,U1)),Cj),2)]
        actual=list(map(ar.decode,record['cofactors_in_source_first3']))
        for actualcoef,owncoef in zip(actual,expected):same(actualcoef,owncoef)
        check([(actual[i],K[i])for i in range(3)],mul(M,K[j+3]))
    # Geometric regularity, independent of the producer's QQ gcd calls.
    assert min(map(sum,q0))==min(map(sum,q1))==2
    f={m:c for m,c in q0.items()if sum(m)==2};g={m:c for m,c in q1.items()if sum(m)==2}
    same(f,ar.decode(prep['initial_quadrics'][0]));same(g,ar.decode(prep['initial_quadrics'][1]))
    symmetric=[[Q()for j in range(4)]for i in range(4)]
    for m,c in f.items():
        idx=[i for i,k in enumerate(m)for _ in range(k)];i,j=idx
        if i==j:symmetric[i][i]+=c
        else:symmetric[i][j]+=c/2;symmetric[j][i]+=c/2
    fdet=ar.det(symmetric);assert fdet==144
    v2mon=(0,0,2,0);assert f[v2mon]==-5 and not g.get(v2mon) and g
    def substitute_one(p,i,image):
        cache={0:ar.one(4)};ans={}
        for m,c in p.items():
            k=m[i]
            if k not in cache:cache[k]=power(image,k)if image else{}
            mm=list(m);mm[i]=0;ans=add(ans,mul({tuple(mm):c},cache[k]))
        return ans
    fr=substitute_one(f,2,add(a,d));gr=substitute_one(g,2,add(a,d))
    same(fr,scale(add(add(scale(mul(a,a),-1),scale(mul(a,d),-3)),mul(d,t)),8));same(gr,scale(mul(a,add(a,t,-1)),40))
    rtests=[substitute_one(fr,0,{}),substitute_one(fr,3,a)];assert all(rtests)
    fs=substitute_one(f,2,scale(a,-1));gs=substitute_one(g,2,scale(a,-1))
    same(fs,scale(mul(d,add(add(scale(a,-1),scale(d,2)),scale(t,2))),4))
    stests=[substitute_one(gs,1,{}),substitute_one(gs,0,scale(add(d,t),2))];assert all(stests)
    same(rtests[0],scale(mul(d,t),8));same(rtests[1],add(scale(mul(a,a),-8),scale(mul(a,d),-16)))
    same(stests[0],add(scale(mul(a,a),10),scale(mul(a,t),-20)));same(stests[1],add(scale(mul(d,d),10),scale(mul(d,t),20)))
    proof={'all_passed':True,'arithmetic':'Own sparse QQ dictionaries plus integer-denominator-cleared convolutions, no producer/GB/CAS.',
      'all7_original_K_source_bindings':True,'Q0_minus2Q2_and_first3_ideal_source_identity':True,
      'all4_three_generator_membership_identities':True,'all12_source_lift_cofactor_formula_equalities':True,'all4_original_K_source_cofactor_identities':True,
      'source_common_multiplier':enc(M),'source_multiplier_origin':-5,'only_units':'D=1-d², P=d²+5, nonzeroQQ constants',
      'initial_f2':enc(f),'initial_g2':enc(g),'f2_symmetric_matrix_det':str(fdet),
      'geometric_pair_coprimality':'f2 has rank4 and is irreducible over an algebraic closure. g2 is not proportional (its v2² coefficient is0 while f2 has−5). Thus gcd(f2,g2)=1 geometrically.',
      'r_zero_nondivisibility_tests':[enc(p)for p in rtests],'s_zero_nondivisibility_tests':[enc(p)for p in stests],
      'restricted_coprimality':'At r=0, g2=40a(a−v6); f2 does not vanish identically on either factor. At s=0, f2=4d(−a+2d+2v6); g2 does not vanish identically on either factor. Thus both restrictions are geometrically coprime.',
      'initial_regular_sequence':'The two quadrics form a homogeneous CM complete intersection of dimension2. The restriction tests show heights(f2,g2,r)=(f2,g2,s)=3, so r,s are nonzerodivisors on that CM ring; sr² is too. The initial sequence has degrees(2,2,3). Filtered regular-sequence lifting gives a1-dimensional complete intersection for Q0,Q1,sr², multiplicity12, CM and Gorenstein. The four positive literal membership certificates identify the actual full seven-function local ideal with it.',
      'actual_core_dimension':1,'actual_core_multiplicity':12,'actual_core_CM_Gorenstein':True,
      'normalized_E4_argument':'In S=QQ[[u5,a,d,v2,v6]], J=(u5,K7). Multiplication by d gives0→S/J→S/dJ→S/(d)→0 since d is regular on S. The first term is CM dimension1; S/(d) is regular dimension8. H0 of both is0, so H0(E4)=0; H1 of the first term injects and is nonzero, so depth(E4)=1. Whole E4 is not thereby1-dimensional or CM.',
      'normalized_E4_H0_zero_depth':1,'no_D0_depth_H0_T0_inference':True,'new_GB_sourcejet_point_prime_or_CAS_jobs':0,'counts':ar.COUNT}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report={'id':ID,'all_passed':True,'all4_membership_all12_source_cofs_all7_actual_core_checked':True,
      'actual_core_CI_2_2_3_dimension1_CM_Gorenstein_multiplicity12':True,'normalized_E4_H0_zero_depth1':True,'D0_not_identified':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(ap.parse_args().output.resolve())
