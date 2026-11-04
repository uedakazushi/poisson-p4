#!/usr/bin/env python3
"""Independent sparse-QQ source audit of five core roots and all seven residuals."""
import hashlib,importlib.util,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    helper=ROOT/'code/ambient_B18_proper_rational_elimination_review_20261104_v63.py'
    spec=importlib.util.spec_from_file_location('own_independent_sparse_QQ',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
    add,mul,scale,power,same=ar.add,ar.mul,ar.scale,ar.power,ar.same
    pb=ROOT/'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63';pr=read(pb/'report.json');assert pr['all_passed']
    paths={Path(__file__),helper,ROOT/'research/notes'/(ID+'.md'),pb/'report.json'}
    for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in pr['evidence_sha256'].items():assert sha(pb/rel)==h,rel;paths.add(pb/rel)
    previous=ROOT/'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63_REVIEW/report.json';assert read(previous)['all_passed'];paths.add(previous)
    c=read(pb/'actual4variable7core-polynomials.json');pre=read(pb/'sequential-unit-preflight.json')
    sp=ROOT/'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63/actual12-flag-divided-polynomials.json'
    src=read(sp);paths.add(sp);H=[ar.decode(r['divided_primitive_polynomial'])for r in src['all12_actual_remaining_and_divided_records']]
    names=['a','d','v2','v3','v4','v5','v6','v7','v8'];Rindices=[0,1,2,6];Rnames=[names[i]for i in Rindices]
    assert c['source9_variable_order']==names and c['actual_core4_variable_order']==Rnames
    source_rows=[0,4,1,2,3];pivots=[4,3,5,7,8];assert c['selected5_source_H_rows']==source_rows and c['selected5_pivot_order']==[names[i]for i in pivots]
    X=[ar.var(4,i)for i in range(4)];a,d,v2,v6=X;D=add(ar.one(4),mul(d,d),-1);P=add(mul(d,d),scale(ar.one(4),5))
    assert D[ar.zero(4)]==1 and P[ar.zero(4)]==5
    roots={}
    for i in pivots:
        r=c['five_actual_rational_roots'][names[i]];assert r['equation_row']==source_rows[pivots.index(i)]
        roots[i]=(ar.decode(r['numerator']),r['D_power'],r['P_power']);assert min(r['D_power'],r['P_power'])>=0
        assert not roots[i][0].get(ar.zero(4),Q())
    def base(p):
        assert all(all(not m[i]for i in range(9)if i not in Rindices)for m in p)
        return {tuple(m[i]for i in Rindices):v for m,v in p.items()}
    def coefficient(p,i):
        ans={}
        for m,v in p.items():
            if m[i]==1:
                mm=list(m);mm[i]=0;ans[tuple(mm)]=v
        return ans
    dp={0:ar.one(4)};pp={0:ar.one(4)};cache={}
    def denpower(cache_,value,k):
        if k not in cache_:cache_[k]=power(value,k)
        return cache_[k]
    def denom(ds,ps):return mul(denpower(dp,D,ds),denpower(pp,P,ps))
    def substitute(p,allowed):
        needed={i for i in pivots if any(m[i]for m in p)};assert needed<=set(allowed)
        md=max((sum(m[i]*roots[i][1]for i in needed)for m in p),default=0)
        mp=max((sum(m[i]*roots[i][2]for i in needed)for m in p),default=0);acc={}
        for m,v in p.items():
            term={tuple(m[i]for i in Rindices):v};dd=ep=0
            for i in needed:
                if m[i]:
                    key=i,m[i]
                    if key not in cache:cache[key]=power(roots[i][0],m[i])
                    term=mul(term,cache[key]);dd+=m[i]*roots[i][1];ep+=m[i]*roots[i][2]
            term=mul(term,denom(md-dd,mp-ep));acc=add(acc,term)
        return acc,md,mp
    # Fixed original two-equation linear block, not a generic rank search.
    lc=[[base(coefficient(H[i],j))for j in [3,4]]for i in [0,4]]
    same(lc[0][0],{});same(lc[0][1],scale(P,2));same(lc[1][0],scale(D,8))
    same(lc[1][1],scale(mul(d,add(add(scale(a,-3),d),v2)),4))
    determinant=add(mul(lc[0][0],lc[1][1]),mul(lc[0][1],lc[1][0]),-1);same(determinant,scale(mul(D,P),-16));assert determinant[ar.zero(4)]==-80
    checks=[];prior=[]
    for k,(row,i)in enumerate(zip(source_rows,pivots)):
        assert max(m[i]for m in H[row])==1
        coef9=coefficient(H[row],i);coef=base(coef9);unit=c['unit_tests'][k]
        assert unit==pre['tests'][k] and unit['row']==row and unit['pivot']==names[i] and unit['degree']==1
        factor=Q(unit['constant_factor']);assert factor
        same(coef,scale(denom(unit['D_power'],unit['P_power']),factor))
        assert coef.get(ar.zero(4))==Q(unit['coefficient_at_origin'])==[10,8,80,16,192][k]
        target=ar.var(9,i);tail=add(H[row],mul(coef9,target),-1)
        assert all(all(not m[j]for j in pivots[k:])for m in tail)
        tn,td,tp=substitute(tail,prior);rn,rd,rp=roots[i]
        # Exact triangular source equality after preceding roots.
        same(add(mul(tn,denom(rd,rp)),mul(mul(coef,rn),denom(td,tp))),{})
        prior.append(i);whole,_,_=substitute(H[row],pivots);same(whole,{})
        checks.append({'source_H':row,'pivot':names[i],'unit_origin':str(coef[ar.zero(4)]),'root_terms':len(rn)})
    residual=[];K=[]
    source_rest=[i for i in range(12)if i not in source_rows];assert len(c['remaining7_actual_core_polynomials'])==7
    for row,rec in zip(source_rest,c['remaining7_actual_core_polynomials']):
        assert rec['source_H_row']==row
        N,md,mp=substitute(H[row],pivots);kpoly=ar.decode(rec['actual_core_polynomial'])
        sd=rec['D_denominator_power'];spow=rec['P_denominator_power'];qq=Q(rec['QQ_unit_content_numerator'],rec['QQ_unit_content_denominator']);assert qq
        same(mul(N,denom(sd,spow)),scale(mul(kpoly,denom(md,mp)),qq));assert not kpoly.get(ar.zero(4),Q())
        K.append(kpoly);residual.append({'source_H':row,'numerator_terms':len(N),'K_terms':len(kpoly),'source_denominator_D_P':[md,mp],'stored_denominator_D_P':[sd,spow],'QQ_unit':str(qq)})
    J=[[p.get(tuple(int(i==j)for i in range(4)),Q())for j in range(4)]for p in K]
    assert J==[[Q(v)for v in row]for row in c['origin_Jacobian7_by4']] and ar.rank(J)==c['origin_rank']==0
    proof={'all_passed':True,'arithmetic':'Own sparse Fraction polynomial arithmetic, imported only the independently reviewed reviewer helper; no producer/SymPy/CAS.',
      'source_H12_binding':True,'first2_actual_linear_block_and_det_origin_minus80':True,'five_triangular_source_unit_identities':checks,
      'all5_whole_H_source_zero_identities':True,'all7_original_rational_residual_crossproducts':residual,'known_denominator_units':'D=1-d², P=d²+5, origins1,5 only',
      'core4_origin_Jacobian_rank':0,'whole_E4':'QQ[[u1..u5,a,d,v2,v6]] / d*(u1..u5,K1..K7)',
      'divided_core':'QQ[[a,d,v2,v6]]/(K1..K7)',
      'why_whole_E4_keeps_d':'All final five roots depend only on base(a,d,v2,v6); translating pivot_i to u_i+root_i is an ambient coordinate automorphism fixing d. The checked triangular unit relations identify the selected H ideal with (u1..u5). Every other H equals its checked unit*K modulo that graph ideal, so the full H ideal is (u1..u5,K7). Multiplication of the entire ideal by the unchanged d yields the whole E4 presentation. No selected H equation is imposed alone in E4.',
      'no_D0_identification_CM_depth_GB_sourcejet_or_membership_search':True,'counts':ar.COUNT}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report={'id':ID,'all_passed':True,'five_roots_all7_residuals_source_verified':True,'whole_E4_d_times_full_ideal_retained':True,'core4_input_accepted':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser(description=__doc__);pa.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(pa.parse_args().output.resolve())
