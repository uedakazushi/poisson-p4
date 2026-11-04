#!/usr/bin/env python3
"""Own ALL22/source-open integer audit and two scheme maps; no producer/query."""
import hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63_REVIEW'
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(exist_ok=False,parents=True)
    helper=ROOT/'code/s61_nonsquare_raw_H_B_K13_d_Euclid_review_20261104_v63.py';spec=importlib.util.spec_from_file_location('own_quotient_table',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
    pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63';r=load(pb/'report.json');assert r['all_passed']
    paths={Path(__file__),helper,ROOT/'research/notes'/(ID+'.md'),pb/'report.json'}
    for rel,h in r['input_sha256'].items():assert sha(ROOT/rel)==h,rel;paths.add(ROOT/rel)
    for rel,h in r['evidence_sha256'].items():assert sha(pb/rel)==h,rel;paths.add(pb/rel)
    mp=load(pb/'monomial-map.json');g=ar.ip(mp['g13']);A=ar.OwnAlgebra(g)
    hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_H4_INPUT_20261104_v63/H4-monic-field-coefficients.json';h=load(hb)
    db=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63/result.json';dresult=load(db)
    assert ar.ip(h['g13'])==g
    a=A.make({1:1});t=A.make(ar.ip(h['e_numerator']),int(h['e_denominator']));s=A.dec(dresult['root_s'])
    literal={'one':(A.o,(0,0,0)),'a':(a,(1,0,0)),'e=t':(t,(0,1,0)),'d=s':(s,(0,0,1))}
    nodes=[];exponents=[];product_count=0
    for line in (pb/'monomial-nodes.jsonl').read_text().splitlines():
        n=json.loads(line);j=n['id'];assert j==len(nodes);v=A.dec(n['value'])
        if n['op']=='literal_accepted_input':
            expected,exp=literal[n['name']];assert v==expected
        else:
            assert n['op']=='multiply_mod_g13';i,k=n['left'],n['right'];assert i<j and k<j
            x,y=nodes[i],nodes[k];numx={i:c for i,c in enumerate(x[0])if c};numy={i:c for i,c in enumerate(y[0])if c};numv={i:c for i,c in enumerate(v[0])if c}
            diff=ar.plus({i:c*v[1]for i,c in ar.times(numx,numy).items()},{i:c*x[1]*y[1]for i,c in numv.items()},-1)
            assert diff==ar.times(g,ar.ip(n['g13_multiple_quotient'])) and A.mul(x,y)==v
            exp=tuple(p+q for p,q in zip(exponents[i],exponents[k]));product_count+=1
        nodes.append(v);exponents.append(exp)
    assert len(nodes)==312 and mp['input_node_ids']==[0,1,2,3]
    mon={tuple(m):i for m,i in mp['monomial_nodes']};assert len(mon)==len(mp['monomial_nodes'])
    for m,i in mon.items():assert exponents[i]==m
    master=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63/certificate.json';source=load(master);paths.add(master)
    Pi=source['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling'];U=source['FULLU_product']
    assert len(Pi)==22 and sum(map(len,Pi))==2963 and len(U)==49
    assert set(mon)=={tuple(m)for p in Pi+[U]for m,v in p}
    ev=load(pb/'all22-source-evaluations.json')['evaluations'];assert len(ev)==23
    values=[]
    for k,(p,record)in enumerate(zip(Pi+[U],ev)):
        assert record['source_index']==(k if k<22 else 'FULLU') and record['source_terms']==len(p)
        den=math.lcm(*(Q(v).denominator*nodes[mon[tuple(m)]][1]for m,v in p));assert den==int(record['common_numeric_denominator'])
        raw={}
        for m,c in p:
            x=nodes[mon[tuple(m)]];c=Q(c);factor=c.numerator*(den//(c.denominator*x[1]))
            raw=ar.plus(raw,{i:v*factor for i,v in enumerate(x[0])if v})
        assert raw==ar.ip(record['raw_QQ_a_linear_sum_numerator']) and max(raw,default=-1)<13
        value=A.make(raw,den);assert value==A.dec(record['value'])
        if k<22:assert not raw and record['zero']
        values.append(value)
    Uvalue=values[-1]
    # Reuse and independently check every one of the seven integer unit certs.
    nb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_NAMED_SATURATION_20261104_v63/progress.json';named=load(nb);paths.add(nb)
    assert named['status']=='clopen_U_unit_cover_retained' and not named['removals']
    un={u['tag']:u for u in named['units']};assert len(un)==7
    inverse={};unitinput={}
    for tag,u in un.items():
        x=A.dec(u['input']);p=ar.ip(u['primitive_numerator']);content=int(u['numerator_content']);S,T=ar.ip(u['S']),ar.ip(u['T']);R=int(u['R']);assert content and R
        assert {i:v*content for i,v in p.items()}=={i:v for i,v in enumerate(x[0])if v}
        assert ar.plus(ar.times(S,p),ar.times(T,g))=={0:R}
        y=A.dec(u['inverse']);assert y==A.make({i:v*x[1]for i,v in S.items()},content*R) and A.mul(x,y)==A.o
        inverse[tag]=y;unitinput[tag]=x
    def sc(x,c):return A.mul(x,A.c(c))
    t2=A.mul(t,t);F=A.add(A.add(sc(a,500),sc(t2,90)),A.add(sc(t,3),A.c(-24)))
    C=A.add(A.add(sc(a,100),sc(t2,18)),A.add(sc(t,-9),A.c(-2)))
    b1=A.add(A.add(sc(a,50),sc(t2,18)),A.add(sc(t,-9),A.c(-2)));b0=A.add(sc(t,60),A.c(-40))
    for tag,x in [('a',a),('F',F),('C',C),('b1',b1)]:assert unitinput[tag]==x
    ia,iF,iC,ib1=[inverse[tag]for tag in ['a','F','C','b1']]
    rB=A.neg(A.mul(b0,ib1));rJ=sc(A.mul(A.add(sc(t,6),A.c(4)),ia),Q(1,5))
    H=[A.uni(p)for p in h['H4_coefficients_d0_to4']]
    def evaluate(p,x):
        v=A.z
        for c in reversed(p):v=A.add(A.mul(v,x),c)
        return v
    assert H[-1]==A.o and evaluate(H,s)==A.z
    lead=sc(A.mul(A.mul(A.mul(A.mul(a,a),F),C),b1),5)
    leadinv=sc(A.mul(A.mul(A.mul(A.mul(ia,ia),iF),iC),ib1),Q(1,5));assert A.mul(lead,leadinv)==A.o
    stored=load(pb/'source-open-units.json');assert A.dec(stored['root_s'])==s and A.dec(stored['leading_unit'])==lead and A.dec(stored['leading_unit_inverse'])==leadinv
    Ui=leadinv
    for (tag,root),rec in zip([('0',A.z),('B',rB),('J',rJ)],stored['three_root_unit_recipes']):
        assert rec['root_tag']==tag and A.dec(rec['root'])==root
        q=[A.z]*4;q[3]=H[4]
        for j in [2,1,0]:q[j]=A.add(H[j+1],A.mul(root,q[j+1]))
        hr=A.add(H[0],A.mul(root,q[0]));assert A.pa(A.pm([A.neg(root),A.o],q),[hr])==H
        assert unitinput['remaining_eval_at_'+tag]==hr==A.dec(rec['H_value']) and A.poly(rec['quotient_H_by_d_minus_root'])==q
        ih=inverse['remaining_eval_at_'+tag];qs=evaluate(q,s);ir=A.neg(A.mul(qs,ih))
        assert A.dec(rec['H_value_inverse'])==ih and A.dec(rec['quotient_at_s'])==qs and A.dec(rec['inverse_s_minus_root'])==ir
        assert A.mul(A.sub(s,root),ir)==A.o;Ui=A.mul(Ui,ir)
    B=A.add(A.mul(b1,s),b0);J=A.sub(sc(A.mul(a,s),5),A.add(sc(t,6),A.c(4)));factors=[a,s,F,B,C,J]
    assert factors==[A.dec(x)for x in stored['FULLU_product_all6_factors']]
    product=A.o
    for x in factors:product=A.mul(product,x)
    assert product==Uvalue==A.dec(stored['FULLU_value']) and Ui==A.dec(stored['FULLU_inverse']) and A.mul(Uvalue,Ui)==A.o
    # G(a,t) unit is a separate consistency check; it is NOT substituted for
    # the accepted global theorem that G is already a unit modulo the actual I.
    gb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_G_BASE_GATE_20261104_v63';gc=load(gb/'certificate.json');gu=load(gb/'G-unit-certificate.json');paths.update([gb/'certificate.json',gb/'G-unit-certificate.json'])
    np=ROOT/'reports/S61_NONSQUARE_RAW_H_B_NORM_PENCIL_20261104_v63/native-certificate.json';Gsource=load(np)['final_common_curve_gcd'];paths.add(np)
    Grec=gc['records'][2];assert {tuple(m):Q(v)for m,v in Grec['source']}=={tuple(m):Q(v)for m,v in Gsource}
    ap=[A.o];tp=[A.o]
    for _ in range(max(m[0]for m,v in Gsource)):ap.append(A.mul(ap[-1],a))
    for _ in range(max(m[1]for m,v in Gsource)):tp.append(A.mul(tp[-1],t))
    Gval=A.z
    for m,c in Gsource:assert m[2]==0;Gval=A.add(Gval,sc(A.mul(ap[m[0]],tp[m[1]]),Q(c)))
    assert Gval==A.dec(Grec['value'])==A.dec(stored['G_value'])
    gp=ar.ip(gc['G_unit_primitive_numerator']);gcontent=int(gc['G_unit_numerator_content']);assert {i:v*gcontent for i,v in gp.items()}=={i:v for i,v in enumerate(Gval[0])if v}
    GS,GT=ar.ip(gu['S']),ar.ip(gu['T']);GR=int(gu['R']);assert GR and ar.plus(ar.times(GS,gp),ar.times(GT,g))=={0:GR}
    Gi=A.make({i:v*Gval[1]for i,v in GS.items()},gcontent*GR);assert Gi==A.dec(stored['G_inverse']) and A.mul(Gval,Gi)==A.o
    dependencies=['S61_NONSQUARE_RAW_H_B_G0_REMOVED_SCHEME_BOUND_20261104_v63_REVIEW',
      'S61_NONSQUARE_RAW_H_B_K13_SOURCE_BRIDGE_20261104_v63_REVIEW','S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63_REVIEW',
      'S61_NONSQUARE_RAW_H_B_MIXED_NORM16_20261104_v63_REVIEW','S61_NONSQUARE_RAW_H_B_MINORS_GCD_20261104_v63_REVIEW']
    for stem in dependencies:
        q=ROOT/'reports'/stem/'report.json';assert load(q)['all_passed'];paths.add(q)
    proof={'all_passed':True,'arithmetic':'Own fixed-table QQ algebra and direct ordinary integer polynomial identities; no author Ring/query.',
      'literal_input_and_all312_monomial_nodes_quotient_identities':True,'product_nodes':product_count,'all_monomial_exponent_labels':True,
      'all22_actual_source_polynomials_2963_terms_zero':True,'FULLU49_term_source_and_all6_original_factors':True,
      'all7_old_integer_Bezout_and_inverse_certificates':True,'three_exact_H4_synthetic_quotients_and_root_inverses':True,'FULLU_explicit_inverse_on_entire_Ka13':True,
      'G_source49_evaluation_and_integer_Bezout':True,
      'upper_map':'Let R=QQ[a,e,d,1/U], I=actualI3(B), identified with the22 normalized generators by the accepted35-minor source review. The accepted G0 exclusion gives G invertible in R/I independently of any evaluation at(a,t,s). Accepted norm/cofactor identities give H,h=Gh1,K=Gk1,M1=Gm1 in I, hence H,h1,k1,m1 in I. The accepted exact T=Ka13 source bridge forces g13(a)=0 and e=t(a). Actual g15,k20 belong to I; their H4 reductions and the accepted whole-Ka13 ideal(H4,gR,kR)=(d-s) force d=s(a). The explicit U(a,t,s) inverse verified here expresses even the localized generator1/U as a polynomial in a. Therefore alpha:Ka13→R/I is a well-defined SURJECTION.',
      'lower_map':'Every one of the22 original source generators was independently evaluated to zero, and U(a,t,s) has its explicit inverse on the entire Ka13. Thus evaluation defines a unital beta:R/I→Ka13.',
      'two_map_isomorphism':'beta∘alpha fixes a, which generates Ka13 overQQ, so the composition is the identity. Surjectivity of alpha and this left inverse make alpha an isomorphism. No field, irreducibility, reducedness, etaleness or length guess is used.',
      'actual_bad_algebra':'R/I3(B) is exactly QQ[a]/g13','exact_QQ_length':13,'actual_B_Fitting_bad_scheme_nonempty':True,
      'scope':'Only the actual7x3 B Fitting bad algebra. This is not cokerB-module length, fullcubic/Poisson rankloss, quartic239 or genericexhaustion. No geometric pointcount/etaleness/reducedness claim.',
      'quotient_arithmetic_operations':A.ops,'new_inverse_GCD_solve_GB_point_prime_jobs':0}
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    report={'id':ID,'all_passed':True,'all22_original_source_and_FULLU_inverse_independently_checked':True,'both_scheme_maps_and_isomorphism_accepted':True,
      'actual_B_bad_algebra_exact_Ka13_length13_nonempty':True,'no_full_cubic_Poisson_rank_loss_inference':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)},'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},
      'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser();pa.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(pa.parse_args().output.resolve())
