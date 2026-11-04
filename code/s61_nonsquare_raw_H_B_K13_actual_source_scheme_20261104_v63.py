#!/usr/bin/env python3
"""One ALL22/source-open Ka13 gate, with literal monomial product certificates."""
import argparse,hashlib,json,math,resource,signal,sys,time
from fractions import Fraction as QQ
from pathlib import Path
from s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63 import Ring,Stop,pa,pm
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_B_K13_ACTUAL_SOURCE_SCHEME_20261104_v63'
sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ip(p):return {int(i):int(v)for i,v in p}
def ser(p):return [[i,str(v)]for i,v in sorted(p.items(),reverse=True)]
def exactdiv(p,g):
    p=dict(p);q={};deg=max(g);lc=g[deg]
    while p and max(p)>=deg:
        k=max(p)-deg;v=p[max(p)];assert v%lc==0;c=v//lc;q[k]=c
        p=pa(p,{i+k:c*v for i,v in g.items()},-1)
    assert not p;return q
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s total cap')));signal.alarm(60)
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    mb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63'
    hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_H4_INPUT_20261104_v63'
    db=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63'
    ub=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_NAMED_SATURATION_20261104_v63'
    gb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_G_BASE_GATE_20261104_v63'
    assert sha(db/'report.json')=='004096adc77ddee31067f3f6bf58e804fc6af2057a7552d43462c88b2bb55f45'
    assert sha(mb/'certificate.json')=='c4de0f35c9e30b771f7aa661b8ea74ce8470967f0e5fec63f546f19bf46de767'
    assert read(db/'report.json')['all_passed']
    d=read(db/'result.json');assert d['gcd_degree']==1
    h=read(hb/'H4-monic-field-coefficients.json');g=ip(h['g13']);ring=Ring(g);one=ring.c(1);zero=ring.c(0)
    a=ring.make({1:1});t=ring.make(ip(h['e_numerator']),int(h['e_denominator']));s=ring.dec(d['root_s'])
    source=read(mb/'certificate.json');ps=source['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling'];U=source['FULLU_product']
    assert len(ps)==22 and sum(map(len,ps))==2963 and len(U)==49
    S=set(tuple(m)for p in ps for m,v in p);SX=S|set(tuple(m)for m,v in U);pairs=set((i,j)for i,j,k in SX)
    degrees=[max(m[j]for m in SX)for j in range(3)];assert degrees==[7,11,7]
    pre=dict(source_generators=22,source_terms=2963,source_union_monomials=len(S),FULLU_source_terms=49,
        degrees_a_e_d=degrees,extended_monomials_including_U=len(SX),extended_ae_pairs=len(pairs),
        power_multiplications_upper=sum(degrees),ae_and_d_products_upper=len(pairs)+len(SX),
        saved_field_coefficient_slots_upper=13*(len(SX)+len(pairs)+sum(degrees)+4),
        method='Every monomial once; certified products modg13. Each22 evaluation is then a degree<=12 QQ[a] linear cancellation, never unreduced t^11*s^7.',
        hard_seconds=60,hard_AS_bytes=2000000000,hard_artifact_bytes=32*1024*1024,hard_explicit_terms=150000,
        no_new_inverse_search_factor_solver_point_prime_GB_or_automatic_next_job=True)
    (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
    vals=[];meta=[];termcount=0;maxbits=0;stage='monomial_nodes';evals=[];nodespath=out/'monomial-nodes.jsonl'
    def guard(extra=0):
        n=sum(p.stat().st_size for p in out.iterdir()if p.is_file())
        if n+extra+2*1024*1024>32*1024*1024:raise Stop('artifact_size_guard',dict(existing=n,extra=extra))
        if termcount>150000:raise Stop('explicit_support_guard',dict(terms=termcount))
    def checkpoint(status,extra=None):
        p=dict(status=status,stage=stage,node_count=len(vals),saved_explicit_terms=termcount,max_coefficient_bits=maxbits,
            source_evaluations_done=len(evals),extra=extra,seconds=time.monotonic()-start)
        (out/'checkpoint.json').write_text(json.dumps(p,indent=2)+'\n')
    def node(value,record):
        nonlocal termcount,maxbits
        j=len(vals);rec=dict(id=j,value=ring.enc(value),**record)
        line=json.dumps(rec,separators=(',',':'))+'\n';guard(len(line.encode()))
        with nodespath.open('a')as f:f.write(line)
        vals.append(value);meta.append(rec);termcount+=len(value[0])+len(record.get('g13_multiple_quotient',[]))
        maxbits=max([maxbits,value[1].bit_length()]+[abs(v).bit_length()for v in value[0].values()]);checkpoint('node_saved');return j
    def mulnode(i,j,tag):
        x,y=vals[i],vals[j];z=ring.mul(x,y)
        diff=pa({k:v*z[1]for k,v in pm(x[0],y[0]).items()},{k:v*x[1]*y[1]for k,v in z[0].items()},-1)
        q=exactdiv(diff,g)
        return node(z,dict(op='multiply_mod_g13',left=i,right=j,tag=tag,g13_multiple_quotient=ser(q),
            identity='den_out*num_left*num_right-den_left*den_right*num_out=g13*quotient'))
    status='pending';unitdata=None;allzero=False;Uunit=False
    try:
        ids=[node(v,dict(op='literal_accepted_input',name=n))for v,n in[(one,'one'),(a,'a'),(t,'e=t'),(s,'d=s')]]
        powers=[]
        for name,base,n in zip(['a','e','d'],ids[1:],degrees):
            p=[ids[0],base]
            for k in range(2,n+1):p.append(mulnode(p[-1],base,name+'_power_'+str(k)))
            powers.append(p)
        ae={}
        for i,j in sorted(pairs):
            ae[i,j]=powers[1][j]if i==0 else powers[0][i]if j==0 else mulnode(powers[0][i],powers[1][j],'ae_'+str((i,j)))
        mon={}
        for i,j,k in sorted(SX):
            mon[i,j,k]=ae[i,j]if k==0 else powers[2][k]if i==j==0 else mulnode(ae[i,j],powers[2][k],'aed_'+str((i,j,k)))
        (out/'monomial-map.json').write_text(json.dumps(dict(g13=ser(g),input_node_ids=ids,
            powers=powers,monomial_nodes=[[list(m),j]for m,j in sorted(mon.items())]),indent=2)+'\n');guard()
        stage='all22_original_source_evaluations'
        for index,p in enumerate(ps+[U]):
            D=math.lcm(*(QQ(v).denominator*vals[mon[tuple(m)]][1]for m,v in p));raw={}
            for m,v in p:
                q=QQ(v);x=vals[mon[tuple(m)]];scale=q.numerator*(D//(q.denominator*x[1]))
                raw=pa(raw,x[0],scale)
            value=ring.make(raw,D)
            rec=dict(source_index=index if index<22 else 'FULLU',source_terms=len(p),common_numeric_denominator=str(D),
                raw_QQ_a_linear_sum_numerator=ser(raw),value=ring.enc(value),zero=not bool(raw),
                identity='sum source_QQ_coefficient*saved_monomial_value = raw_sum/common_numeric_denominator; raw_degree<=12')
            evals.append(rec);checkpoint('source_sum_checked',dict(index=index,zero=not bool(raw)))
            if index<22:assert not raw
            else:Uvalue=value
        allzero=all(r['zero']for r in evals[:22])
        (out/'all22-source-evaluations.json').write_text(json.dumps(dict(all22_zero=allzero,source=str((mb/'certificate.json').relative_to(ROOT)),evaluations=evals),indent=2)+'\n');guard()
        stage='FULLU_and_G_unit_reuse'
        named=read(ub/'progress.json');assert named['status']=='clopen_U_unit_cover_retained'and not named['removals']
        units={x['tag']:x for x in named['units']};checked=[]
        def oldinv(tag):
            r=units[tag];p=ip(r['primitive_numerator']);S0,T0=ip(r['S']),ip(r['T']);R=int(r['R'])
            assert pa(pm(S0,p),pm(T0,g))=={0:R}
            x=ring.dec(r['input']);y=ring.dec(r['inverse']);assert ring.mul(x,y)==one
            checked.append(dict(tag=tag,source_unit_index=named['units'].index(r),all_integer_Bezout_coefficients_and_inverse_checked=True));return y
        ia,iF,iC,ib1=[oldinv(tag)for tag in ['a','F','C','b1']]
        t2=ring.mul(t,t);F=ring.add(ring.add(ring.scale(a,500),ring.scale(t2,90)),ring.add(ring.scale(t,3),ring.c(-24)))
        C=ring.add(ring.add(ring.scale(a,100),ring.scale(t2,18)),ring.add(ring.scale(t,-9),ring.c(-2)))
        b1=ring.add(ring.add(ring.scale(a,50),ring.scale(t2,18)),ring.add(ring.scale(t,-9),ring.c(-2)));b0=ring.add(ring.scale(t,60),ring.c(-40))
        params={n:ring.dec(v)for n,v in named['parameters'].items()}
        assert all(params[n]==v for n,v in [('a',a),('t',t),('F',F),('C',C),('b1',b1),('b0',b0)])
        rB=ring.neg(ring.mul(b0,ib1));rJ=ring.scale(ring.mul(ring.add(ring.scale(t,6),ring.c(4)),ia),1,5)
        assert rB==params['rB']and rJ==params['rJ']
        H=[]
        for p in h['H4_coefficients_d0_to4']:
            D=math.lcm(*(QQ(v).denominator for i,v in p));H.append(ring.make({int(i):int(QQ(v)*D)for i,v in p},D))
        assert H[4]==one and ring.eval(H,s)==zero
        lead=params['U_leading_unit'];leadinv=ring.scale(ring.mul(ring.mul(ring.mul(ring.mul(ia,ia),iF),iC),ib1),1,5);assert ring.mul(lead,leadinv)==one
        Ui=leadinv;rootcerts=[]
        for tag,r in [('0',zero),('B',rB),('J',rJ)]:
            ih=oldinv('remaining_eval_at_'+tag);q=[zero]*4;q[3]=H[4]
            for j in [2,1,0]:q[j]=ring.add(H[j+1],ring.mul(r,q[j+1]))
            hr=ring.add(H[0],ring.mul(r,q[0]));assert hr==ring.dec(units['remaining_eval_at_'+tag]['input'])
            qs=ring.eval(q,s);inv=ring.neg(ring.mul(qs,ih));assert ring.mul(ring.sub(s,r),inv)==one
            Ui=ring.mul(Ui,inv)
            rootcerts.append(dict(root_tag=tag,root=ring.enc(r),H_value=ring.enc(hr),H_value_inverse=ring.enc(ih),
                quotient_H_by_d_minus_root=[ring.enc(v)for v in q],quotient_at_s=ring.enc(qs),inverse_s_minus_root=ring.enc(inv)))
        B=ring.add(ring.mul(b1,s),b0);J=ring.sub(ring.scale(ring.mul(a,s),5),ring.add(ring.scale(t,6),ring.c(4)))
        factors=[a,s,F,B,C,J];product=one
        for x in factors:product=ring.mul(product,x)
        assert product==Uvalue and ring.mul(Uvalue,Ui)==one;Uunit=True
        Gcert=read(gb/'certificate.json');Gc=read(gb/'G-unit-certificate.json');p=ip(Gcert['G_unit_primitive_numerator'])
        assert pa(pm(ip(Gc['S']),p),pm(ip(Gc['T']),g))=={0:int(Gc['R'])}
        Gvalue=ring.dec(Gcert['records'][2]['value']);Gi=ring.make({i:v*Gvalue[1]for i,v in ip(Gc['S']).items()},int(Gcert['G_unit_numerator_content'])*int(Gc['R']))
        assert ring.mul(Gvalue,Gi)==one
        unitdata=dict(all7_named_certificates_reused_and_checked=checked,source_named_progress=str((ub/'progress.json').relative_to(ROOT)),
            root_s=ring.enc(s),FULLU_value=ring.enc(Uvalue),FULLU_inverse=ring.enc(Ui),FULLU_product_all6_factors=[ring.enc(v)for v in factors],
            leading_unit=ring.enc(lead),leading_unit_inverse=ring.enc(leadinv),three_root_unit_recipes=rootcerts,
            G_value=ring.enc(Gvalue),G_inverse=ring.enc(Gi),G_integer_Bezout_reused=True,no_new_inverse_search=True)
        payload=json.dumps(unitdata,indent=2)+'\n';guard(len(payload.encode()));(out/'source-open-units.json').write_text(payload)
        status='all22_zero_FULLU_unit_actual_Bbad_EXACT_Ka13';checkpoint(status)
    except Stop as e:status=e.status;checkpoint(status,e.data)
    complete=status=='all22_zero_FULLU_unit_actual_Bbad_EXACT_Ka13'
    proof=dict(status=status,R='QQ[a,e,d,1/(adFBCJ)]',I='actual I_3(B7x3)=(ALL22 normalized source minor generators)',
        upper_map='Accepted ALLQ/L exclusion: G unitinR/I. Norm memberships give H,h1,k1,m1 inI. Accepted Ka13 bridge gives g13=0,e=t; accepted d-Euclid givesd=s. Thusalpha:Ka13->R/I is a SURJECTION.',
        lower_map='All22 original source polynomials evaluateEXACTzero at(a,t,s), and FULLU has explicitinverse onwholeKa13. Thusbeta:R/I->Ka13 is a well-defined unitalmap.',
        maps_fix_a='beta(alpha(a))=a, andKa13 generatedbya overQQ. Thusbeta*alpha=id; alpha isbothsurjectiveandinjective. No lengthguess or field assumption.',
        actual_Fitting_bad_algebra_EXACT_Ka13_given_accepted_gates=complete,exact_QQ_length=13 if complete else None,
        actual_Bbad_nonempty_as_scheme_given_accepted_gates=complete,
        no_reducedness_etaleness_irreducibility_or_geometric_point_count_claim=True,
        no_actual_Poisson_full_cubic_rank_loss_quartic239_or_generic_exhaustion_claim=True,
        not_coker_B_module_length=True,independent_review_of_new_source_gate_pending=True)
    (out/'actual-scheme-isomorphism.json').write_text(json.dumps(proof,indent=2)+'\n')
    extra_stems=['S61_NONSQUARE_RAW_H_B_G0_REMOVED_SCHEME_BOUND_20261104_v63_REVIEW','S61_NONSQUARE_RAW_H_B_K13_SOURCE_BRIDGE_20261104_v63_REVIEW','S61_NONSQUARE_RAW_H_B_K13_D_EUCLID_20261104_v63_REVIEW','S61_NONSQUARE_RAW_H_B_K13_B21_ROOT_TEST_20261104_v63_REVIEW','S61_NONSQUARE_RAW_H_B_MIXED_NORM16_20261104_v63_REVIEW','S61_NONSQUARE_RAW_H_B_MINORS_GCD_20261104_v63_REVIEW']
    paths=[Path(__file__),ROOT/'code/s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63.py',ROOT/f'research/notes/{ID}.md',mb/'report.json',mb/'certificate.json',hb/'report.json',hb/'H4-monic-field-coefficients.json',db/'report.json',db/'result.json',ub/'report.json',ub/'progress.json',gb/'report.json',gb/'certificate.json',gb/'G-unit-certificate.json',ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_E_RECONSTRUCTION_20261104_v63/projection-isomorphism.json']+[ROOT/'reports'/stem/'report.json'for stem in extra_stems]
    report=dict(id=ID,all_passed=complete,status=status,all22_original_source_polynomials_zero=allzero,
        FULLU_explicit_unit_on_whole_Ka13=Uunit,actual_Fitting_bad_algebra_EXACT_Ka13_given_accepted_gates=complete,
        exact_actual_Fitting_bad_QQ_length=13 if complete else None,actual_Bbad_nonempty_given_accepted_gates=complete,
        no_field_irreducibility_or_etale_assumption=True,no_Poisson_fullrankloss_or_exhaustion_inference=True,
        monomial_nodes=len(vals),explicit_saved_terms=termcount,max_monomial_coefficient_bits=maxbits,field_integer_operations=ring.operations,
        no_new_inverse_search_norm_point_prime_GB_or_automatic_next_job=True,independent_review_pending=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
    assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:main(a.output)
    except Exception as e:
        if a.output.exists():(a.output/'failure.json').write_text(json.dumps(dict(all_passed=False,error=str(e),nodes_preserved=(a.output/'monomial-nodes.jsonl').exists(),no_auto_retry=True),indent=2)+'\n')
        raise
