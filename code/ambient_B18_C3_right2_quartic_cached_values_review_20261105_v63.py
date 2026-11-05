#!/usr/bin/env python3
"""Independent dual-number six-root replay and coherent RIGHT2 fixed rank."""
import hashlib,importlib.util,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63'
ID=PRI+'_REVIEW'
PIN='13b27dd23843107e9fe3ad97658fe58ef61d582bedbf6a61a1ecf6bd03d6e7e6'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def module(name,p,h):
    assert sha(p)==h;s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
    def pin(p,h=None):
        nonlocal raw
        p=pathlib.Path(p);s=sha(p)
        if h is not None:assert s==h,(str(p),s,h)
        k=str(p.relative_to(ROOT))
        if k not in pins:raw+=p.stat().st_size
        pins[k]=s
    p=ROOT/'reports'/PRI/'report.json';pin(p,PIN);pr=read(p);assert pr['all_passed']
    for p,h in pr['input_sha256'].items():pin(ROOT/p,h)
    for p,h in pr['evidence_sha256'].items():pin(ROOT/'reports'/PRI/p,h)
    pin(ROOT/'code/ambient_B18_C3_right2_quartic_cached_values_20261105_v63.py',pr['code_sha256'])
    arithmetic=ROOT/'code/ambient_B18_C3_alpha4_cached_values_review_20261105_v63.py'
    E=module('own_dualnumber',arithmetic,'7187f5c0c19f9d91d961563071e863291114a38f8ed605808036c70df16d4788');pin(arithmetic)
    h=ROOT/'code/ambient_B18_C3_alpha4_cached_source_review_20261105_v63.py'
    own=module('own_reader',h,'430cc27fc33b263cdcaec0fd616b87e88824bef702a3ff49fbc7136eaf428313');pin(h)
    c=read(ROOT/'reports/AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63/compiled-six-even-roots-and-provided-joint105-plane3-cache.json')
    phase1=read(ROOT/'reports'/PRI/'phase1-six-root-existing-cache-readbacks.json')
    phase2=read(ROOT/'reports'/PRI/'phase2-actual-RIGHT2-quartics-two-factors-and-rank.json')
    sp=ROOT/'reports/AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63_REVIEW/report.json'
    assert sha(sp)=='431c456552d3bea3845bfca334e2528aa4608fae749caa2f17bcc692509279a2' and read(sp)['all_passed']
    mathp=ROOT/'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63_REVIEW/report.json';mathr=read(mathp)
    assert sha(mathp)=='05d8c4eb681961627e8d3d6210856052b52a6c6945efebfc055e2621b0b975af'
    for p,hsh in mathr['evidence_sha256'].items():
        ep=mathp.parent/p;pin(ep,hsh);proof=read(ep)
    forms=proof['independently_multiplied_HB_polynomials_and_direct_NF4']
    c2=read(ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
    c2={r['physical_row']:dict(r['delta']) for r in c2['physical12_coefficients']}
    records=c['operator_records_originalID_op_args'];assert len(records)==7095
    bound=sum(map(own.nc,[c,phase1,phase2,proof,pr]))+2*7095+2940+735+162*7+22000
    assert bound<75000<100000 and raw+200000<32<<20
    save(out/'pre-arithmetic-source-and-bound.json',{'input_hashes_verified':True,'owned_numeric_upper':bound,'hard_numeric':100000,'hard_bits':2048,'requested_values_not_yet_replayed':True})
    R=([Q(0)]*3,[Q(1),Q(0),Q(0)]);S=([Q(0),Q(1,5),Q(0)],[Q(1),Q(0),Q(0)])
    rp=[E.one(),R];ss=[E.one(),S]
    for _ in range(2):rp.append(E.multiply(rp[-1],R));ss.append(E.multiply(ss[-1],S))
    mons=c['provided_monomials10'];monvalues=[E.multiply(rp[a],ss[b]) for a,b in mons]
    qpath=ROOT/'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase1-common105-and-114-source-equations.json'
    qm=own.field(qpath,'base27_monomials');qv=dict(own.field(qpath,'all105_cached_composite_values'))
    qi=[qm.index(m) for m in [[],[7],[8],[7,7],[7,8],[8,8]]];lowmons=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
    phi={};overlaps=0
    for i,vals in c['provided105_originalID_r_s_coefficients10']:
        x=E.zero();lower=[Q(0)]*6
        for (nr,ns),q,mv in zip(mons,map(E.parse,vals),monvalues):
            if q:x=E.plus(x,E.scale(mv,q))
            if nr+ns<=2:
                for j in range(ns+1):
                    k=lowmons.index((nr+ns-j,j));lower[k]=E.add(lower[k],E.times(q,E.constant(math.comb(ns,j),5**j)))
        assert lower==[E.parse(qv[i][k]) for k in qi]
        assert E.flatten(x)[:5]==[lower[j] for j in [0,1,2,4,5]];overlaps+=6;phi[i]=x
    assert overlaps==630 and set(phi)==set(qv);del qv,qm,rp,ss,monvalues
    roots=[r[k] for r in c['root_pairs_physical3'] for k in ['aA_node','aB_node']]
    uses=[None]*(max(r[0] for r in records)+1);cache=[None]*len(uses)
    for r in records:uses[r[0]]=0
    for r in records:
        if r[1] not in ['c','input']:
            for i in r[2:]:uses[i]+=1
    for i in roots:uses[i]+=1
    live=peak=0;units=[]
    try:
        for r in records:
            i,op=r[:2]
            if op=='input':x=phi[i]
            elif op=='c':x=([E.constant(r[2],r[3]),Q(0),Q(0)],[Q(0)]*3)
            elif op=='a':x=E.plus(cache[r[2]],cache[r[3]])
            elif op=='n':x=E.negative(cache[r[2]])
            elif op=='m':x=E.multiply(cache[r[2]],cache[r[3]])
            elif op=='i':
                x=E.invert(cache[r[2]]);assert E.multiply(x,cache[r[2]])==E.multiply(cache[r[2]],x)==E.one();units.append([i,r[2],True])
            else:raise AssertionError(op)
            cache[i]=x;live+=1;peak=max(peak,live)
            if op not in ['input','c']:
                for j in r[2:]:
                    uses[j]-=1
                    if not uses[j]:cache[j]=None;live-=1
            if not uses[i]:cache[i]=None;live-=1
    except BaseException as ex:
        save(out/'ROOT_PREFIX_STOP.json',{'error':repr(ex),'completed_node':i,'roots':[[n,list(map(str,E.flatten(cache[n])))] for n in roots if cache[n] is not None],'no_cap_raise':True});raise
    assert peak==162 and units==[[a,b,True] for a,b,z in c['both_original_inverse_centers']]
    rootdata=[[i,list(map(str,E.flatten(cache[i])))] for i in roots];assert rootdata==phase1['all6_even_roots_B6']
    lower=dict(c['expected6_lower_B6_coefficients_order_1_r_ell_rell_ell2']);lam={};deltas=[]
    for meta in c['root_pairs_physical3']:
        a,b=cache[meta['aA_node']],cache[meta['aB_node']]
        assert E.flatten(a)[:5]==list(map(E.parse,lower[meta['aA_node']]))
        assert E.flatten(b)[:5]==list(map(E.parse,lower[meta['aB_node']]))
        d=E.flatten(E.plus(a,E.negative(b)));assert d[:5]==[Q(0)]*5
        lam[meta['physical_row']]=d[5];deltas.append([meta['physical_row'],list(map(str,d))])
    assert E.premax<=79
    assert lam=={18:Q(29,30),60:Q(53,15),67:Q(17,180)}
    # The source weights are our independently constructed source/NF4
    # evidence, then evaluated in the reverse term order from the producer.
    matrix=[];termproof=[]
    for row in forms:
        values=[];logs=[]
        for form in row:
            q=Q(0);log=[]
            for key,w in reversed(form):
                kind,p,j=key;v=lam[p] if kind=='lambda' else E.parse(c2[p].get(j,'0'))
                assert kind in ['lambda','cache'];weight=E.parse(w);prod=E.times(weight,v);q=E.add(q,prod)
                log.append([key,str(weight),str(v),str(prod)])
            values.append(q);logs.append(log)
        matrix.append(values);termproof.append(logs)
    assert [[str(v) for v in row] for row in matrix]==phase2['actual_quartic_matrix_rows_A4_B4_columns_HB0_HB1']
    for t in range(2):
        for j in range(2):
            assert sorted(termproof[t][j])==sorted(phase2['all_termwise_affine_products'][t][j])
    a0,a1=matrix[0];b0,b1=matrix[1]
    C0=E.add(b0,E.times(Q(1,4),a0));C1=E.add(b1,E.times(Q(1,4),a1));F=E.add(E.times(Q(10),a0),a1)
    det=E.add(E.times(a0,b1),-E.times(a1,b0))
    assert E.add(C1,E.times(Q(10),C0))==0 and det==-E.times(C0,F)
    assert [str(C0),str(F)]==phase2['two_affine_factors'] and str(det)==phase2['determinant']
    assert not det and any(v for row in matrix for v in row)
    cert=phase2['rank_certificate'];i,j=cert['pivot_row_column']
    column=list(map(E.parse,cert['pivot_column']));factors=list(map(E.parse,cert['all_two_source_column_factors']))
    assert column==[matrix[t][j] for t in range(2)] and column[i]
    assert all(E.times(column[t],factors[k])==matrix[t][k] for t in range(2) for k in range(2))
    assert phase2['actual_RIGHT2_rank']==1 and not phase2['actual_C3_length19_proved']
    remaining_dual=[-Q(121,48),Q(1)]
    assert all(sum(remaining_dual[t]*matrix[t][j] for t in range(2))==0 for j in range(2))
    assert E.premax<=675 and E.valmax<=56
    evidence={'own_original_scalar_composition':'QQ[ell]/ell3 dual numbers A+rB; closed inverse Ainv-rB*Ainv²',
      'all7095_source_ops_630_input_overlaps_30_prior_root_values':True,'both_inverse_products':units,
      'all6_even_roots_B6':rootdata,'three_actual_lambda_deltas':deltas,
      'all72_affine_term_products_from_independent_source_NF_weights_match':True,
      'actual_RIGHT2_matrix':[[str(v) for v in row] for row in matrix],
      'two_factors':[str(C0),str(F)],'determinant':str(det),'rank':1,
      'all4_fixed_column_factorizations_verified':True,'column_factors':list(map(str,factors)),
      'known_subset_remaining_quartic_dual':list(map(str,remaining_dual)),
      'remaining_dual_not_asserted_to_survive_all_actual_quartic_relations':True,
      'C3_length_still_19_or20':True,'no_fullquarticrank_or_length20_inference':True,
      'no_new_source_solve_physical_gauge_IFT_jet_or_parameter':True,'no_author_producer_imports':True,
      'different_inverse_and_reversed_summand_order':True,'peak_live_nodes':peak,
      'max_before_QQ_bits':E.premax,'max_stored_QQ_bits':E.valmax,'owned_numeric_upper':bound,
      'hard_numeric':100000,'RAW_input':raw,'hard_RAW':32<<20,'hard_bits':2048,'hard_seconds':60,'hard_address_bytes':2<<30}
    ep=out/'independent-sixroot-dualnumber-and-all-four-quartic-source-folds.json';save(ep,evidence)
    result={'id':ID,'all_passed':True,'frozen':True,'primary_report_sha256':PIN,'input_sha256':pins,
      'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
      'actual_RIGHT2_rank1_accepted':True,'C3_length_19_or20_unresolved':True,'elapsed_seconds':time.monotonic()-start}
    save(out/'report.json',result);signal.alarm(0)
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':result['elapsed_seconds']}))
if __name__=='__main__':main()
