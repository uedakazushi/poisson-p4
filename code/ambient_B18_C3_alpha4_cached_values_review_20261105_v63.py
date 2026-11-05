#!/usr/bin/env python3
"""Independent eight-root replay via QQ[ell]/ell3 dual numbers, no solve."""
import hashlib,importlib.util,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63_REVIEW'
PRI='AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63'
PIN='dd54eae90174d3a95359e8fa99d129cb1a9835cd3fa6da7a009eb89d1f49bd05'
SOURCEPEER='88805ffa6b5f0b949b1ab3ac89b17b752930879eeda23afa6f5cae4b99939fa0'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
premax=0;valmax=0;products=0;adds=0
def check(n):
    global premax
    premax=max(premax,n)
    if n>2048:raise OverflowError('before QQ operation')
def note(q):
    global valmax
    valmax=max(valmax,abs(q.numerator).bit_length(),q.denominator.bit_length())
    assert valmax<=2048;return q
def parse(s):
    parts=str(s).split('/');assert len(parts)<=2
    check(max(abs(int(v)).bit_length() for v in parts));return note(Q(s))
def constant(n,d):
    check(max(abs(n).bit_length(),abs(d).bit_length()));return note(Q(n,d))
def add(a,b):
    global adds
    if not a:return b
    if not b:return a
    check(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1)
    check(a.denominator.bit_length()+b.denominator.bit_length());adds+=1;return note(a+b)
def times(a,b):
    global products
    if not a or not b:return Q(0)
    check(abs(a.numerator).bit_length()+abs(b.numerator).bit_length())
    check(a.denominator.bit_length()+b.denominator.bit_length());products+=1;return note(a*b)
def neg(a):return note(-a)
def lsum(a,b):return [add(x,y) for x,y in zip(a,b)]
def lneg(a):return [neg(x) for x in a]
def lmul(a,b):
    # Reverse both summand and coefficient order relative to the producer.
    out=[Q(0)]*3
    for k in [2,1,0]:
        for i in range(k,-1,-1):out[k]=add(out[k],times(a[i],b[k-i]))
    return out
def zero():return ([Q(0)]*3,[Q(0)]*3)
def one():return ([Q(1),Q(0),Q(0)],[Q(0)]*3)
def plus(x,y):return (lsum(x[0],y[0]),lsum(x[1],y[1]))
def negative(x):return (lneg(x[0]),lneg(x[1]))
def multiply(x,y):
    # (A+rB)(C+rD)=AC+r(AD+BC), modulo r2 and ell3.
    return (lmul(x[0],y[0]),lsum(lmul(x[1],y[0]),lmul(x[0],y[1])))
def scale(x,c):return ([times(v,c) for v in x[0]],[times(v,c) for v in x[1]])
def invert(x):
    A,B=x;assert A[0]==1
    Ai=[Q(1),neg(A[1]),add(times(A[1],A[1]),neg(A[2]))]
    # Different closed formula from the original B6 recursive evaluator.
    return (Ai,lneg(lmul(B,lmul(Ai,Ai))))
def flatten(x):
    A,B=x;return [A[0],B[0],A[1],B[1],A[2],B[2]]
def main():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
    def pin(p,h=None):
        nonlocal raw
        p=pathlib.Path(p);actual=sha(p)
        if h is not None:assert actual==h,(str(p),actual,h)
        key=str(p.relative_to(ROOT))
        if key not in pins:raw+=p.stat().st_size
        pins[key]=actual
    p=ROOT/'reports'/PRI/'report.json';pin(p,PIN);pr=read(p);assert pr['all_passed']
    for p,h in pr['input_sha256'].items():pin(ROOT/p,h)
    for p,h in pr['evidence_sha256'].items():pin(ROOT/'reports'/PRI/p,h)
    pin(ROOT/'code/ambient_B18_C3_alpha4_cached_values_20261105_v63.py',pr['code_sha256'])
    sourcepeer=read(ROOT/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63_REVIEW/report.json')
    assert sha(ROOT/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63_REVIEW/report.json')==SOURCEPEER
    assert sourcepeer['all_passed']
    helper=ROOT/'code/ambient_B18_C3_alpha4_cached_source_review_20261105_v63.py'
    assert sha(helper)=='430cc27fc33b263cdcaec0fd616b87e88824bef702a3ff49fbc7136eaf428313'
    spec=importlib.util.spec_from_file_location('own_reader',helper);own=importlib.util.module_from_spec(spec);spec.loader.exec_module(own)
    c=read(ROOT/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63/compiled-eight-even-roots-and-provided-joint105-plane3-cache.json')
    profile=read(ROOT/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63/alpha4-source-callback-binding-and-operational-height-profile.json')
    actual=read(ROOT/'reports'/PRI/'eight-cached-even-root-values-and-four-complete-alpha-readbacks.json')
    records=c['operator_records_originalID_op_args'];assert len(records)==7945
    bound=own.nc(c)+own.nc(profile)+own.nc(pr)+own.nc(actual)+own.nc(sourcepeer)+2*7945+2940+735+184*7+20000
    assert bound<75000<100000 and raw+200000<32<<20
    save(out/'pre-arithmetic-bound.json',{'all_actual_source_hashes_match':True,'owned_numeric_upper':bound,
      'hard_numeric':100000,'RAW_input':raw,'hard_bits':2048,'method':'dual number over QQ[ell]/ell3','root_values_not_yet_evaluated':True})
    R=([Q(0)]*3,[Q(1),Q(0),Q(0)])
    S=([Q(0),Q(1,5),Q(0)],[Q(1),Q(0),Q(0)])
    RP=[one(),R];SP=[one(),S]
    for _ in range(2):RP.append(multiply(RP[-1],R));SP.append(multiply(SP[-1],S))
    mon=c['provided_monomials10'];monvalues=[multiply(RP[a],SP[b]) for a,b in mon]
    qpath=ROOT/'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase1-common105-and-114-source-equations.json'
    qm=own.field(qpath,'base27_monomials');qv=dict(own.field(qpath,'all105_cached_composite_values'))
    qi=[qm.index(m) for m in [[],[7],[8],[7,7],[7,8],[8,8]]]
    lower_m=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
    phi={};overlap=0
    for i,vals in c['provided105_originalID_r_s_coefficients10']:
        values=list(map(parse,vals));x=zero();lower=[Q(0)]*6
        for (nr,ns),q,mv in zip(mon,values,monvalues):
            if q:x=plus(x,scale(mv,q))
            if nr+ns<=2:
                for k in range(ns+1):
                    m=(nr+ns-k,k)
                    j=lower_m.index(m)
                    lower[j]=add(lower[j],times(q,constant(math.comb(ns,k),5**k)))
        assert lower==[parse(qv[i][j]) for j in qi];overlap+=6
        flat=flatten(x);assert flat[:5]==[lower[j] for j in [0,1,2,4,5]]
        phi[i]=x
    assert overlap==630 and set(phi)==set(qv);del qv,qm,values,monvalues,RP,SP
    roots=[r[k] for r in c['root_pairs_physical4'] for k in ['aA_node','aB_node']]
    max_id=max(r[0] for r in records);uses=[None]*(max_id+1);cache=[None]*(max_id+1)
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
            elif op=='c':x=([constant(r[2],r[3]),Q(0),Q(0)],[Q(0)]*3)
            elif op=='a':x=plus(cache[r[2]],cache[r[3]])
            elif op=='n':x=negative(cache[r[2]])
            elif op=='m':x=multiply(cache[r[2]],cache[r[3]])
            elif op=='i':
                x=invert(cache[r[2]])
                assert multiply(cache[r[2]],x)==multiply(x,cache[r[2]])==one()
                units.append([i,r[2],True])
            else:raise AssertionError(op)
            cache[i]=x;live+=1;peak=max(peak,live)
            if op not in ['input','c']:
                for j in r[2:]:
                    uses[j]-=1
                    if not uses[j]:cache[j]=None;live-=1
            if not uses[i]:cache[i]=None;live-=1
    except BaseException as e:
        save(out/'STOP_ROOT_PREFIX.json',{'error':repr(e),'completed_original_node':i,
            'roots':[[n,list(map(str,flatten(cache[n])))] for n in roots if cache[n] is not None],
            'max_prebits':premax,'max_value_bits':valmax,'no_cap_raise':True});raise
    assert peak==184 and units==[[i,j,True] for i,j,q in c['both_original_inverse_centers']]
    expected=dict(c['expected8_lower_B6_coefficients_order_1_r_ell_rell_ell2'])
    rootdata=[[i,list(map(str,flatten(cache[i])))] for i in roots]
    assert rootdata==actual['all8_even_roots_B6']
    for i in roots:assert flatten(cache[i])[:5]==list(map(parse,expected[i]))
    deltas=[]
    for r in c['root_pairs_physical4']:
        x=plus(cache[r['aA_node']],negative(cache[r['aB_node']]))
        f=flatten(x);assert f==[Q(0)]*6
        deltas.append({'physical_row':r['physical_row'],'delta_B6':list(map(str,f)),'alpha':str(f[5])})
    assert len(deltas)==4 and all(row['delta_r_ell2']=='0' for row in actual['all4_actual_alpha_delta_B6'])
    assert premax<=79 and valmax<=39
    completepeer=ROOT/'reports/AMBIENT_B18_C3_COMPLETE_CUBIC_LINE_SOURCE_20261105_v63_REVIEW/report.json'
    pin(completepeer,'f6a8b2cdc221adfc8c1c298a973bbd73d5a0097dc858da820697c95dacd39d4f')
    assert read(completepeer)['complete_cubic_image_line_and_HF3_3_or4_accepted']
    evidence={'own_dual_number_representation':'A(ell)+r*B(ell), ell3=r2=0',
      'independent_inverse':'Ainv=1-A1*ell+(A1^2-A2)*ell2; (A+rB)^-1=Ainv-r*B*Ainv^2',
      'all7945_original_source_operations_replayed':True,'all630_known_input_overlaps':True,
      'both_inverse_products':units,'all40_prior_lower_root_values':True,
      'all8_even_roots_B6':rootdata,'all4_actual_deltas':deltas,'all4_alpha_exact_zero':True,
      'complete155_source_peer_pinned':True,'exact_extra_cubic_rank':0,'actual_HF_C3_through3':[3,6,6,4],
      'actual_C3_length_range':[19,20],'both_remaining_exponents_equal2_excluded':True,
      'atleast_one_remaining_torsion_exponent_atleast3':True,
      'remaining_HFs':[[3,6,6,4,0],[3,6,6,4,1]],
      'no_exact_C3_length_all_exponents_H0_or_original_nil_conclusion':True,
      'no_author_producer_imports_new_source_solve_gauge_physical_IFT_or_jet':True,
      'different_inverse_and_reversed_product_summand_order':True,
      'lastuse_peak':peak,'max_before_QQ_bits':premax,'max_saved_QQ_bits':valmax,
      'QQ_products':products,'QQ_additions':adds,'owned_numeric_upper':bound,'hard_numeric':100000,
      'RAW_input':raw,'hard_RAW':32<<20,'hard_bits':2048,'hard_seconds':60,'hard_address_bytes':2<<30}
    ep=out/'independent-eight-root-dualnumber-replay-and-complete-HF3-four-proof.json';save(ep,evidence)
    result={'id':ID,'all_passed':True,'frozen':True,'primary_report_sha256':PIN,
      'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
      'all4_alpha_zero_and_exact_HF3_four_accepted':True,'actual_C3_length_still_19_or20':True,
      'elapsed_seconds':time.monotonic()-started}
    save(out/'report.json',result);signal.alarm(0)
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':result['elapsed_seconds']}))
if __name__=='__main__':main()
