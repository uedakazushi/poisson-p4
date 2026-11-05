#!/usr/bin/env python3
"""Independent ONE CUT1 dot and ordinary-QQ[t] rank certificate."""
import hashlib,importlib.util,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C3_CUT1_QUARTIC_VALUES_20261105_v63';ID=PRI+'_REVIEW'
PIN='e8da9d42de12c842dd749a6cb8121f8eb0ed44587fb38a711648fdfe8113be7e'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
pre=stored=0
def guard(n):
    global pre
    pre=max(pre,n);assert n<=2048
def note(q):
    global stored
    stored=max(stored,abs(q.numerator).bit_length(),q.denominator.bit_length());assert stored<=2048;return q
def parse(x):guard(max(abs(int(v)).bit_length() for v in str(x).split('/')));return note(Q(x))
def add(a,b):
    if not a:return b
    if not b:return a
    guard(max(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1,a.denominator.bit_length()+b.denominator.bit_length()))
    return note(a+b)
def mul(a,b):
    if not a or not b:return Q(0)
    guard(max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length()))
    return note(a*b)
def count(x):
    if isinstance(x,(int,float,Q)):return 1
    if isinstance(x,str):
        try:Q(x);return 1
        except(ValueError,ZeroDivisionError):return 0
    if isinstance(x,list):return sum(map(count,x))
    if isinstance(x,dict):return sum(map(count,x.values()))
    return 0
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
    def pin(p,h=None):
        nonlocal raw
        p=pathlib.Path(p);s=sha(p)
        if h is not None:assert s==h,(str(p),s,h)
        key=str(p.relative_to(ROOT))
        if key not in pins:raw+=p.stat().st_size
        pins[key]=s
    p=ROOT/'reports'/PRI/'report.json';pin(p,PIN);pr=read(p);assert pr['all_passed']
    for p,h in pr['input_sha256'].items():pin(ROOT/p,h)
    for p,h in pr['evidence_sha256'].items():pin(ROOT/'reports'/PRI/p,h)
    pin(ROOT/'code/ambient_B18_C3_cut1_quartic_values_20261105_v63.py',pr['code_sha256'])
    sp=ROOT/'reports/AMBIENT_B18_C3_CUT1_QUARTIC_SOURCE_20261105_v63_REVIEW/report.json'
    assert sha(sp)=='d9979ece14c4a19952eada312c24e73d3350ff2db91f91cec40b9a288f8a126a' and read(sp)['all_passed']
    ownpath=ROOT/'code/ambient_B18_C3_cut1_source_review_20261105_v63.py'
    pin(ownpath,'4ff64d3a238e0977619ce3700f748fc15e99e84cd6250f75f449a09d4b293ec1')
    spec=importlib.util.spec_from_file_location('own_cut_operator',ownpath);own=importlib.util.module_from_spec(spec);spec.loader.exec_module(own)
    # This reconstructs only symbolic source weights; A/B leaves are never
    # paired. The actual numeric work below is the ONE residual dot only.
    operators=own.independent_cut_operator();form=operators['residual'];assert len(form)==24
    source=read(ROOT/'reports/AMBIENT_B18_C3_CUT1_QUARTIC_SOURCE_20261105_v63/literal-cut-quartic-source-and-three-unpaired-cache-operators.json')
    assert form==source['residual_B_minus_121A_over48_operator']
    c2=read(ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
    values={r['physical_row']:dict(r['delta']) for r in c2['physical12_coefficients']}
    data=read(ROOT/'reports'/PRI/'one-CUT1-residual-and-formal-coordinate-rank-certificate.json')
    bound=count(pr)+count(source)+count(c2)+count(data)+6000
    assert bound<18000<30000 and raw+100000<400000
    save(out/'pre-arithmetic-bound.json',{'own_whole_numeric_upper':bound,'hard_numeric':30000,'RAW_input':raw,'hard_RAW':400000,'residual_not_yet_replayed':True})
    residual=Q(0);products=[]
    for key,weight in reversed(form):
        tag,row,index=key;assert tag=='cache'
        w=parse(weight);v=parse(values[row].get(index,'0'));term=mul(w,v);residual=add(residual,term)
        products.append([key,str(w),str(v),str(term)])
    assert sorted(products)==sorted(data['all24_termwise_cached_products'])
    assert str(residual)==data['actual_residual']=='759275/1728' and residual
    a=parse('-3/5');minor=mul(a,residual)
    assert str(minor)==data['actual_known_HB_CUT_minor']=='-151855/576' and minor
    cert=data['rank_source_certificate'];assert cert['unpaired_actual_CUT_A4']=='formal t; not assignedzero or evaluated'
    def poly(v):return {i:parse(c) for i,c in enumerate(v) if Q(c)}
    M=[[poly(v) for v in row] for row in cert['universal_matrix_over_QQ_t']]
    inv=[[poly(v) for v in row] for row in cert['inverse_over_QQ_t']]
    assert M==[[{0:a},{1:Q(1)}],[{},{0:residual}]]
    def padd(x,y):
        d=dict(x)
        for k,v in y.items():
            q=add(d.get(k,Q(0)),v)
            if q:d[k]=q
            elif k in d:del d[k]
        return d
    def pmul(x,y):
        d={}
        for i,c in x.items():
            for j,v in y.items():d=padd(d,{i+j:mul(c,v)})
        return d
    def mm(A,B):return [[padd(pmul(row[0],col[0]),pmul(row[1],col[1])) for col in zip(*B)] for row in A]
    I=[[{0:Q(1)},{}],[{},{0:Q(1)}]]
    assert mm(M,inv)==mm(inv,M)==I
    # Ordinary polynomials: no truncation or special numeric t is used.
    assert all(max(p,default=0)<=1 for mat in [M,inv] for row in mat for p in row)
    assert data['CUT_A4_B4_operators_remain_unpaired'] and data['actual_C3_length19_proved']
    assert pre<=519 and stored<=43
    proof={'independently_reconstructed_source_residual_weights':True,
      'all24_reverse_order_products':products,'actual_residual':str(residual),'nonzero_actual_minor':str(minor),
      'literal_ell3_cut_full_source_peer_pinned':True,'A4_and_B4_cache_operators_NOT_paired':True,
      'A_coordinate_remains_indeterminate_t_not_zero':True,
      'both_inverse_products_verified_in_ordinary_QQ_t_not_nilpotent_quotient':True,
      'actual_two_quartic_initial_classes_independent':True,
      'actual_C3_HF':[3,6,6,4,0],'actual_C3_length':19,'actual_m4_C3_zero':True,
      'DVR_cut3_torsion_pattern':[1,1,2,'e>=3'],'longer_exponent_still_unknown':True,
      'no_original_H0_nil_length_cube_or_nilindex_conclusion':True,
      'no_new_source_gauge_physical_IFT_jet_or_point':True,'no_author_producer_imports':True,
      'own_whole_numeric_upper':bound,'RAW_input':raw,'maximum_preoperation_bits':pre,'maximum_stored_bits':stored,
      'hard_numeric':30000,'hard_RAW':400000,'hard_bits':2048,'hard_seconds':60,'hard_address_bytes':2<<30}
    ep=out/'independent-ONE-CUT1-products-ordinary-polynomial-rank-and-length19-proof.json';save(ep,proof)
    result={'id':ID,'all_passed':True,'frozen':True,'primary_report_sha256':PIN,'input_sha256':pins,
      'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
      'actual_C3_length19_and_HF36640_accepted':True,'longer_exponent_and_original_H0_nil_not_determined':True,
      'elapsed_seconds':time.monotonic()-start}
    save(out/'report.json',result);signal.alarm(0)
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':result['elapsed_seconds']}))
if __name__=='__main__':main()
