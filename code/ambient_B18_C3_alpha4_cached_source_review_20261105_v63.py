#!/usr/bin/env python3
"""Independent alpha4 source/custody audit, never composing requested values."""
import ast,hashlib,json,math,pathlib,re,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
SOURCE='AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63'
ID=SOURCE+'_REVIEW'
PIN='3c4909ff3908664e2986679aabeb99c592ecb7debc107b8faed4dbbdb9e86082'
MON=[(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)]
MI={m:i for i,m in enumerate(MON)}
PAIRS=[[(i,j) for i,a in enumerate(MON) for j,b in enumerate(MON)
        if tuple(x+y for x,y in zip(a,b))==m] for m in MON]
def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def field(p,k):
    text=pathlib.Path(p).read_text();depth=0;i=0;decoder=json.JSONDecoder()
    while i<len(text):
        ch=text[i]
        if ch=='"':
            start=i;i+=1
            while True:
                if text[i]=='\\':i+=2
                elif text[i]=='"':i+=1;break
                else:i+=1
            if depth==1:
                j=i
                while text[j].isspace():j+=1
                if text[j]==':' and decoder.raw_decode(text,start)[0]==k:
                    j+=1
                    while text[j].isspace():j+=1
                    return decoder.raw_decode(text,j)[0]
        else:
            if ch in '{[':depth+=1
            elif ch in '}]':depth-=1
            i+=1
    raise KeyError(k)
def nc(x):
    if isinstance(x,bool) or x is None:return 0
    if isinstance(x,(int,float,Q)):return 1
    if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
    if isinstance(x,(list,tuple)):return sum(map(nc,x))
    if isinstance(x,dict):return sum(map(nc,x.values()))
    return 0
def ceilabs(q):return (abs(q.numerator)+q.denominator-1)//q.denominator
def encode(x):return json.dumps(x,separators=(',',':'),sort_keys=True)

class Majorants:
    """Exact origin; otherwise absolute integer bound and denominator divisor.

    Also counts BEFORE-operation numerator/denominator bit bounds. It never
    represents an actual positive-order requested scalar coefficient.
    """
    def __init__(self):self.pre=0;self.node_coarse=0;self.point_bits=0
    def limbs(self,x):
        if isinstance(x,Q):return abs(x.numerator).bit_length(),x.denominator.bit_length(),not x
        U,D=x;return (U*D).bit_length(),D.bit_length(),not U
    def check_add(self,x,y):
        nx,dx,zx=self.limbs(x);ny,dy,zy=self.limbs(y)
        if zx or zy:return
        self.pre=max(self.pre,max(nx+dy,ny+dx)+1,dx+dy)
    def check_mul(self,x,y):
        nx,dx,zx=self.limbs(x);ny,dy,zy=self.limbs(y)
        if zx or zy:return
        self.pre=max(self.pre,nx+ny,dx+dy)
    def add(self,x,y):
        self.check_add(x,y)
        if isinstance(x,Q):assert isinstance(y,Q);return x+y
        return (x[0]+y[0],math.lcm(x[1],y[1]))
    def mul(self,x,y):
        self.check_mul(x,y)
        if isinstance(x,Q) and isinstance(y,Q):return x*y
        if isinstance(y,Q):x,y=y,x
        if isinstance(x,Q):
            if not x:return (0,1)
            # Preserve the primary's harmless denominator divisor even
            # when U=0, so its complete metadata hash is comparable.
            return (ceilabs(x)*y[0],x.denominator*y[1])
        if not x[0] or not y[0]:return (0,1)
        return (x[0]*y[0],x[1]*y[1])
    def vector_add(self,a,b):return [self.add(x,y) for x,y in zip(a,b)]
    def conv(self,a,b):
        out=[Q(0)]+[(0,1)]*5
        for m,terms in enumerate(PAIRS):
            for i,j in terms:out[m]=self.add(out[m],self.mul(a[i],b[j]))
        return out
    def inverse(self,a):
        assert a[0]==1;v=[Q(1)]+[(0,1)]*5
        for m,terms in enumerate(PAIRS):
            if m:
                for i,j in terms:
                    if i:v[m]=self.add(v[m],self.mul(a[i],v[j]))
        return v
    def stored(self,v):
        self.point_bits=max(self.point_bits,abs(v[0].numerator).bit_length(),v[0].denominator.bit_length())
        for U,D in v[1:]:self.node_coarse=max(self.node_coarse,U.bit_length()+D.bit_length())
    def serialize(self,v):return [[str(v[0]),'1']]+[[str(U),str(D)] for U,D in v[1:]]

def input_envelopes(jets,mons):
    provided={}
    for i,values in jets:
        a=list(map(Q,values));v=[Q(a[0])]+[(0,1)]*5
        for (nr,ns),q in zip(mons,a):
            for ld in range(ns+1):
                mon=(nr+ns-ld,ld)
                if mon in MI and mon!=(0,0):
                    k=MI[mon];t=q*math.comb(ns,ld)/5**ld
                    U,D=v[k];v[k]=(U+ceilabs(t),math.lcm(D,t.denominator))
        provided[i]=v
    return provided

def main():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw_read=0
    def pin(p,h=None):
        nonlocal raw_read
        p=pathlib.Path(p);actual=sha(p)
        if h is not None:assert actual==h,(str(p),actual,h)
        k=str(p.relative_to(ROOT))
        if k not in pins:raw_read+=p.stat().st_size
        pins[k]=actual
    sourcep=ROOT/'reports'/SOURCE/'report.json';pin(sourcep,PIN);pr=read(sourcep);assert pr['all_passed']
    for p,h in pr['input_sha256'].items():pin(ROOT/p,h)
    for p,h in pr['evidence_sha256'].items():pin(ROOT/'reports'/SOURCE/p,h)
    pin(ROOT/'code/ambient_B18_C3_alpha4_cached_source_20261105_v63.py',pr['code_sha256'])
    pin(ROOT/'research/notes'/f'{SOURCE}.md',pr['note_sha256'])
    c=read(ROOT/'reports'/SOURCE/'compiled-eight-even-roots-and-provided-joint105-plane3-cache.json')
    profile=read(ROOT/'reports'/SOURCE/'alpha4-source-callback-binding-and-operational-height-profile.json')
    records=c['operator_records_originalID_op_args'];assert len(records)==7945
    originalpath=ROOT/'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63/compiled-original164-DAG-and-provided-pure-ell-inputs.json'
    original=field(originalpath,'operator_records_originalID_op_args')
    # Merge the sorted lists; no second18584-key table is allocated.
    k=0
    for r in original:
        if k<len(records) and r[0]==records[k][0]:assert r==records[k];k+=1
    assert k==len(records);original_numeric=nc(original);del original
    meta=field(originalpath,'physical12_coefficient_roots')
    assert c['root_pairs_physical4']==[r for r in meta if r['physical_row'] in [19,28,29,72]]
    roots=[r[name] for r in c['root_pairs_physical4'] for name in ['aA_node','aB_node']]
    assert roots==[20117,40895,20132,40910,20135,40913,20159,40937]
    needed=set(roots)
    for r in reversed(records):
        i,op=r[:2]
        if i in needed and op not in ['input','c']:needed.update(r[2:])
    assert needed=={r[0] for r in records}
    seen=set();fixtures=[]
    for r in records:
        i,op=r[:2]
        assert i not in seen
        if op not in ['input','c']:assert set(r[2:])<=seen
        elif op=='c':assert len(r)==4 and r[3]>0;fixtures.append([i,r[2],r[3],str(Q(r[2],r[3]))])
        seen.add(i)
    del seen,needed
    jointpath=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    jets=field(jointpath,'all105_unknown_jets');mons=field(jointpath,'monomials')
    unknown=field(jointpath,'all105_unknown_original_node_order')
    equations=field(jointpath,'all105_equation_original_node_order')
    assert field(jointpath,'all105_equations_zero_through_degree3')
    assert jets==c['provided105_originalID_r_s_coefficients10'] and mons==c['provided_monomials10']
    assert [i for i,v in jets]==unknown
    assert equations==field(originalpath,'equation105_order')
    assert set(unknown)=={i for i,v in field(originalpath,'provided105_originalID_pure_ell_coefficients')}
    q20path=ROOT/'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase1-common105-and-114-source-equations.json'
    qm=field(q20path,'base27_monomials');qvals=dict(field(q20path,'all105_cached_composite_values'))
    qi=[qm.index(m) for m in [[],[7],[8],[7,7],[7,8],[8,8]]]
    assert set(qvals)==set(unknown)
    conversion_pre=0;input_max=0
    def exact_add(a,b):
        nonlocal conversion_pre
        if not a:return b
        if not b:return a
        conversion_pre=max(conversion_pre,max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1,a.denominator.bit_length()+b.denominator.bit_length())
        return a+b
    def exact_mul(a,b):
        nonlocal conversion_pre
        if not a or not b:return Q(0)
        conversion_pre=max(conversion_pre,abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length())
        return a*b
    for i,vals in jets:
        a=list(map(Q,vals));input_max=max(input_max,max(max(abs(v.numerator).bit_length(),v.denominator.bit_length()) for v in a))
        lower=[a[0],exact_add(a[1],a[2]),exact_mul(a[2],Q(1,5)),
          exact_add(exact_add(a[3],a[4]),a[5]),exact_mul(exact_add(a[4],exact_mul(Q(2),a[5])),Q(1,5)),exact_mul(a[5],Q(1,25))]
        assert lower==[Q(qvals[i][j]) for j in qi]
        # Input transformation, not any requested scalar-DAG readback.
        exact_mul(exact_add(a[8],exact_mul(Q(3),a[9])),Q(1,25))
    del qvals,qm
    # Check the algebra quotient/substitution independently on all6x6 basis
    # products and all10 input monomials. No requested root is composed.
    assert MON==[tuple(m) for m in c['B6_monomials']]
    for m in mons:
        assert len(m)==2 and sum(m)<=3
    for a in MON:
        for b in MON:
            x=tuple(i+j for i,j in zip(a,b));assert (x in MI)==(x[0]<2 and x[1]<3)
    assert all(sum(m)<=3 for m in MON) and all(sum(m)<4 for m in MON)
    # Replay metadata bounds only. Positive-order requested values are never
    # stored or calculated. Origins are fixed prior point data.
    provided=input_envelopes(jets,mons);env=Majorants();cache={};uses=Counter()
    for r in records:
        if r[1] not in ['input','c']:uses.update(r[2:])
    uses.update(roots);peak=0;hasher=hashlib.sha256();inverse_pre=0;inverses=[]
    max_u=max_d=0
    for r in records:
        i,op=r[:2]
        if op=='input':v=provided[i]
        elif op=='c':v=[Q(r[2],r[3])]+[(0,1)]*5
        elif op=='a':v=env.vector_add(cache[r[2]],cache[r[3]])
        elif op=='n':v=[-cache[r[2]][0]]+cache[r[2]][1:]
        elif op=='m':v=env.conv(cache[r[2]],cache[r[3]])
        elif op=='i':
            v=env.inverse(cache[r[2]]);inverses.append([i,r[2],'1'])
            check=Majorants();check.conv(v,cache[r[2]]);check.conv(cache[r[2]],v)
            inverse_pre=max(inverse_pre,check.pre)
        else:raise AssertionError(op)
        env.stored(v)
        for U,D in v[1:]:max_u=max(max_u,U.bit_length());max_d=max(max_d,D.bit_length())
        hasher.update((encode([i,env.serialize(v)])+'\n').encode())
        cache[i]=v;peak=max(peak,len(cache))
        if op not in ['input','c']:
            for j in r[2:]:
                uses[j]-=1
                if not uses[j]:del cache[j]
        if not uses[i]:del cache[i]
    assert peak==184 and inverses==c['both_original_inverse_centers']
    assert hasher.hexdigest()==profile['every_node_majorant_sha256']
    assert env.node_coarse==profile['max_coefficient_height_bound']==39
    assert max_u==profile['max_positive_abs_ceil_bits'] and max_d==profile['max_denominator_divisor_bits']
    expected_bounds=[[i,env.serialize(cache[i])] for i in roots]
    assert expected_bounds==profile['root10_final_majorants']
    diffenv=Majorants()
    for m in c['root_pairs_physical4']:diffenv.vector_add(cache[m['aA_node']],cache[m['aB_node']])
    whole_pre=max(env.pre,conversion_pre,inverse_pre,diffenv.pre,input_max)
    assert whole_pre<=79<2048
    prelim={'all_original7945_records_and_minimal_closure_match':True,
      'all630_input_overlaps_and_same105_unknown_equations_match':True,
      'all45_constants_exact_numerator_denominator':fixtures,
      'all_node_absolute_denominator_majorants_hash_match':True,'lastuse_peak':peak,
      'node_coefficient_bound':env.node_coarse,'original_DAG_operation_prebound':env.pre,
      'input_conversion_operation_prebound':conversion_pre,'inverse_diagnostic_prebound':inverse_pre,
      'final_delta_subtraction_prebound':diffenv.pre,'whole_operational_prebound':whole_pre,
      'no_positive_order_requested_root_or_alpha_value_computed':True,
      'source_static_numeric_upper':125000,'global_hard_numeric':150000,
      'old_original_record_numeric_count':original_numeric,'RAW_hash_read':raw_read}
    save(out/'source-majorant-and-frame-checkpoint.json',prelim)
    save(out/'phase.json',{'status':'SOURCE_CORE_PASS_pending_prepared_custody_and_lower_target_bindings','no_value_GO_yet':True})
    print(json.dumps({'core_passed':True,'whole_prebound':whole_pre,'inverse_pre':inverse_pre,'DAG_pre':env.pre,'conversion_pre':conversion_pre,'delta_pre':diffenv.pre,'elapsed':time.monotonic()-started}))
    signal.alarm(0)
if __name__=='__main__':main()
