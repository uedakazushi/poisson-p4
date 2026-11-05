#!/usr/bin/env python3
"""Prepared ONE existing-cache root38648 + coherent HEAD23 job. UNRUN."""
import argparse,hashlib,json,re,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63'
CACHE='AMBIENT_B18_C_CONNECTING_HEAD23_CACHE_SOURCE_20261105_v63'
HEAD='AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63'
CUSTODY='AMBIENT_B18_C_CONNECTING_HEAD23_VALUE_CUSTODY_20261105_v63'
MON=[(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)];MI={m:i for i,m in enumerate(MON)}
PRELIMIT=83;MAXPRE=0;MAXVALUE=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def size(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def before(h):
 global MAXPRE
 MAXPRE=max(MAXPRE,h)
 if h>PRELIMIT or h>2048:raise ArithmeticError(('QQ_PREOP_STOP',h,PRELIMIT))
def check(q):
 global MAXVALUE
 MAXVALUE=max(MAXVALUE,size(q));assert size(q)<=2048
 return q
def parse(s):
 if isinstance(s,int):before(abs(s).bit_length());return check(Q(s))
 assert isinstance(s,str)and re.fullmatch(r'-?\d+(?:/\d+)?',s)
 a,b=(s.split('/')+[1])[:2]if '/'not in s else s.split('/')
 a,b=int(a),int(b);assert b>0;before(max(abs(a).bit_length(),b.bit_length()));return check(Q(a,b))
def add(a,b):
 before(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1)
 before(a.denominator.bit_length()+b.denominator.bit_length());return check(a+b)
def neg(a):before(size(a));return check(-a)
def mul(a,b):
 before(max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length()));return check(a*b)
def zero():return[Q(0)]*6
def conv(a,b):
 v=zero()
 for i,x in enumerate(MON):
  for j,y in enumerate(MON):
   k=MI.get((x[0]+y[0],x[1]+y[1]))
   if k is not None:v[k]=add(v[k],mul(a[i],b[j]))
 return v
def inv(x):
 assert x[0]==1;v=zero();v[0]=Q(1)
 for j,mon in enumerate(MON[1:],1):
  for k,other in enumerate(MON[1:],1):
   t=MI.get((mon[0]-other[0],mon[1]-other[1]))
   if t is not None and t<j:v[j]=add(v[j],neg(mul(x[k],v[t])))
 assert conv(x,v)==[Q(1)]+[Q(0)]*5==conv(v,x)
 return v
def numeric(x):
 if isinstance(x,(bool,int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple,set)):return sum(numeric(v)for v in x)
 if isinstance(x,dict):return sum(numeric(k)+numeric(v)for k,v in x.items())
 return 0
def encode(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return{str(k):encode(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return[encode(v)for v in x]
 return x
def self_test():
 global PRELIMIT
 assert parse('0')==0 and parse('-7/12')==Q(-7,12)and parse(5)==5
 x=list(map(Q,[1,2,-1,3,4,-2]));assert inv(x)
 old=PRELIMIT;PRELIMIT=3
 try:mul(Q(15),Q(17))
 except ArithmeticError:pass
 else:raise AssertionError('before-operation guard did not stop')
 PRELIMIT=old;print('PURE_GUARD_AND_B6_INVERSE_FIXTURES_PASS_NO_ACTUAL_SOURCE_VALUES')
def main(source_peer_path,head_peer_path):
 global PRELIMIT,MAXPRE,MAXVALUE
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20))
 pins={};RAW=0
 def read(path):
  nonlocal RAW
  p=ROOT/path if not Path(path).is_absolute()else Path(path);raw=p.read_bytes();RAW+=len(raw);assert RAW<=2000000
  pins[str(p.relative_to(ROOT))]=hashlib.sha256(raw).hexdigest();return json.loads(raw)
 custody=read(f'reports/{CUSTODY}/prepared-worker-code-source-custody.json');assert custody['prepared_worker_sha256']==sha(Path(__file__))
 report=read(f'reports/{CACHE}/report.json');headreport=read(f'reports/{HEAD}/report.json')
 assert report['all_passed']and headreport['all_passed']and not report['actual_root_or_edge_values_computed']
 peer=read(source_peer_path);hp=read(head_peer_path);assert peer['all_passed']and hp['all_passed']
 assert sha(ROOT/f'reports/{CACHE}/report.json')in peer['input_sha256'].values()
 assert sha(ROOT/f'reports/{HEAD}/report.json')in hp['input_sha256'].values()
 assert sha(Path(__file__))in peer['input_sha256'].values()
 compiled_path=f'reports/{CACHE}/compiled-single-bB63-root-and-provided-joint105-cache.json'
 profile_path=f'reports/{CACHE}/single-root-head23-coherent-positive-height-and-whole-budget.json'
 head_path=f'reports/{HEAD}/full-coherent-head23-two-socle-cache-polynomials-and-low-two-row-source.json'
 C=read(compiled_path);P=read(profile_path);H=read(head_path)
 assert pins[compiled_path]==report['evidence_sha256'][Path(compiled_path).name]
 assert pins[profile_path]==report['evidence_sha256'][Path(profile_path).name]
 assert pins[head_path]==headreport['evidence_sha256'][Path(head_path).name]
 assert C['root_original_ID']==38648 and C['B6_monomials']==[list(m)for m in MON]
 assert P['original_records']==2194 and P['lastuse_peak']==120 and P['whole_before_operation_bound']==207
 assert P['root_before_operation_bound']==83 and P['post_before_operation_bound']==207
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);saved_numeric=0;owned_peak=0
 base=numeric([C,P,H,report,headreport,custody,peer,hp])+2000
 assert base<30000
 def own(*live):
  nonlocal owned_peak
  owned_peak=max(owned_peak,base+saved_numeric+sum(numeric(x)for x in live));assert owned_peak<=40000,('OWNED_NUMERIC_STOP',owned_peak)
 def save(name,data,*live):
  nonlocal RAW,saved_numeric
  own(data,*live);raw=(json.dumps(encode(data),separators=(',',':'))+'\n').encode();RAW+=len(raw);assert RAW<=2000000
  (out/name).write_bytes(raw);saved_numeric+=numeric(data)
  (out/'phase.json').write_text(json.dumps({'phase':name,'seconds':time.monotonic()-start,'owned_numeric_peak':owned_peak,'RAW_committed':RAW})+'\n')
 try:
  PRELIMIT=83;MAXPRE=MAXVALUE=0;given={};conversion=[]
  for i,arr in C['provided105_originalID_r_s_coefficients10']:
   q=list(map(parse,arr));v=[q[0],add(q[1],q[2]),mul(q[2],Q(1,5)),mul(add(q[4],mul(Q(2),q[5])),Q(1,5)),mul(q[5],Q(1,25)),mul(add(q[8],mul(Q(3),q[9])),Q(1,25))]
   given[i]=v;conversion.append([i,v])
  assert len(given)==105;records=C['operator_records_originalID_op_args'];uses=Counter()
  for r in records:
   if r[1]not in['c','input']:uses.update(r[2:])
  uses[38648]+=1;cache={};peak=0;unit_checks=[]
  for r in records:
   i,op=r[:2]
   if op=='input':v=given[i]
   elif op=='c':
    v=zero();a,b=r[2:];assert b;before(max(abs(a).bit_length(),abs(b).bit_length()));v[0]=check(Q(a,b))
   elif op=='a':v=[add(a,b)for a,b in zip(cache[r[2]],cache[r[3]])]
   elif op=='n':v=[neg(a)for a in cache[r[2]]]
   elif op=='m':v=conv(cache[r[2]],cache[r[3]])
   elif op=='i':v=inv(cache[r[2]]);unit_checks.append([i,r[2],'both_products_identity'])
   else:raise AssertionError(op)
   cache[i]=v;peak=max(peak,len(cache));own(cache,given,uses,conversion)
   if op not in['c','input']:
    for j in r[2:]:
     uses[j]-=1
     if not uses[j]:del cache[j]
   if not uses[i]:del cache[i]
  assert peak==120 and unit_checks==[[21179,21178,'both_products_identity']]
  result=cache[38648];assert result[:5]==list(map(parse,C['expected_root_lower_B6_order_1_r_ell_rell_ell2']))
  source_root_pre,source_root_bits=MAXPRE,MAXVALUE
  phase1={'root_original_ID':38648,'same_given105_count':105,'all5_lower_root_overlaps':True,'root_B6_coefficients':result,'requested_bB63_r_ell2':result[5],
    'all105_converted_input_fingerprints_sha256':hashlib.sha256(json.dumps(encode(conversion),separators=(',',':')).encode()).hexdigest(),
    'unit_inverse_checks':unit_checks,'source_nodes':2194,'lastuse_peak':peak,'max_root_precheck_bits':MAXPRE,'max_root_committed_bits':MAXVALUE}
  save('phase1-single-root-existing-cache-readback.json',phase1,given,cache,uses,conversion)
  lam=result[5];del given,cache,uses,conversion,result
  PRELIMIT=207;inputs={tuple(k):parse(v)for k,v,loc in H['all_input_values_already_cached']};inputs[('NEW_bB63_lambda12',)]=lam
  expressions=H['head1_genuine_high_A_B_UNPAIRED']+H['head1_canonical_socle_coordinates_rg_k2_UNPAIRED']+[p for pair in H['five_head_low_rows_rH_ellH_UNPAIRED']for p in pair]
  values=[];termwise=[]
  for j,poly in enumerate(expressions):
   val=Q(0);terms=[]
   for mon,w in poly:
    term=parse(w)
    for label in mon:term=mul(term,inputs[tuple(label)])
    val=add(val,term);terms.append([mon,w,term])
   values.append(val);termwise.append({'output':j,'terms':terms,'value':val});own(inputs,values,termwise)
  A,B,rg,k2=values[:4];low=[values[4+2*j:6+2*j]for j in range(5)]
  assert rg==mul(Q(5,3),A)and k2==add(B,mul(Q(-121,48),A))
  assert all(low[j]==[0,0]for j in[1,2,3])
  chi=[add(a,mul(Q(-10),b))for a,b in low];assert chi[:4]==[0,0,0,0]and chi[4]
  witness=mul(low[0][1],k2);positive=bool(witness)
  final={'head23_genuine_high_A_B':[A,B],'head23_canonical_J2_k2_socle_coordinates':[rg,k2],
    'fivehead_low_rH_ellH_rows':low,'actual_cyclic_top_chi':chi,
    'all14_output_termwise_products':termwise,'head23_lowtwo_literal_zero_verified':True,
    'simple_target_h0_h1_minor':witness,'simple_target_rank2_witness_nonzero':positive,
    'actual_ker_M7_to_G4_module_type_if_positive':'J2(r) plus k'if positive else None,
    'negative_scope':'Ifwitnesszero,thissubsetonlynegative; ker(M→G)typeisNOTselected.',
    'remaining_other_head_higher_coords':'notcomputed/formal; no arbitrary targetcoordinatezero-fill',
    'no_original_H0_Ktorsion_exponents_or_nilcube_claim':True,
    'max_root_precheck_bits':source_root_pre,'max_root_committed_bits':source_root_bits,
    'max_all_precheck_bits':MAXPRE,'max_all_committed_bits':MAXVALUE,'owned_numeric_peak':owned_peak,'RAW_committed_before_final':RAW}
  save('phase2-full-coherent-head23-and-low-two-values.json',final,inputs,values)
  rep={'id':STEM,'all_passed':True,'status':'ACTUAL_CACHED_SINGLE_HEAD23_AND_LOWTWO_VALUES_PENDING_INDEPENDENT_PEER',
    'head23_values_computed':True,'all_other_head_higher_values_computed':False,
    'simple_target_rank2_witness_nonzero':positive,'kernel_M_type_if_witness_positive':'J2 plus k'if positive else None,
    'input_sha256':pins,'code_sha256':sha(Path(__file__)),'evidence_sha256':{p.name:sha(p)for p in out.glob('phase[12]-*.json')},
    'seconds':time.monotonic()-start,'whole_numeric_peak':owned_peak,'max_QQ_precheck_bits':MAXPRE,'RAW_committed':RAW}
  save('report.json',rep);signal.alarm(0);print(json.dumps(rep))
 except BaseException as exc:
  failure={'status':'STOP_PRESERVE_NO_AUTOMATIC_RETRY','error_type':type(exc).__name__,'error':str(exc),'seconds':time.monotonic()-start,
    'phase1_exists':(out/'phase1-single-root-existing-cache-readback.json').exists(),'max_precheck':MAXPRE,'owned_peak':owned_peak,'RAW':RAW}
  (out/'failure.json').write_text(json.dumps(failure,indent=2)+'\n');raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-peer-report');p.add_argument('--head-source-peer-report');p.add_argument('--self-test',action='store_true');a=p.parse_args()
 if a.self_test:self_test()
 else:
  assert a.source_peer_report and a.head_source_peer_report,'BOTH SOURCE peers and parentONEGO required; preparedworkerUNRUN'
  main(a.source_peer_report,a.head_source_peer_report)
