#!/usr/bin/env python3
"""Independent binomial transport/dictionary B6/geometric inverse VALUE peer."""
import hashlib,json,math,pathlib,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
import ambient_B18_C_connecting_head23_cache_source_review_20261105_v63 as POS
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63';ID=PRI+'_REVIEW'
MON=[(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)];SMON=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(3,0),(2,1),(1,2),(0,3)]
PRE=0;STORED=0;LIMIT=83
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def check(q):
 global STORED
 STORED=max(STORED,abs(q.numerator).bit_length(),q.denominator.bit_length());assert STORED<=103
 return q
def calc(a,b,mul=False):
 global PRE
 na,da=abs(a.numerator).bit_length(),a.denominator.bit_length();nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
 h=max(na+nb,da+db) if mul else max(max(na+db,nb+da)+1,da+db)
 PRE=max(PRE,h);assert h<=LIMIT<=207<2048
 return check(a*b if mul else a+b)
def parse(x):
 q=Q(x);assert max(abs(q.numerator).bit_length(),q.denominator.bit_length())<=LIMIT;return check(q)
def plus(a,b):
 d=dict(a)
 for k,v in b.items():
  d[k]=calc(d.get(k,Q(0)),v)
  if not d[k]:del d[k]
 return d
def scale(a,c):return {k:calc(v,c,True) for k,v in a.items() if v and c}
def times(a,b):
 d={}
 for (r,e),x in a.items():
  for (s,f),y in b.items():
   if r+s<2 and e+f<3:
    k=(r+s,e+f);d[k]=calc(d.get(k,Q(0)),calc(x,y,True))
 return {k:v for k,v in d.items() if v}
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,list):return [enc(v) for v in x]
 if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
 return x
def main():
 global LIMIT
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None,read=True):
  nonlocal RAW
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h,(str(p),s,h)
  k=str(p.relative_to(ROOT))
  if k not in pins:RAW+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text()) if read else None
 pr=pin(ROOT/'reports'/PRI/'report.json','7d51e4a12f69f4c760cc547c4f83cb844bdd9dd67074c0ed78276e6fd64c8e53')
 for p,h in pr['input_sha256'].items():pin(ROOT/p,h,False)
 pin(ROOT/'code/ambient_B18_C_connecting_head23_cached_values_20261105_v63.py',pr['code_sha256'],False)
 for name,h in pr['evidence_sha256'].items():pin(ROOT/'reports'/PRI/name,h,False)
 def get(s):return ROOT/next(p for p in pr['input_sha256'] if p.endswith(s))
 C=pin(get('compiled-single-bB63-root-and-provided-joint105-cache.json'))
 P=pin(get('single-root-head23-coherent-positive-height-and-whole-budget.json'))
 H=pin(get('full-coherent-head23-two-socle-cache-polynomials-and-low-two-row-source.json'))
 cp=pin(ROOT/'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHE_SOURCE_20261105_v63_REVIEW/report.json','bd4980ca33e99428dbab56797a10a5e6e3102d727c74061e4d2d0c76638f0b9a')
 hp=pin(ROOT/'reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63_REVIEW/report.json','25cb1e164eb54b85f6835a51781df9cfa0f5adae5279193e22c1a48e34d2b281')
 assert cp['all_passed'] and hp['all_passed']
 records=C['operator_records_originalID_op_args'];index={r[0]:r for r in records};given=C['provided105_originalID_r_s_coefficients10']
 # New inverse method's positive SOURCE gate: unit517-closure only, no
 # signed root/composition. The geometric nilpotent series is finite inB6.
 POS.HEIGHT=0;pb={}
 for i,arr in given:
  q=list(map(Q,arr));v=POS.zero();v[0]=(q[0],1)
  v[1]=POS.add(POS.bound(q[1]),POS.bound(q[2]));v[2]=POS.mul((Q(1,5),1),POS.bound(q[2]))
  v[3]=POS.mul((Q(1,5),1),POS.add(POS.bound(q[4]),POS.mul((Q(2),1),POS.bound(q[5]))));v[4]=POS.mul((Q(1,25),1),POS.bound(q[5]))
  v[5]=POS.mul((Q(1,25),1),POS.add(POS.bound(q[8]),POS.mul((Q(3),1),POS.bound(q[9]))));pb[i]=v
 want=set();todo=[21178]
 while todo:
  i=todo.pop()
  if i not in want:
   want.add(i);r=index[i]
   if r[1] not in ['input','c']:todo.extend(r[2:])
 assert len(want)==517;maj={}
 for i in sorted(want):
  r=index[i];op=r[1]
  if op=='input':v=pb[i]
  elif op=='c':v=POS.zero();v[0]=(Q(r[2],r[3]),1)
  elif op=='a':v=[POS.add(x,y) for x,y in zip(maj[r[2]],maj[r[3]])]
  elif op=='n':v=list(maj[r[2]]);v[0]=(-v[0][0],1)
  else:assert op=='m';v=POS.conv(maj[r[2]],maj[r[3]])
  maj[i]=v
 x=maj[21178];assert x[0][0]==1;nil=list(x);nil[0]=(Q(0),1)
 power=POS.zero();power[0]=(Q(1),1);invbound=list(power)
 for _ in range(3):power=POS.conv(power,nil);invbound=[POS.add(a,b) for a,b in zip(invbound,power)]
 POS.conv(x,invbound);POS.conv(invbound,x)
 assert POS.HEIGHT<=41 and 2*POS.HEIGHT+1<=83
 owned=count([pr,C,P,H,cp,hp,pins])+2*2194+105*20+120*20+15000
 assert owned<60000 and RAW+1000000<2000000
 save(out/'positive-before-independent-geometric-B6-values.json',{'root_original_ID':38648,'same2194_provided105':True,'different_unit_method_positive517closure':True,'geometric_inverse_and_both_diagnostic_H':POS.HEIGHT,'root_prelimit':83,'post_prelimit':207,'own_whole_owned_upper':owned,'RAW_plus_output_reserve':RAW+1000000,'no_signed_requested_root_yet':True})
 del pb,maj,want,nil,invbound,power
 # Generic binomial substitution ofALL10 given ordinary r,s coefficients.
 phi={};conversion=[]
 for i,arr in given:
  v={}
  for (nr,ns),literal in zip(SMON,arr):
   q=parse(literal)
   for k in range(ns+1):
    r,e=nr+ns-k,k
    if r<2 and e<3:
     c=Q(math.comb(ns,k),5**k);value=calc(q,c,True)
     v=plus(v,{(r,e):value})
  phi[i]=v;conversion.append([i,[v.get(m,Q(0)) for m in MON]])
 uses=Counter()
 for r in records:
  if r[1] not in ['input','c']:uses.update(r[2:])
 uses[38648]+=1;live={};peak=0;unitproof=[]
 for r in records:
  i,op=r[:2]
  if op=='input':v=phi[i]
  elif op=='c':q=check(Q(r[2],r[3]));v={(0,0):q} if q else {}
  elif op=='a':v=plus(live[r[2]],live[r[3]])
  elif op=='n':v=scale(live[r[2]],Q(-1))
  elif op=='m':v=times(live[r[2]],live[r[3]])
  else:
   assert op=='i';x=live[r[2]];assert x.get((0,0))==1
   t=dict(x);del t[0,0];p={(0,0):Q(1)};v=dict(p)
   for power in range(1,4):p=times(p,t);v=plus(v,scale(p,Q((-1)**power)))
   assert not times(p,t) and times(x,v)=={(0,0):Q(1)}==times(v,x)
   unitproof.append([i,r[2],[str(x.get(m,Q(0))) for m in MON],[str(v.get(m,Q(0))) for m in MON]])
  live[i]=v;peak=max(peak,len(live));assert peak<=120
  if op not in ['input','c']:
   for j in r[2:]:
    uses[j]-=1
    if not uses[j]:del live[j]
  if not uses[i]:del live[i]
 assert peak==120 and len(unitproof)==1 and unitproof[0][:2]==[21179,21178]
 result=[live[38648].get(m,Q(0)) for m in MON];rootpre=PRE;rootstored=STORED
 phase1=pin(ROOT/'reports'/PRI/'phase1-single-root-existing-cache-readback.json',pr['evidence_sha256']['phase1-single-root-existing-cache-readback.json'])
 assert result==list(map(Q,phase1['root_B6_coefficients']))
 assert result[:5]==list(map(Q,C['expected_root_lower_B6_order_1_r_ell_rell_ell2']))
 fingerprint=hashlib.sha256(json.dumps(enc(conversion),separators=(',',':')).encode()).hexdigest()
 assert fingerprint==phase1['all105_converted_input_fingerprints_sha256']
 assert result[5]==Q(2039,360)
 lam=result[5];del live,uses,phi,conversion,index
 LIMIT=207;inputs={tuple(k):parse(q) for k,q,loc in H['all_input_values_already_cached']};inputs['NEW_bB63_lambda12',]=lam
 polys=H['head1_genuine_high_A_B_UNPAIRED']+H['head1_canonical_socle_coordinates_rg_k2_UNPAIRED']+[p for pair in H['five_head_low_rows_rH_ellH_UNPAIRED'] for p in pair]
 values=[];termwise=[]
 for j,poly in enumerate(polys):
  terms=[None]*len(poly);value=Q(0)
  for k in reversed(range(len(poly))):
   labels,w=poly[k];q=parse(w)
   for label in reversed(labels):q=calc(q,inputs[tuple(label)],True)
   value=calc(value,q);terms[k]=[labels,w,q]
  values.append(value);termwise.append({'output':j,'terms':terms,'value':value})
 A,B,rg,k2=values[:4];low=[values[4+2*j:6+2*j] for j in range(5)]
 assert rg==calc(Q(5,3),A,True) and k2==calc(B,calc(Q(-121,48),A,True))
 chi=[calc(a,calc(Q(-10),b,True)) for a,b in low]
 assert chi[:4]==[0,0,0,0] and chi[4] and all(low[j]==[0,0] for j in [1,2,3])
 minor=calc(low[0][1],k2,True);assert minor==Q(1331,115200) and k2==Q(1331,9216) and low[0][1]==Q(2,25)
 phase2=pin(ROOT/'reports'/PRI/'phase2-full-coherent-head23-and-low-two-values.json',pr['evidence_sha256']['phase2-full-coherent-head23-and-low-two-values.json'])
 assert enc(termwise)==phase2['all14_output_termwise_products']
 assert [str(A),str(B)]==phase2['head23_genuine_high_A_B'] and [str(rg),str(k2)]==phase2['head23_canonical_J2_k2_socle_coordinates']
 assert enc(low)==phase2['fivehead_low_rH_ellH_rows'] and str(minor)==phase2['simple_target_h0_h1_minor']
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'independent_dictionary_B6_and_binomial_s_r_ell5_transport':True,'all2194_original_records_and105branch_inputs':True,
  'all105_conversion_fingerprint_and5lower_root_overlaps':True,'unit_inverse_by_finite_geometric_series':unitproof,'both_inverse_products_and_nil4':True,
  'actual_root_B6':list(map(str,result)),'actual_bB63_r_ell2':str(lam),
  'all14_reversed_fold_values_and_all_original_termwise_products_match':True,
  'head23_high_A_B':[str(A),str(B)],'canonical_rg_k2':[str(rg),str(k2)],'fivehead_low_rows':enc(low),'cyclic_chi':enc(chi),
  'simple_target_minor':str(minor),'actual_kernel_M7_to_G4_type':'J2(r) plus k',
  'type_proof':'The exact surjectionM7=J2²+k³→G4=J2+k² andpositive simpletargeth0/h1minor forceimageSocMdimension>=3. Asemisimplekernel3wouldgiveimageSocMdimension2. HencekerM3hasr-rank1,typeJ2+k. Thisis the honestquotientkernel, notKbar/IKbar.',
  'no_other_head_high_original_H0_K_exponent_or_nil_inference':True,
  'source_prelimits_root_post':[83,207],'actual_root_prebits':rootpre,'actual_root_stored_bits':rootstored,'actual_all_prebits':PRE,'actual_all_stored_bits':STORED,'owned_upper':owned}
 ep=out/'independent-geometric-cache-root-and-all14-head23-values.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'independent_SOURCE_bound_VALUE':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'actual_one_root_all14_outputs_positive_minor_PASS':True,'actual_kernel_M_type':'J2 plus k','seconds':time.monotonic()-started}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'prebits':PRE,'seconds':report['seconds']}))
if __name__=='__main__':main()
