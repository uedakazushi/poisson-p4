#!/usr/bin/env python3
"""Independent fixed known-cubic source/minor/factorization verifier."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def pairs(v):return {int(i):F(c)for i,c in v if F(c)}
def plus(v,w,c=F(1)):
 z=dict(v)
 for i,a in w.items():
  z[i]=z.get(i,F(0))+c*a
  if not z[i]:del z[i]
 return z
def sp(v):return [[i,c]for i,c in sorted(v.items())]
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json',h=None):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p)
  if h:assert sha(p)==h
  return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63';rep=load(primary,h='504583d643c70c463bdf0692acf89010847b3c7ce08f88ec15eefd52d4155cea');assert rep['all_passed']
 for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 for n,h in rep['evidence_sha256'].items():assert sha(ROOT/'reports'/primary/n)==h
 assert load('AMBIENT_B18_CPRIME_CUBIC_ELL3_SOURCE_20261105_v63_REVIEW',h='daf51e0a6992ec7a5e12c1c0e730758aece2a5af53499f8e4776ece263741085')['all_passed']
 assert load('AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63_REVIEW',h='5faf10001ed74c2b9e107d7fa1f4fe873600423affa75305a75a8f1d3fad5bb3')['actual_Cprime_HF_through2']==[3,3,3]
 e=load(primary,'phase2-complete-streamed-cubic-upper-certificate.json',h='256fe65cd0f6b47b62fea8baea845c614c28ac5a24dd78a71eb957397e24b1ae');phase=load(primary,'phase1-streamed-fixed-pivots.json')
 val=load('AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63','phase4-all-twenty-actual-quadratic-values.json');L=[pairs(v)for v in load('AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63','phase1-actual-fullnine-first-V-columns.json')['all44_actual_linear_columns27_sparse']];q4=[pairs(v)for v in load('AMBIENT_B18_CPRIME_PRODUCT_ANN_QUOTIENT_SOURCE_PREFLIGHT_20261105_v63')['literal_fixed_pi4_rows_by135_sparse']]
 v=list(map(F,val['all20_actual_quadratic_values4'][18]));q=[plus(q4[i],q4[1],-v[i]/v[1])for i in [0,2,3]];E=[44,133,134];m2=list(itertools.combinations_with_replacement(range(9),2));i2={m:i for i,m in enumerate(m2)};m3=list(itertools.combinations_with_replacement(range(9),3));i3={m:i for i,m in enumerate(m3)};P={}
 for j,c in pairs(val['source_constant_kernel20'][18]).items():
  for t in range(3):
   for k,mon in enumerate([(7,7),(7,8),(8,8)]):P=plus(P,{45*t+i2[mon]:c*F(val['all44_actual_V2_pureplane_nine_values'][j][3*t+k])})
 proj=[]
 for t in range(3):
  for mon in m3:
   d=45*t+i2[mon[:2]];proj.append({9*j+mon[2]:row[d]for j,row in enumerate(q)if d in row})
 assert q==[pairs(v)for v in e['exact_NF3_rows_by135_sparse']]and P==pairs(e['actual_kernel18_quadratic_representative135_sparse'])and E==e['unit_embedding_rows135']and proj==[pairs(v)for v in e['all495_source_projection_columns27_sparse']]
 def project(raw):
  z={}
  for i,c in raw.items():z=plus(z,proj[i],c)
  return z
 def stream():
  for j,col in enumerate(L):
   for mon in m2:
    raw={}
    for i,c in col.items():t,x=divmod(i,9);raw=plus(raw,{165*t+i3[tuple(sorted(mon+(x,)))]:c})
    yield ['first44_times_S2',j,list(mon)],raw,project(raw)
  for x in range(9):
   raw={}
   for i,c in P.items():t,j=divmod(i,45);raw=plus(raw,{165*t+i3[tuple(sorted(m2[j]+(x,)))]:c})
   yield ['actual_kernel18_quad_times_S1',x],raw,project(raw)
  for t in range(3):
   raw={165*t+i3[(8,8,8)]:F(1)};yield ['ell_cube',t],raw,project(raw)
  for j,row in enumerate(E):
   t,k=divmod(row,45)
   for x in range(9):yield ['consistency',j,x],{},plus({9*j+x:F(1)},proj[165*t+i3[tuple(sorted(m2[k]+(x,)))]],-1)
 cp=e['fixed_columns'];rr=e['fixed_rows'];n=len(cp);assert n==len(rr)==26 and set(range(27))-set(rr)=={25};picked={};sourcehash=hashlib.sha256();count=0
 for j,(label,raw,col)in enumerate(stream()):
  sourcehash.update(json.dumps(enc([label,sp(raw),sp(col)]),separators=(',',':')).encode()+b'\n');count+=1
  if j in cp:picked[j]=(label,col)
 assert count==2019 and sourcehash.hexdigest()==e['all2019_source_and_projection_stream_sha256']==phase['source_projection_stream_sha256']
 basis=[picked[j][1]for j in cp];assert basis==[pairs(v)for v in e['selected_original_compressed_columns']]==[pairs(v)for v in phase['selected_original_compressed_columns']]and [picked[j][0]for j in cp]==e['fixed_labels']
 A=[[col.get(i,F(0))for col in basis]for i in rr];iv=[[F(c)for c in row]for row in e['fixed_inverse']];assert A==[list(map(F,row))for row in e['fixed_minor']]
 assert all(sum((A[i][k]*iv[k][j]for k in range(n)),F(0))==F(i==j)and sum((iv[i][k]*A[k][j]for k in range(n)),F(0))==F(i==j)for i in range(n)for j in range(n))
 # Fixed integer Bareiss, no row/source pivot search and no inverse solve.
 den=13365000000;C=[[int(c*den)for c in row]for row in A];assert all(F(C[i][j],den)==A[i][j]for i in range(n)for j in range(n));prev=1;steps=[]
 for k in range(n-1):
  pivot=C[k][k];assert pivot;divisions=0
  for i in range(k+1,n):
   for j in range(k+1,n):a=C[i][j]*pivot-C[i][k]*C[k][j];quot,rem=divmod(a,prev);assert not rem;C[i][j]=quot;divisions+=1
  for i in range(k+1,n):C[i][k]=0
  steps.append([k,str(pivot),str(prev),divisions]);prev=pivot
 det=F(C[-1][-1],den**n);assert det==F(e['fixed_minor_determinant'])==F(-61,18000)
 ih=[{i:c for i,c in enumerate(row)if c}for row in iv];fh=hashlib.sha256()
 for label,raw,col in stream():
  co=[sum((c*col.get(rr[j],F(0))for j,c in row.items()),F(0))for row in ih];back={}
  for p,c in zip(basis,co):back=plus(back,p,c)
  assert back==col;fh.update(json.dumps(enc(co),separators=(',',':')).encode()+b'\n')
 assert fh.hexdigest()==e['all2019_fixed_factorization_sha256']
 dual={25:F(1)}
 for j,i in enumerate(rr):
  c=sum((basis[k].get(25,F(0))*iv[k][j]for k in range(n)),F(0))
  if c:dual[i]=-c
 assert all(sum((c*col.get(i,F(0))for i,c in dual.items()),F(0))==0 for label,raw,col in stream())
 pull={}
 for i,col in enumerate(proj):
  c=sum((a*col.get(j,F(0))for j,a in dual.items()),F(0))
  if c:pull[i]=c
 assert not rep['positive_Cprime_ell3_length9_and_long_exponent2']and e['known_cubic_quotient_dimension']==1
 save('independent-known-cubic-fixed-certificate-and-open-dual.json',dict(all2019_source_projection_sha256=sourcehash.hexdigest(),all2019_factorization_sha256=fh.hexdigest(),fixed_rank=26,fixed_minor_determinant=det,fixed_integer_Bareiss_steps=steps,both_inverse_products=True,all_selected_source_labels_and_matrix_match=True,known_cubic_quotient_dimension=1,known_quotient_dual27=dual,known_quotient_dual495=pull,actual_cubic_class_or_long_exponent_NOT_decided=True,full_cubic_values_or_higherjets_NOT_evaluated=True))
 for p in[Path(__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 raw=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir()if p.is_file());assert raw<32<<20 and time.monotonic()-st<60
 save('report.json',dict(id=ID,all_passed=True,frozen=True,known_cubic_rank=26,known_cubic_quotient_dimension=1,all_source_fixed_minor_inverse_factorization_gates_verified=True,upper9_and_long_exponent2_not_proved=True,no_actual_HF3_or_e_gt2_inference=True,no_new_source_jet_or_rank_search=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},actual_RAW=raw,seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
 print(json.dumps(dict(all_passed=True,known_rank=26,known_cokernel=1,long_exponent_open=True,seconds=time.monotonic()-st),indent=2))
if __name__=='__main__':main()
