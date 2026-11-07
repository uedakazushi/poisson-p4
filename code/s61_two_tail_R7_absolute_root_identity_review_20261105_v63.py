"""Independent sparse root-power and failed-coefficient R7 absolute certificate."""
from pathlib import Path
import hashlib,json,time,resource,signal
R=Path(__file__).resolve().parents[1];P=32003;Z=(0,0);ONE=(1,0)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out=R/'reports/S61_TWO_TAIL_R7_ABSOLUTE_ROOT_IDENTITY_REVIEW_20261105_v63'
out.mkdir(exist_ok=False);start=time.monotonic();pins={};products=0;bits=0
resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
state={'phase':'SOURCE','mask':None,'root':None,'candidate':None}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;return json.loads(q.read_text())
def bm(a,b):
 global products,bits
 assert 0<=a<P and 0<=b<P;z=a*b
 products+=1;bits=max(bits,z.bit_length());assert products<=100000 and bits<=30
 return z%P
def power(a,n):
 z=1
 while n:
  if n%2:z=bm(z,a)
  a=bm(a,a);n//=2
 return z
def add(a,b):return((a[0]+b[0])%P,(a[1]+b[1])%P)
def neg(a):return((-a[0])%P,(-a[1])%P)
def km(a,b):
 # Ordinary four cross products, not producer Karatsuba.
 ac=bm(a[0],b[0]);bd=bm(a[1],b[1])
 return((ac-bm(q0,bd))%P,(bm(a[0],b[1])+bm(a[1],b[0])-bm(q1,bd))%P)
def scale(a,n):return(bm(a[0],n),bm(a[1],n))
def insert(d,k,v):
 z=add(d.get(k,Z),v)
 if z==Z:d.pop(k,None)
 else:d[k]=z
def smul(a,b):
 z={}
 for e,x in reversed(list(a.items())):
  for f,y in reversed(list(b.items())):
   if e+f<7:insert(z,e+f,km(x,y))
 return z
def allocation(poly,b):
 if not b:return dict(poly)
 z={}
 for (a,e),x in poly.items():
  insert(z,(a,e),scale(x,55))
  if e<6:insert(z,(a,e+1),scale(x,66))
 return z
try:
 pr=load('reports/S61_TWO_TAIL_R7_ABSOLUTE_REUSE_ENTRY_20261105_v63/report.json','93be721c249fa0340f8e66027c1d77a13ad3ab9a9b9e55afee50b867ada69b13')
 assert pr['all_passed'] and pr['R7_geometrically_irreducible_PRIMARY']
 name='four-cached-jets-and14-products-certificate.json'
 c=load('reports/S61_TWO_TAIL_R7_ABSOLUTE_REUSE_ENTRY_20261105_v63/'+name,pr['evidence_sha256'][name])
 q0,q1=c['q'];roots=[[tuple(x)for x in row]for row in c['roots']]
 failures=c['failures'];assert len(roots)==4 and all(len(x)==7 for x in roots)
 del c
 f=load('reports/S61_TWO_TAIL_H_EXACT_QQ_ENTRY_20261105_v63/canonical-factor-candidates.json','4dea15b612a158ff9fedfac45b493c1ae3c1c0218f6b196842b242577a2bd979')
 native=next(x for x in f['factors']if x['total_degree']==7)
 terms=[(tuple(k),int(v))for k,v in native['terms']];del f,native
 (out/'phase.json').write_text(json.dumps(dict(SOURCE_prior=True,method='naiveFp2/sparse direct powers, fresh reverse-root subset products, one exact failed coefficient per each14candidate; no producer helpers or seed/Hensel replay',fieldbits30=True,productcap100000=True,owned12000=True,RAW2097152=True))+'\n')
 assert all(P%d for d in range(2,179))
 disc=(bm(q1,q1)-bm(4,q0))%P;assert disc and power(disc,(P-1)//2)==P-1
 alpha=140625%P;ia=power(alpha,P-2);assert bm(alpha,ia)==1
 target={}
 for (a,e,d),v in terms:insert(target,(a,e),(bm(v%P,ia),0))
 assert {e:x for (a,e),x in target.items()if a==4}=={0:(55,0),1:(66,0)}
 assert len(set(row[0]for row in roots))==4
 residuals=[]
 for j,row in enumerate(roots):
  state.update(phase='root_powers',root=j)
  series={e:x for e,x in enumerate(row)if x!=Z};powers=[{0:ONE}]
  for a in range(1,5):powers.append(smul(powers[-1],series))
  residual={}
  for (a,e),v in target.items():
   for f,x in powers[a].items():
    if e+f<7:insert(residual,e+f,km(v,x))
  assert not residual;residuals.append(True)
 subsets={}
 for mask in range(16):
  state.update(phase='reverse_sparse_subsets',mask=mask);poly={(0,0):ONE}
  for j in reversed(range(4)):
   if not(mask>>j&1):continue
   z={}
   for (a,e),v in poly.items():
    insert(z,(a+1,e),v)
    for f,x in enumerate(roots[j]):
     if e+f<7 and x!=Z:insert(z,(a,e+f),neg(km(v,x)))
   poly=z
  subsets[mask]=poly
 assert allocation(subsets[15],1)==target
 expected={(m,b)for m in range(1,15,2)for b in (0,1)}
 assert {(x['mask'],x['LC_allocation'])for x in failures}==expected and len(failures)==14
 witnessed=[]
 for record in reversed(failures):
  mask=record['mask'];b=record['LC_allocation'];a=record['a'];e=record['e']
  state.update(phase='independent_failure_witness',candidate=[mask,b,a,e])
  G=allocation(subsets[mask],b);H=allocation(subsets[15^mask],1-b);value=Z
  for (i,j),x in G.items():
   y=H.get((a-i,e-j))
   if y is not None:value=add(value,km(x,y))
  diff=add(value,neg(target.get((a,e),Z)))
  assert diff!=Z and diff==tuple(record['difference'])
  witnessed.append(dict(mask=mask,LC_allocation=b,a=a,e=e,nonzero_difference=diff))
 ep=out/'independent-sparse-roots-subsets-and14-nonzero-coefficients.json'
 ep.write_text(json.dumps(dict(all28_R7_root_coefficients_zero=True,full_L_times_four_root_product=True,all14_candidate_failure_witnesses=witnessed),indent=2)+'\n')
 r=dict(id=out.name,all_passed=True,input_sha256=pins,independent_naive_field_sparse_finite_certificate=True,
  all28_root_coefficients_zero=True,all14_geometric_factor_candidates_excluded=True,
  R7_geometrically_irreducible_proved=True,H_has_exactly_two_geometric_irreducible_factors_R7_S4=True,
  theorem='Fp2 is a field by prime/discriminant certificate. Four distinct roots plus unitLC55 give unique simple-root Hensel factors. Primitive a-polynomial Rbar has e-degree6, LC55+66e: every geometric factor normalizes leading1orLC, is a seed subset product overFp2 and hasedegree<=6. All14 possibilities excluded by independent nonzerofull-product coefficient. Preserveda4 andgeometricglobala-content then implychar0/globalR7 absoluteirreducibility.',
  API_factor_irreducibility_contract_not_used=True,no_badFitting_full57_239_or_newPoisson_component_claim=True,
  field_products=products,field_temporary_bits=bits,hard_field_bits=30,owned_upper=12000,
  code_sha256=sha(Path(__file__)),evidence_sha256={ep.name:sha(ep),'phase.json':sha(out/'phase.json')},seconds=time.monotonic()-start)
 p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n')
 assert sum((R/k).stat().st_size for k in pins)+sum(p.stat().st_size for p in out.iterdir())<2097152
 print(sha(p),products)
except BaseException as e:
 (out/'failure-prefix.json').write_text(json.dumps(dict(error=repr(e),state=state,no_retry=True),indent=2)+'\n')
 raise
