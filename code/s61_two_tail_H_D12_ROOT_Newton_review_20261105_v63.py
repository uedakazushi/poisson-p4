"""Independent Newton/Faddeev trace determinant over the SAME quartic ring."""
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1];P=32003
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports/S61_TWO_TAIL_H_D12_ROOT_NEWTON_REVIEW_20261105_v63';out.mkdir(exist_ok=False);pins={};ops=0
 def read(n,h=None):
  p=ROOT/n;s=sha(p)
  if h:assert s==h,n
  pins[n]=s;return json.loads(p.read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 cust=read('reports/S61_TWO_TAIL_H_D12_AXIS_CUSTODY_20261105_v63/report.json','4d2044ed5e3220730ea0564819c52d029531b02592d1d4e909e81b1334a74472')
 pri=read('reports/S61_TWO_TAIL_H_D12_AXIS_QUOTIENT_ENTRY_20261105_v63/report.json','23646500355809a05d19af1774fa33ce75a8ec56f0e31151f52fa56ed76063c8');assert cust['all_passed']and pri['all_passed']
 m=read(cust['matrix_input'],cust['runtime_input_sha256'][cust['matrix_input']]);assert len(m['matrix'])==12 and sum(len(c)for row in m['matrix']for c in row)==726
 assert max(abs(int(v)).bit_length()for row in m['matrix']for cell in row for k,v in cell)<=110
 assert all(0<=k<=8 for row in m['matrix']for cell in row for k,v in cell)
 owned=6000+2*(1453+3*576+12*576+48);assert owned==26282 and owned<28000
 raw=sum((ROOT/p).stat().st_size for p in pins)+Path(__file__).stat().st_size+600000;assert raw<1048576
 dump('pre-arithmetic-Newton-SOURCE.json',dict(input_sha256=pins,fixed_same_p=P,fixed_same_axis=[0,0],ring='Fp[a]/monicHaxis dimension4; notassumedfield',algorithm='Faddeev-LeVerrier/Newton trace recursion; integers1..12 areunits since p32003',ring_products=20736,field_product_upper=600000,hard_field_products=750000,field_temp_bits=30,ROM_input_bits=110,ROM_LCD_bits=519,owned_with_failure_copy=owned,hard_owned=28000,RAW_prior=raw,soft_seconds=45,hard_seconds=60))
 h=[v%P for v in cust['H_axis_integer']];assert len(h)==5 and h[-1]==7157 and h[0]==30819
 il=pow(h[-1],P-2,P);h=[v*il%P for v in h]
 state=dict(stage='ring',matrix=[],B=[],T=[],B_stages=[],characteristic_coefficients=[],cursor={})
 def mul(a,b):
  nonlocal ops
  v=a*b;assert v.bit_length()<=30;ops+=1;assert ops<=750000
  if ops%4096==0:assert time.monotonic()-start<45
  return v%P
 def rm(a,b):
  # Reverse both coefficient orders from the subset producer.
  c=[0]*7
  for i in range(3,-1,-1):
   for j in range(3,-1,-1):c[i+j]=(c[i+j]+mul(a[i],b[j]))%P
  for k in range(6,3,-1):
   q=c[k]
   for j in range(3,-1,-1):c[k-4+j]=(c[k-4+j]-mul(q,h[j]))%P
  return c[:4]
 try:
  powers=[[1,0,0,0]]
  for _ in range(8):powers.append(rm(powers[-1],[0,1,0,0]))
  A=[]
  for row in m['matrix']:
   rr=[]
   for cell in row:
    v=[0]*4
    for k,c in reversed(cell):
     for j in range(4):v[j]=(v[j]+mul(int(c)%P,powers[k][j]))%P
    rr.append(v)
   A.append(rr)
  state['matrix']=A;B=[[[int(i==j),0,0,0]for j in range(12)]for i in range(12)];state['B']=B
  for step in range(1,13):
   state['stage']='Newton_trace';state['step']=step
   T=[[[0]*4 for j in range(12)]for i in range(12)];state['T']=T
   for i in range(12):
    for j in range(12):
     for k in range(11,-1,-1):
      state['cursor']=[step,i,j,k];v=rm(A[i][k],B[k][j]);T[i][j]=[(a+b)%P for a,b in zip(T[i][j],v)]
   iv=pow(step,P-2,P);ops+=30
   c=[mul((-sum(T[i][i][j]for i in range(12)))%P,iv)for j in range(4)]
   state['characteristic_coefficients'].append(c)
   for i in range(12):T[i][i]=[(a+b)%P for a,b in zip(T[i][i],c)]
   B=T;state['B']=B;state['B_stages'].append(B)
  assert not any(v for row in B for cell in row for v in cell),'Cayley-Hamilton residual'
  det=state['characteristic_coefficients'][-1] # n=12 is even.
  lcd=int(m['LCD_product']);assert lcd%P==30563;iv=pow(lcd%P,P-2,P);ops+=30
  rem=[mul(v,iv)for v in rm(det,powers[5])]
  assert det==pri['detC_remainder']and rem==pri['D12_remainder']and any(rem)
  dump('all-twelve-Newton-stages-and-Cayley-Hamilton-zero.json',state)
  rep=dict(id=out.name,all_passed=True,independent_method='12 Newton/Faddeev matrixproducts/traces; no subsetdeterminants and no producercode import',all12_Cayley_Hamilton_recursion_stages_saved=True,Cayley_Hamilton_matrix_zero=True,
   determinant_C12_remainder=det,D12_remainder=rem,H_nondivides_D12_in_BASE6_independently_proved=True,
   proof='H coprimeU6 and D12 representedusingonlya/Bdenominators. HypotheticallocalizedHdivisibility yieldsHdividinga/BclearedpolynomialN byUFDcancellation. LegalfixedaxisretainsnonzeroHaxisdegree4,aunit,Baxis−40andLCDunit; ourNONzeroD12quotient contradicts divisibility.',
   no_field_or_etale_assumption=True,no_actual_bad_Fitting_or_rank_drop_claim=True,no_extra_axis_prime_or_source_oracle=True,field_products_with_inverse_priors=ops,source_owned_with_failure_copy=owned,RAW_source_prior=raw,
   input_sha256=pins,code_sha256=sha(Path(__file__)),evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
  dump('report.json',rep);print(json.dumps(rep))
 except BaseException as e:dump('full-Newton-failure-prefix.json',dict(error=repr(e),state=state,input_sha256=pins,no_retry=True));raise
if __name__=='__main__':main()
