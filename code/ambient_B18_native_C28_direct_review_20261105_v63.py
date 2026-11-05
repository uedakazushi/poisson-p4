"""Independent direct rational-root evaluation with a precomputed trace table."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,resource,signal,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import count

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));o=ROOT/'reports/AMBIENT_B18_NATIVE_C28_DIRECT_REVIEW_20261105_v63';o.mkdir(exist_ok=False)
 primary=ROOT/'reports/AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_VALUES_20261105_v63/report.json';pr=json.loads(primary.read_text());assert pr['all_passed'];ep=ROOT/'reports/AMBIENT_B18_ELL2_NATIVE_CANONICAL33_INPUT_GATE_V2_20261105_v63/canonical33-to-B17-C68-G8-source-input-and-response-budget.json';d=json.loads(ep.read_text());pp=ROOT/'reports/AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_SOURCE_V2_20261105_v63/native-P7-C28-canonical-third-source-and-positive34-response-input.json';pd=json.loads(pp.read_text());sel=pd['P7_raw33_selector_indices'];bt={(i,j):k for i,j,k in pd['P7_unital_multiplication_entries']}
 kap=[{4*j+m:Q(v[4*i+m])for j,i in enumerate(sel)for m in range(4)if Q(v[4*i+m])}for v in d['actual_kappa9_C132_cached_UNEVALUATED']];nu=[{j:Q(v[i])for j,i in enumerate(sel)if Q(v[i])}for v in d['actual_nu_A_B33_cached_UNEVALUATED']]
 # Positive support paths first, exact numeric budget before monic/G8 arithmetic.
 routes={}
 for u in range(28):
  i,m=divmod(u,4)
  for v in range(28):
   j,n=divmod(v,4);b=bt.get((i,j))
   if b is None:continue
   paths=[(b,())]
   for side in range(2):
    if m&n&(1<<side):paths=[(bt[x,h],steps+((side,h),))for x,steps in paths for h in nu[side]if (x,h)in bt]
   if paths:routes[u,v]=[(4*x+(m^n),steps)for x,steps in paths]
 bound=count(d)+8000+3*sum(len(v)for v in routes.values())+2*len(routes)+3*28*(9+40+8+6);assert bound<60000
 operands=[q for v in kap for q in v.values()]+[q for v in nu for q in v.values()];DD=__import__('math').lcm(*(q.denominator for q in operands));UU=max(abs(q.numerator)*(DD//q.denominator)for q in operands);coarse=6*max(DD.bit_length(),UU.bit_length())+21+100;assert coarse<=313
 (o/'pre-arithmetic-source.json').write_text(json.dumps(dict(owned_upper=60000,source_derived_owned=bound,committed_prior313=True,temporary_prior940=True,RAW_upper4000000=True,positive_table_paths=sum(map(len,routes.values())),algorithm='materialized trace-basis multiplication; direct five rational roots with geometric P inverse powers plus three original core preparations; no expanded G8 program or producer power cache'),indent=2)+'\n')
 def ck(x):assert max(abs(x.numerator).bit_length(),x.denominator.bit_length())<=313;return x
 table={}
 for key,paths in routes.items():
  row={}
  for k,steps in paths:
   x=Q(1)
   for side,h in steps:x=ck(x*nu[side][h])
   row[k]=ck(row.get(k,Q(0))+x)
  table[key]={k:x for k,x in row.items()if x}
 def add(A,B,q=Q(1)):
  C=A.copy()
  for k,x in B.items():C[k]=ck(C.get(k,Q(0))+ck(q*x))
  return {k:x for k,x in C.items()if x}
 def mul(A,B):
  C={}
  for i,x in A.items():
   for j,y in B.items():
    if (i,j)not in table:continue
    for k,v in table[i,j].items():C[k]=ck(C.get(k,Q(0))+ck(ck(x*y)*v))
  assert time.monotonic()-start<45;return {k:x for k,x in C.items()if x}
 one={0:Q(1)};core=[0,1,2,6];powcache={(i,0):one for i in core}
 def power(i,n):
  if (i,n)not in powcache:powcache[i,n]=mul(power(i,n-1),kap[i])
  return powcache[i,n]
 dsq=mul(kap[1],kap[1]);invP={};v=one
 for k in range(3):invP=add(invP,v,Q((-1)**k,5**(k+1)));v=mul(v,dsq)
 IP=[one,invP]
 for _ in range(2):IP.append(mul(IP[-1],invP))
 Kp=ROOT/'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json';Pp=ROOT/'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json';K=json.loads(Kp.read_text());prep=json.loads(Pp.read_text())
 def evaluate(rows):
  R={}
  for ex,n,den in rows:
   v=one
   for i,k in reversed(list(zip(core,ex))):v=mul(v,power(i,k))
   R=add(R,v,Q(n,den))
  return R
 result=[]
 for name in K['selected5_pivot_order']:
  f=K['five_actual_rational_roots'][name];assert f['D_power']==0
  result.append(add(kap[K['source9_variable_order'].index(name)],mul(evaluate(f['numerator']),IP[f['P_power']]),Q(-1)))
 result += [evaluate(rows)for rows in prep['generators3']]
 target=primary.parent/'actual32-G8-r2ell-and-two-monic-cached-coefficients.json';assert sha(target)==pr['evidence_sha256'][target.name];expect=json.loads(target.read_text())
 dense=[[str(R.get(24+i,Q(0)))for i in range(4)]for R in result];assert dense==expect['actual_G8_r2ell_eight_by4']
 assert [str(v.get(6,Q(0)))for v in nu]==expect['actual_nu_A_B_r2ell']
 rr=dict(all_passed=True,independent_original_rational_roots_and_preparation_all32_coefficients=True,all32_G8_and2_monic_third_coefficients_verified=True,no_producer_expanded_program_or_same_multiplication_algorithm=True,source_owned_upper=60000,source_derived_bound=bound,committed_prior313=True,no_honest_quotient_upper_HF_or_new_physical_jets=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [primary,target,ep,pp,Kp,Pp]},code_sha256=sha(Path(__file__)),seconds=time.monotonic()-start)
 (o/'report.json').write_text(json.dumps(rr,indent=2)+'\n');print(sha(o/'report.json'))
if __name__=='__main__':main()
