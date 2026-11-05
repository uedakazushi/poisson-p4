"""Bounded cached G8 values in C68. --prepare SOURCE first; --run needs exact gate SHA."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,resource,signal,time,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import count,fields

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(a):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));stem='AMBIENT_B18_NATIVE_CANONICAL68_'+('VALUE_GATE'if a.prepare else'VALUES')+'_20261105_v63';o=ROOT/'reports'/stem;o.mkdir(exist_ok=False)
 source=ROOT/'reports/AMBIENT_B18_ELL2_NATIVE_CANONICAL33_INPUT_GATE_V2_20261105_v63/report.json';sr=json.loads(source.read_text());assert sr['all_passed']
 peer=ROOT/'reports/AMBIENT_B18_NATIVE_CANONICAL68_ROOT_SOURCE_PEER_20261105_v63/report.json';assert json.loads(peer.read_text())['all_passed']
 ep=source.parent/'canonical33-to-B17-C68-G8-source-input-and-response-budget.json';assert sha(ep)==sr['evidence_sha256'][ep.name];d=json.loads(ep.read_text());sel=d['B17_raw33_selector_indices'];BT={(i,j):k for i,j,k in d['B17_multiplication_nonzero_entries']};program=d['G8_formal_kappa9_polynomials_mod_radical7'];maxpow=[max(ex[i]for poly in program for ex,q in poly)for i in range(9)]
 pins={str(p.relative_to(ROOT)):sha(p)for p in [source,peer,ep,Path(__file__)]};clock=0;products=0;peak=0;state={};powers={};results=[]
 def put(n,x):
  def clean(v):
   if isinstance(v,Q):return str(v)
   if isinstance(v,dict):return {str(k):clean(w)for k,w in v.items()}
   if isinstance(v,(list,tuple)):return [clean(w)for w in v]
   return v
  raw=(json.dumps(clean(x),separators=(',',':'))+'\n').encode();assert sum(p.stat().st_size for p in o.iterdir())+len(raw)<4000000;(o/n).write_bytes(raw)
 # Parse only existing input values. No G8 response has been computed.
 kap=[{4*j+m:Q(v[4*i+m])for j,i in enumerate(sel)for m in range(4)if Q(v[4*i+m])}for v in d['actual_kappa9_C132_cached_UNEVALUATED']]
 nu=[{j:Q(v[i])for j,i in enumerate(sel)if Q(v[i])}for v in d['actual_nu_A_B33_cached_UNEVALUATED']]
 def routes(u,v):
  i,m=divmod(u,4);j,n=divmod(v,4);k=BT.get((i,j))
  if k is None:return []
  masks=m^n;rows=[(k,())]
  for side in range(2):
   if (m&(1<<side))and(n&(1<<side)):
    nxt=[]
    for b,steps in rows:
     for h in nu[side]:
      z=BT.get((b,h))
      if z is not None:nxt.append((z,steps+((side,h),)))
    rows=nxt
  return [(4*b+masks,steps)for b,steps in rows]
 # SOURCE-only positive propagation, discarding no cancelling path.
 def supports(A,B):
  nonlocal clock
  C=set()
  for u in A:
   for v in B:
    paths=routes(u,v);clock+=1
    for k,steps in paths:clock+=len(steps)+2;C.add(k)
  return C
 Sp={};one={0}
 for i in range(9):
  Sp[i,0]=one;Sp[i,1]=set(kap[i])
  for n in range(2,maxpow[i]+1):Sp[i,n]=supports(Sp[i,n-1],Sp[i,1])
 for poly in program:
  for ex,q in poly:
   S=one
   for i,n in enumerate(ex):
    if n:S=supports(S,Sp[i,n])
   clock+=2*len(S)
 predicted=clock;assert predicted<2000000
 # Every sparse entry is charged as key+numerator+denominator. Include all
 # preflight support sets until released plus full expression/input metadata.
 ownprior=count(d)+4000+3*(9*68+2*17+sum(maxpow)*68+8*68+4*68)+sum((n+1)*68 for n in maxpow)
 assert ownprior<25258
 gate=dict(all_passed=True,SOURCE_no_new_G8_values=True,code_sha256=sha(Path(__file__)),input_sha256=pins,positive_operation_count_upper=predicted,hard_products2000000=2000000,whole_numeric_upper=25258,source_derived_owned_upper=ownprior,committed_prior313_hard384=True,temporary_prior940_hard1153=True,RAW_upper4000000=4000000,soft45_hard60_RAM2GB=True,only_existing_canonical_cache_and_formal_unit_program=True,no_new_physical_or_gauge_solve=True)
 if a.prepare:put('report.json',gate);print(sha(o/'report.json'));return
 assert a.gate_sha;gp=ROOT/'reports/AMBIENT_B18_NATIVE_CANONICAL68_VALUE_GATE_20261105_v63/report.json';assert sha(gp)==a.gate_sha and json.loads(gp.read_text())==gate
 del Sp;put('phase-before-values.json',dict(ONE_GO=True,gate_sha256=a.gate_sha,gate=gate))
 def own(extra=0):
  nonlocal peak
  n=count(d)+4000+3*(sum(map(len,kap))+sum(map(len,nu))+sum(map(len,powers.values()))+sum(map(len,results))+extra)
  peak=max(peak,n);assert n<=25258
 def check(q):assert abs(q.numerator).bit_length()<=313 and q.denominator.bit_length()<=313;return q
 def times(A,B):
  nonlocal products
  C={};state['partial_product']=C
  for u,x in A.items():
   for v,y in B.items():
    state['pair']=[u,v];paths=routes(u,v);products+=1;assert products<=predicted
    if not paths:continue
    t=check(x*y)
    for k,steps in paths:
     w=t
     for side,h in steps:w=check(w*nu[side][h]);products+=1
     z=check(C.get(k,Q(0))+w);products+=2;own(len(A)+len(B)+len(C)+4)
     if z:C[k]=z
     else:C.pop(k,None)
  assert time.monotonic()-start<45;return C
 try:
  for i in range(9):
   powers[i,0]={0:Q(1)};powers[i,1]=kap[i]
   for n in range(2,maxpow[i]+1):state.update(phase='powers',i=i,n=n);powers[i,n]=times(powers[i,n-1],kap[i])
  for pi,poly in enumerate(program):
   R={};state['current_G8']=R
   for ti,(ex,q)in enumerate(poly):
    state.update(phase='G8',pi=pi,ti=ti);V={0:Q(1)}
    for i,n in enumerate(ex):
     if n:state['factor_i']=i;V=times(V,powers[i,n])
    for k,v in V.items():
     z=check(R.get(k,Q(0))+check(v*Q(q)));products+=2;assert products<=predicted;own(len(R)+len(V)+4)
     if z:R[k]=z
     else:R.pop(k,None)
   results.append(R);own()
  # Both old quotients are literal monomial projections of this B17 ring.
  b11=d['old_B11_projection_raw33_indices'];loc={j:i for i,j in enumerate(sel)}
  for R,expected in zip(results,d['old_B11_expected_G8_C44']):assert [str(R.get(4*loc[j]+m,Q(0)))for j in b11 for m in range(4)]==expected
  cpath=ROOT/'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json';old=fields(cpath,['actualG8'])['actualG8'];pins[str(cpath.relative_to(ROOT))]=sha(cpath)
  for R,expected in zip(results,old):assert [str(R.get(4*loc[j]+m,Q(0)))for j in d['old_C16_projection_raw33_indices']for m in range(4)]==expected
  full=[[str(R.get(i,Q(0)))for i in range(68)]for R in results];second=[[row[4*j+m]for j in range(7,17)for m in range(4)]for row in full];nusecond=[[str(v.get(j,Q(0)))for j in range(7,17)]for v in nu]
  data=dict(G8_order=d['source_G8_order'],B17_basis=d['basis_B17'],C68_index='4*B17+mask',full_G8_C68=full,G8_second10_free4=second,monic_second10=nusecond,all352_B11_and128_C16_overlap_coefficients_verified=True,no_honest_quotient_or_upper_HF_values=True)
  assert count(d)+count(data)+4000+3*(sum(map(len,powers.values()))+sum(map(len,results))+sum(map(len,kap)))<25258;put('actual-G8-C68-and-monic-second-values.json',data)
  put('report.json',dict(all_passed=True,all_G8_scalar_second320=True,all_monic_second20=True,full_G8_coefficients544=True,all480_prior_overlaps_exact=True,positive_operation_upper=predicted,actual_operations=products,peak_owned=peak,source_prior313=True,source_future_unrelated_honest_frontier_still_opaque=True,input_sha256=pins,code_sha256=sha(Path(__file__)),evidence_sha256={p.name:sha(p)for p in o.iterdir()},seconds=time.monotonic()-start));print(sha(o/'report.json'))
 except Exception as e:
  put('failure-prefix.json',dict(error=repr(e),state=state,powers={str(k):v for k,v in powers.items()},results=results,input_sha256=pins,operation_count=products,no_auto_retry=True));raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--gate-sha');main(p.parse_args())
