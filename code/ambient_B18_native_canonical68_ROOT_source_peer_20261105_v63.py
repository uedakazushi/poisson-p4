"""Independent B17 ring and formal G8 recipe replay, no new response values."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,itertools,json,sys,signal,resource
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import count,fields

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));o=ROOT/'reports/AMBIENT_B18_NATIVE_CANONICAL68_ROOT_SOURCE_PEER_20261105_v63';o.mkdir(exist_ok=False)
 p=ROOT/'reports/AMBIENT_B18_ELL2_NATIVE_CANONICAL33_INPUT_GATE_V2_20261105_v63/report.json';r=json.loads(p.read_text());assert r['all_passed']
 pins={str(p.relative_to(ROOT)):sha(p)}
 for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h;pins[n]=h
 ep=p.parent/'canonical33-to-B17-C68-G8-source-input-and-response-budget.json';assert sha(ep)==r['evidence_sha256'][ep.name];d=json.loads(ep.read_text());pins[str(ep.relative_to(ROOT))]=sha(ep)
 proof=dict(phase='before all formal recipe arithmetic',formal_integer_bound_bits=256,source_owned_upper=60000,no_new_functionresponse=True);(o/'pre-arithmetic.json').write_text(json.dumps(proof,indent=2)+'\n')
 B=list(map(tuple,d['basis_B17']));raw=list(map(tuple,d['raw33_basis']));ids=d['B17_raw33_selector_indices'];assert [raw[i]for i in ids]==B
 T=json.loads((ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json').read_text())
 for i in range(33):
  for j in range(33):
   v={int(k):Q(x)for k,x in T['table_columns_for_each_left_basis'][i][j]};lhs=tuple(v.get(k,Q(0))for k in ids)
   mon=tuple(x+y for x,y in zip(raw[i],raw[j]));rhs=tuple(Q(int(raw[i]in B and raw[j]in B and mon==m))for m in B);assert lhs==rhs
 assert all(Q(v[0])==0 for v in d['actual_nu_A_B33_cached_UNEVALUATED']) and all(Q(v[0])==0 for v in d['actual_kappa9_C132_cached_UNEVALUATED'])
 # Independent formal inverse: compute (1/5-d²/25+d⁴/125-d⁶/625)^p
 # by polynomial multiplication, not producer binomial-coefficient expansion.
 z=(0,)*9
 def mono(i,n=1):m=[0]*9;m[i]=n;return tuple(m)
 def times(a,b):
  c={}
  for m,x in a.items():
   for n,y in b.items():
    e=tuple(v+w for v,w in zip(m,n))
    if sum(e)<=6:
     c[e]=c.get(e,Q(0))+x*y;assert max(abs(c[e].numerator).bit_length(),c[e].denominator.bit_length())<=256
  return {m:x for m,x in c.items()if x}
 ip={mono(1,2*k):Q((-1)**k,5**(k+1)) for k in range(4)};powers=[{z:Q(1)}]
 for _ in range(3):powers.append(times(powers[-1],ip))
 K=json.loads((ROOT/'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json').read_text());prep=json.loads((ROOT/'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json').read_text());program=[];core=[0,1,2,6]
 for name in K['selected5_pivot_order']:
  src=K['five_actual_rational_roots'][name];assert src['D_power']==0
  num={}
  for ex,n,den in src['numerator']:
   e=[0]*9
   for i,x in zip(core,ex):e[i]=x
   num[tuple(e)]=num.get(tuple(e),Q(0))+Q(n,den)
  p0={mono(K['source9_variable_order'].index(name)):Q(1)}
  for m,x in times(num,powers[src['P_power']]).items():p0[m]=p0.get(m,Q(0))-x
  program.append({m:x for m,x in p0.items()if x})
 for source in prep['generators3']:
  poly={}
  for ex,n,den in source:
   e=[0]*9
   for i,x in zip(core,ex):e[i]=x
   if sum(e)<=6:poly[tuple(e)]=poly.get(tuple(e),Q(0))+Q(n,den)
  program.append({m:x for m,x in poly.items()if x})
 assert program==[{tuple(m):Q(v)for m,v in row}for row in d['G8_formal_kappa9_polynomials_mod_radical7']]
 weights=[2*sum(m)+mask.bit_count()for m in B for mask in range(4)];assert max(weights)==6
 assert r['coarse_future_committed_bits']==313 and r['coarse_future_temporary_bits']==940 and r['future_numeric_upper']==25258
 rr=dict(all_passed=True,source_only=True,all1089_ring_products_independently_projected=True,eight_G8_formal_recipes_independently_expanded_by_geometric_power=True,
 exact_radical7_and_fourterm_inverse_scope='B17 m3=0; basevariables weight2 and trace weight1, monic carry nu in base maximal ideal has weight>=2. Maximum basis weight6. Scalar residue of each kappa zero => kappa^7=0; d^8=0. G8 degree<=6 truncation and P inversefourterms exact. At most2 nonzero nu carries in any surviving product.',
 no_new_native_response_evaluated=True,positive_future313_940_25258_bound_accepted=True,input_sha256=pins,code_sha256=sha(Path(__file__)))
 (o/'report.json').write_text(json.dumps(rr,indent=2)+'\n');print(sha(o/'report.json'))
if __name__=='__main__':main()
