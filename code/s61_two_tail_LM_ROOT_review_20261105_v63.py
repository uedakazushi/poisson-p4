"""Independent streamed LM, complete localization, and UFD nondivision review."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import fields

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 o=ROOT/'reports/S61_NONSQUARE_PREPOWER_TWO_TAIL_LM_ROOT_REVIEW_20261105_v63';o.mkdir(exist_ok=False)
 src=ROOT/'reports/S61_NONSQUARE_PREPOWER_TWO_TAIL_LM_SOURCE_20261105_v63/report.json';r=json.loads(src.read_text());assert r['all_passed']
 for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h
 for n,h in r['evidence_sha256'].items():assert sha(src.parent/n)==h
 (o/'pre-arithmetic.json').write_text(json.dumps(dict(SOURCE_ONLY=True,unit_product_prior_bits=128,owned_upper=12000,RAW_upper4000000=4000000),indent=2)+'\n')
 unit=ROOT/'reports/S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63/unit-inventory.json';u=json.loads(unit.read_text());assert u['ring']=='QQ[a,e,d,1/(a*d*F*B*C*J)]' and u['factor_order']==['a','d','F','B','C','J']
 fs=[{tuple(m):Q(x)for m,x in f}for f in u['named_factors']]
 normalized=[sorted((m,x/f[max(f)])for m,x in f.items())for f in fs];assert len({str(v)for v in normalized})==6
 # Literal prime-form coefficients: variables a,d; F,C have constant leading a;
 # J has constant leading e. B primitive in d, since e=2/3 gives b1=50a.
 a,d,F,B,C,J=fs;assert a=={(1,0,0):Q(1)} and d=={(0,0,1):Q(1)}
 assert max(m[0]for m in F)==1 and {(m[1],m[2]):x for m,x in F.items()if m[0]}=={(0,0):Q(500)}
 assert max(m[0]for m in C)==1 and {(m[1],m[2]):x for m,x in C.items()if m[0]}=={(0,0):Q(100)}
 assert max(m[1]for m in J)==1 and {(m[0],m[2]):x for m,x in J.items()if m[1]}=={(0,0):Q(-6)}
 assert B=={(1,0,1):Q(50),(0,2,1):Q(18),(0,1,1):Q(-9),(0,0,1):Q(-2),(0,1,0):Q(60),(0,0,0):Q(-40)}
 assert Q(18)*Q(2,3)**2-Q(9)*Q(2,3)-2==0
 product={(0,0,0):Q(1)}
 for f in reversed(fs):
  tmp={}
  for m,x in product.items():
   for n,y in f.items():
    e=tuple(i+j for i,j in zip(m,n));tmp[e]=tmp.get(e,Q(0))+x*y;assert max(abs(tmp[e].numerator).bit_length(),tmp[e].denominator.bit_length())<=128
  product={m:x for m,x in tmp.items()if x}
 p=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63/certificate.json';v=fields(p,['FULLU_product']);assert product=={tuple(m):Q(x)for m,x in v['FULLU_product']} and len(product)==49
 stream=[]
 for stem,magic,n in [('S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_ENTRY_20261105_v63','T0C4QV001',4500),('S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_ENTRY_20261105_v63','Q1B4QV001',2317)]:
  p=ROOT/'reports'/stem/'quotient-integer-common-D.txt';mn=[100,100,100];lm=(-1,-1,-1);count=0
  with p.open()as f:
   assert f.readline().strip()==magic;D=int(f.readline());assert D.bit_length()==197 and int(f.readline())==n
   for line in f:
    k,c=map(int,line.split());assert c;aa,x=divmod(k,1048576);ee,dd=divmod(x,1024);m=(aa,ee,dd);lm=max(lm,m);mn=[min(i,j)for i,j in zip(mn,m)];count+=1
  assert count==n;stream.append((lm,mn,D))
 assert stream[0][0]==(16,10,18) and stream[1][0]==(12,11,14) and stream[1][1]==[0,0,0] and stream[0][2]==stream[1][2]
 cp=json.loads((src.parent/'N2-coprime-to-ALL-inverted-factors.json').read_text())
 for kind,h in cp['saved_valuation_dependencies'].items():
  assert h in r['input_sha256'].values(),kind
 # Exact valuation dependencies are separate accepted primary/peer receipts;
 # this proof uses additivity in a UFD and the already checked B4 identity.
 rr=dict(all_passed=True,source_only=True,streamed4500_2317_LMs_independently_verified=True,complete_six_factor_product49_rebuilt_in_reverse_order=True,all_six_prime_and_nonassociate_forms_verified=True,
 mathematical_reason='N2 has no a,d factors by literal minima; F,J by v(Q1)=0 and Q1=B^4 N2/D; C by vC(T1)=4 and T1=a^3 C^4 B^4 N2/D; B by exact vB(Q1)=4 and full B4 identity. Thus N2 is coprime to all six inverted primes. If N2 divides s*N0, UFD cancellation forces N2|N0; lex leading monomials contradict this. Nondivision also overC by same primitive/linear prime forms.',
 exact_scope='(Q2) properly contained in (Q0,Q2) in QQ[a,e,d,1/(a*d*F*B*C*J)] ONLY; no further D12/T/unknown localization, strictopen/radical/Fitting/emptycomplement claim',input_sha256={str(src.relative_to(ROOT)):sha(src),str(unit.relative_to(ROOT)):sha(unit)},code_sha256=sha(Path(__file__)))
 (o/'report.json').write_text(json.dumps(rr,indent=2)+'\n');print(sha(o/'report.json'))
if __name__=='__main__':main()
