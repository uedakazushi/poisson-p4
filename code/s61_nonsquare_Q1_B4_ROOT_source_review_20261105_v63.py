"""Independent fixed B4 quotient source/code/budget audit; no actual quotient."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B4_ROOT_SOURCE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_WORKER_CUSTODY_20261105_v63';r=json.loads((base/'report.json').read_text());assert r['all_passed'] and r['value_worker_NOT_run'];pins={str((base/'report.json').relative_to(ROOT)):sha(base/'report.json')}
 for paths,where in [(r['input_sha256'],ROOT),(r['evidence_sha256'],base)]:
  for n,h in paths.items():
   p=where/n;assert sha(p)==h,(n,h);pins[str(p.relative_to(ROOT))]=h
 m=json.loads((base/'fixedB-fourpass-lifetime-height-cursors-and-guards.json').read_text())
 native=(base/'native-input.txt').read_text().splitlines();old=(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt').read_text().splitlines()
 assert native[0]=='Q1B4Q001 4811 4 0' and old[0]=='T1C4QN001' and old[2]=='4811'
 assert native[1]==old[1] and native[2:]==old[3:] and len(native)==4813
 assert int(native[1]).bit_length()==197
 # Source bounds only; no arithmetic on an actual quotient or B-Hasse value.
 Ns=[]
 for j in range(5):
  A,E,D,T=16-j,33-2*j,18-j,50-3*j
  Ns.append(sum(a+e+d<=T for a in range(A+1) for e in range(E+1) for d in range(D+1)))
 assert Ns==[10013,8400,6970,5712,4615]==m['stage_support_bounds']
 L=(1<<301)-1
 for j,st in enumerate(m['exact_L1_bounds']):
  factor=sum((F(129,50)**h)/50 for h in range(16-j));q=F(L)*factor;Q=(q.numerator+q.denominator-1)//q.denominator
  assert str(L)==st['input_L1_upper'] and factor==F(int(st['factor_num']),int(st['factor_den']))
  assert str(Q)==st['quotient_L1_upper'] and (L+179*Q).bit_length()==st['all_partial_remainder_bits']
  assert Q.bit_length()==st['quotient_bits'] and (60*Q).bit_length()==st['scalarproduct_bits'];L=Q
 verify=((1<<301)-1)+179**4*L;assert verify.bit_length()==386
 assert max(2*(Ns[j]+Ns[j+1])+16000 for j in range(4))==52826
 assert 2*(Ns[0]+Ns[-1]+225)+16000==45706
 products=8*sum(Ns[1:])+225*Ns[-1]+5400;assert products==1249351<1300000
 # Independent small polynomial fixture uses tuple monomials, no packed code.
 B={(1,0,1):50,(0,2,1):18,(0,1,1):-9,(0,0,1):-2,(0,1,0):60,(0,0,0):-40}
 target={(2,1,2):7,(0,3,0):-11,(0,0,0):5}
 p={}
 for x,a in B.items():
  for y,b in target.items():
   z=tuple(i+j for i,j in zip(x,y));p[z]=p.get(z,0)+a*b
 p={x:a for x,a in p.items() if a};q={}
 while p:
  lead=max(p);assert lead[0] and lead[2] and p[lead]%50==0
  key=(lead[0]-1,lead[1],lead[2]-1);c=p[lead]//50;assert key not in q;q[key]=c
  for x,a in B.items():
   z=tuple(i+j for i,j in zip(key,x));p[z]=p.get(z,0)-c*a
   if not p[z]:del p[z]
 assert q==target
 cpp=ROOT/'code/s61_nonsquare_Q1_B4_quotient_20261105_v63.cpp';driver=ROOT/'code/s61_nonsquare_prePOWER_Q1_B4_quotient_entry_20261105_v63.py'
 c=cpp.read_text();d=driver.read_text()
 required=['while(!R.empty())','fmpz_divisible_ui(lead->second.v,50)','for(int j=0;j<6;++j)','R.swap(Q)','F.clear();nf=0;F.swap(N)','std::ifstream again(argv[2])','for(auto&[qkey,qc]:Q)for(auto&[key,c]:F)','if(!R.empty())throw Stop("FULL_B4_multiplyback_NONZERO")','if(fmpz_bits(z.v)>518)','Q1B4QP001','if(v>OWN)']
 assert all(x in c for x in required)
 assert 'for n,h in pins.items():assert sha(ROOT/n)==h,n' in d
 assert c.index('F.clear();nf=0;F.swap(N)')<c.index('std::ifstream again(argv[2])')
 assert 'std::vector' not in c and 'fmpz_inv' not in c
 raw=json.loads((base/'RAW-closure.json').read_text());assert raw['total_upper']==7174350<8000000
 evidence=dict(all_passed=True,SOURCE_ONLY_no_actual_quotient=True,input_BODY4811_D197_identical=True,
 independent_degree_boxes=Ns,intermediate_integer_bits_max363=True,full_multiplyback_bits=386,hard518=True,
 whole_owned=52826,final_multiplyback_owned=45706,all_products_upper=products,RAW_upper=raw['total_upper'],
 independent_small_tuple_lex_fixture_PASS=True,
 source_reasoning='Exact Gauss divisibility from independently known vB=4. Lex LM50ad decreases strictly; all six subtraction terms retained, quotient coefficients integral. Partial R stays in true old-stage box from R=I-B Qremoved. Four empty R imply original I=B^4 Qfinal. Full multiplyback separately rereads immutable original SHA-bound input. Packed bases injective since exponents <1024. No original coefficient vector retained next to R/Q. F/N two small powers released before original reread.',
 limitation='Code counters are telemetry. Proven degree/term bounds imply product upper; no quotient value or new factor is evaluated here.',input_sha256=pins,code_sha256=sha(Path(__file__)))
 (out/'report.json').write_text(json.dumps(evidence,indent=2)+'\n');print(sha(out/'report.json'))
if __name__=='__main__':main()
