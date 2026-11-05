"""Independent squarefreeness via derivative gcd at the SAME fixed H point."""
from pathlib import Path
import hashlib,json,math,time
ROOT=Path(__file__).resolve().parents[1];P=32003
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();src=ROOT/'reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/factor-sparse.txt';assert sha(src)=='241488fafa9b45e5353f0b707f75902bc5bc32e85ee093b9f2e0a5f07c7d571c'
 with src.open()as f:assert int(f.readline())==94;H=[tuple(map(int,l.split()))for l in f]
 assert len(H)==94 and max(abs(c).bit_length()for a,e,d,c in H)<=40
 h6={(e,d):c for a,e,d,c in H if a==6};assert h6=={(0,3):1406250000,(0,4):-439453125,(1,4):1160156250}
 h00={d:c for a,e,d,c in H if a==0 and e==0};assert h00=={0:-1025280,1:-12816}
 assert -1025280*5-12816*16==-5331456
 out=ROOT/'reports/S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63';out.mkdir(exist_ok=False)
 pin={str(src.relative_to(ROOT)):sha(src)}
 (out/'pre-arithmetic-SOURCE.json').write_text(json.dumps(dict(input_sha256=pin,same_fixed_e=0,same_fixed_d=1,same_fixed_p=P,no_new_parameter_prime_or_source_query=True,field_max_degree=6,field_product_bound=10000,field_temp_bits=30,input_conversion_bits=41,source_numeric_upper=4096,RAW_upper=262144,all_inputs_existing=True),indent=2)+'\n')
 f=[0]*7
 for a,e,d,c in H:
  if e==0:f[a]=(f[a]+c%P)%P
 assert f[-1]==18248
 f=[(c*pow(f[-1],P-2,P))%P for c in f]
 derivative=[(i*f[i])%P for i in range(1,7)]
 def trim(x):
  while x and not x[-1]:x.pop()
  return x
 a=f[:];b=trim(derivative);chain=[];products=0
 while b:
  assert time.monotonic()-t<45
  r=a[:];q=[0]*max(0,len(a)-len(b)+1);iv=pow(b[-1],P-2,P);products+=30
  while r and len(r)>=len(b):
   k=len(r)-len(b);v=r[-1]*iv%P;q[k]=v;products+=1
   for j,c in enumerate(b):r[k+j]=(r[k+j]-v*c)%P;products+=1
   trim(r)
  chain.append(dict(a=a,b=b,quotient=trim(q),remainder=r));a,b=b,r
 assert len(a)==1 and a[0] and products<10000
 rep=dict(id=out.name,all_passed=True,independent_method='derivative Euclid gcd, no Rabin/Frobenius producer import',same_fixed_e=0,same_fixed_d=1,same_fixed_p=P,specialized_monic=f,derivative=derivative,complete_Euclid_chain=chain,gcd_is_one=True,
 H_primitive_as_polynomial_in_a=True,content_proof='leadinga6=17578125*d^3*L, L=80-25d+66ed irreducible (primitive linearine); a0 at e0=-1025280-12816d neither divisiblebyd norL since valueat d16/5 is-5331456/5',
 exact_H_squarefree_over_QQ_a_e_d=True,squarefree_proof='a-primitiveH has no a-freefactor. Repeatedchar0positive-a factor wouldspecialize withsamepositive-degree since totaldegreea6 andLCsurvive32003, givingrepeatedfactorof f. But independentgcd(f,fprime)=1. Thus no repeatedHfactor.',
 irreducibility_proved=False,no_Fitting_full_rank_or_exhaustion_claim=True,source_numeric_upper=4096,field_products_including_inverse_prior=products,input_sha256=pin,code_sha256=sha(Path(__file__)),evidence_sha256={'pre-arithmetic-SOURCE.json':sha(out/'pre-arithmetic-SOURCE.json')},seconds=time.monotonic()-t)
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');print(sha(p))
if __name__=='__main__':main()
