"""Independent closed multinomial B^4 multiplyback, no quotient division."""
from pathlib import Path
import hashlib,itertools,json,math,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B4_ROOT_MULTINOMIAL_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
 primary=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_ENTRY_20261105_v63/report.json';p=json.loads(primary.read_text());assert p['all_passed']
 inp=ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt';qp=primary.parent/'quotient-integer-common-D.txt'
 assert sha(qp)==p['evidence_sha256'][qp.name]
 L=(1<<301)-1
 from fractions import Fraction as Q
 for j in range(4):
  f=sum(Q(129,50)**h/50 for h in range(16-j));v=L*f;L=(v.numerator+v.denominator-1)//v.denominator
 bound=((1<<301)-1)+179**4*L;assert bound.bit_length()==386
 proof=dict(SOURCE_BEFORE_MULTINOMIAL_ARITHMETIC=True,algorithm='closed6-multinomial compositions of4, not repeated multiplication or lexdivision. Convolve supplied quotient with literalB4; stream-subtract original4811 then require empty result.',quotient_L1_upper=str(L),all_partial_accumulation_bound=str(bound),bits=386,hard518=True,
 coefficient_key_cells=2,owned_upper=45706,product_upper=225*4615+4811,RAW_proposed=8000000,soft45_hard60_RAM2GB=True,
 pins={str(primary.relative_to(ROOT)):sha(primary),str(inp.relative_to(ROOT)):sha(inp),str(qp.relative_to(ROOT)):sha(qp)})
 (out/'pre-arithmetic-source.json').write_text(json.dumps(proof,indent=2)+'\n')
 def own(q,f,res):assert 2*(len(q)+len(f)+len(res))+16000<=45706
 def checked(n):assert abs(n)<=bound and abs(n).bit_length()<=518;return n
 with qp.open()as f:
  assert f.readline().strip()=='Q1B4QV001';D=int(f.readline());n=int(f.readline());q=[];last=-1;Lq=0
  for _ in range(n):
   key,c=map(int,f.readline().split());assert key>last and c;last=key
   a,x=divmod(key,1048576);e,d=divmod(x,1024);assert 0<=a<=12 and 0<=e<=25 and 0<=d<=14 and a+e+d<=38
   q.append((key,c));Lq+=abs(c)
  assert not f.read().strip() and Lq<=L and n<=4615
 # all126 compositions (h0,...,h5) of4, closed multinomial coefficient4!/prod h!.
 terms=[((1,0,1),50),((0,2,1),18),((0,1,1),-9),((0,0,1),-2),((0,1,0),60),((0,0,0),-40)]
 F={};compositions=0
 for h in itertools.product(range(5),repeat=6):
  if sum(h)!=4:continue
  compositions+=1;m=tuple(sum(h[j]*terms[j][0][i] for j in range(6)) for i in range(3));v=math.factorial(4)//math.prod(math.factorial(x) for x in h)
  for j,x in enumerate(h):v*=terms[j][1]**x
  key=m[0]*1048576+m[1]*1024+m[2];F[key]=F.get(key,0)+v
 F={k:v for k,v in F.items() if v};assert compositions==126 and len(F)<=225 and sum(abs(x) for x in F.values())<=179**4
 res={};products=0
 for k,c in q:
  for z,v in F.items():
   key=k+z;w=checked(res.get(key,0)+checked(c*v));products+=1
   a,x=divmod(key,1048576);e,d=divmod(x,1024);assert a<=16 and e<=33 and d<=18 and a+e+d<=50
   if w:res[key]=w
   else:res.pop(key,None)
  own(q,F,res)
  assert time.monotonic()-start<45
 with inp.open()as f:
  assert f.readline().strip()=='T1C4QN001' and int(f.readline())==D and int(f.readline())==4811
  for _ in range(4811):
   key,c=map(int,f.readline().split());w=checked(res.get(key,0)-c)
   if w:res[key]=w
   else:res.pop(key,None)
  assert not f.read().strip()
 assert not res
 r=dict(all_passed=True,all_original4811_coefficients_full_multiplyback=True,independent_closed_multinomial126=True,no_division_or_native_factor_engine=True,quotient_terms=n,B4_terms=len(F),convolution_products=products,original_D197=D,source_prior_bits=386,owned_upper=45706,
 input_sha256=proof['pins'],code_sha256=sha(Path(__file__)),seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()
