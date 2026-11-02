from pathlib import Path
import json,time,functools,sys
import sympy as s,numpy as np
from sympy.polys.rings import ring
sys.path.insert(0,'code');from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank,determinant_mod
start=time.time();R,a,e=ring('a,e',s.QQ);parse=lambda pp:R.from_dict({tuple(m):s.QQ(c)for m,c in pp});sparse=lambda pp:[[list(m),str(c)]for m,c in pp.items()]
cr=json.loads(Path('work/binary_hyperplane_Bfree_infinity_cubic_product_delta1_20261102.json').read_text());qr=json.loads(Path('work/generic_Bfree_infinity_quadratic_syzygy_delta1_20261102.json').read_text());M=[{j:parse(p)for j,p in row['numerators']}for row in cr['rows']];dens=[parse(row['denominator'])for row in cr['rows']]
def readpoly(line):
 z=line.split();assert len(z)==1+3*int(z[0]);return R.from_dict({(int(z[2+3*i]),int(z[3+3*i])):s.QQ(z[1+3*i])for i in range(int(z[0]))})
def addto(out,k,p):
 z=out.get(k,R.zero)+p
 if z:out[k]=z
 elif k in out:del out[k]
rows=[{j:parse(p)for j,p in rr}for rr in qr['obstruction_numerators']];rowden=[R.from_expr(s.sympify(qr['common_denominator']))]*42
with Path('/tmp/hyperplane-Bfree-infinity-reduced-left-kernels-20261102.out').open()as f:
 for i in range(116):
  assert f.readline().strip()==f'KERNEL {i}';anchor=readpoly(f.readline());L=[readpoly(f.readline())for _ in range(147)];used=[j for j,z in enumerate(L)if z];den=R.one
  for j in used:den=den.lcm(dens[j])
  out={}
  for j in used:
   c=L[j]*den.exquo(dens[j])
   for k,v in M[j].items():addto(out,k,c*v)
  g=den
  for z in out.values():g=g.gcd(z)
  out={k:v.exquo(g)for k,v in out.items()};den=den.exquo(g);rows.append(out);rowden.append(den)
  if i%10==0:print('LEFT',i,'used',len(used),'nnz',len(out),'terms',sum(len(v)for v in out.values()),'seconds',time.time()-start,flush=True)
 assert not f.read().strip()
prime=32003
@functools.lru_cache(None)
def ev(v):return sum(int(c.numerator)*pow(int(c.denominator),-1,prime)*pow(2,m[1],prime)for m,c in v.items())%prime
num=np.array([[ev(row.get(j,R.zero))for j in range(106)]for row in rows],dtype=np.int64);kern=[]
with Path('/tmp/hyperplane-Bfree-infinity-reduced-kernels-20261102.out').open()as f:
 for i in range(75):
  assert f.readline().strip()==f'KERNEL {i}';anchor=readpoly(f.readline());assert ev(anchor);kern.append([ev(readpoly(f.readline()))for _ in range(106)])
 assert not f.read().strip()
prod=num@np.array(kern,dtype=np.int64).T%prime;rank,pc=elimination_rank(prod,prime);_,pr=elimination_rank(prod[:,pc].T,prime);minor=determinant_mod(prod[np.ix_(pr,pc)],prime)
for d in rowden:assert ev(d)
print('NEXT',rank,'MINOR',minor,'TERMS',sum(len(v)for row in rows for v in row.values()),'SECONDS',time.time()-start,flush=True)
rep=dict(all_passed=True,scope='Exact 158x106 first-order left obstruction before multiplying the certified 106x75 right kernel; numeric product rank is a lower bound only.',shape=[158,106],rows=[[[j,sparse(v)]for j,v in row.items()]for row in rows],denominators=list(map(sparse,rowden)),numeric_point=[1,2],prime=prime,numeric_product_rank=rank,numeric_product_pivot_rows=pr,numeric_product_pivot_cols=pc,numeric_product_minor=minor,seconds=time.time()-start)
p=Path(__file__).with_suffix('.json');assert not p.exists();p.write_text(json.dumps(rep,indent=2)+'\n');print('PASS')
