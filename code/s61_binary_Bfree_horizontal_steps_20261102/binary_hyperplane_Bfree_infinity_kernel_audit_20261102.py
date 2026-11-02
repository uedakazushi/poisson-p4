from pathlib import Path
import json,time,hashlib
from fractions import Fraction
p=32003;start=time.time();inp=Path('/tmp/hyperplane-Bfree-infinity-reduced-kernel-20261102.in');lines=inp.read_text().splitlines();nr,nc,n,k=map(int,lines[0].split());free=list(map(int,lines[3].split()));assert(nr,nc,n,k)==(147,106,31,75)
def readpoly(line):
 z=line.split();assert len(z)==1+3*int(z[0]);return{(int(z[2+3*i]),int(z[3+3*i])):Fraction(z[1+3*i])for i in range(int(z[0]))}
def evalp(q):
 v=0
 for(i,j),c in q.items():assert c.denominator%p;v=(v+c.numerator*pow(c.denominator,-1,p)*pow(2,j,p))%p
 return v
out=Path('/tmp/hyperplane-Bfree-infinity-reduced-kernels-20261102.out');counts=[];denvalues=[];maxdegrees=[]
with out.open()as f:
 for j in range(75):
  assert f.readline().strip()==f'KERNEL {j}';den=readpoly(f.readline());dv=evalp(den);assert dv;denvalues.append(dv);vectors=[readpoly(f.readline())for _ in range(106)];assert vectors[free[j]]==den
  for h in range(75):
   if h!=j:assert not vectors[free[h]]
  for q in vectors:evalp(q)
  counts.append(sum(map(len,vectors)));maxdegrees.append(max((sum(m)for q in vectors for m in q),default=0))
 assert not f.read().strip()
r=dict(all_passed=True,scope='Audit of the exact reduced leading-kernel output. Full 147x75 polynomial identities were checked by the native producer; this audit independently checks all free anchors, p-integrality, and nonzero denominator evaluations.',rank_lower=31,kernel_columns=75,free_anchor_diagonal=True,denominator_point_values=denvalues,max_numerator_degrees=maxdegrees,total_numerator_terms=sum(counts),maximum_column_terms=max(counts),input_hashes={str(q):hashlib.sha256(q.read_bytes()).hexdigest()for q in [inp,out]},seconds=time.time()-start)
Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:r[k]for k in ['all_passed','kernel_columns','total_numerator_terms','maximum_column_terms','seconds']}))
