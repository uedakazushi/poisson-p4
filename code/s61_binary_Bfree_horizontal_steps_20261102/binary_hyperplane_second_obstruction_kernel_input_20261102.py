from pathlib import Path
import json,time,sys,functools
import sympy as s,numpy as np
sys.path.insert(0,'code');from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank,determinant_mod
st=time.time();p=32003;src=Path('/tmp/hyperplane-Bfree-infinity-second-dual-second-product-20261102.out');lines=src.read_text().splitlines();nr,nc=map(int,lines[0].split());assert(nr,nc)==(117,34);matrix=[];cost=[];scales=[];values=[]
@functools.lru_cache(None)
def scalar(v):
 q=s.QQ(v);assert int(q.denominator)%p;return int(q.numerator)*pow(int(q.denominator),-1,p)%p
for i in range(nr):
 scales.append(lines[1+i*(nc+1)]);row=lines[2+i*(nc+1):2+i*(nc+1)+nc];matrix.append(row);zz=[];vv=[]
 for ln in row:
  t=ln.split();n=int(t[0]);zz.append(n);vv.append(sum(scalar(t[1+3*j])*pow(2,int(t[3+3*j]),p)for j in range(n))%p)
 cost.append(zz);values.append(vv)
a=np.array(values,dtype=np.int64);cost=np.array(cost,dtype=np.int64);rowcost=cost.sum(axis=1);colcost=cost.sum(axis=0);rr=sorted(range(nr),key=lambda i:(rowcost[i],i));cc=sorted(range(nc),key=lambda j:(colcost[j],j));rank,pc0=elimination_rank(a[np.ix_(rr,cc)],p);pc=[cc[j]for j in pc0];_,pr0=elimination_rank(a[np.ix_(rr,pc)].T,p);pr=[rr[i]for i in pr0];assert rank==18;fr=[j for j in range(nc)if j not in pc];minor=determinant_mod(a[np.ix_(pr,pc)],p)
out=Path('/tmp/hyperplane-Bfree-infinity-third-kernel-20261102.in')
with out.open('w')as f:
 f.write(f'{nr} {nc} {rank} {len(fr)}\n'+' '.join(map(str,pr))+'\n'+' '.join(map(str,pc))+'\n'+' '.join(map(str,fr))+'\n')
 for row in matrix:f.write('\n'.join(row)+'\n')
rep=dict(all_passed=True,point=[1,2],prime=p,rank=rank,minor=minor,rows=pr,cols=pc,free=fr,zero_rows=[i for i,v in enumerate(rowcost)if not v],zero_cols=[j for j,v in enumerate(colcost)if not v],total_terms=int(cost.sum()),pivot_matrix_terms=int(cost[np.ix_(pr,pc)].sum()),augmented_terms=int(cost[pr,:].sum()),input=str(out),seconds=time.time()-st)
Path(__file__).with_suffix('.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
