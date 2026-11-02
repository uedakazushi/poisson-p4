from pathlib import Path
import json,time
import sympy as s
from sympy.polys.rings import ring
start=time.time();R,a,e=ring('a,e',s.QQ);light=[0,1,2,3,4,6,7]
parse=lambda p:R.from_dict({tuple(m):s.QQ(c)for m,c in p})
sparse=lambda p:[[list(m),str(c)]for m,c in p.items()]
info=json.loads(Path('work/binary_hyperplane_Bfree_infinity_cubic_reduction_20261102.json').read_text());cc=json.loads(Path('work/binary_hyperplane_Bfree_infinity_corrected_cubic_delta1_20261102.json').read_text());qq=json.loads(Path('work/generic_Bfree_infinity_quadratic_syzygy_delta1_20261102.json').read_text());tr=json.loads(Path('work/binary_hyperplane_infinity_traced_projection_20261102.json').read_text())['projections']['4'];bp=json.loads(Path('work/binary_hyperplane_infinity_constant_quadratic_projection_20261102.json').read_text());Q1=[{j:parse(p)for j,p in row}for row in qq['Q1_projected_numerators']];DQ=R.from_expr(s.sympify(qq['common_denominator']));ms=[tuple(m)for m in tr['monomials']];ix={m:i for i,m in enumerate(ms)};proj=[{j:s.QQ(c)for j,c in row}for row in bp['projection_sparse']];lifts=[{j:s.QQ(c)for j,c in row}for row in tr['relation_lifts']];reps3=[tuple(m)for m in info['cubic_free_monomials']]
def addto(out,k,p):
 z=out.get(k,R.zero)+p
 if z:out[k]=z
 elif k in out:del out[k]
def prod_project(poly,i):
 raw={};res={}
 for m,c in poly.items():
  n=tuple(v+int(k==i)for k,v in enumerate(m));j=ix[n];addto(raw,j,c)
  for col,v in proj[j].items():addto(res,col,c*v)
 return raw,res
rows=[];leading=[]
for j,rr in enumerate(cc['rows']):
 assert rr['source']==info['cubic_sources'][j];c0={reps3[k]:parse(p)for k,p in info['cubic_rows_sparse'][j]};c1={tuple(m):parse(p)for m,p in rr['numerators']};den=parse(rr['denominator']);common=den.lcm(DQ);mulc=common.exquo(den);mulq=common.exquo(DQ)
 for i in light:
  raw0,out0=prod_project(c0,i);_,num=prod_project(c1,i);num={k:v*mulc for k,v in num.items()};alpha={}
  for pivot,lift in zip(tr['pivot_columns'],lifts):
   scale=raw0.get(pivot,R.zero)
   if scale:
    for k,c in lift.items():addto(alpha,k,scale*c)
  for k,c in alpha.items():
   for col,v in Q1[k].items():addto(num,col,-mulq*c*v)
  h=R.zero
  for v in out0.values():h=h.gcd(v)
  if not h:h=R.one
  normalized={k:v.exquo(h)for k,v in out0.items()};dd=common*h;g=dd
  for v in num.values():g=g.gcd(v)
  dd=dd.exquo(g);num={k:v.exquo(g)for k,v in num.items()}
  leading.append(normalized);rows.append(dict(source=rr['source'],light=i,leading_content=sparse(h),denominator=sparse(dd),numerators=[[k,sparse(v)]for k,v in num.items()]))
 print('PRODUCT',rr['source'],'terms',sum(len(v)for row in rows[-7:]for _,p in row['numerators']for v in [p]),'sec',time.time()-start,flush=True)
# Compare every leading coefficient with the actual input of the certified left/right kernels.
with Path('/tmp/hyperplane-Bfree-infinity-reduced-kernel-20261102.in').open()as f:
 assert tuple(map(int,f.readline().split()))==(147,106,31,75)
 for _ in range(3):f.readline()
 for row in leading:
  for k in range(106):
   t=f.readline().split();n=int(t[0]);p=R.from_dict({(int(t[2+3*j]),int(t[3+3*j])):s.QQ(t[1+3*j])for j in range(n)});assert p==row.get(k,R.zero)
 assert not f.read().strip()
rep=dict(all_passed=True,scope='First-delta rows of the exact 147 leading cubic-product generators, after traced elimination of the canonical quadratic rows. All leading rows match the certified kernel input.',rows=rows,seconds=time.time()-start)
p=Path(__file__).with_suffix('.json');assert not p.exists();p.write_text(json.dumps(rep,indent=2)+'\n');print('PASS',time.time()-start,flush=True)
