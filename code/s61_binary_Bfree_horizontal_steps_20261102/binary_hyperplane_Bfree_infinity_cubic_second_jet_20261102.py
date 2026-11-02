from pathlib import Path
import json,time,functools
import sympy as s
from sympy.polys.rings import ring
st=time.time();R,a,e=ring('a,e',s.QQ);light=[0,1,2,3,4,6,7];tw=[1,0,0,0,1,0,1,1];info=json.loads(Path('work/binary_hyperplane_Bfree_infinity_cubic_reduction_20261102.json').read_text());inits=json.loads(Path('work/binary_hyperplane_Bfree_infinity_initial_forms_20261102.json').read_text());tr=json.loads(Path('work/binary_hyperplane_infinity_traced_projection_20261102.json').read_text())['projections']['3'];ms=[tuple(m)for m in tr['monomials']];ix={m:i for i,m in enumerate(ms)};lifts=[{i:s.QQ(c)for i,c in row}for row in tr['relation_lifts']];contents={i:R.from_expr(s.sympify(c)) for i,c in zip(info['cubic_sources'],info['cubic_contents'])};leading_content={r['source']:R.from_expr(s.sympify(r['content']))for r in inits['rows']};N={i:[R.from_expr(s.sympify(c))for c in row]for i,row in zip(light,info['normal_num'])}
def addto(out,k,p):
 z=out.get(k,R.zero)+p
 if z:out[k]=z
 elif k in out:del out[k]
@functools.lru_cache(None)
def image(mon):
 out={(0,)*8:R.one}
 for i in range(8):
  for _ in range(mon[i]):
   linear={tuple(int(k==5)for k in range(8)):R.one} if i==5 else {tuple(int(k==j)for k in range(8)):N[i][j]for j in light if N[i][j]};new={}
   for m,c in out.items():
    for n,d in linear.items():addto(new,tuple(x+y for x,y in zip(m,n)),c*d)
   out=new
 return out
def sparse(p):return[[list(m),str(c)]for m,c in p.items()]
rows=[]
with Path('/tmp/poisson-hyperplane-Bfree-cubic-20261102/rows').open()as f:
 for _ in range(int(f.readline())):
  idx,w,ct=map(int,f.readline().split());rr=[]
  for _ in range(ct):
   m=tuple(map(int,f.readline().split()));z=f.readline().split()
   if w!=3:continue
   for k in range(int(z[0])):rr.append((m,tuple(map(int,z[2+4*k:5+4*k])),s.QQ(z[1+4*k])))
  if w!=3:continue
  degree=max(ex[2]-sum(i*j for i,j in zip(m,tw))for m,ex,c in rr);parts=[{},{},{}]
  for m,ex,c in rr:
   delta=degree-ex[2]+sum(i*j for i,j in zip(m,tw))
   if delta<3:parts[delta].setdefault(m,{})[ex[:2]]=c
  transformed=[]
  for pp in parts:
   out={}
   for mon,coeff in pp.items():
    coeff=R.from_dict(coeff)
    for m,c in image(mon).items():addto(out,ix[m],coeff*c)
   transformed.append(out)
  g=leading_content[idx];F0={i:c.exquo(g)for i,c in transformed[0].items()};content=contents[idx] or R.one;alpha={}
  for row,piv in zip(lifts,tr['pivot_columns']):
   c=F0.get(piv,R.zero)
   for k,scale in row.items():addto(alpha,k,c*scale)
  reduced=dict(F0)
  # Reconstruct the polynomial leading remainder from the traced original relations.
  for row,piv in zip(tr['relations'],tr['pivot_columns']):
   c=F0.get(piv,R.zero)
   for k,scale in row:addto(reduced,k,-c*s.QQ(scale))
  expected=info['cubic_rows_sparse'][info['cubic_sources'].index(idx)];evec={tr['free_columns'][i]:R.from_dict({tuple(m):s.QQ(c)for m,c in p})for i,p in expected};assert {i:c.exquo(content)for i,c in reduced.items()}==evec
  rows.append(dict(source=idx,delta_shift=degree,raw_leading_content=sparse(g),primitive_cubic_content=sparse(content),first_numerator=[[[*ms[i]],sparse(c)]for i,c in transformed[1].items()],second_numerator=[[[*ms[i]],sparse(c)]for i,c in transformed[2].items()],correction_numerator=[[i,sparse(c)]for i,c in alpha.items()]));print('SECOND',idx,'terms',sum(len(c)for c in transformed[1].values()),'alpha',sum(len(c)for c in alpha.values()),'sec',time.time()-st,flush=True)
rep=dict(all_passed=True,scope='Exact through-second-order cubic data in constant quadratic normal coordinates; actual first-order correction requires the separately reconstructed Q1.',normal_common_denominator=info['D'],formula='C1_j = first_numerator/(raw_leading_content*primitive_cubic_content) - sum(correction_numerator[k]*Q1_i*y_shift)/primitive_cubic_content',quadratic_product_labels=tr['row_labels'],rows=rows,seconds=time.time()-st)
p=Path(__file__).with_suffix('.json');assert not p.exists();p.write_text(json.dumps(rep,indent=2)+'\n');print('PASS',time.time()-st)
