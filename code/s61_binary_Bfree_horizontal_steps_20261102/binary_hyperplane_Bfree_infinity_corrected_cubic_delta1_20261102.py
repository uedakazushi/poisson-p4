from pathlib import Path
import json,time,sympy as s
from sympy.polys.rings import ring
start=time.time();R,a,e=ring('a,e',s.QQ);data=json.loads(Path('work/binary_hyperplane_Bfree_infinity_cubic_first_jet_20261102.json').read_text());qd=json.loads(Path('work/generic_Bfree_infinity_quadratic_delta1_20261102.json').read_text());denq=R.from_expr(s.sympify(qd['common_denominator']));parse=lambda arr:R.from_dict({tuple(m):s.QQ(c)for m,c in arr});qn=[{tuple(m):parse(c)for m,c in row}for row in qd['delta1_numerators']];labels=data['quadratic_product_labels']
def addto(r,m,p):
 c=r.get(m,R.zero)+p
 if c:r[m]=c
 elif m in r:del r[m]
def sparse(p):return[[list(m),str(c)]for m,c in p.items()]
rows=[]
for row in data['rows']:
 g=parse(row['raw_leading_content']);c=parse(row['primitive_cubic_content']);num={tuple(m):parse(p)*denq for m,p in row['first_numerator']}
 for index,poly in row['correction_numerator']:
  coeff=parse(poly)*g;i,shift=labels[index]
  for m,p in qn[i].items():addto(num,tuple(x+y for x,y in zip(m,shift)),-coeff*p)
 den=g*c*denq;gcd=den
 for p in num.values():gcd=gcd.gcd(p)
 den=den.exquo(gcd);num={m:p.exquo(gcd)for m,p in num.items()};assert den
 rows.append(dict(source=row['source'],denominator=sparse(den),numerators=[[list(m),sparse(p)]for m,p in num.items()]));print('CORRECTED',row['source'],'terms',sum(len(p)for p in num.values()),'den_degree',max(map(sum,den)),'sec',time.time()-start,flush=True)
p=Path(__file__).with_suffix('.json');assert not p.exists();p.write_text(json.dumps(dict(all_passed=True,scope='First-order corrected cubic generators whose constant term is the saved primitive remainder modulo constant quadrics.',rows=rows,seconds=time.time()-start),indent=2)+'\n');print('PASS',time.time()-start)
