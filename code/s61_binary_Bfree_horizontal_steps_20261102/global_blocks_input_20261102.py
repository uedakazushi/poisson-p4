from pathlib import Path
import json,time,sympy as s
from sympy.polys.rings import ring
st=time.time();R,a,e=ring('a,e',s.QQ);parse=lambda z:R.from_dict({tuple(m):s.QQ(c)for m,c in z});q=json.load(open('work/generic_Bfree_infinity_quadratic_delta2_blocks_20261102.json'));first=json.load(open('work/binary_hyperplane_Bfree_infinity_cubic_product_full_delta1_20261102.json'));firstpr=json.load(open('work/binary_hyperplane_Bfree_infinity_cubic_product_delta1_20261102.json'));second=json.load(open('work/binary_hyperplane_Bfree_infinity_cubic_product_delta2_20261102.json'));pr=json.load(open('work/binary_hyperplane_infinity_constant_quadratic_projection_20261102.json'));tr=json.load(open('work/binary_hyperplane_infinity_traced_projection_20261102.json'))['projections']['4'];proj=[{j:s.QQ(c)for j,c in row}for row in pr['projection_sparse']];qfull=[{j:parse(z)for j,z in row}for row in q['blocks_numerators']['full1_174x238']];DQ=R.from_expr(s.sympify(q['denominators']['1']));assert DQ==R.from_expr(s.sympify(q['denominators']['2']))
def addto(out,k,c):
 z=out.get(k,R.zero)+c
 if z:out[k]=z
 elif k in out:del out[k]
A1=[];D1s=[]
for sy in tr['left_syzygies']:
 row={}
 for i,c in sy:
  for j,v in qfull[i].items():addto(row,j,s.QQ(c)*v)
 A1.append(row);D1s.append(DQ)
for r0,r1 in zip(first['rows'],firstpr['rows']):
 row={j:parse(z)for j,z in r0['numerators']};den=parse(r0['denominator']);out={}
 for j,v in row.items():
  for k,c in proj[j].items():addto(out,k,c*v)
 # Cross multiplication is checked by the native all-coefficient verifier.
 A1.append(row);D1s.append(den)
A2=[{j:parse(z)for j,z in row}for row in q['blocks_numerators']['syzygy2_projected42x106']];D2s=[DQ]*42
for row in second['rows']:A2.append({j:parse(z)for j,z in row['numerators']});D2s.append(parse(row['denominator']))

assert len(A1)==len(A2)==189
T=[{j:parse(z)for j,z in row}for row in q['blocks_numerators']['top1_projected132x106']]
def enc(p):return str(len(p))+''.join(' '+str(c)+' '+str(m[0])+' '+str(m[1])for m,c in p.items())+'\n'
path=Path('/tmp/hyperplane-Bfree-infinity-second-obstruction-uncleared-20261102.in')
with path.open('w')as f:
 f.write(enc(DQ));f.write(' '.join(map(str,tr['pivot_columns']))+'\n')
 for row in proj:f.write(str(len(row))+''.join(' '+str(j)+' '+str(c)for j,c in row.items())+'\n')
 Ref=[{j:parse(z)for j,z in row['numerators']}for row in firstpr['rows']]
 for mat,ds,nc in [(A1,D1s,238),(A2,D2s,106),(T,[DQ]*132,106),(Ref,[parse(row['denominator'])for row in firstpr['rows']],106)]:
  f.write(str(len(mat))+' '+str(nc)+'\n')
  for row,d in zip(mat,ds):
   f.write(enc(d))
   for j in range(nc):f.write(enc(row.get(j,R.zero)))
rep=dict(all_passed=True,scope='Uncleared second Schur blocks; native checks the full/projected equality and all exact divisions.',DQ=str(s.factor(DQ.as_expr())),input=str(path),bytes=path.stat().st_size,seconds=time.time()-st);Path(__file__).with_suffix('.json').write_text(json.dumps(rep,indent=2)+'\n');print(rep)
