from pathlib import Path
import hashlib,json,time
import sympy as s
from sympy.polys.rings import ring
st=time.time();path=Path('/tmp/poisson-hyperplane-Bfree-cubic-20261102/rows');frame_path=Path('reports/S61_BINARY_HYPERPLANE_INFINITY_QUADRATIC_NORMAL_FORM_20261102.json');cert_path=Path('code/certificates/s61_binary_hyperplane_infinity_quadrics_20261102.json')
frame=json.loads(frame_path.read_text());cert=json.loads(cert_path.read_text());assert frame['all_passed']
for name,sha in frame['input_sha256'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha
R,a,e=ring('a,e',s.QQ);K=R.to_field();parse=lambda v:K.from_expr(s.sympify(v));tw=[1,0,0,0,1,0,1,1];ids=frame['quadric_source_ids'];saved={r['source']:r for r in cert['rows']};rr={};extracted={};metadata=[]
with path.open() as f:
 for _ in range(int(f.readline())):
  idx,wt,count=map(int,f.readline().split());terms=[]
  for _ in range(count):
   mon=tuple(map(int,f.readline().split()));words=f.readline().split()
   if idx not in ids:continue
   assert sum(mon)==2 and mon[5]==0
   n=int(words[0]);assert len(words)==1+4*n
   for k in range(n):terms.append((mon,tuple(map(int,words[2+4*k:5+4*k])),s.QQ(words[1+4*k])))
  if idx not in ids:continue
  assert wt==2
  D=max(ex[2]-sum(b*t for b,t in zip(mon,tw)) for mon,ex,c in terms);assert D==saved[idx]['d_shift']
  levels=[];g=parse(saved[idx]['content'])
  for degree in [0,1]:
   raw={}
   for mon,ex,c in terms:
    if D-ex[2]+sum(b*t for b,t in zip(mon,tw))==degree:
     z=raw.setdefault(mon,{});z[ex[:2]]=z.get(ex[:2],s.QQ.zero)+c
   levels.append({mon:K(R.from_dict(co))/g for mon,co in raw.items() if co})
  expected={tuple(m):parse(c) for m,c in saved[idx]['terms']}
  assert levels[0]==expected
  rr[idx]=levels;metadata.append(dict(source=idx,d_shift=D,content=str(g.as_expr()),original_core_terms=len(terms),delta0_normal_terms=len(levels[0]),delta1_normal_terms=len(levels[1])))
 assert not f.read().strip()
assert set(rr)==set(ids)
Ti=[[parse(v) for v in row] for row in frame['coordinate_inverse8']];G=[[parse(v) for v in row] for row in frame['generator_matrix']]
cache={}
for mon in {m for lv in rr.values() for pp in lv for m in pp}:
 ij=[i for i,n in enumerate(mon) for _ in range(n)];assert len(ij)==2
 poly={}
 for j,c in enumerate(Ti[ij[0]]):
  if not c:continue
  for k,d in enumerate(Ti[ij[1]]):
   if not d:continue
   mm=[0]*8;mm[j]+=1;mm[k]+=1;mm=tuple(mm);poly[mm]=poly.get(mm,K.zero)+c*d
 cache[mon]={m:c for m,c in poly.items() if c}
def change(poly):
 out={}
 for mon,c in poly.items():
  for nn,d in cache[mon].items():out[nn]=out.get(nn,K.zero)+c*d
 return {m:c for m,c in out.items() if c}
new=[]
for level in [0,1]:
 old=[change(rr[idx][level]) for idx in ids];res=[]
 for row in G:
  out={}
  for c,poly in zip(row,old):
   if c:
    for mon,d in poly.items():out[mon]=out.get(mon,K.zero)+c*d
  res.append({m:c for m,c in out.items() if c})
 new.append(res)
y=s.symbols('L W C P E H r s');
for i,g in enumerate(frame['normal_generators']):
 expected={m:parse(str(c)) for m,c in s.Poly(s.sympify(g,locals={str(z):z for z in y}),*y).terms()};assert new[0][i]==expected,(i,new[0][i],expected)
common=R.one
for pp in new[1]:
 for c in pp.values():common=common.lcm(c.denom)
nums=[{m:c.numer*common.exquo(c.denom) for m,c in pp.items()} for pp in new[1]]
p=32003
for value in [common]+[c for pp in nums for c in pp.values()]:
 for co in value.values():assert int(co.denominator)%p
cv=sum((c*2**ex[1] for ex,c in common.items()),s.QQ.zero);assert cv
cvp=int(cv.numerator)*pow(int(cv.denominator),-1,p)%p;assert cvp
out={'scope':'Exact delta0 and delta1 of the six supplied B-free quadrics after fixed weighted d-infinity normalization and constant normal-form coordinate/generator changes; no higher-order claim.','all_passed':True,'source_rows':ids,'old_d_infinity_weights':tw,'new_variable_order':frame['new8_order'],'rows':metadata,'delta0':[[[list(m),str(c.as_expr())] for m,c in pp.items()] for pp in new[0]],'delta1':[[[list(m),str(c.as_expr())] for m,c in pp.items()] for pp in new[1]],'common_denominator':str(common.as_expr()),'common_denominator_factorization':str(s.factor(common.as_expr())),'delta1_numerators':[[[list(m),[[list(ex),str(c)] for ex,c in co.items()]] for m,co in pp.items()] for pp in nums],'prime':p,'good_point':[1,2],'common_denominator_point_residue':cvp,'delta1_normal_term_counts':[len(pp) for pp in nums],'delta1_numerator_core_terms':sum(len(co) for pp in nums for co in pp.values()),'checks':{'all_original_delta0_equal_saved_input':True,'all_transformed_delta0_equal_constant_generators':True,'all_input_and_output_coefficient_denominators_p_units':True,'common_denominator_nonzero_modp_at_point':True,'coordinate_and_generator_changes_independent_of_delta':True},'input_sha256':{str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in [Path(__file__),path,frame_path,cert_path]},'seconds':time.time()-st}
# Check p-goodness of the raw rational input too, without relying on its transformed cancellation.
with path.open() as f:
 for _ in range(int(f.readline())):
  idx,wt,count=map(int,f.readline().split())
  for _ in range(count):
   f.readline();z=f.readline().split()
   if idx in ids:
    for k in range(int(z[0])):assert int(s.QQ(z[1+4*k]).denominator)%p
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print({k:out[k] for k in ['all_passed','common_denominator_factorization','common_denominator_point_residue','delta1_normal_term_counts','delta1_numerator_core_terms','seconds']})
