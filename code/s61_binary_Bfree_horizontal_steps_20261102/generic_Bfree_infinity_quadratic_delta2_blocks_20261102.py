"""Exact delta0..2 quadratic blocks; old files are read-only dependencies."""
from pathlib import Path
import hashlib,json,time
import sympy as s
from sympy.polys.rings import ring
st=time.time();raw_path=Path('/tmp/poisson-hyperplane-Bfree-cubic-20261102/rows');frame_path=Path('reports/S61_BINARY_HYPERPLANE_INFINITY_QUADRATIC_NORMAL_FORM_20261102.json');cert_path=Path('code/certificates/s61_binary_hyperplane_infinity_quadrics_20261102.json');old_path=Path('work/generic_Bfree_infinity_quadratic_delta1_20261102.json');trace_path=Path('work/binary_hyperplane_infinity_traced_projection_20261102.json');projection_path=Path('work/binary_hyperplane_infinity_constant_quadratic_projection_20261102.json');old_syzygy_path=Path('work/generic_Bfree_infinity_quadratic_syzygy_delta1_20261102.json')
frame=json.loads(frame_path.read_text());cert=json.loads(cert_path.read_text());old=json.loads(old_path.read_text());tr=json.loads(trace_path.read_text())['projections']['4'];B=json.loads(projection_path.read_text());old_syzygy=json.loads(old_syzygy_path.read_text());assert frame['all_passed'] and old['all_passed'] and B['all_passed'] and old_syzygy['all_passed']
for name,sha in frame['input_sha256'].items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha
R,a,e=ring('a,e',s.QQ);K=R.to_field();parse=lambda v:K.from_expr(s.sympify(v));tw=[1,0,0,0,1,0,1,1];ids=frame['quadric_source_ids'];saved={r['source']:r for r in cert['rows']};raw={};metadata=[];prime=32003
with raw_path.open() as f:
 for _ in range(int(f.readline())):
  idx,wt,count=map(int,f.readline().split());terms=[]
  for _ in range(count):
   mon=tuple(map(int,f.readline().split()));words=f.readline().split()
   if idx not in ids:continue
   assert len(mon)==8 and sum(mon)==2 and mon[5]==0
   n=int(words[0]);assert len(words)==1+4*n
   for k in range(n):
    co=s.QQ(words[1+4*k]);assert int(co.denominator)%prime
    terms.append((mon,tuple(map(int,words[2+4*k:5+4*k])),co))
  if idx not in ids:continue
  assert wt==2
  D=max(ex[2]-sum(b*t for b,t in zip(mon,tw)) for mon,ex,c in terms);assert D==saved[idx]['d_shift'];g=parse(saved[idx]['content']);levels=[]
  for degree in range(3):
   temp={}
   for mon,ex,c in terms:
    if D-ex[2]+sum(b*t for b,t in zip(mon,tw))==degree:
     z=temp.setdefault(mon,{});z[ex[:2]]=z.get(ex[:2],s.QQ.zero)+c
   levels.append({mon:K(R.from_dict(co))/g for mon,co in temp.items() if co})
  assert levels[0]=={tuple(m):parse(c) for m,c in saved[idx]['terms']}
  raw[idx]=levels;metadata.append(dict(source=idx,d_shift=D,content=str(g.as_expr()),normal_terms_by_delta_degree=[len(p) for p in levels]))
 assert not f.read().strip()
assert set(raw)==set(ids)
Ti=[[parse(v) for v in row] for row in frame['coordinate_inverse8']];G=[[parse(v) for v in row] for row in frame['generator_matrix']]
def add(out,k,c):
 if c:
  v=out.get(k,0)+c
  if v:out[k]=v
  elif k in out:del out[k]
cache={}
for mon in {m for lv in raw.values() for pp in lv for m in pp}:
 ij=[i for i,n in enumerate(mon) for _ in range(n)];poly={}
 for j,c in enumerate(Ti[ij[0]]):
  if c:
   for k,d in enumerate(Ti[ij[1]]):
    if d:
     mm=[0]*8;mm[j]+=1;mm[k]+=1;add(poly,tuple(mm),c*d)
 cache[mon]=poly
def change(poly):
 out={}
 for mon,c in poly.items():
  for nn,d in cache[mon].items():add(out,nn,c*d)
 return out
new=[]
for level in range(3):
 rr=[change(raw[idx][level]) for idx in ids];res=[]
 for row in G:
  out={}
  for c,poly in zip(row,rr):
   if c:
    for mon,d in poly.items():add(out,mon,c*d)
  res.append(out)
 new.append(res)
for level in [0,1]:
 assert new[level]==[{tuple(m):parse(c) for m,c in pp} for pp in old['delta'+str(level)]]
assert new[0]==[{tuple(m):parse(c) for m,c in pp} for pp in B['quadrics']]
den={};nums={}
for level in [1,2]:
 common=R.one
 for pp in new[level]:
  for c in pp.values():common=common.lcm(c.denom)
 den[level]=common;nums[level]=[{m:c.numer*common.exquo(c.denom) for m,c in pp.items()} for pp in new[level]]
 for pp in nums[level]:
  for co in pp.values():
   for c in co.values():assert int(c.denominator)%prime
 cv=sum((c*2**ex[1] for ex,c in common.items()),s.QQ.zero)
 assert int(cv.denominator)%prime and int(cv.numerator)%prime
assert den[1]==R.from_expr(s.sympify(old['common_denominator']))
mons=[tuple(v) for v in B['quartic_monomials']];assert mons==[tuple(v) for v in tr['monomials']];assert B['free_columns']==tr['free_columns'];index={m:i for i,m in enumerate(mons)}
labels=tr['row_labels'];pr=[[(int(j),s.QQ(c)) for j,c in row] for row in B['projection_sparse']];syzygies=[[(int(j),s.QQ(c)) for j,c in row] for row in tr['left_syzygies']];lifts=[[(int(j),s.QQ(c)) for j,c in row] for row in tr['relation_lifts']];assert len(labels)==174 and len(syzygies)==42 and len(lifts)==132
full={};projected={}
def project(row):
 out={}
 for j,c in row.items():
  for k,d in pr[j]:add(out,k,c*d)
 return out
for level in [0,1,2]:
 full[level]=[];projected[level]=[]
 coeffs=new[0] if level==0 else nums[level]
 for iq,mon in labels:
  row={}
  for nn,c in coeffs[iq].items():add(row,index[tuple(x+y for x,y in zip(mon,nn))],c)
  full[level].append(row);projected[level].append(project(row))
assert all(not r for r in projected[0])
def combine(transform,rows):
 out=[]
 for rel in transform:
  v={}
  for i,c in rel:
   for j,d in rows[i].items():add(v,j,c*d)
  out.append(v)
 return out
assert all(not r for r in combine(syzygies,full[0]))
assert combine(lifts,full[0])==[{int(j):parse(c) for j,c in row} for row in tr['relations']]
# RREF pivots must form the identity, so its lifted top block may be used for Schur correction.
rref=combine(lifts,full[0]);piv=tr['pivot_columns']
assert all(rref[i].get(j,0)==int(i==k) for i in range(132) for k,j in enumerate(piv))
top1full=combine(lifts,full[1]);top1proj=combine(lifts,projected[1]);syzygy2=combine(syzygies,projected[2])
def decode(rows):return [{int(j):R.from_dict({tuple(ex):s.QQ(c) for ex,c in pp}) for j,pp in row} for row in rows]
assert projected[1]==decode(old_syzygy['Q1_projected_numerators'])
assert combine(syzygies,projected[1])==decode(old_syzygy['obstruction_numerators'])
assert [project(r) for r in top1full]==top1proj
def encode(rows):return [[[j,[[list(ex),str(c)] for ex,c in poly.items()]] for j,poly in sorted(row.items())] for row in rows]
def count(rows):return {'nonzero_rows':sum(bool(row) for row in rows),'nonzero_entries':sum(len(row) for row in rows),'core_terms':sum(len(c) for row in rows for c in row.values())}
blocks={'full1_174x238':full[1],'full2_174x238':full[2],'projected1_174x106':projected[1],'projected2_174x106':projected[2],'top1_full132x238':top1full,'top1_projected132x106':top1proj,'syzygy2_projected42x106':syzygy2}
paths=[Path(__file__),raw_path,frame_path,cert_path,old_path,trace_path,projection_path,old_syzygy_path]
out={'all_passed':True,'scope':'Exact delta0..2 of supplied six quadrics and their174 weighted-four products, with fixed constant projection/top lifts/left syzygies. No multiplication by the large right kernels, no second Schur correction inferred.','normal_variable_order':frame['new8_order'],'projection_variable_aliases':{'r':'free1','s':'free2'},'rows':metadata,'quartic_monomials':[list(m) for m in mons],'quotient_free_columns':B['free_columns'],'rref_pivot_columns':piv,'row_labels':labels,'denominators':{str(k):str(v.as_expr()) for k,v in den.items()},'denominator_factorizations':{str(k):str(s.factor(v.as_expr())) for k,v in den.items()},'delta2_quadrics':[[[list(m),str(c.as_expr())] for m,c in pp.items()] for pp in new[2]],'delta2_quadrics_numerators':[[[list(m),[[list(ex),str(c)] for ex,c in pp.items()]] for m,pp in row.items()] for row in nums[2]],'blocks_numerators':{k:encode(v) for k,v in blocks.items()},'block_delta_degree':{k:2 if ('2_' in k or k.startswith('full2') or k.startswith('projected2')) else 1 for k in blocks},'block_stats':{k:count(v) for k,v in blocks.items()},'checks':{'all_raw_delta0_equal_input':True,'delta0_delta1_match_previous_exact_report':True,'all174_constant_rows_project_zero':True,'all42_constant_syzygies_zero':True,'all132_lifts_equal_saved_rref':True,'rref_pivot_block_identity':True,'projected_delta1_and_syzygies_match_previous_report':True,'all_coefficient_denominators_and_point_good_mod32003':True},'input_sha256':{str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in paths},'seconds':time.time()-st}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print('denominators',out['denominator_factorizations']);print('stats',out['block_stats']);print('PASS',out['seconds'])
