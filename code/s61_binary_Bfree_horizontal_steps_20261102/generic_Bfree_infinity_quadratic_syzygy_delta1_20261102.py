from pathlib import Path
import json,time,hashlib
import sympy as s
from sympy.polys.rings import ring
st=time.time();p1=Path('work/generic_Bfree_infinity_quadratic_delta1_20261102.json');pt=Path('work/binary_hyperplane_infinity_traced_projection_20261102.json');pb=Path('work/binary_hyperplane_infinity_constant_quadratic_projection_20261102.json');q1=json.loads(p1.read_text());jt=json.loads(pt.read_text());B=json.loads(pb.read_text());tr=jt['projections']['4'];assert q1['all_passed'] and jt['all_passed'] and B['all_passed'];assert q1['new_variable_order'][:6]==B['normal_variables'][:6];assert q1['new_variable_order'][6:]==['r','s'] and B['normal_variables'][6:]==['free1','free2'];assert tr['monomials']==B['quartic_monomials'];assert tr['free_columns']==B['free_columns'];assert tr['rank']==B['rank']==132
R,a,e=ring('a,e',s.QQ);mons=[tuple(v) for v in B['quartic_monomials']];index={v:i for i,v in enumerate(mons)};pr=[[(int(j),s.QQ(c)) for j,c in row] for row in B['projection_sparse']];labels=tr['row_labels'];syzygies=[[(int(j),s.QQ(c)) for j,c in row] for row in tr['left_syzygies']];assert len(labels)==174 and len(syzygies)==42
Q0=[{tuple(m):s.QQ(c) for m,c in pp} for pp in B['quadrics']];Q1=[{tuple(m):R.from_dict({tuple(ex):s.QQ(c) for ex,c in pp}) for m,pp in qq} for qq in q1['delta1_numerators']]
def add(out,k,c):
 if c:
  v=out.get(k,0)+c
  if v:out[k]=v
  elif k in out:del out[k]
raw0=[];proj1=[]
for iq,mon in labels:
 r0={};r1={}
 for n,c in Q0[iq].items():add(r0,index[tuple(x+y for x,y in zip(mon,n))],c)
 raw0.append(r0)
 z={}
 for j,c in r0.items():
  for k,d in pr[j]:add(z,k,c*d)
 assert not z
 for n,c in Q1[iq].items():
  j=index[tuple(x+y for x,y in zip(mon,n))]
  for k,d in pr[j]:add(r1,k,c*d)
 proj1.append(r1)
ob=[]
for rel in syzygies:
 z0={};z1={}
 for row,c in rel:
  for j,d in raw0[row].items():add(z0,j,c*d)
  for j,d in proj1[row].items():add(z1,j,c*d)
 assert not z0
 ob.append(z1)
def encode(rows):return [[[j,[[list(ex),str(c)] for ex,c in poly.items()]] for j,poly in sorted(row.items())] for row in rows]
p=32003;ptval=(1,2)
for row in ob:
 for pp in row.values():
  for co in pp.values():assert int(co.denominator)%p
out={'all_passed':True,'scope':'Exact first-delta obstruction from the 42 constant quadratic left syzygies, projected to the constant 106-dimensional quartic quotient; no right75-kernel multiplication or higher-delta claim.','new_variable_order':B['normal_variables'],'quotient_monomials':[list(mons[i]) for i in B['free_columns']],'common_denominator':q1['common_denominator'],'common_denominator_factorization':q1['common_denominator_factorization'],'Q1_projected_shape':[174,106],'Q1_projected_numerators':encode(proj1),'obstruction_shape':[42,106],'obstruction_numerators':encode(ob),'nonzero_obstruction_rows':sum(bool(row) for row in ob),'nonzero_obstruction_entries':sum(len(row) for row in ob),'obstruction_total_core_terms':sum(len(v) for row in ob for v in row.values()),'Q1_projected_nonzero_entries':sum(len(row) for row in proj1),'Q1_projected_total_core_terms':sum(len(v) for row in proj1 for v in row.values()),'row_labels':labels,'left_syzygies':tr['left_syzygies'],'checks':{'same_238_monomials_and_106_free_columns':True,'all174_constant_rows_project_to_zero':True,'all42_constant_left_syzygies_zero':True,'no_new_denominator_from_constant_products':True},'input_sha256':{str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in [Path(__file__),p1,pt,pb]},'seconds':time.time()-st}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print({k:out[k] for k in ['nonzero_obstruction_rows','nonzero_obstruction_entries','obstruction_total_core_terms','Q1_projected_nonzero_entries','Q1_projected_total_core_terms','seconds']})
