#!/usr/bin/env python3
"""Independent full9 first640 callback reader; definitions only on import.

Inputs are previously accepted source caches. No source IFT or gauge is
solved here. Invocation belongs to the frozen mu(T2) postcheck.
"""
from fractions import Fraction as Q

def sparse(p):return {int(i):Q(c) for i,c in p if Q(c)}
def const(c):return {0:Q(c)} if c else {}
def add(*ps):
 out={}
 for p in ps:
  for i,c in p.items():
   v=out.get(i,Q(0))+c
   if v:out[i]=v
   else:out.pop(i,None)
 return out
def scale(p,c):return {i:c*v for i,v in p.items() if c*v}
def times(a,b):
 a0=a.get(0,Q(0));b0=b.get(0,Q(0));out=const(a0*b0)
 for i in (set(a)|set(b))-{0}:
  c=a0*b.get(i,Q(0))+b0*a.get(i,Q(0))
  if c:out[i]=c
 return out
def full_to_new(p):
 assert len(p)==10
 out={i:Q(c) for i,c in enumerate(p[:8]) if Q(c)}
 r=Q(p[8])+Q(p[9]);ell=Q(p[9])/5
 if r:out[8]=r
 if ell:out[9]=ell
 return out
def encode(p):return [[i,str(c)] for i,c in sorted(p.items())]

def first640_full10(C,B,J,Q4,source_values,old_first,X0,boundary_ids=()):
 """Dictionaries indexed row56 then(0,zA0..6,rA,ell)."""
 ops=B['operator_scalar_nodes'];roots=[i for col in B['B48_by640_sparse_columns'] for row,i in col]
 seen=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i);n=ops[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
  else:assert n[0] in ['c','actual_C0_scalar_function']
 atoms=sorted({ops[i][1] for i in seen if ops[i][0]=='actual_C0_scalar_function'})
 assert len(roots)==2808 and len(seen)==1302 and len(atoms)==79
 vals={int(i):full_to_new(p) for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA'] for i,p in J[key].items()}
 basis=[tuple(m) for m in source_values['basis56']];ix={m:i for i,m in enumerate(basis)}
 def first(p):
  p=sparse(p);v=const(p.get(0,Q(0)))
  for k in range(9):
   c=p.get(ix[tuple(int(i==k) for i in range(9))],Q(0))
   if c:v[k+1]=c
  return v
 vals.update({i:first(p) for i,p in source_values['all164_original_source_roots']})
 for j,i in enumerate(C['actual_C0_basis8_scalar_nodes']):vals[i]={j+1:Q(1)} if j<7 else {}
 words=J['source2_monomials'];prefix=[]
 for record in C['all35_regular_coordinate_divisions']:
  f=record['original_pure_CA_source_function'];coeff=J['both57_source2']['CA']['pure_function2jets'][str(f)]
  for j,q in enumerate(record['all7_regular_coordinate_divisions']):
   raw=[Q(0)]*10
   for word,c in zip(words,coeff):
    c=Q(c)
    if not c or j not in word or any(k>j and k<7 for k in word):continue
    rest=list(word);rest.remove(j);assert len(rest)<=1
    raw[0 if not rest else rest[0]+1]+=c
   jet=full_to_new(raw);assert jet.get(0,Q(0))==Q(q['origin']);vals[q['node']]=jet
   prefix.append({'node':q['node'],'source':f,'prefix_index':j,'full10':encode(jet)})
 assert len(prefix)==35
 mons=[tuple(m) for m in Q4['monomials_bivariate_through4']];adapter=[[Q(c) for c in row] for row in J['YB_to_YA_linear']]
 gradients=[];source4_checks=0
 for family in ['CB','CA']:
  data=Q4['families'][family];g={m:Q(c) for m,c in zip(mons,data['g_all15_coefficients']) if Q(c)}
  assert all(sum(m)>=3 for m in g)
  for num,lin in zip(data['four_numerator_all15_coefficients'],data['quotient_firsts4_by2']):
   a,b=map(Q,lin);product={}
   assert all(not Q(c) for m,c in zip(mons,num) if sum(m)<=3)
   for m,c in g.items():
    for j,d in enumerate([a,b]):
     n=list(m);n[j]+=1;n=tuple(n)
     if sum(n)<=4:product[n]=product.get(n,Q(0))+c*d
   assert all(product.get(m,Q(0))==Q(c) for m,c in zip(mons,num));source4_checks+=15
   grad=[Q(0)]*9
   if family=='CB':grad=[a*x+b*y for x,y in zip(adapter[7],adapter[8])]
   else:grad[7:]=[a,b]
   gradients.append(grad)
 assert gradients==[list(map(Q,p)) for p in Q4['all8_regular_function_first_gradients_in_FULL9_CA_F9']]
 ids=Q4['original_regular_function_nodes_unmodified']['old_CB']+Q4['original_regular_function_nodes_unmodified']['new_CA']
 for i,p in zip(ids,gradients):vals[i]=full_to_new([0]+p)
 tail=next(x for x in B['accepted_actual_regular_function_identities'] if x['kind']=='accepted_actual_CB_f63_tail')
 assert tail['function']==J['both57_source2']['CB']['f63'] and tail['function']!=tail['gB']
 assert B['old_multiplier_last_is_f63_not_new_gA'];vals[tail['gB']]={}
 dag=C['DAG_scalar_nodes'];memo=dict(vals)
 # Exact negative aliases preserve the bound function as a boundary.
 changed=True
 while changed:
  changed=False
  for i,p in list(memo.items()):
   if dag[i][0]=='n':
    j=dag[i][1];v=scale(p,-1)
    if j in memo:assert memo[j]==v
    else:memo[j]=v;changed=True
 def original(i):
  if i in memo:return memo[i]
  n=dag[i];op=n[0]
  if op=='c':p=const(Q(n[1],n[2]))
  elif op=='n':p=scale(original(n[1]),-1)
  elif op=='a':p=add(original(n[1]),original(n[2]))
  elif op=='m':
   a=original(n[1]);p=times(a,original(n[2])) if a else {}
  elif op=='i':
   a=original(n[1]);c=a.get(0,Q(0));assert c,('NONUNIT',i)
   p={j:-q/c**2 for j,q in a.items() if j};p[0]=1/c
  else:raise AssertionError(('UNBOUND_REGULAR_SOURCE',i,n))
  memo[i]=p;assert sum(map(len,memo.values()))<150000;return p
 atom_values={i:original(i) for i in atoms}
 raw=J['both57_source2']['CB']['pure_function2jets'][str(tail['function'])]
 f=[Q(raw[words.index([j])]) for j in range(9)];adapted=[sum((f[i]*adapter[i][j] for i in range(9)),Q(0)) for j in range(9)]
 assert original(tail['function'])==full_to_new([Q(raw[0])]+adapted)
 typed={}
 for i in sorted(seen):
  n=ops[i]
  if n[0]=='c':p=const(Q(n[1],n[2]))
  elif n[0]=='actual_C0_scalar_function':p=atom_values[n[1]]
  elif n[0]=='sum':p=add(*(scale(typed[j],Q(a,b)) for j,a,b in n[1]))
  elif n[0]=='product':
   p=const(Q(n[1],n[2]))
   for j in n[3]:p=times(p,typed[j])
  else:raise AssertionError(n)
  typed[i]=p
 output=[];P81=Q4['all640_exact_D0_firstjet81_columns']
 for j,col in enumerate(B['B48_by640_sparse_columns']):
  values={r:typed[i] for r,i in col if typed[i]}
  assert len(P81[j])==81
  for r in range(1,9):
   p=full_to_new([0]+P81[j][9*r:9*(r+1)])
   if p:values[47+r]=p
  if 'columns_sparse' in X0:assert all(values.get(i,{}).get(0,Q(0))==X0['columns_sparse'][j].get(i,Q(0)) for i in range(56))
  else:assert all(values.get(i,{}).get(0,Q(0))==Q(X0['full_X0_56_by768'][i][j]) for i in range(56))
  old={int(i):sparse(p) for i,p in old_first['source_columns640_rows56_sparse_dual9'][j]}
  assert {i:{k:c for k,c in p.items() if k<9} for i,p in values.items() if any(k<9 for k in p)}==old
  output.append(values)
 return output,{'all79_atoms_1302typed_and2808_B_roots_rebuilt':True,'all35_ordered_prefixes_full9_no_factor2':True,'all120_source4_numerator_coefficients_verified':source4_checks==120,'old_CB_f63_distinct_from_zero7_gB_plane':True,'all640_old8_projections_and_origins_exact':True,'source_prefix35':prefix,'all79_atom_values_full10':[[i,encode(p)] for i,p in sorted(atom_values.items())],'source_boundary_values_full10':[[i,encode(original(i))] for i in boundary_ids],'source_nonzero_full10_coefficients':sum(len(p) for col in output for p in col.values())}
