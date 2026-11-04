"""Independent reusable source3/B11 readers; no implicit solve or worker main."""
from fractions import Fraction as Q
from collections import Counter
from ambient_B18_T_gdiff_second_univariate5_review_20261104_v63 import Symbolic,closure,numeric_eval,series_mul,rec

def ordinary3(D,roots):
 nodes=D['C0']['DAG_scalar_nodes'];joint=D['joint'];common={}
 for i,p in joint['all105_unknown_jets']:
  q=[Q(0)]*4
  for m,c in zip(joint['monomials'],p):q[sum(m)]+=Q(c)
  common[i]=tuple(q)
 y=[i for i,p in D['copy2']['actual_YB_last2_from_saved105']];eq=joint['all105_equation_original_node_order'];v,units=numeric_eval(nodes,eq+y,common,4)
 assert all(v[i]==(0,1 if j in [55,56] else 0,0,0) for j,i in enumerate(eq))
 copied={};labels=dict(D['copy1']['explicit57_source_copy_bindings'])
 for i,p in D['fourth']['families']['CB']['all57_source4_jets']:
  q=(Q(0),)*4
  for m,c in zip(D['fourth']['monomials_bivariate_through4'],p):
   if sum(m)>3 or not Q(c):continue
   t=(Q(c),0,0,0)
   for k,n in enumerate(m):
    for _ in range(n):t=series_mul(t,v[y[k]],4)
   q=tuple(a+b for a,b in zip(q,t))
  copied[labels[i]]=q
 ring=Symbolic([(i,) for i in range(4)]);bind={i:({k:c for k,c in enumerate(p) if c},{}) for i,p in (common|copied).items()};opaque={}
 for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,p in D['FJ'][key].items():opaque[int(i)]=(Q(p[0]),Q(p[8])+Q(p[9]))
 for row in D['atom1']['all35_ordered_CA_prefix_firsts']:
  p=rec(row['first9_bar']);opaque[row['node']]=(p.get(0,Q(0)),p.get(8,Q(0)))
 for i,p in D['atom1']['all8_g_quotient_overrides']:
  p=rec(p);opaque[i]=(p.get(0,Q(0)),p.get(8,Q(0)))
 for i,(c,a) in opaque.items():bind[i]=({k:c for k,c in enumerate([c,a]) if c},{(i,(7,7),2):Q(1),(i,(7,7,7),3):Q(1)})
 ga=[Q(0)]*4
 for m,c in zip(joint['monomials'],joint['g_CA_all10_coefficients']):ga[sum(m)]+=Q(c)
 for j,i in enumerate(D['C0']['actual_C0_basis8_scalar_nodes']):bind[i]=({k:c for k,c in enumerate(ga) if c},{}) if j==7 else ring.c(0)
 copyeq=D['copy1']['all57_copied_equation_indices'];want=sorted(set(roots)|set(copyeq));order=closure(nodes,want,set(bind));uses=Counter()
 for i in order:
  if i in bind:continue
  n=nodes[i];uses.update(n[1:3] if n[0] in ['a','m'] else n[1:] if n[0] in ['n','i'] else [])
 uses.update(want);values={};peak=0
 for i in order:
  n=nodes[i];op=n[0]
  if i in bind:p=bind[i]
  elif op=='c':p=ring.c(Q(n[1],n[2]))
  elif op=='a':p=ring.add(values[n[1]],values[n[2]])
  elif op=='n':p=ring.scale(values[n[1]],-1)
  elif op=='m':p=ring.mul(values[n[1]],values[n[2]])
  elif op=='i':p=ring.inv(values[n[1]]);units+=1
  else:raise AssertionError((i,n))
  values[i]=p;peak=max(peak,sum(len(a)+len(b) for a,b in values.values()));assert peak<150000
  if i not in bind:
   for j in n[1:3] if op in ['a','m'] else n[1:] if op in ['n','i'] else []:
    uses[j]-=1
    if not uses[j]:values.pop(j,None)
 assert all(values[i]==ring.c(0) for i in copyeq)
 return {i:values[i] for i in roots},{'source105_and_copy57_equations_checked':True,'unit_checks':units,'QQ_peak':peak}

def typed_B11(D,roots):
 def ex(m):return tuple(m.count(i) for i in range(8))
 bm=[(0,)*8,ex((7,)),ex((7,7)),ex((0,)),ex((1,)),ex((4,)),ex((5,)),ex((0,7)),ex((1,7)),ex((4,7)),ex((5,7))]
 ring=Symbolic(bm);ops=D['chain']['operator_scalar_nodes'];cv={i:rec(v) for i,v,g in D['cache']['all108_composite_values_and_opaque_second_gradients']}
 assert all(not g for i,v,g in D['cache']['all108_composite_values_and_opaque_second_gradients'])
 atoms={i:rec(p) for i,p in D['atom1']['all79_atom_values']};lin={0:0,8:1,1:3,2:4,5:5,6:6};needed=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);n=ops[i]
  if n[0]=='sum':todo.extend(j for j,a,b in n[1])
  elif n[0]=='product':todo.extend(n[3])
 values={}
 for i in sorted(needed):
  n=ops[i]
  if n[0]=='c':p=ring.c(Q(n[1],n[2]))
  elif n[0]=='actual_C0_scalar_function':
   j=n[1]
   if j in cv:p=(cv[j],{})
   else:p=({b:atoms[j][k] for k,b in lin.items() if atoms[j].get(k)},{(j,tuple(k for k,c in enumerate(bm[b]) for _ in range(c)),b):Q(1) for b in [2,7,8,9,10]})
  elif n[0]=='sum':
   p=ring.c(0)
   for j,a,b in n[1]:p=ring.add(p,values[j],Q(a,b))
  elif n[0]=='product':
   p=ring.c(Q(n[1],n[2]))
   for j in n[3]:p=ring.mul(p,values[j])
  else:raise AssertionError(n)
  values[i]=p
 return {i:values[i] for i in roots},cv

def coefficient(p,index):
 return p[0].get(index,Q(0)),{(i,tag):c for (i,tag,out),c in p[1].items() if out==index}
