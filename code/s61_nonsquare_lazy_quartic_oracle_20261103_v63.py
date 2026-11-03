#!/usr/bin/env python3
"""Exact current-chart weighted quartic DAG and lazy coefficient API.

No parameter polynomial matrix is expanded. Coefficients remain arithmetic
circuits over explicitly encoded QQ(a,e,d) atoms with known-unit denominators.
"""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

class Oracle:
 def __init__(self,variable_weights=(1,1,1,1,1,2,1,1)):
  self.R,self.a,self.e,self.d=ring('a,e,d',s.QQ);self.K=self.R.to_field();self.variable_weights=list(variable_weights)
  self.known=[self.a,500*self.a+90*self.e**2+3*self.e-24,(50*self.a+18*self.e**2-9*self.e-2)*self.d+60*self.e-40,100*self.a+18*self.e**2-9*self.e-2]
  self.atoms=[];self.atom_index={};self.nodes=[['z']];self.node_index={('z',):0};self.weights=[-1];self.degrees=[-1];self.coefficient_cache={}
  self.atom(self.K.zero);self.atom(self.K.one)
 def atom(self,c):
  c=self.K(c);key=str(c)
  if key in self.atom_index:return self.atom_index[key]
  den=c.denom
  for f in self.known:
   while den.rem(f)==0:den=den.exquo(f)
  assert den.is_ground and den
  i=len(self.atoms);self.atoms.append(c);self.atom_index[key]=i;return i
 def intern(self,key,node,w,d):
  if key in self.node_index:return self.node_index[key]
  assert len(self.nodes)<100000
  i=len(self.nodes);self.nodes.append(node);self.node_index[key]=i;self.weights.append(w);self.degrees.append(d);return i
 def var(self,i):return self.intern(('v',i),['v',i],self.variable_weights[i],1)
 def const(self,c):
  c=self.K(c)
  if not c:return 0
  j=self.atom(c);return self.intern(('c',j),['c',j],0,0)
 def lin(self,terms):
  combined={}
  for c,n in terms:
   if n and c:combined[n]=combined.get(n,self.K.zero)+self.K(c)
  combined={n:c for n,c in combined.items()if c}
  if not combined:return 0
  if len(combined)==1:
   n,c=next(iter(combined.items()))
   if c==self.K.one:return n
  ws={self.weights[n]for n in combined};assert len(ws)==1
  ds={self.degrees[n]for n in combined};degree=next(iter(ds))if len(ds)==1 else None
  terms=tuple((self.atom(c),n)for n,c in sorted(combined.items()));return self.intern(('l',terms),['l',[list(t)for t in terms]],next(iter(ws)),degree)
 def mul(self,i,j):
  if not i or not j:return 0
  w=self.weights[i]+self.weights[j]
  if w>4:return 0
  if j<i:i,j=j,i
  d=None if self.degrees[i]is None or self.degrees[j]is None else self.degrees[i]+self.degrees[j]
  return self.intern(('m',i,j),['m',i,j],w,d)
 def coefficient(self,node,monomial):
  """Return a weight0 DAG node for one coefficient; do not expand core atoms."""
  monomial=tuple(monomial);assert len(monomial)==len(self.variable_weights)and min(monomial)>=0
  key=(node,monomial)
  if key in self.coefficient_cache:return self.coefficient_cache[key]
  if not node or sum(w*m for w,m in zip(self.variable_weights,monomial))!=self.weights[node]:ans=0
  else:
   nd=self.nodes[node]
   if nd[0]=='v':ans=self.const(1)if monomial==tuple(int(i==nd[1])for i in range(len(monomial)))else 0
   elif nd[0]=='c':ans=node if not any(monomial)else 0
   elif nd[0]=='l':ans=self.lin([(self.atoms[c],self.coefficient(n,monomial))for c,n in nd[1]])
   else:
    assert nd[0]=='m';terms=[]
    for left in itertools.product(*(range(m+1)for m in monomial)):
     if sum(w*m for w,m in zip(self.variable_weights,left))!=self.weights[nd[1]]:continue
     right=tuple(x-y for x,y in zip(monomial,left));i=self.coefficient(nd[1],left);j=self.coefficient(nd[2],right)
     if i and j:terms.append((1,self.mul(i,j)))
    ans=self.lin(terms)
  self.coefficient_cache[key]=ans;return ans
 def core_value(self,node):
  """Explicit optional core expansion of a coefficient node, not used in build."""
  cache={0:self.K.zero}
  def go(i):
   if i in cache:return cache[i]
   nd=self.nodes[i];assert self.weights[i]==0
   if nd[0]=='c':v=self.atoms[nd[1]]
   elif nd[0]=='m':v=go(nd[1])*go(nd[2])
   else:assert nd[0]=='l';v=sum((self.atoms[c]*go(n)for c,n in nd[1]),self.K.zero)
   cache[i]=v;return v
  return go(node)
 def encoded(self):
  enc=lambda p:[[list(ex),str(c)]for ex,c in p.items()]
  return {'variable_weights':self.variable_weights,'atoms':[[enc(c.numer),enc(c.denom)]for c in self.atoms],'nodes':self.nodes,'node_weights':self.weights,'node_ordinary_degrees':self.degrees}
 @classmethod
 def load(cls,path):
  data=json.loads(Path(path).read_text());record=data['arithmetic'];o=cls(record['variable_weights']);o.atoms=[];o.atom_index={}
  load=lambda ts:o.R.from_dict({tuple(ex):s.QQ(c)for ex,c in ts})
  for num,den in record['atoms']:o.atom(o.K(load(num))/o.K(load(den)))
  o.nodes=record['nodes'];o.weights=record['node_weights'];o.degrees=record['node_ordinary_degrees'];o.node_index={('z',):0}
  for i,nd in enumerate(o.nodes[1:],1):
   key=(nd[0],tuple(tuple(t)for t in nd[1]))if nd[0]=='l'else tuple(nd);o.node_index[key]=i
  return o,data

def fixture():
 o=Oracle((1,1,2));u,v,h=[o.var(i)for i in range(3)]
 f=o.lin([(1,o.mul(u,u)),(1,o.mul(v,v)),(3,h)]);g=o.lin([(1,o.mul(u,v)),(-1,h)]);root=o.mul(f,g)
 U,V,H=s.symbols('U V H');expected=s.Poly((U*U+V*V+3*H)*(U*V-H),U,V,H)
 checked=0
 for mon in itertools.product(range(5),range(5),range(3)):
  if mon[0]+mon[1]+2*mon[2]!=4:continue
  assert o.core_value(o.coefficient(root,mon))==o.K(expected.coeff_monomial(U**mon[0]*V**mon[1]*H**mon[2]));checked+=1
 return checked

def build(output):
 from s61_model import S61
 from s61_flag_transverse_obstruction_20261101 import specialize
 start=time.monotonic();output.mkdir(parents=True,exist_ok=False);o=Oracle();R,K=o.R,o.K;a,e,d=o.a,o.e,o.d
 chartpath=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json';oldpath=ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json';chart=json.loads(chartpath.read_text());old=json.loads(oldpath.read_text())['linear']
 for path in [chartpath,oldpath]:
  for name,want in json.loads(path.read_text()).get('input_sha256',{}).items():assert sha(ROOT/name)==want,name
 localweights=chart['local_coordinate_weights'];rowweights=chart['equation_weights'];selected=chart['selected34_rows'];assert len(selected)==34 and max(rowweights[i]for i in selected)<=4
 frac=lambda x:K.from_expr(s.sympify(x));T=[[frac(row[j])for j in chart['normal34']]for row in chart['T34']];U=[[frac(c)for c in row]for row in chart['U34']];assert len(T)==len(U)==45
 for i in [19,23]:assert not any(T[i])and not any(U[i])
 for i in range(45):
  for j,c in enumerate(T[i]):
   if c:assert localweights[i]==o.variable_weights[j]
  for j,c in enumerate(U[i]):
   if c:assert localweights[i]==rowweights[selected[j]]
 beta={int(k[1:])-1:s.sympify(v)for k,v in chart['beta'].items()};_,F=specialize(S61(),beta,True);assert len(F)==135
 fixed,cols=old['fixed'],old['cols'];J=[{}for _ in F];B=[[]for _ in F]
 for i,poly in enumerate(F):
  for mon,value in poly.items():
   cc=K(R.from_expr(s.sympify(value)))
   if len(mon)==2:
    x,y=mon
    assert not(x==fixed and y==fixed)
    if x!=fixed and y==fixed:j=cols.index(x);J[i][j]=J[i].get(j,K.zero)+cc
    if y!=fixed and x==fixed:j=cols.index(y);J[i][j]=J[i].get(j,K.zero)+cc
   elif len(mon)==1 and mon[0]==fixed:J[i][44]=J[i].get(44,K.zero)+cc
   if fixed not in mon:
    ij=[cols.index(j)for j in mon]+[44]*(2-len(mon));assert len(ij)==2;B[i].append((ij[0],ij[1],cc))
  J[i]={j:c for j,c in J[i].items()if c}
 for i,row in enumerate(J):
  for k in range(8):assert sum((c*T[j][k]for j,c in row.items()),K.zero)==0
 P=[]
 for i,row in enumerate(J):P.append([sum((c*U[j][k]for j,c in row.items()),K.zero)for k in range(34)])
 for r,i in enumerate(selected):assert P[i]==[K.one if r==k else K.zero for k in range(34)]
 variables=[o.var(i)for i in range(8)];E=[o.lin(list(zip(row,variables)))for row in T]
 def bilinear(left,right):
  result=[]
  for i,row in enumerate(B):
   terms=[]
   if rowweights[i]<=4:
    for j,k,c in row:
     if j==k:terms.append((c,o.mul(left[j],right[k])))
     else:terms.extend([(c/2,o.mul(left[j],right[k])),(c/2,o.mul(left[k],right[j]))])
   result.append(o.lin(terms))
  return result
 def lift(w):return[o.lin(list(zip(row,w)))for row in U]
 EE=bilinear(E,E);w2=[o.lin([(-1,EE[i])])for i in selected];Z2=lift(w2)
 EZ2=bilinear(E,Z2);w3=[o.lin([(-2,EZ2[i])])for i in selected];Z3=lift(w3)
 Z2Z2=bilinear(Z2,Z2);EZ3=bilinear(E,Z3);b4=[o.lin([(1,x),(2,y)])for x,y in zip(Z2Z2,EZ3)];w4=[o.lin([(-1,b4[i])])for i in selected];Z4=lift(w4)
 for degree,local in [(1,E),(2,Z2),(3,Z3),(4,Z4)]:
  for j,node in enumerate(local):
   if node:assert o.weights[node]==localweights[j]and o.degrees[node]==degree and degree<=localweights[j]
 def project(b):return[o.lin([(1,b[i])]+[(-c,b[selected[r]])for r,c in enumerate(P[i])if c])for i in range(135)]
 f2=project(EE);f3=project([o.lin([(2,n)])for n in EZ2]);f4=project(b4);roots=[o.lin([(1,x),(1,y),(1,z)])for x,y,z in zip(f2,f3,f4)]
 assert all(not roots[i]for i in selected)
 for i,node in enumerate(roots):
  if node:assert o.weights[node]==rowweights[i]<=4
 genuine=[i for i,w in enumerate(rowweights)if w==4 and i not in selected];assert len(genuine)==33
 for i,row in enumerate(J):
  for j,c in row.items():o.atom(c)
 for row in B:
  for j,k,c in row:o.atom(c)
 for row in T+U+P:
  for c in row:o.atom(c)
 evidence={'source_original_equations':135,'fixed_original_coordinate':fixed,'other_original_coordinate_order':cols,'selected34_rows':selected,'local_coordinate_weights':localweights,'equation_weights':rowweights,'normal_variable_names':['v0','v1','v2','v3','v4','v5','v6','v7'],'T_atom_ids':[[o.atom(c)for c in row]for row in T],'U_atom_ids':[[o.atom(c)for c in row]for row in U],'J_sparse_atom_ids':[[[j,o.atom(c)]for j,c in row.items()]for row in J],'B_sparse_atom_ids':[[[j,k,o.atom(c)]for j,k,c in row]for row in B],'P_JU_atom_ids':[[o.atom(c)for c in row]for row in P],'linear_graph_roots':E,'implicit_w2_roots':w2,'implicit_w3_roots':w3,'implicit_w4_roots':w4,'local_Z2_roots':Z2,'local_Z3_roots':Z3,'local_Z4_roots':Z4,'residual_by_ordinary_degree_roots':[f2,f3,f4],'residual_weight_le4_roots':roots,'genuine_quartic_original_rows':genuine,'genuine_quartic_roots':[roots[i]for i in genuine],'rows_above_weight4_not_asserted_to_vanish':[i for i,w in enumerate(rowweights)if w>4],'normal_graph_formula':'x=T*v+U*(w2+w3+w4); w2=-Bsel(E,E),w3=-2Bsel(E,U*w2),w4=-Bsel(U*w2,U*w2)-2Bsel(E,U*w3)','arithmetic':o.encoded()}
 text=json.dumps(evidence,separators=(',',':'))+'\n';assert len(text)<8*1024*1024;(output/'oracle.json').write_text(text)
 fixture_count=fixture();o2,data=Oracle.load(output/'oracle.json');assert o2.encoded()==evidence['arithmetic']
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),chartpath,oldpath,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py']
 report=dict(id=ID,all_exact_source_and_graph_checks_passed=True,original_equation_count=135,all_JT_zero_and_selected_JU_identity_exact=True,all_core_atom_denominators_only_a_F_B_C=True,all_nuisance_zero_weight_coordinates_identically_zero=True,implicit_selected_weight_le4_roots_identically_zero=True,all33_genuine_quartic_roots_present=True,positive_weights_make_degree4_graph_exact=True,no_higher_source_row_vanishing_claim=True,lazy_coefficient_API_fixture_monomials=fixture_count,oracle_roundtrip_exact=True,arithmetic_nodes=len(o.nodes),core_atoms=len(o.atoms),proof_packet_bytes=len(text),no_expanded_parametric_quartic_matrix=True,no_new_actual_point_or_field_or_GB=True,no_rank64_full_quartic_or_ambient_claim=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},oracle_sha256=sha(output/'oracle.json'),seconds=time.monotonic()-start,frozen=True,independent_review_pending=True)
 (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_exact_source_and_graph_checks_passed','arithmetic_nodes','core_atoms','proof_packet_bytes','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);build(p.parse_args().output.resolve())
