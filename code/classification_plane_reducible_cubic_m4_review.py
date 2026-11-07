#!/usr/bin/env python3
"""Separate-author Fraction review: all24 flat connections and m4 scalar bounds.

Imports only previous independent standard-library Fraction arithmetic, never
the producer or its helpers. The concurrent upper31 is explicitly NOT closed.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,itertools as it,json,resource,time
from classification_projected_veronese_zero_review import Ring,plus,times,scale,derivative,monomials,sparse_rref,det
from classification_plane_rank_zero_fullrank_review import PAIRS,QUADS,anti,wedge,tensor_add,strict,projective
ROOT=Path(__file__).resolve().parents[1]
def transpose(M):return [list(row)for row in zip(*M)]
def decode(p):return {tuple(e):F(c)for e,c in p}
def det3(M):
 return plus(*(scale(times(*(M[i][j]for i,j in enumerate(perm))),(-1)**sum(perm[i]>perm[j]for i in range(3)for j in range(i+1,3)))for perm in it.permutations(range(3))))
def matbasis():
 out=[[[F(a==i and b==j)for b in range(3)]for a in range(3)]for i,j in it.product(range(3),repeat=2)if i!=j]
 out += [[[F(1),F(0),F(0)],[F(0),F(-1),F(0)],[F(0),F(0),F(0)]],[[F(0),F(0),F(0)],[F(0),F(1),F(0)],[F(0),F(0),F(-1)]]]
 return out
def numerator(E,V):return [plus(times(E[2],V[1]),scale(times(E[1],V[2]),-1)),plus(times(E[0],V[2]),scale(times(E[2],V[0]),-1)),plus(times(E[1],V[0]),scale(times(E[0],V[1]),-1))]
def divides(poly,g):
 rem=dict(poly);lead=max(g)
 while rem:
  mon=max(rem)
  if any(a<b for a,b in zip(mon,lead)):return False
  term={tuple(a-b for a,b in zip(mon,lead)):rem[mon]/g[lead]};rem=plus(rem,scale(times(term,g),-1))
 return True
def flat_connection_review(saved):
 R=Ring([f'p{i}'for i in range(24)]+['x','y','z']);p=R.parse;E=[R.var(n)for n in ['x','y','z']];B=matbasis();V=[[plus(*(scale(times(R.var(f'p{8*j+k}'),E[l]),B[k][i][l])for k in range(8)for l in range(3)))for i in range(3)]for j in range(3)];N=[numerator(E,v)for v in V];mons=monomials(3,3);records=[]
 def curv(f,W):
  out=[]
  for i,j in it.combinations(range(3),2):
   H=plus(times(W[1][i],W[2][j]),scale(times(W[2][i],W[1][j]),-1));A=scale(plus(times(W[0][i],W[1][j]),scale(times(W[1][i],W[0][j]),-1)),2);Bv=scale(plus(times(W[0][i],W[2][j]),scale(times(W[2][i],W[0][j]),-1)),-2)
   out.append([plus(times(f,plus(derivative(W[k][j],24+i),scale(derivative(W[k][i],24+j),-1))),scale(times(derivative(f,24+i),W[k][j]),-1),times(derivative(f,24+j),W[k][i]),[H,A,Bv][k])for k in range(3)])
  return out
 for case in saved['all_four_types']:
  f=p(case['cubic']);value=lambda v:plus(*(times(v[i],derivative(f,24+i))for i in range(3)));J=[plus(value(V[0]),det3([E,V[1],V[2]])),plus(value(V[1]),scale(det3([E,V[0],V[1]]),2)),plus(value(V[2]),scale(det3([E,V[0],V[2]]),-2))];co=[]
  for expr in J:
   for mon in mons:co.append({e[:24]+(0,0,0):c for e,c in expr.items()if e[24:]==mon})
  assert co==[decode(v)for v in case['full30_flatness_polynomials']]
  raw=curv(f,N)
  for at,(i,j)in enumerate(it.combinations(range(3),2)):
   k=next(k for k in range(3)if k not in[i,j]);assert raw[at]==[scale(times(E[k],v),(-1)**k)for v in J]
  T=[[decode(v)for v in row]for row in case['nonlog_flat_witness_matrix']];Q=plus(times(T[0][0],T[1][1]),scale(times(T[0][1],T[1][0]),-1));L=p('z'if case['type']!='concurrent'else'x-y');assert times(Q,L)==f;adj=[[T[1][1],scale(T[0][1],-1)],[scale(T[1][0],-1),T[0][0]]];W=[[],[],[]]
  for k in range(3):
   A=[[plus(*(times(adj[i][l],derivative(T[l][j],24+k))for l in range(2)))for j in range(2)]for i in range(2)];W[0].append(scale(times(L,plus(A[0][0],scale(A[1][1],-1))),F(1,2)));W[1].append(times(L,A[0][1]));W[2].append(times(L,A[1][0]))
  assert W==[[decode(v)for v in row]for row in case['nonlog_flat_witness_numerators']];assert all(not v for row in curv(f,W)for v in row);assert all(not plus(*(times(E[i],v[i])for i in range(3)))for v in W)
  nonlog=[]
  for factor in case['components']:
   g=p(factor)
   if any(not divides(plus(times(derivative(g,24+i),form[j]),scale(times(derivative(g,24+j),form[i]),-1)),g)for form in W for i,j in it.combinations(range(3),2)):nonlog.append(factor)
  assert nonlog==case['nonlog_components']and nonlog;records.append(dict(type=case['type'],all24_parameters_and30_flatness_equations_exact=True,all_nine_raw_curvature_components_exact=True,genuine_nonlog_flat_witness_exact=True,nonlog_components=nonlog))
 return records
def matrix_cert(M):
 basis,_,rows=sparse_rref(M);cols=sorted(basis);value=det([[M[i][j]for j in cols]for i in rows]);assert value;return dict(rank=len(cols),rows=rows,columns=cols,minor_det=str(value))
def kernel(M,n):
 basis,_,_=sparse_rref(M);free=[j for j in range(n)if j not in basis];out=[]
 for j in free:
  v=[F(0)]*n;v[j]=F(1)
  for p,row in basis.items():v[p]=-row.get(j,F(0))
  assert all(sum(a*b for a,b in zip(row,v))==0 for row in M);out.append(v)
 return out
def scalar_review(saved):
 R=Ring(('x','y','z','u','v'));p=R.parse;X=[R.var(n)for n in R.names];base=X[:3];mons=monomials(5,4);bmons=monomials(3,2);B=matbasis();fields=[[plus(*(scale(base[j],M[i][j])for j in range(3)))for i in range(3)]for M in B];nums=[numerator(base,V)for V in fields];A=[[n[i].get(e+(0,0),F(0))for n in nums]for i in range(3)for e in bmons];records=[]
 first=[]
 for normal in [X[3],X[4]]:
  for ij in it.combinations(range(3),2):
   for a in base:first.append(anti({ij:times(normal,a)}))
 for normalmon,j in[(times(X[3],X[4]),3),(times(X[4],X[4]),3),(times(X[3],X[3]),4),(times(X[3],X[4]),4)]:
  for i in range(3):first.append(anti({(i,j):normalmon}))
 second=[anti({ij:normal})for ij in it.combinations(range(3),2)for normal in[times(X[3],X[3]),times(X[3],X[4]),times(X[4],X[4])]];assert len(first)==30 and len(second)==9
 def jvec(P):
  J=projective(strict(P),R);return[J[q].get(e,F(0))for q in QUADS for e in mons]
 def actual(P):
  pf=[]
  for i in range(5):j,k,l,m=[r for r in range(5)if r!=i];pf.append(scale(plus(times(P[j,k],P[l,m]),scale(times(P[j,l],P[k,m]),-1),times(P[j,m],P[k,l])),(-1)**i))
  return plus(*(times(X[i],pf[i])for i in range(5)))
 for case in saved['cases']:
  f=p(case['cubic']);L=p(case['line_component']);Q={}
  # All cubics here have the displayed line factor; exact division is linear.
  possibilities={'conic_transverse':'x*y-z*z','conic_tangent':'x*z-y*y','triangle':'x*y','concurrent':'x*y'};Q=p(possibilities[case['type']]);assert times(Q,L)==f;target=[plus(scale(derivative(f,i),F(1,6)),scale(times(Q,derivative(L,i)),-F(1,2)))for i in range(3)];rhs=[v.get(e+(0,0),F(0))for v in target for e in bmons];K=kernel([row+[-v]for row,v in zip(A,rhs)],9);assert len(K)==1 and K[0][-1];coeff=[3*a/K[0][-1]for a in K[0][:-1]];M=[[sum(coeff[k]*B[k][i][j]for k in range(8))for j in range(3)]for i in range(3)];assert M==[[F(v)for v in row]for row in case['scalar_normal_matrix']];V=[plus(*(scale(base[j],M[i][j])for j in range(3)))for i in range(3)]+[{},{}];En=[{},{},{},X[3],X[4]];Pf=anti({(0,1):derivative(f,2),(0,2):scale(derivative(f,1),-1),(1,2):derivative(f,0)});P0=tensor_add(Pf,wedge(En,V));assert not any(jvec(P0))and not actual(P0)
  J1=transpose([jvec(tensor_add(P0,P))for P in first]);cross=[actual(tensor_add(P0,P))for P in first];d3=[e for e in monomials(5,5)if e[3]+e[4]==3];assert len(d3)==24;D3=[[v.get(e,F(0))for v in cross]for e in d3];necessary=J1+D3;K1=kernel(necessary,30);assert J1==[[F(v)for v in row]for row in case['full350_weight1_matrix']]and D3==[[F(v)for v in row]for row in case['all24_D_normal_degree3_linear_rows']];assert K1==[[F(v)for v in row]for row in case['full_weight1_and_order4_necessary_kernel']]
  J2=transpose([jvec(tensor_add(P0,P))for P in second]);assert J2==[[F(v)for v in row]for row in case['full350_weight2_matrix']];c1=matrix_cert(necessary);c2=matrix_cert(J2);d=len(K1)+9-c2['rank'];aut={'conic_transverse':1,'conic_tangent':2,'triangle':2,'concurrent':3}[case['type']];image=24+d-aut-10;assert image==case['coarse_whole_projective_image_dimension_upper'];assert len(K1)==case['weight1_kernel_dimension']and c2['rank']==case['weight2_linear_rank']
  normal4=[]
  for terms in it.combinations_with_replacement([3,4],4):
   ex=[0]*5
   for i in terms:ex[i]+=1
   normal4.append(ex)
  def shape(poly):
   out=[]
   for ex in normal4:
    values=[]
    for i in range(3):e=list(ex);e[i]+=1;values.append(poly.get(tuple(e),F(0)))
    out.extend([values[0]+values[1],values[2]]if case['type']=='concurrent'else values[:2])
   return out
  L2=J2+transpose([shape(actual(tensor_add(P0,P)))for P in second]);assert L2==[[F(v)for v in row]for row in case['full350_plus10_leading_shape_weight2_linear_matrix']];cshape=matrix_cert(L2);assert cshape['rank']==case['shape_weight2_rank']and cshape['rank']<9;assert case['exact_proper_quadratic_compatibility']is None
  records.append(dict(type=case['type'],P0_scalar_connection_and_Jacobi_D_zero_exact=True,all350_by30_and24_order_rows_exact=True,complete_kernel30_exact=True,all350_by9_and10_leading_shape_rows_exact=True,weight1_kernel_dimension=len(K1),weight2_rank=c2['rank'],weight2_affine_fiber=9-c2['rank'],full1_rank_certificate=c1,full2_rank_certificate=c2,shape_rank_certificate=cshape,connected_stabilizer_dimension=aut+10,whole_projective_image_upper_bound=image,G24_accepted=image<=24,quadratic_P1_compatibility_still_unresolved=True))
 return records
def main():
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(119,120));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);flatpath=ROOT/'reports/CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_CONNECTION_20261107.json';scalarpath=ROOT/'reports/CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SCALAR_V2_20261107.json';flat=json.loads(flatpath.read_text());saved=json.loads(scalarpath.read_text());connections=flat_connection_review(flat);records=scalar_review(saved);assert[case['whole_projective_image_upper_bound']for case in records]==[24,23,23,31]
 paths=[Path(__file__),flatpath,scalarpath,ROOT/'code/classification_projected_veronese_zero_review.py',ROOT/'code/classification_plane_rank_zero_fullrank_review.py',ROOT/'code/classification_unimodular_logarithmic_review.py',ROOT/'research/classification/PLANE_REDUCIBLE_CUBIC_CONNECTION_ENTRANCE.md',ROOT/'research/classification/PLANE_REDUCIBLE_CUBIC_M4_SCALAR_BOUND.md'];result=dict(task_id='G',all_independent_checks_passed=True,whole_connection_sources=connections,scalar_m4_sources=records,accepted_whole_bounds=[24,23,23],concurrent_bound31_not_closed=True,geometry_scope='m4 squarefree binary leading, both I,J nonzero; parallel invariant/UFD/regularity/Omega1(1) and group-fiber bridge reviewed separately in note',complements=['concurrent m4','special I0 or J0','non-squarefree leading','m2 and m3','nonreduced cubic R','unresolved P1 quadratic compatibility'],method='Separate-author standard-library Fraction, no producer or producer-helper/CAS import.',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},checker_bytes=Path(__file__).read_text(),seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);out=ROOT/'reports/CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SCALAR_REVIEW_20261107.json';assert not out.exists();out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,accepted_bounds=[24,23,23],concurrent_open=True,seconds=result['seconds'],maximum_RSS_KiB=result['maximum_RSS_KiB'])))
if __name__=='__main__':main()
