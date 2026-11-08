#!/usr/bin/env python3
"""Separate Fraction full-source finite-shear/factor review, no producer/CAS.

Only the frozen q4 arithmetic primitives are reused. Both existing complete K1
spaces, all six actual finite shears, all nine P2 variables, every original
Jacobi row and the signed Euler-Pfaffian are retained.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
import ast,hashlib,json,resource,signal,time
import classification_degenerate_quadric_q4_nil31_I_review as q
ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/plane_reduced_m3_whole_source_separate_review.json'
OUT='reports/CLASSIFICATION_PLANE_REDUCED_M3_WHOLE_SOURCE_SEPARATE_REVIEW_20261108.json'
ORIGINAL='reports/CLASSIFICATION_PLANE_REDUCED_M3_NINE_FULLHIGHER_20261108.json'
PRODUCER='reports/CLASSIFICATION_PLANE_REDUCED_M3_WHOLE_SHEAR_SOURCE_20261108.json'
PAIRS=list(combinations(range(5),2));TRIPLES=list(combinations(range(5),3));QUADS=list(combinations(range(5),4))
def atom(n,c=1,i=None):
 e=[0]*n
 if i is not None:e[i]=1
 return {tuple(e):F(c)}if c else{}
def mono(n,ex,c=1):return{tuple(ex)+(0,)*(n-len(ex)):F(c)}if c else{}
def parse(text,n,names):
 def walk(node):
  if isinstance(node,ast.Constant)and type(node.value)is int:return atom(n,node.value)
  if isinstance(node,ast.Name):return atom(n,1,names.index(node.id))
  if isinstance(node,ast.UnaryOp):return q.scale(walk(node.operand),-1 if isinstance(node.op,ast.USub)else 1)
  if isinstance(node,ast.BinOp):
   A=walk(node.left)
   if isinstance(node.op,ast.Pow):
    assert isinstance(node.right,ast.Constant)and type(node.right.value)is int and node.right.value>=0
    B=atom(n)
    for _ in range(node.right.value):B=q.mul(B,A)
    return B
   B=walk(node.right)
   if isinstance(node.op,ast.Add):return q.add(A,B)
   if isinstance(node.op,ast.Sub):return q.add(A,q.scale(B,-1))
   if isinstance(node.op,ast.Mult):return q.mul(A,B)
   if isinstance(node.op,ast.Div):assert set(B)=={(0,)*n};return q.scale(A,1/B[(0,)*n])
  raise ValueError(ast.dump(node))
 return walk(ast.parse(str(text),mode='eval').body)
def mons(d):
 out=[]
 for entries in combinations_with_replacement(range(5),d):
  ex=[0]*5
  for i in entries:ex[i]+=1
  out.append(tuple(ex))
 return out
COLS=[(ij,e)for ij in PAIRS for e in mons(2)];ROWS=[(ij,e)for ij in QUADS for e in mons(4)]
def skew(P):return q.antisymmetric(P)
def plus(*Ps):return skew({ij:q.add(*(P.get(ij,{})for P in Ps))for ij in PAIRS})
def wedge(A,B):return skew({(i,j):q.add(q.mul(A[i],B[j]),q.scale(q.mul(A[j],B[i]),-1))for i,j in PAIRS})
def fromvector(vector,n):return skew({ij:{e+(0,)*(n-5):F(c)for (kl,e),c in zip(COLS,vector)if ij==kl and c}for ij in PAIRS})
def vector(P,n):return[P.get(ij,{}).get(e+(0,)*(n-5),F(0))for ij,e in COLS]
def polyvector(P,n):return[{ex[5:]:c for ex,c in P.get(ij,{}).items()if ex[:5]==e}for ij,e in COLS]
def lift(P,n):return{ij:{e+(0,)*(n-len(e)):c for e,c in p.items()}for ij,p in P.items()}
def scaleP(P,p):return{ij:q.mul(pol,p)for ij,pol in P.items()}
def mulmat(A,B):return[[sum((a*b for a,b in zip(row,col)),F())for col in zip(*B)]for row in A]
def eye(n):return[[F(i==j)for j in range(n)]for i in range(n)]
def mat(rows):return[[F(v)for v in row]for row in rows]
def parseP(saved,n,names):return skew({ast.literal_eval(ij):parse(v,n,names)for ij,v in saved.items()})
def compose(P,images,n):
 out={}
 for ex,c in P.items():
  term=atom(n,c)
  for i,k in enumerate(ex):
   for _ in range(k):term=q.mul(term,images[i])
  out=q.add(out,term)
 return out
def push(P,g,n):
 x=[atom(n,1,i)for i in range(n)];images=x[:]
 for i in range(3):images[i]=q.add(x[i],q.scale(q.mul(g[2*i],x[3]),-1),q.scale(q.mul(g[2*i+1],x[4]),-1))
 M=[[atom(n,int(i==j))for j in range(5)]for i in range(5)]
 for i in range(3):M[i][3],M[i][4]=g[2*i:2*i+2]
 Psub={ij:compose(pol,images,n)for ij,pol in P.items()}
 return skew({(i,j):q.add(*(q.mul(q.mul(M[i][k],Psub.get((k,l),{})),M[j][l])for k in range(5)for l in range(5)if M[i][k]and M[j][l]))for i,j in PAIRS})
def jac(P,n):
 strict={ijk:q.add(*(q.mul(P.get((i,l),{}),q.der(P.get((j,k),{}),l))for i,j,k in[ijk,(ijk[1],ijk[2],ijk[0]),(ijk[2],ijk[0],ijk[1])]for l in range(5)))for ijk in TRIPLES}
 return{quad:q.add(*(q.scale(q.mul(atom(n,1,quad[k]),strict[quad[:k]+quad[k+1:]]),(-1)**k)for k in range(4)))for quad in QUADS}
def jvec(J):return[{ex[5:]:c for ex,c in J[quad].items()if ex[:5]==e}for quad,e in ROWS]
def actual(P,n):
 out={}
 for k in range(5):
  i,j,l,m=[z for z in range(5)if z!=k];pf=q.add(q.mul(P.get((i,j),{}),P.get((l,m),{})),q.scale(q.mul(P.get((i,l),{}),P.get((j,m),{})),-1),q.mul(P.get((i,m),{}),P.get((j,l),{})))
  out=q.add(out,q.scale(q.mul(atom(n,1,k),pf),(-1)**k))
 return out

def source_bases():
 B1=[];n=5
 for normal in(3,4):
  for ij in combinations(range(3),2):
   for i in range(3):ex=[0]*5;ex[normal]+=1;ex[i]+=1;B1.append(skew({ij:mono(n,ex)}))
 for normal,indices in[(3,(3,4)),(3,(4,4)),(4,(3,3)),(4,(3,4))]:
  for i in range(3):ex=[0]*5;ex[indices[0]]+=1;ex[indices[1]]+=1;B1.append(skew({(i,normal):mono(n,ex)}))
 B2=[]
 for ij in combinations(range(3),2):
  for indices in[(3,3),(3,4),(4,4)]:ex=[0]*5;ex[indices[0]]+=1;ex[indices[1]]+=1;B2.append(skew({ij:mono(n,ex)}))
 x=[atom(n,1,i)for i in range(5)];G=[]
 for normal in(3,4):
  for i in range(3):V=[{}for _ in range(5)];V[i]=x[normal];G.append(wedge(x,V))
 raw=[i for i,(ij,e)in enumerate(COLS)if e[3]+e[4]-sum(t>=3 for t in ij)==1];assert len(raw)==36
 E=[[vector(P,n)[i]for P in B1+G]for i in raw];EI=q.inverse(E);assert mulmat(E,EI)==eye(36)
 return B1,B2,G,raw,EI

def encpoly(P):return{','.join(map(str,e)):str(c)for e,c in P.items()}
def main():
 start=time.monotonic();signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(3072*1024**2,)*2);resource.setrlimit(resource.RLIMIT_CPU,(175,180));dest=ROOT/OUT;assert not dest.exists();job=json.loads((ROOT/JOB).read_text())
 for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
 original=json.loads((ROOT/ORIGINAL).read_text());producer=json.loads((ROOT/PRODUCER).read_text());B1,B2,G,raw,EI=source_bases();out=[]
 for saved in producer['cases']:
  name=saved['name'];old=next(z for z in original['cases']if z['name']==name);K=q.transpose(mat(old['constrained_weight1_certificate']['complete_kernel']));dim=len(K[0]);assert dim in(7,8);free=saved['kernel_free_coordinate_indices'];assert [[K[i][j]for j in range(dim)]for i in free]==eye(dim);C=mat(saved['global_first_source_matrix']);CI=mat(saved['global_first_source_inverse']);assert mulmat(C,CI)==mulmat(CI,C)==eye(dim);d=q.determinant(C);assert d==F(saved['global_constant_determinant'])and d
  residual=saved['residual_kernel_indices'];nr=len(residual);n=5+6+nr+9;names=['a','b','c','u','v']+['s'+str(i)for i in range(6)]+['r'+str(i)for i in range(nr)]+['q'+str(i)for i in range(9)];g=[atom(n,1,5+i)for i in range(6)];rpars=[atom(n,1,11+i)for i in range(nr)];qpars=[atom(n,1,11+nr+i)for i in range(9)];P0=fromvector(old['full150_P0'],n);P1=plus(*(scaleP(lift(B1[i],n),q.add(*(q.scale(rpars[j],K[i][col])for j,col in enumerate(residual))))for i in range(30)));assert P1==parseP(saved['normalized_P1'],n,names)
  Pnorm=plus(P0,P1,*(scaleP(lift(B,n),param)for B,param in zip(B2,qpars)));Ppush=push(Pnorm,g,n);Pback=push(Ppush,[q.scale(p,-1)for p in g],n);assert Pback==Pnorm
  # Project all weight1 coefficients by the independent raw36/Euler6 inverse.
  pv=polyvector(Ppush,n);w1=[pv[i]for i in raw];coeff=[{(0,)*5+e:c for e,c in q.add(*(q.scale(pol,z)for pol,z in zip(w1,row))).items()}for row in EI];change=[g[i]for i in range(6)]+rpars;alpha=[q.add(*(q.scale(p,c)for p,c in zip(change,row)))for row in C];expected=[q.add(*(q.scale(p,c)for p,c in zip(alpha,row)))for row in K];assert coeff[:30]==expected
  canonical=plus(Ppush,*(scaleP(lift(B,n),q.scale(co,-1))for B,co in zip(G,coeff[30:])));cv=polyvector(canonical,n)
  translations=[]
  for B,qp in zip(B2,qpars):
   label=next(i for i,v in enumerate(vector(B,5))if v);pol={ (0,)*5+e:c for e,c in cv[label].items()};diff=q.add(pol,q.scale(qp,-1));assert all(not any(e[11+nr:])for e in diff);translations.append(encpoly(diff))
  assert all((e[3]+e[4]-sum(t>=3 for t in ij))in(0,1,2)for ij in PAIRS for e in canonical[ij]);assert actual(Ppush,n)==compose(actual(Pnorm,n),[q.add(atom(n,1,i),q.scale(q.mul(g[2*i],atom(n,1,3)),-1),q.scale(q.mul(g[2*i+1],atom(n,1,4)),-1))if i<3 else atom(n,1,i)for i in range(n)],n);assert actual(canonical,n)==actual(Ppush,n)
  # Original full normalized equations, before any P2 variables are eliminated.
  JV=jvec(jac(Pnorm,n));A=mat(old['complete350x9']);rr=saved['P2_pivot_rows'];cc=saved['P2_pivot_columns'];minor=[[A[i][j]for j in cc]for i in rr];assert q.determinant(minor)==F(saved['P2_minor']);II=q.inverse(minor)
  for i in rr:assert JV[i]==q.add(*(q.scale({e[5:]:c for e,c in qp.items()},z)for qp,z in zip(qpars,A[i])))
  sol={int(key[1:]):parse(val,n,names)for key,val in saved['P2_solution'].items()};assert set(sol)==set(cc)and all(not p for p in sol.values());free2=[j for j in range(9)if j not in cc];assert [[A[i][j]for j in free2]for i in rr]==[[F(0)for j in free2]for i in rr]
  Pfinal=plus(P0,P1,*(scaleP(lift(B,n),sol.get(i,qpars[i]))for i,B in enumerate(B2)));assert Pfinal==parseP(saved['normalized_complete_source'],n,names);JF=jvec(jac(Pfinal,n));assert JF==[{e[5:]:c for e,c in parse(poly,n,names).items()}for poly in saved['all350_after_P2_substitution']];D=actual(Pfinal,n);assert D==parse(saved['actual_divisor'],n,names)==parse(saved['actual_divisor_factorization'],n,names)
  if name=='M':factor='r0*(b*u+c*v)*(3*a*u**2+6*b*u*v+3*c*v**2+r0*v**3)/2';assert all(not p for p in JF)and not free2
  else:factor='-c*(2*a*r1*u**3-2*b*r0*v**3+(q1+r0*r1)*u**2*v**2)';assert free2==[1];nonzero={str(c)for pol in JF for c in pol.values()};assert nonzero=={'-2','2','4'};assert all(not pol or set(pol)=={tuple([0]*6+[1,1]+[0]*9)}for pol in JF)
  assert D==parse(factor,n,names)
  out.append({'name':name,'all_complete_K1_directions_retained':dim,'all_six_finite_shears_and_global_inverse_checked':True,'raw36_Euler6_projection_and_constant_fullsource_inverse':True,'full_source_change_matrix':saved['global_first_source_matrix'],'full_source_inverse':saved['global_first_source_inverse'],'constant_determinant':str(d),'whole_actual_P1_P2_pushforward_and_inverse':True,'all9_P2_triangular_translations':translations,'P2_linear_identity_on_all9_original_variables':True,'all_parameter_Euler_Pfaffian_push_and_regauge_checked':True,'all350_original_normalized_equations_after_fullP2_elimination_checked':True,'all9_P2_pivot_minor':str(q.determinant(minor)),'P2_free_variables':[str(i)for i in free2],'remaining_exact_Poisson_compatibility':'none'if name=='M'else'r0*r1=0','whole_actual_divisor_factor':factor,'irreducible_quintic_whole_source_empty':True})
 result={'task_id':'E_reduced_nonnormal','all_independent_checks_passed':True,'no_producer_or_CAS_imports':True,'alternate_Fraction_engine':'classification_degenerate_quadric_q4_nil31_I_review.py','cases':out,'pure_scope':'Both complete prior constrained K1 sources globally equal finite-shear translates of displayed full9P2 normalized sources; original Jacobi and full signed Pfaffian imply proper factors on all support. Prior root-cover completeness read separately, no extension to concurrent/nonreduced/otherm.','inputs':job['exact_inputs_and_hashes'],'reviewer_source_bytes':Path(__file__).read_text(),'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 signal.alarm(0);dest.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in('all_independent_checks_passed','seconds','maximum_RSS_KiB')}))
if __name__=='__main__':main()
