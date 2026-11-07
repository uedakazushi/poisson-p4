#!/usr/bin/env python3
"""Independent full quadric-cone source/Ore/rank certificate review.

Only standard-library Fraction sparse arithmetic helpers are reused. Neither
producer nor a CAS is imported. Generic kernels are not needed for the bound.
"""
import argparse,hashlib,itertools as it,json,resource,time
from pathlib import Path
from classification_cubic_scroll_abelian_review import Q,add,scale,mul,diff,parse,monomials,rref,determinant

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2));TRIPLES=list(it.combinations(range(5),3));QUADS=list(it.combinations(range(5),4));P3=list(it.combinations(range(3),2))
def var(n,i):return {tuple(int(i==j)for j in range(n)):Q(1)}
def constant(n,c):return {(0,)*n:Q(c)}if c else{}
def zero(d):return [[{}for j in range(d)]for i in range(d)]
def setentry(P,i,j,p):P[i][j]=p;P[j][i]=scale(p,-1)
def plusmatrix(A,B):return [[add(a,b)for a,b in zip(row,col)]for row,col in zip(A,B)]
def scalematrix(A,c):return [[scale(p,c)for p in row]for row in A]
def strict(P,offset=0):
 return {tri:add(*(mul(P[i][l],diff(P[j][k],offset+l))for i,j,k in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(len(P))))for tri in it.combinations(range(len(P)),3)}
def projective(P,offset=0):
 n=len(next(iter(next(p for row in P for p in row if p))));J=strict(P,offset)
 return {quad:add(*(scale(mul(var(n,offset+quad[k]),J[quad[:k]+quad[k+1:]]),(-1)**k)for k in range(4)))for quad in QUADS}
def rank(A):return len(rref(A))
def matmul(A,B):return [[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]
def transpose(A):return list(map(list,zip(*A)))
def dot(row,v):return sum(c*d for c,d in zip(row,v))
def restrict(poly,index,value=1):
 ans={}
 for ex,c in poly.items():
  ee=list(ex);c*=value**ee[index];ee[index]=0;ee=tuple(ee);ans[ee]=ans.get(ee,Q(0))+c
 return {ex:c for ex,c in ans.items()if c}
def substitute(poly,values):
 n=len(next(iter(poly)))if poly else len(next(iter(next(iter(values.values())))))
 out={}
 for ex,c in poly.items():
  p=constant(n,c)
  for j,k in enumerate(ex):
   value=values.get(j,var(n,j))
   for _ in range(k):p=mul(p,value)
  out=add(out,p)
 return out
def coefficient_space(polys,indices):
 result={}
 for ex,c in polys.items():
  key=tuple(ex[i]for i in indices);ee=list(ex)
  for i in indices:ee[i]=0
  result[key]=add(result.get(key,{}),{tuple(ee):c})
 return result
def rowspace_polys(polys):
 keys=sorted(set().union(*(p.keys()for p in polys)));R=rref([[p.get(k,Q(0))for k in keys]for p in polys])
 return {tuple((keys[j],c)for j,c in sorted(row.items()))for row in R.values()}

def review_source(saved):
 xs=[var(5,i)for i in range(5)];q=add(mul(xs[0],xs[2]),scale(mul(xs[1],xs[1]),-1));w=xs[4];mons2=monomials(5,2);mons3=monomials(5,3);mons4=monomials(5,4)
 def vec(P):return [P[i][j].get(ex,Q(0))for i,j in PAIRS for ex in mons2]
 raw=[]
 for i,j in PAIRS:
  for k in range(5):P=zero(5);setentry(P,i,j,mul(w,xs[k]));raw.append(P)
 for i,j in PAIRS:P=zero(5);setentry(P,i,j,q);raw.append(P)
 source=[vec(P)for P in raw];assert rank(source)==60
 vertical=[]
 for a,b in it.product(range(5),repeat=2):
  P=zero(5)
  for i,j in PAIRS:setentry(P,i,j,add(mul(xs[i],xs[b])if j==a else{},scale(mul(xs[j],xs[b]),-1)if i==a else{}))
  vertical.append(vec(P))
 assert rank(vertical)==24 and rank(vertical+source)==79
 E=[[Q((j==a and k==i)-(i==a and k==j))for i,j in PAIRS for k in range(5)]for a in range(5)];assert rank(E)==5
 membership=[add(*(mul(P[4][i],diff(q,i))for i in range(5)))for P in raw]
 ideals=[mul(mul(w,w),v)for v in xs]+[mul(w,q)]
 C=[[p.get(ex,Q(0))for p in membership]+[-p.get(ex,Q(0))for p in ideals]for ex in mons3]
 assert C==[[Q(v)for v in row]for row in saved['full_conormal_matrix']]and rank(C)==17
 K=[[Q(v)for v in row]for row in saved['complete_conormal_kernel_including_membership']]
 assert len(K)==49 and rank(K)==49 and all(dot(row,v)==0 for row in C for v in K)
 constraints=[[Q(v)for v in row]for row in saved['all11_eliminated_raw_source_constraints']]
 assert len(constraints)==rank(constraints)==11 and all(dot(row,v[:60])==0 for row in constraints for v in K)
 image=[[v[50+j]for j in range(10)]for v in K];assert rank(image)==7
 forbidden=[PAIRS.index((i,4))for i in range(3)]
 assert all(v[j]==0 for v in image for j in forbidden)
 cases=[]
 for case in saved['cases']:
  Avec=list(map(Q,case['constant_A_upper']));A=zero(5)
  for c,(i,j)in zip(Avec,PAIRS):setentry(A,i,j,scale(q,c))
  constantJ=projective(A);mixed=[]
  for P in raw[:50]:
   J=projective(plusmatrix(A,P));Jp=projective(P)
   mixed.append({quad:add(J[quad],scale(constantJ[quad],-1),scale(Jp[quad],-1))for quad in QUADS})
  M=[];rhs=[]
  for quad in QUADS:
   for ex in mons4:
    if ex[4]:continue
    M.append([J[quad].get(ex,Q(0))for J in mixed]);rhs.append(-constantJ[quad].get(ex,Q(0)))
  M+=[row[:50]for row in constraints];rhs+=[-dot(row[50:],Avec)for row in constraints]
  point=list(map(Q,case['raw_L_particular']));kernel=[[Q(v)for v in row]for row in case['raw_L_homogeneous_kernel']]
  assert rank(M)==20 and len(kernel)==rank(kernel)==30
  assert all(dot(row,point)==b for row,b in zip(M,rhs))and all(dot(row,v)==0 for row in M for v in kernel)
  assert all(dot(row,v)==0 for row in M for v in E)
  qr=list(map(Q,case['quotient_particular']));qk=[[Q(v)for v in row]for row in case['complete_quotient_kernel']];cuts=case['Euler_gauge_coordinate_indices']
  assert len(qk)==rank(qk)==25 and determinant([[v[j]for j in cuts]for v in E])==Q(case['Euler_gauge_minor'])!=0
  assert all(qr[j]==0 for j in cuts)and all(v[j]==0 for v in qk for j in cuts)
  assert all(dot(row,qr)==b for row,b in zip(M,rhs))and all(dot(row,v)==0 for row in M for v in qk)
  cases.append(dict(name=case['name'],all_modH_full350_and_conormal_rows_reconstructed=True,raw_L_rank=20,complete_raw_kernel=30,Euler_redundancy=5,complete_projective_kernel=25))
 # Universal inverse-shear polynomial identity, with no chosen C/b pivot.
 names=['l','b0','b1','b2','c01','c02','c12'];l=parse('l',names);bs=[parse('b'+str(i),names)for i in range(3)];cs=[parse(v,names)for v in ['c01','c02','c12']];AA=zero(5)
 for p,(i,j)in zip(cs,[(0,1),(0,2),(1,2)]):setentry(AA,i,j,p)
 for i,p in enumerate(bs):setentry(AA,i,3,p)
 setentry(AA,4,3,l);G=[[l if i==j else{}for j in range(5)]for i in range(5)];GI=[[l if i==j else{}for j in range(5)]for i in range(5)]
 for i in range(3):G[i][4]=scale(bs[i],-1);GI[i][4]=bs[i]
 def pm(A,B):return [[add(*(mul(a,b)for a,b in zip(row,col)))for col in zip(*B)]for row in A]
 assert pm(G,GI)==[[mul(l,l)if i==j else{}for j in range(5)]for i in range(5)]
 target=[[dict(p)for p in row]for row in AA]
 for i in range(3):target[i][3]={};target[3][i]={}
 assert pm(pm(G,AA),transpose(G))==[[mul(mul(l,l),p)for p in row]for row in target]
 return dict(raw60_independent=True,vertical_rank=24,source_vertical_rank=79,projective_source_dimension=55,full35_by66_conormal_rank=17,complete_raw_conormal_kernel=49,allowed_A_exactly7=True,all11_conditions=True,universal_inverse_shear_verified=True,all_two_modH_sources=cases),E

def wedge(A,B):return [[add(mul(A[i],B[j]),scale(mul(B[i],A[j]),-1))for j in range(len(A))]for i in range(len(A))]
def lie(Y,B,offset):
 return [[add(*(mul(Y[k],diff(B[i][j],offset+k))for k in range(3)),
              *(scale(mul(B[k][j],diff(Y[i],offset+k)),-1)for k in range(3)),
              *(scale(mul(B[i][k],diff(Y[j],offset+k)),-1)for k in range(3)))for j in range(3)]for i in range(3)]
def matrix_add(*AA):return [[add(*(A[i][j]for A in AA))for j in range(len(AA[0]))]for i in range(len(AA[0]))]

def review_bridge(case,ore,E):
 packet=ore['complete_affine25_bridge'];names=packet['parameter_order']+['x0','x1','x2','z','w'];assert len(names)==30
 pars=[var(30,i)for i in range(25)];uu=[var(30,25+i)for i in range(3)];z,w=var(30,28),var(30,29);q=add(mul(uu[0],uu[2]),scale(mul(uu[1],uu[1]),-1))
 C=zero(3);ij=(0,2)if case['name']=='Delta_nonzero'else(0,1);setentry(C,*ij,constant(30,1));X=[add(*(mul(C[i][j],diff(q,25+j))for j in range(3)))for i in range(3)]
 DD=zero(3)
 for value,pair in zip(pars[:3],P3):setentry(DD,*pair,value)
 aa=pars[3:6];bb=pars[6:9];cc=pars[9];MM=[pars[10+3*i:13+3*i]for i in range(3)];dd=pars[19:22];ee=pars[22:25]
 L=zero(5)
 for i,j in P3:setentry(L,i,j,mul(w,DD[i][j]))
 for i in range(3):
  setentry(L,4,i,add(scale(X[i],-1),mul(w,aa[i])))
  setentry(L,3,i,add(*(mul(MM[i][j],uu[j])for j in range(3)),mul(z,dd[i]),mul(w,ee[i])))
 setentry(L,4,3,add(*(mul(bb[i],uu[i])for i in range(3)),mul(cc,z)))
 point=[Q(0)]*50;bridge=[[Q(0)]*25 for j in range(50)]
 for at,(i,j)in enumerate(PAIRS):
  for ex,value in L[i][j].items():
   spatial=ex[25:];assert sum(spatial)==1;row=5*at+spatial.index(1);parameters=ex[:25];assert sum(parameters)<=1
   if sum(parameters):bridge[row][parameters.index(1)]+=value
   else:point[row]+=value
 assert point==list(map(Q,packet['rawL_particular']))and bridge==[[Q(v)for v in row]for row in packet['rawL_parameter_matrix']]
 Erows=transpose(E);assert Erows==[[Q(v)for v in row]for row in packet['Euler5_matrix']]
 full=[row+v for row,v in zip(bridge,Erows)];recovery=[[Q(v)for v in row]for row in packet['all30_recovery_matrix']]
 assert matmul(recovery,full)==[[Q(i==j)for j in range(30)]for i in range(30)]
 oldK=transpose([[Q(v)for v in row]for row in case['raw_L_homogeneous_kernel']]);difference=[[Q(v)-c]for v,c in zip(case['raw_L_particular'],point)]
 assert matmul(full,matmul(recovery,oldK))==oldK and matmul(full,matmul(recovery,difference))==difference
 A=zero(5)
 for i,j in P3:setentry(A,i,j,C[i][j])
 setentry(A,4,3,constant(30,1));P=[[add(mul(w,L[i][j]),mul(q,A[i][j]))for j in range(5)]for i in range(5)]
 coords=[*uu,z,w];aff=[[restrict(add(P[i][j],scale(mul(coords[i],P[4][j]),-1),mul(coords[j],P[4][i])),29)for j in range(4)]for i in range(4)]
 B=matrix_add([[mul(q,p)for p in row]for row in C],wedge(uu,X),DD,scalematrix(wedge(uu,aa),-1))
 T=[add(aa[i],scale(dd[i],-1),scale(mul(cc,uu[i]),-1),scale(X[i],-1))for i in range(3)]
 bu=add(*(mul(bb[i],uu[i])for i in range(3)))
 Y=[scale(add(mul(uu[i],q),mul(uu[i],bu),*(mul(MM[i][j],uu[j])for j in range(3)),ee[i]),-1)for i in range(3)]
 expected=zero(4)
 for i,j in P3:setentry(expected,i,j,B[i][j])
 for i in range(3):setentry(expected,i,3,add(Y[i],mul(z,T[i])))
 assert aff==expected
 Jaff=strict(aff,25);PJ=projective(P,25)
 for tri,value in Jaff.items():assert add(value,restrict(PJ[(*tri,4)],29))=={}
 assert add(*(scale(mul(coords[i],PJ[tuple(j for j in range(5)if j!=i)]),(-1)**i)for i in range(5)))=={}
 BJ=strict(B,25);LT=lie(T,B,25);FY=matrix_add(lie(Y,B,25),scalematrix(wedge(Y,T),-1))
 assert Jaff[(0,1,2)]==BJ[(0,1,2)]
 assert all(add(Jaff[(i,j,3)],FY[i][j],mul(z,LT[i][j]))=={}for i,j in P3)
 return dict(all25_parameters_and_Euler5_recovered_exactly=True,explicit30_by50_left_inverse_verified=True,old30_kernel_and_particular_recovered=True,universal_projective350_to_affine_identity=True,universal_full_Ore_Jacobi_equivalence=True)

def review_BT(which):
 # First recover the entire B-Poisson support using only constant QQ rows.
 names=['D0','D1','D2','a0','a1','a2','x0','x1','x2'];n=9;u=[var(n,6+i)for i in range(3)];C=zero(3);ij=(0,2)if which==0 else(0,1);setentry(C,*ij,constant(n,1));q=add(mul(u[0],u[2]),scale(mul(u[1],u[1]),-1));X=[add(*(mul(C[i][j],diff(q,6+j))for j in range(3)))for i in range(3)];D=zero(3)
 for v,pair in zip([var(n,i)for i in range(3)],P3):setentry(D,*pair,v)
 B=matrix_add([[mul(q,p)for p in row]for row in C],wedge(u,X),D,scalematrix(wedge(u,[var(n,3+i)for i in range(3)]),-1));BJ=strict(B,6)[(0,1,2)];spatial=coefficient_space(BJ,[6,7,8]);rows=[]
 for ex,p in spatial.items():
  if sum(ex)>=1:
   assert all(sum(e[:6])==1 for e in p);rows.append([p.get(tuple(int(j==k)for k in range(n)),Q(0))for j in range(6)])
 assert rank(rows)==3
 allowed=[1,3,5]if which==0 else[0,3,4]
 assert all(row[j]==0 for row in rows for j in allowed)
 forbidden=[j for j in range(6)if j not in allowed]
 assert substitute(BJ,{j:{}for j in forbidden})=={}
 # All translation equations and their remaining support are checked over QQ.
 names=['p0','p1','delta','c','t0','t1','t2','x0','x1','x2'];n=10;p0,p1,delta,c=[var(n,i)for i in range(4)];u=[var(n,7+i)for i in range(3)];ts=[var(n,4+i)for i in range(3)];C=zero(3);setentry(C,*ij,constant(n,1));q=add(mul(u[0],u[2]),scale(mul(u[1],u[1]),-1));X=[add(*(mul(C[i][j],diff(q,7+j))for j in range(3)))for i in range(3)]
 av=[p0,{},p1]if which==0 else[p0,p1,{}]
 B=matrix_add([[mul(q,p)for p in row]for row in C],wedge(u,X),[[mul(delta,p)for p in row]for row in C],scalematrix(wedge(u,av),-1));T=[add(ts[i],scale(mul(c,u[i]),-1),scale(X[i],-1))for i in range(3)];LT=lie(T,B,7)
 linear=[]
 for i,j in P3:
  for ex,value in coefficient_space(LT[i][j],[7,8,9]).items():
   if sum(ex)==1:linear.append(value)
 Mt=[]
 for value in linear:
  row=[]
  for j in range(3):
   derivative=diff(value,4+j);assert not derivative or set(derivative)=={(0,)*n};row.append(derivative.get((0,)*n,Q(0)))
  Mt.append(row)
 assert rank(Mt)==(3 if which==0 else 2)
 if which==0:
  solution=[mul(add(c,constant(n,1)),p0),{},mul(add(constant(n,1),scale(c,-1)),p1)]
  conditions=[mul(c,add(delta,scale(mul(p0,p1),-1)))]
 else:
  solution=[ts[0],add(p1,scale(mul(c,p0),Q(-1,2))),{}]
  conditions=[mul(c,p1),add(scale(mul(c,delta),-4),mul(c,mul(p0,p0)),scale(mul(p0,p1),-2),scale(mul(p1,ts[0]),2))]
 sub={4+i:v for i,v in enumerate(solution)};residual=[]
 for i,j in P3:residual+=list(coefficient_space(substitute(LT[i][j],sub),[7,8,9]).values())
 assert rowspace_polys([p for p in residual if p])==rowspace_polys(conditions)
 # Translation solution is forced modulo the same constant-row compatibilities.
 t_equations=[add(ts[i],scale(solution[i],-1))for i in range(3)if ts[i]!=solution[i]]
 all_linear_monomials=sorted(set().union(*(p.keys()for p in linear+t_equations+conditions)))
 # Derive translation pivots with constant QQ elimination, keeping polynomial RHS.
 basis={};compat=[]
 for value,row in zip(linear,Mt):
  rhs=scale(substitute(value,{4+i:{}for i in range(3)}),-1);row={j:v for j,v in enumerate(row)if v}
  for pivot in sorted(basis):
   if pivot in row:
    old,b=basis[pivot];factor=row[pivot];rhs=add(rhs,scale(b,-factor));row={j:row.get(j,Q(0))-factor*old.get(j,Q(0))for j in set(row)|set(old)};row={j:v for j,v in row.items()if v}
  if row:
   pivot=min(row);factor=row[pivot];basis[pivot]=({j:v/factor for j,v in row.items()},scale(rhs,1/factor))
  elif rhs:compat.append(rhs)
 assert set(basis)==({0,1,2}if which==0 else{1,2})
 for pivot,(row,rhs)in basis.items():
  assert add(*(scale(solution[j],v)for j,v in row.items()),scale(rhs,-1))=={}
 # Substituting the forced pivots leaves exactly the displayed full residual
 # ideal; hence no c=0/p1=0 branch is lost.
 return dict(full_B_linear_rank=3,full_B_support_dimension=3,full_B_all_remaining_Jacobi_zero=True,full_T_constant_pivot_rank=rank(Mt),full_BT_residual_ideal_checked_by_QQ_rowspace=True,c_zero_complement_retained=True,BT_dimension_upper=3,fixed_nonzero_c_dimension_upper=2)

def polynomial_exact_division(A,B):
 assert B;result={};left=dict(A);db=max(e[0]for e in B);cb=B[(db,)]
 while left and max(e[0]for e in left)>=db:
  da=max(e[0]for e in left);value=left[(da,)]/cb;term={(da-db,):value};result=add(result,term);left=add(left,scale(mul(term,B),-1))
 assert not left;return result
def polynomial_determinant(matrix):
 n=len(matrix);assert n;A=[[dict(p)for p in row]for row in matrix];previous=constant(1,1);sign=1
 for k in range(n-1):
  at=next((i for i in range(k,n)if A[i][k]),None)
  if at is None:return{}
  if at!=k:A[k],A[at]=A[at],A[k];sign=-sign
  pivot=A[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):A[i][j]=polynomial_exact_division(add(mul(A[i][j],pivot),scale(mul(A[i][k],A[k][j]),-1)),previous)
   A[i][k]={}
  previous=pivot
 return scale(A[-1][-1],sign)
def evaluate_univariate(p,c):return sum(value*c**ex[0]for ex,value in p.items())

def review_principal_blocks(which,ore):
 n=4;c=var(n,0);u=[var(n,i+1)for i in range(3)];q=add(mul(u[0],u[2]),scale(mul(u[1],u[1]),-1));C=zero(3);ij=(0,2)if which==0 else(0,1);setentry(C,*ij,constant(n,1));X=[add(*(mul(C[i][j],diff(q,1+j))for j in range(3)))for i in range(3)]
 B=matrix_add([[mul(q,p)for p in row]for row in C],wedge(u,X));T=[scale(add(X[i],mul(c,u[i])),-1)for i in range(3)]
 radial=[[scale(mul(u[i],u[j]),-1)for i in range(3)]for j in range(3)];linear=[];constants=[]
 for i,j in it.product(range(3),repeat=2):Y=[{}for k in range(3)];Y[i]=scale(u[j],-1);linear.append(Y)
 for i in range(3):Y=[{}for k in range(3)];Y[i]=constant(n,-1);constants.append(Y)
 matrices=[];certificates=[];roots={Q(0)}
 for basis,degree,saved in zip([radial,linear,constants],[3,2,1],ore['all_principal_blocks']):
  mons=monomials(3,degree);columns=[]
  for Y in basis:
   F=matrix_add(lie(Y,B,1),scalematrix(wedge(Y,T),-1));column=[]
   for i,j in P3:
    coefficients=coefficient_space(F[i][j],[1,2,3])
    for ex in mons:
     p=coefficients.get(ex,{});assert all(e[1:]==(0,0,0)for e in p);column.append({(e[0],):v for e,v in p.items()})
   columns.append(column)
  matrix=transpose(columns);assert matrix==[[parse(v,['c'])for v in row]for row in saved['all_matrix_entries']]
  rows=saved['generic_minor_rows'];cols=saved['generic_minor_columns'];assert len(rows)==len(cols)==saved['generic_rank']
  minor=polynomial_determinant([[matrix[i][j]for j in cols]for i in rows]);assert minor==parse(saved['generic_minor'],['c'])and minor
  leftover=minor;allroots=list(map(Q,saved['all_exceptional_roots']));roots.update(allroots)
  for value in allroots:
   factor={(1,):Q(1),(0,):-value};factor={e:c for e,c in factor.items()if c}
   while max(e[0]for e in leftover)>0 and evaluate_univariate(leftover,value)==0:leftover=polynomial_exact_division(leftover,factor)
  assert set(leftover)=={(0,)}and leftover[(0,)]!=0
  matrices.append(matrix);certificates.append(dict(shape=saved['shape'],rank_lower=len(cols),exact_polynomial_minor=saved['generic_minor'],all_complex_roots_exhausted_by_exact_division=list(map(str,allroots))))
 assert roots=={Q(e['c'])for e in ore['all_exceptional_c_including_zero']}
 exceptions=[]
 for case in ore['all_exceptional_c_including_zero']:
  value=Q(case['c']);ranks=[]
  for matrix,cert,expected_rank in zip(matrices,case['rank_minors'],case['exact_principal_ranks']):
   rows,cols=cert['rows'],cert['columns'];assert len(rows)==len(cols)==expected_rank
   A=[[evaluate_univariate(matrix[i][j],value)for j in cols]for i in rows];det=determinant(A);assert det==Q(cert['minor'])and det!=0;ranks.append(expected_rank)
  fiber=sum(len(matrix[0])-lower for matrix,lower in zip(matrices,ranks));base=3 if value==0 else 2
  assert fiber==case['Y_affine_fiber_dimension_upper']and base==case['BT_base_dimension_upper']
  exceptions.append(dict(c=str(value),all_exact_minors_verified=True,rank_lower_bounds=ranks,fiber_upper=fiber,base_upper=base,normalized_source_upper=base+fiber))
 generic_fiber=sum(len(matrix[0])-c['rank_lower']for matrix,c in zip(matrices,certificates));normalized=max([3+generic_fiber]+[e['normalized_source_upper']for e in exceptions]);assert normalized==(5 if which==0 else 6)
 orbit=7 if which==0 else 6;fixed=orbit+normalized;assert fixed==12
 return dict(all_three_full_diagonal_matrices_reconstructed=True,all_generic_polynomial_minors_Bareiss_exact=True,all_generic_kernel_recomputations_unnecessary=True,generic_rank_lower_certificates=certificates,all_exceptional_c_including_zero=exceptions,generic_Y_fiber_upper=generic_fiber,whole_normalized_source_upper=normalized,restored_A_orbit_dimension=orbit,whole_fixed_surface_affine_source_upper=fixed,PGL_saturated_projective_dimension_upper=24+fixed-12-1)

def main():
 resource.setrlimit(resource.RLIMIT_CPU,(179,180));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();started=time.monotonic()
 manifest_path=ROOT/'research/classification/inputs/quadric_cone_zero_packet_manifest_20261107.json';manifest=json.loads(manifest_path.read_text())
 for item in manifest['files']:assert hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()==item['sha256']
 assert hashlib.sha256((ROOT/manifest['helper']).read_bytes()).hexdigest()==manifest['helper_sha256']
 source_path=ROOT/'reports/CLASSIFICATION_QUADRIC_CONE_ZERO_SOURCE_20261107.json';ore_path=ROOT/'reports/CLASSIFICATION_QUADRIC_CONE_ZERO_ORE_20261107.json';source=json.loads(source_path.read_text());ore=json.loads(ore_path.read_text())
 src,E=review_source(source);cases=[]
 for which,(s_case,o_case)in enumerate(zip(source['cases'],ore['cases'])):
  assert s_case['name']==o_case['name']
  cases.append(dict(name=s_case['name'],whole25_Ore_bridge=review_bridge(s_case,o_case,E),whole_BT_support=review_BT(which),whole_Y_rank_bound=review_principal_blocks(which,o_case)))
 assert all(c['whole_Y_rank_bound']['PGL_saturated_projective_dimension_upper']==23 for c in cases)
 paths=[manifest_path,source_path,ore_path,ROOT/manifest['helper'],ROOT/'research/classification/QUADRIC_CONE_ZERO_SOURCE.md',ROOT/'research/classification/QUADRIC_CONE_ORE_NONDOMINANCE.md']
 result=dict(task_id='G',all_passed=True,no_producer_import_or_CAS=True,independent_implementation='Fraction sparse polynomials, all characteristiczero exact matrices; generic minors by univariate fractionfree Bareiss, all complexroots exhausted by exact linear division',whole_source=src,all_two_normalized_cases=cases,whole_quadric_cone_rankzero_PGL_image_upper=23,
  full_A_orbit_cover_review='Allowed A has lambda,b0,b1,b2,C01,C02,C12. Lambda0 or C0 makes D|H=-lambda*q^2*(C12*u0-C02*u1+C01*u2)=0. Otherwise universal normal shear removes allb without pivot; vertexscale normalizeslambda. Dual q-norm Delta=4C01*C12-C02^2 has exactly the nonisotropic and nonzeroisotropic SO3 strata, by orthogonal/Witt bases overC; tensor scalar normalizes Delta to-1. All complements are retained.',
  triangular_fiber_logic_review='For each fixed wholeBT value, homogeneousYfiber is solved by spatial degrees3(b3),2(M9),1(e3). The exact diagonalblocks depend only on c; previous solutions/lowercouplings alter affine RHS only. Sum of nullity bounds therefore bounds every affinefiber, even inconsistent ones. All minor roots and c0 are treated separately, with fixednonzeroc base<=2.',
  restoration_review='AutS dimension12 and tensor scalar1 act freely on PGL5 x nonzero fixedS sources through the PGL factor plus scalar, without assuming any individualPoissonstabilizer trivial. Normalizedbounds5/6 plusAorbit7/6 givefixedsource<=12, so image<=24+12-12-1=23.',
  accepted_external_geometry_boundary='PrimeCI square primary/global conormal abelian condition and nonplane multiplicitytwo firstjet input are the earlier cited rankzero surface results; no clean/node/vertical assumption added.',
  no_E_coverage_or_newcomponent_or_globalG_completion=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},reviewer_bytes=Path(__file__).read_text(),reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=round(time.monotonic()-started,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert not args.output.exists();args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,PGL_upper=23,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))

if __name__=='__main__':main()
