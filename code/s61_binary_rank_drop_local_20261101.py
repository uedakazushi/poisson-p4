#!/usr/bin/env python3
"""Exact local dimension at a dominant binary rank-drop point, not saturation."""
import sys,json,time,hashlib,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import sympy as s,numpy as np
from itertools import combinations_with_replacement
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
parser=argparse.ArgumentParser(description='Exact local dimension3 at a dominant binary rank-drop point.');parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
if args.output.exists():raise FileExistsError(args.output)
start=time.time();p=32003;m=S61()
family_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json';r=json.loads(family_path.read_text());assert r['all_passed'];b={int(k[1:])-1:s.Rational(v) for k,v in r['sample'].items()};two,F=specialize(m,b,True);fixed=two.index(29);cols=[i for i in range(45) if i!=fixed];q0=s.zeros(45,1);q0[fixed]=1
J=s.zeros(135,45)
for row,f in enumerate(F):
 for mon,c in f.items():
  if len(mon)==2:
   i,j=mon
   if i!=fixed:J[row,cols.index(i)]+=c*q0[j]
   if j!=fixed:J[row,cols.index(j)]+=c*q0[i]
  elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
T=s.Matrix.hstack(*J.nullspace());n=T.cols;assert n==11 and T[-1,:]==s.zeros(1,11)
_,piv=J.to_DM().rref();free=[i for i in range(45) if i not in piv];assert T[free,:]==s.eye(n)
z,uu,bf=primitive_beta(m);bf=bf.subs({s.Symbol(k):s.Rational(v) for k,v in r['sample'].items()});a,c,d,h,k=s.symbols('a c d h k');v1=s.Matrix([a,1,0,c]);v2=s.Matrix([-c,0,1,d]);origin={a:0,c:0,d:0,h:0,k:0}
sec=s.Matrix([(v1.T*bf.subs(dict(zip(z,v1)))*v2)[0],(v1.T*bf.subs(dict(zip(z,v2)))*v2)[0]])
assert sec.subs(origin)==s.zeros(2,1)
curve_jac=sec.jacobian([a,c,d]).subs(origin);assert curve_jac.rank()==2
direction=curve_jac.nullspace()[0]
x,y,zz,w=X[1:];ell1=x-a*y+c*zz;ell2=w-c*y-d*zz;fam=(ell1*ell2+h*ell1**2+k*ell2**2)*(v1*v2.T-v2*v1.T)
ds=[]
for v in [a,c,d,h,k]:
 P=X[0]**2*J0;P[1:,1:]+=fam.diff(v).subs(origin);qs=s.Matrix([m.normalize(P)['parameters'][i] for i in two]);qs-=q0*qs[fixed];ds.append(qs[cols,:].col_join(s.zeros(1,1)))
family=s.Matrix.hstack(s.Matrix.hstack(*ds[:3])*direction,ds[3],ds[4]);coordinates=family[free,:];assert T*coordinates==family and J*family==s.zeros(135,3)
_,params=coordinates.T.to_DM().rref();params=list(params);normal=[i for i in range(n) if i not in params];assert params==[0,1,3] and len(normal)==8
parameter_minor=coordinates[params,:].det();assert parameter_minor
def mod(v):v=s.Rational(v);return int(v.p)*pow(int(v.q),-1,p)%p
def arr(M):return np.array([[mod(v) for v in row] for row in M.tolist()],dtype=np.int64)
def rref(A):
 A=A.copy()%p;row=0;piv=[]
 for col in range(A.shape[1]):
  choices=np.flatnonzero(A[row:,col])
  if len(choices)==0:continue
  j=row+int(choices[0]);A[[row,j]]=A[[j,row]];A[row]=A[row]*pow(int(A[row,col]),-1,p)%p
  mask=np.flatnonzero(A[:,col]);mask=mask[mask!=row]
  if len(mask):A[mask]=(A[mask]-A[mask,col,None]*A[row,None,:])%p
  piv.append(col);row+=1
  if row==A.shape[0]:break
 return A,piv
j=arr(J);t=arr(T);_,C=rref(j);_,R=rref(j[:,C].T);rank=len(C);rr,_=rref(np.concatenate([j[np.ix_(R,C)],np.eye(rank,dtype=np.int64)],axis=1));U=np.zeros((45,rank),dtype=np.int64);U[C,:]=rr[:,rank:];JU=j@U%p
assert len(C)==len(R)==34
assert np.array_equal(j[R,:]@U%p,np.eye(34,dtype=np.int64))
assert np.count_nonzero(j@t%p)==0
H=np.zeros((135,45,45),dtype=np.int64)
for row,f in enumerate(F):
 for mon,c in f.items():
  if fixed in mon:continue
  i,l=[cols.index(v) for v in mon]+[44]*(2-len(mon));v=mod(c)
  if i==l:H[row,i,l]+=v
  else:H[row,i,l]+=v*pow(2,-1,p)%p;H[row,l,i]+=v*pow(2,-1,p)%p
H%=p
# Establish the exact QQ quadratic rank before lifting modular row cancellations.
monsQQ=list(combinations_with_replacement(range(n),2));indQQ={mon:i for i,mon in enumerate(monsQQ)};cacheQQ={}
def multQQ(i,j):
 key=tuple(sorted((i,j)))
 if key not in cacheQQ:
  result=[s.S.Zero]*len(monsQQ)
  for a,v in enumerate(T[i,:]):
   if not v:continue
   for b,w in enumerate(T[j,:]):
    if w:result[indQQ[tuple(sorted((a,b)))]]+=v*w
  cacheQQ[key]=result
 return cacheQQ[key]
Qexact=[]
for f in F:
 row=[s.S.Zero]*len(monsQQ)
 for mon,c in f.items():
  if fixed in mon:continue
  i,jj=[cols.index(v) for v in mon]+[44]*(2-len(mon))
  for k,v in enumerate(multQQ(i,jj)):row[k]+=c*v
 Qexact.append(row)
co=s.Matrix.hstack(*J.T.nullspace()).T
Aexact=(co.to_DM()*s.Matrix(Qexact).to_DM()).to_Matrix();_,qpiv=Aexact.T.to_DM().rref();assert len(qpiv)==14
exact_basis=Aexact[list(qpiv),:]
m2=list(combinations_with_replacement(range(n),2));m3=list(combinations_with_replacement(range(n),3));i3={mon:i for i,mon in enumerate(m3)}
Q=np.zeros((135,len(m2)),dtype=np.int64)
for row in range(135):
 q=t.T@(H[row]@t%p)%p
 for i,(a,c) in enumerate(m2):Q[row,i]=q[a,c]*(1 if a==c else 2)%p
w2=-U@Q[R,:]%p;A2=(Q-JU@Q[R,:])%p;assert len(rref(A2)[1])==14
C3=np.zeros((135,len(m3)),dtype=np.int64)
for row in range(135):
 cross=2*t.T@(H[row]@w2%p)%p
 for i in range(n):
  for k,mon in enumerate(m2):C3[row,i3[tuple(sorted((i,)+mon))]]+=cross[i,k]
C3%=p;A3=(C3-JU@C3[R,:])%p
rr,pc=rref(np.concatenate([A2,A3],axis=1));qrank=sum(k<len(m2) for k in pc);q=rr[:qrank,:len(m2)];c=rr[qrank:len(pc),len(m2):]
assert qrank==14 and len(c)==54
nmons={d:list(combinations_with_replacement(normal,d)) for d in range(6)};nidx={d:{v:i for i,v in enumerate(ms)} for d,ms in nmons.items()}
nq=np.array([[row[m2.index(mon)] for mon in nmons[2]] for row in q]);nc=np.array([[row[m3.index(mon)] for mon in nmons[3]] for row in c]);assert len(rref(nq)[1])==14
results={};matrices={};row_indices={}
for degree in [3,4,5]:
 mat=[]
 for d,gens in [(2,nq),(3,nc)]:
  for row in gens:
   for monmul in nmons[degree-d]:
    out=np.zeros(len(nmons[degree]),dtype=np.int64)
    for mon,v in zip(nmons[d],row):out[nidx[degree][tuple(sorted(mon+monmul))]]+=v
    mat.append(out%p)
 mat=np.array(mat);rank=len(rref(mat)[1]);results[degree]=(rank,len(nmons[degree]));matrices[degree]=mat
 print('degree',degree,results[degree],flush=True)
 if rank==len(nmons[degree]):break
assert results=={3:(86,120),4:(328,330),5:(792,792)}
_,rows5=rref(matrices[5].T);assert len(rows5)==792
minor=matrices[5][rows5,:].copy();det=1
for col in range(792):
 ii=col+int(np.flatnonzero(minor[col:,col])[0])
 if ii!=col:minor[[col,ii]]=minor[[ii,col]];det=-det
 value=int(minor[col,col]);det=det*value%p;minor[col]=minor[col]*pow(value,-1,p)%p
 mask=np.flatnonzero(minor[col+1:,col])+col+1
 if len(mask):minor[mask]=(minor[mask]-minor[mask,col,None]*minor[col,None,:])%p
assert det
enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
paths=[Path(__file__),family_path,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
result=dict(id='S61_BINARY_RANK_DROP_LOCAL_20261101',all_passed=True,prime=p,coefficient_field='QQ; modular nonzero minors at32003',
 sample=r['sample'],chart_fixed_q=fixed,full_jacobian_rank=34,tangent_dimension=n,tangent=enc(T),
 exact_quadratic_rank=14,exact_quadratic_basis=enc(exact_basis),quadratic_monomials=m2,
 actual_binary_family_tangent=enc(coordinates),parameter_minor=str(parameter_minor),linear_parameters=params,normal_variables=normal,
 jacobian_pivot_rows=R,jacobian_pivot_columns=C,residual_quadratic_mod_prime=A2.tolist(),residual_cubic_mod_prime=A3.tolist(),
 normal_initial_quadrics_mod_prime=nq.tolist(),normal_initial_cubics_mod_prime=nc.tolist(),
 normal_initial_ranks_by_degree=results,quintic_monomials=nmons[5],quintic_minor_rows=rows5,quintic_minor_determinant_mod_prime=det,
 normal_hilbert_function_upper_bound=[1,8,22,34,2],normal_colength_upper_bound=67,local_dimension=3,
 proof_note='research/notes/S61_BINARY_RANK_DROP_LOCAL_20261101.md',
 limitations=['A specified point on the dominant rank-drop incidence; no universal multiplicity statement.',
              'The actual transverse length, t-nilpotence, and affine-branch exclusion remain open.'],
 input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.time()-start)
with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({k:result[k] for k in ['all_passed','local_dimension','normal_colength_upper_bound','quintic_minor_determinant_mod_prime','seconds']}))
