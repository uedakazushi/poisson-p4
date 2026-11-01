#!/usr/bin/env python3
"""Rank65 binary formal ring across a dominant rank-drop stratum."""
import sys,json,time,hashlib,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import sympy as s,numpy as np
from itertools import combinations_with_replacement
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
parser=argparse.ArgumentParser(description='Rank65 formal extension across a dominant binary rank-drop stratum.');parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
if args.output.exists():raise FileExistsError(args.output)
start=time.time();p=32003;m=S61()
family_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json';local_path=ROOT/'reports/S61_BINARY_RANK_DROP_LOCAL_20261101.json';r=json.loads(family_path.read_text());local=json.loads(local_path.read_text());assert r['all_passed'] and local['all_passed']
assert local['normal_colength_upper_bound']==67 and local['quintic_minor_determinant_mod_prime']!=0;b={int(k[1:])-1:s.Rational(v) for k,v in r['sample'].items()};two,F=specialize(m,b,True);fixed=two.index(29);cols=[i for i in range(45) if i!=fixed];q0=s.zeros(45,1);q0[fixed]=1
J=s.zeros(135,45)
for row,f in enumerate(F):
 for mon,c in f.items():
  if len(mon)==2:
   i,j=mon
   if i!=fixed:J[row,cols.index(i)]+=c*q0[j]
   if j!=fixed:J[row,cols.index(j)]+=c*q0[i]
  elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
T=s.Matrix(local['tangent']).applyfunc(s.Rational);n=T.cols;assert n==11 and J*T==s.zeros(135,11) and T[-1,:]==s.zeros(1,11)
_,piv=J.to_DM().rref();free=[i for i in range(45) if i not in piv];assert T[free,:]==s.eye(n)
z,uu,bf=primitive_beta(m);bf=bf.subs({s.Symbol(k):s.Rational(v) for k,v in r['sample'].items()});a,c,d,h,k=s.symbols('a c d h k');v1=s.Matrix([a,1,0,c]);v2=s.Matrix([-c,0,1,d]);origin={a:0,c:0,d:0,h:0,k:0}
sec=s.Matrix([(v1.T*bf.subs(dict(zip(z,v1)))*v2)[0],(v1.T*bf.subs(dict(zip(z,v2)))*v2)[0]]);direction=sec.jacobian([a,c,d]).subs(origin).nullspace()[0]
x,y,zz,w=X[1:];ell1=x-a*y+c*zz;ell2=w-c*y-d*zz;fam=(ell1*ell2+h*ell1**2+k*ell2**2)*(v1*v2.T-v2*v1.T)
ds=[]
for v in [a,c,d,h,k]:
 P=X[0]**2*J0;P[1:,1:]+=fam.diff(v).subs(origin);qs=s.Matrix([m.normalize(P)['parameters'][i] for i in two]);qs-=q0*qs[fixed];ds.append(qs[cols,:].col_join(s.zeros(1,1)))
family=s.Matrix.hstack(s.Matrix.hstack(*ds[:3])*direction,ds[3],ds[4]);coordinates=family[free,:];assert T*coordinates==family and J*family==s.zeros(135,3)
_,params=coordinates.T.to_DM().rref();params=list(params);normal=[i for i in range(n) if i not in params];assert params==local['linear_parameters']==[0,1,3] and len(normal)==8
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
assert len(C)==len(R)==34 and np.array_equal(j[R,:]@U%p,np.eye(34,dtype=np.int64))
H=np.zeros((135,45,45),dtype=np.int64)
for row,f in enumerate(F):
 for mon,c in f.items():
  if fixed in mon:continue
  i,l=[cols.index(v) for v in mon]+[44]*(2-len(mon));v=mod(c)
  if i==l:H[row,i,l]+=v
  else:H[row,i,l]+=v*pow(2,-1,p)%p;H[row,l,i]+=v*pow(2,-1,p)%p
H%=p;t=t[:,normal];n=len(normal)
m2=list(combinations_with_replacement(range(n),2));m3=list(combinations_with_replacement(range(n),3));m4=list(combinations_with_replacement(range(n),4));i3={mon:i for i,mon in enumerate(m3)};i4={mon:i for i,mon in enumerate(m4)}
Q=np.zeros((135,len(m2)),dtype=np.int64)
for row in range(135):
 q=t.T@(H[row]@t%p)%p
 for k,(a,c) in enumerate(m2):Q[row,k]=q[a,c]*(1 if a==c else 2)%p
w2=-U@Q[R,:]%p;A2=(Q-JU@Q[R,:])%p
C3=np.zeros((135,len(m3)),dtype=np.int64)
for row in range(135):
 cross=2*t.T@(H[row]@w2%p)%p
 for i in range(n):
  for k,mon in enumerate(m2):C3[row,i3[tuple(sorted((i,)+mon))]]+=cross[i,k]
C3%=p;w3=-U@C3[R,:]%p;A3=(C3-JU@C3[R,:])%p
C4=np.zeros((135,len(m4)),dtype=np.int64)
for row in range(135):
 cross=2*t.T@(H[row]@w3%p)%p;square=w2.T@(H[row]@w2%p)%p
 for i in range(n):
  for k,mon in enumerate(m3):C4[row,i4[tuple(sorted((i,)+mon))]]+=cross[i,k]
 for i,mi in enumerate(m2):
  for k,mk in enumerate(m2):C4[row,i4[tuple(sorted(mi+mk))]]+=square[i,k]
C4%=p;w4=-U@C4[R,:]%p;A4=(C4-JU@C4[R,:])%p
base=np.concatenate([A2,A3,A4],axis=1);ncols=base.shape[1];rows=list(base)
for row2,row3 in zip(A2,A3):
 for a in range(n):
  out=np.zeros(ncols,dtype=np.int64)
  for mon,v in zip(m2,row2):out[len(m2)+i3[tuple(sorted((a,)+mon))]]+=v
  for mon,v in zip(m3,row3):out[len(m2)+len(m3)+i4[tuple(sorted((a,)+mon))]]+=v
  if np.any(out):rows.append(out%p)
 for ma in m2:
  out=np.zeros(ncols,dtype=np.int64)
  for mon,v in zip(m2,row2):out[len(m2)+len(m3)+i4[tuple(sorted(ma+mon))]]+=v
  if np.any(out):rows.append(out%p)
rels=np.array(rows);rr,piv=rref(rels);rk=len(piv);length=1+n+ncols-rk
print('actual modular length',length,'rank',rk,'shape',rels.shape,'time',time.time()-start,flush=True)
trow=np.concatenate([w2[-1,:],w3[-1,:],w4[-1,:]])%p
rankt=len(rref(np.concatenate([rr[:rk,:],trow[None,:]],axis=0))[1]);print('t extra rank',rankt-rk,flush=True)
tmult=[trow]
for a in range(n):
 out=np.zeros(ncols,dtype=np.int64)
 for mon,v in zip(m2,w2[-1,:]):out[len(m2)+i3[tuple(sorted((a,)+mon))]]+=v
 for mon,v in zip(m3,w3[-1,:]):out[len(m2)+len(m3)+i4[tuple(sorted((a,)+mon))]]+=v
 tmult.append(out%p)
for ma in m2:
 out=np.zeros(ncols,dtype=np.int64)
 for mon,v in zip(m2,w2[-1,:]):out[len(m2)+len(m3)+i4[tuple(sorted(ma+mon))]]+=v
 tmult.append(out%p)
rankideal=len(rref(np.concatenate([rr[:rk,:],np.array(tmult)],axis=0))[1]);print('t ideal rank',rankideal-rk,'infinity length',length-(rankideal-rk),flush=True)
length3=1+n+len(m2)-len(rref(rels[:,:len(m2)])[1]);length4=1+n+len(m2)+len(m3)-len(rref(rels[:,:len(m2)+len(m3)])[1]);print('Hilbert',[1,n,length3-1-n,length4-length3,length-length4],flush=True)

assert length==65 and rk==430 and rankt==431 and rankideal==431
assert [1,n,length3-1-n,length4-length3,length-length4]==[1,8,22,34,0]
# Save nonzero minors of the actual nonhomogeneous relation matrix, and of
# that matrix with t appended. Their characteristic-zero use is in the note.
_,selected_rows=rref(rels.T);assert len(selected_rows)==430
selected_columns=piv
minor430=rels[np.ix_(selected_rows,selected_columns)].copy()
def determinant_mod(M):
 M=M.copy();det=1
 for col in range(len(M)):
  ii=col+int(np.flatnonzero(M[col:,col])[0])
  if ii!=col:M[[col,ii]]=M[[ii,col]];det=-det
  value=int(M[col,col]);det=det*value%p;M[col]=M[col]*pow(value,-1,p)%p
  mask=np.flatnonzero(M[col+1:,col])+col+1
  if len(mask):M[mask]=(M[mask]-M[mask,col,None]*M[col,None,:])%p
 return det%p
det430=determinant_mod(minor430);assert det430
aug=np.concatenate([rels[selected_rows,:],trow[None,:]],axis=0)
_,augcols=rref(aug);assert len(augcols)==431
det431=determinant_mod(aug[:,augcols]);assert det431
# A regular binary point on the SAME fixed-beta threefold: Q1=xw+x^2+w^2.
P=X[0]**2*J0;xx,yy,zz,ww=X[1:];P[2,3]+=xx*ww+xx*xx+ww*ww;P[3,2]-=xx*ww+xx*xx+ww*ww
q1=s.Matrix([m.normalize(P)['parameters'][i] for i in two]);assert q1[fixed]==1
J1=s.zeros(135,45)
for row,f in enumerate(F):
 value=0
 for mon,c in f.items():
  if len(mon)==2:
   i,jj=mon;value+=c*q1[i]*q1[jj]
   if i!=fixed:J1[row,cols.index(i)]+=c*q1[jj]
   if jj!=fixed:J1[row,cols.index(jj)]+=c*q1[i]
  elif len(mon)==1:J1[row,44]+=c*q1[mon[0]]
 assert value==0
T1=s.Matrix.hstack(*J1.nullspace());assert T1.shape==(45,10) and T1[-1,:]==s.zeros(1,10)
mons1=list(combinations_with_replacement(range(10),2));idx1={mon:i for i,mon in enumerate(mons1)};cache={}
def product1(i,j):
 key=tuple(sorted((i,j)))
 if key not in cache:
  vals=[s.S.Zero]*55
  for a,v in enumerate(T1[i,:]):
   if not v:continue
   for bb,w in enumerate(T1[j,:]):
    if w:vals[idx1[tuple(sorted((a,bb)))]]+=v*w
  cache[key]=vals
 return cache[key]
Q1=[]
for f in F:
 row=[s.S.Zero]*55
 for mon,c in f.items():
  if fixed in mon:continue
  i,jj=[cols.index(v) for v in mon]+[44]*(2-len(mon))
  for k,v in enumerate(product1(i,jj)):row[k]+=c*v
 Q1.append(row)
co1=s.Matrix.hstack(*J1.T.nullspace()).T;A1=(co1.to_DM()*s.Matrix(Q1).to_DM()).to_Matrix();_,goodrows=A1.T.to_DM().rref();assert len(goodrows)==5
hess=[]
for row in A1[list(goodrows),:].tolist():
 M=s.zeros(10)
 for c,(i,jj) in zip(row,mons1):M[i,jj]+=c;M[jj,i]+=c
 hess+=M.tolist()
assert len(s.Matrix(hess).to_DM().rref()[1])==7
enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
paths=[Path(__file__),family_path,local_path,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',ROOT/'research/notes/S61_BINARY_UNIFORM_CUBIC_20261101.md',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
result=dict(id='S61_BINARY_RANK_DROP_FORMAL_RING_20261101',all_passed=True,prime=p,
 sample=r['sample'],linear_parameters=params,normal_variables=normal,jacobian_pivot_rows=R,jacobian_pivot_columns=C,
 residual_quadratic_mod_prime=A2.tolist(),residual_cubic_mod_prime=A3.tolist(),residual_quartic_mod_prime=A4.tolist(),
 homogenizer_quadratic_mod_prime=w2[-1,:].tolist(),homogenizer_cubic_mod_prime=w3[-1,:].tolist(),homogenizer_quartic_mod_prime=w4[-1,:].tolist(),
 relation_matrix_shape=list(rels.shape),relation_minor_rows=selected_rows,relation_minor_columns=selected_columns,relation_minor_mod_prime=det430,
 augmented_homogenizer_minor_columns=augcols,augmented_homogenizer_minor_mod_prime=det431,
 relation_rank_mod_prime=430,augmented_homogenizer_rank_mod_prime=431,
 same_beta_regular_point_q=list(map(str,q1)),same_beta_regular_full_jacobian_rank=35,same_beta_regular_quadratic_rank=5,same_beta_regular_common_radical_dimension=3,
 same_fixed_beta_generic_normal_length_lower_bound=65,transverse_length=65,transverse_hilbert_function=[1,8,22,34],
 completed_ring_finite_free_over='C[[v0,v1,v3]]',completed_ring_rank=65,reduced_completed_ring='C[[v0,v1,v3]]',nilradical_rank=64,nilradical_fourth_power_zero=True,
 homogenizer_nonzero=True,homogenizer_square_zero=True,homogenizer_annihilates_nilradical=True,completed_infinity_quotient_rank=64,local_t_saturation_empty=True,
 proof_note='research/notes/S61_BINARY_RANK_DROP_FORMAL_RING_20261101.md',
 limitations=['No freeness assertion for each individual associated-graded piece: embedding dimension jumps.',
              'Only the specified point and a justified open of the rank-drop stratum; not all its endpoints.',
              'No full affine-fiber exhaustion.'],
 input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.time()-start)
with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({k:result[k] for k in ['all_passed','relation_minor_mod_prime','augmented_homogenizer_minor_mod_prime','completed_ring_rank','homogenizer_nonzero','seconds']}))
