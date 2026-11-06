#!/usr/bin/env python3
"""G: finite source-bound span of translates of ONE normal coordinate.

Derives the45closed quadratic secondjet representation directly from original
S61 reconstruction. Only an exact nonzero rank certificate is sought modulo p;
no modular unit/zero is promoted to characteristic zero. All g are integral
symplectic matrices; the rational plane frame and all formula inputs are pinned.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import numpy as np
import sympy as s
from s61_model import S61,BASIS,X
from weight1_structure import J
ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(4),2))
MON=list(itertools.combinations_with_replacement(range(4),2))
P=32003

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(M):return [[str(c)for c in row]for row in M.tolist()]
def mod(M):
 return np.array([[int(s.numer(c)%P)*pow(int(s.denom(c)%P),-1,P)%P for c in row]for row in M.tolist()],dtype=np.int64)
def inverse(A):
 a=np.concatenate((A.copy()%P,np.eye(len(A),dtype=np.int64)),axis=1);det=1
 for j in range(len(A)):
  nz=np.flatnonzero(a[j:,j])
  if not len(nz):return None,0
  k=j+int(nz[0])
  if k!=j:a[[j,k]]=a[[k,j]];det=-det
  pivot=int(a[j,j]);det=det*pivot%P;a[j]=a[j]*pow(pivot,-1,P)%P
  for i in range(len(A)):
   if i!=j and a[i,j]:a[i]=(a[i]-a[i,j]*a[j])%P
 return a[:,len(A):],det%P

def representation_basis(model):
 weights=model.reduced['parameter_weights'];qids=[i for i,w in enumerate(weights)if w==2];qi={u:i for i,u in enumerate(qids)}
 raw=s.zeros(150,45)
 for i,u in zip(model.certificate['raw_free_indices'],range(61)):
  if u in qi:raw[i,qi[u]]=1
 for rec in model.data['eliminated_raw_variables']:
  for word,num,den in rec['polynomial']:
   if len(word)==1 and word[0]in qi:raw[rec['raw_index'],qi[word[0]]]+=s.Rational(num,den)
 aff=s.zeros(60,45)
 for ri,((i,j),mon)in enumerate(BASIS):
  ex=s.Poly(mon,*X).monoms()[0];word=tuple(k for k,n in enumerate(ex[1:])for _ in range(n))
  if i:
   if len(word)==2:aff[10*PAIRS.index((i-1,j-1))+MON.index(word),:]+=raw[ri,:]
  elif len(word)==1:
   for h in range(4):
    if h==j-1:continue
    pair=tuple(sorted((h,j-1)));sign=-1 if h<j-1 else 1
    mon2=tuple(sorted((h,word[0])));aff[10*PAIRS.index(pair)+MON.index(mon2),:]+=sign*raw[ri,:]
 dual=s.zeros(6)
 for k,(i,j)in enumerate(PAIRS):
  e=s.zeros(4);e[i,j]=1;e[j,i]=-1;b=-J*e*J
  for z,(a,c)in enumerate(PAIRS):dual[z,k]=b[a,c]
 M=s.kronecker_product(dual,s.eye(10))*aff
 # Exterior derivative on homogeneous quadratic2forms:16coefficients, rank15.
 triples=list(itertools.combinations(range(4),3));D=s.zeros(16,60)
 for col,((i,j),(a,b))in enumerate(itertools.product(PAIRS,MON)):
  for k in range(4):
   if k in(i,j):continue
   multiplicity=(a==k)+(b==k)
   if not multiplicity:continue
   remaining=b if a==k else a;tr=(k,i,j);sgn=(-1)**sum(tr[x]>tr[y]for x in range(3)for y in range(x+1,3))
   D[4*triples.index(tuple(sorted(tr)))+remaining,col]+=sgn*multiplicity
 assert D.rank()==15 and D*M==s.zeros(16,45)
 rows=list(M.T.rref()[1]);assert len(rows)==45
 L=s.zeros(45,60);small=M[rows,:].inv()
 for j,k in enumerate(rows):L[:,k]=small[:,j]
 assert L*M==s.eye(45)
 return qids,M,L,D

def finite_action(g,M,L):
 exterior=s.zeros(6);sym=s.zeros(10)
 for i,(a,b)in enumerate(PAIRS):
  for j,(c,d)in enumerate(PAIRS):exterior[i,j]=g[c,a]*g[d,b]-g[d,a]*g[c,b]
 for i,(c,d)in enumerate(MON):
  for j,(a,b)in enumerate(MON):sym[i,j]=g[a,c]*g[b,d]+(g[a,d]*g[b,c]if c!=d else 0)
 T=np.kron(mod(exterior),mod(sym))%P
 acted=T@M%P;A=L@acted%P
 assert np.array_equal(M@A%P,acted)
 return A

def groups():
 #64fixed four-transvection words. No adaptive beta/prime/group search.
 vv=[s.eye(4)[:,i]for i in range(4)]+[s.eye(4)[:,i]+s.eye(4)[:,j]for i,j in PAIRS]
 out=[]
 for k in range(64):
  g=s.eye(4);word=[]
  for j in range(4):
   ix=(k*(2*j+1)+j*j+3*j)%10;t=((k//(j+1)+2*j)%5)-2
   v=vv[ix];g=g*(s.eye(4)+t*J*v*v.T);word.append([ix,t])
  assert g.T*J*g==J and g.det()==1;out.append((g,word))
 return out

def main(output):
 start=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(1500000000,1500000000))
 model=S61();qids,M,L,D=representation_basis(model);Mm,Lm=mod(M),mod(L)
 framepath=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz'
 import gzip
 frame=json.loads(gzip.decompress(framepath.read_bytes()));W=s.Matrix(frame['embedding']);assert W.shape==(45,3)and W[:3,:]==s.eye(3)
 Wm=mod(W);basis={};records=[];selected=[];selectedrows=[];det=1
 for k,(g,word)in enumerate(groups()):
  A=finite_action(g,Mm,Lm);Wg=A@Wm%P;Wi,minor=inverse(Wg[:3])
  rec={'group_index':k,'integral_g':enc(g),'transvection_word':word,'first3_plane_chart_det_modp':minor};records.append(rec)
  if not minor:rec['chart_skipped']=True;continue
  row=(A[3]-Wg[3]@Wi@A[:3])%P;assert not np.any(row@Wm%P)
  rec['pulledback_n0_row45_modp']=row.tolist();v=row[3:].copy();pivot=None
  for j,old in sorted(basis.items()):v=(v-v[j]*old)%P
  nz=np.flatnonzero(v)
  if len(nz):
   pivot=int(nz[0]);pv=int(v[pivot]);basis[pivot]=v*pow(pv,-1,P)%P;det=det*pv%P;selected.append(k);selectedrows.append(row[3:].tolist())
  rec['new_pivot_normal_column']=pivot;rec['cumulative_rank']=len(basis)
  if len(basis)==42:break
 rank=len(basis);positive=rank==42
 if positive:
  _,check=inverse(np.array(selectedrows,dtype=np.int64));assert check and check==det%P
 paths=['code/s61_model.py','code/weight1_structure.py','legacy/v02/data/symplectic_slice_certificate.json','legacy/v02/data/symplectic_slice_equations.json','legacy/v02/data/symplectic_slice_reduced_equations.json',str(framepath.relative_to(ROOT)),'reports/S61_KNOWN_FIBER_AUDIT_20261031.json','reports/TORIC_UNIVERSAL_THREE_PLANE_20261101.json','research/classification/C0_S61_SLICE_SCALAR_AUDIT_20261106.md']
 out={'task_id':'G','positive_rank42_certificate':positive,'rank_modp':rank,'prime':P,'group_bound':64,'groups_evaluated':len(records),'selected_groups':selected,'selected42x42_rows_modp':selectedrows,'selected_determinant_modp':det%P if positive else None,'records':records,'q_to_u_indices_zero_based':qids,'closed_quadratic_beta2_basis60x45':enc(M),'basis_left_inverse45x60':enc(L),'closed_differential16x60':enc(D),'plane_W45x3_QQ':enc(W),'coefficient_field':'Generic C(W16); modular nonzero is a positive QQ minor witness only','span_formula':'a_g=row3(A_g)-row3(A_gW)*((A_gW)first3rows)^-1*(A_g)first3rows','normal_columns':'trailing42 since Wfirst3rows=I3','necessary_hypotheses':['knownfive affine span is3plane generically and equivariant','same fullS61 frame and finitegauge action','all base chart denominators and inverse minors retained'],'unit_or_exhaustion_computed':False,'input_sha256':{p:sha(ROOT/p)for p in paths},'code_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start}
 output.parent.mkdir(parents=True,exist_ok=True)
 with output.open('x')as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:out[k]for k in ['positive_rank42_certificate','rank_modp','groups_evaluated','selected_determinant_modp','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
