"""Reusable complete rational graded source frame for projectivequadratic Poisson tensors."""
import itertools,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,ROWS,jacobi

def graded_frame(P,source_directions,xx,weights,c):
 start=time.monotonic()
 def vector(Q):return s.Matrix([s.Poly(s.expand(Q[ij]),*xx).coeff_monomial(ex)for ij,ex in COLUMNS])
 def sparse(Q):return{ij:dict(s.Poly(s.expand(Q[ij]),*xx).terms())for ij in PAIRS if Q[ij]!=0}
 ps=sparse(P);assert jacobi(ps)==s.zeros(350,1);columns=[]
 for ij,ex in COLUMNS:
  variation={key:dict(val)for key,val in ps.items()};variation.setdefault(ij,{})[ex]=variation.get(ij,{}).get(ex,0)+1;columns.append(jacobi(variation))
 D=s.Matrix.hstack(*columns);E=s.Matrix(xx);V=[];O=[]
 for i,j in itertools.product(range(5),repeat=2):
  A=s.zeros(5);A[i,j]=1;field=A*E;V.append(vector(E*field.T-field*E.T));O.append(vector((A*P+P*A.T-xx[j]*P.diff(xx[i])).applyfunc(s.expand)))
 G=s.Matrix.hstack(*V,*O,vector(P));F=G.row_join(s.Matrix.hstack(*[vector(v)for v in source_directions]));assert D*F==s.zeros(350,F.cols)
 dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ)
 cw=[sum(a*w for a,w in zip(ex,weights))-weights[i]-weights[j]for(i,j),ex in COLUMNS];rw=[sum(a*w for a,w in zip(ex,weights))-sum(weights[i]for i in quad)-2*c for quad,ex in ROWS]
 assert all(cw[i]==c for i,a in enumerate(vector(P))if a!=0)
 assert all(rw[i]==cw[j]-c for i in range(350)for j in range(150)if D[i,j]!=0)
 def colgroups(A):
  g={}
  for j in range(A.cols):
   ws={cw[i]for i,a in enumerate(A[:,j])if a!=0};assert len(ws)<=1
   if ws:g.setdefault(next(iter(ws)),[]).append(j)
  return g
 fg=colgroups(F);gg=colgroups(G);blocks=[];normals=[];normalweights=[];implicit=[];selected=[];rD=rF=rG=0
 for wt in sorted(set(cw)):
  cs=[i for i,w in enumerate(cw)if w==wt];rs=[i for i,w in enumerate(rw)if w==wt-c];BB=D[rs,cs];cols=list(dm(BB).rref()[1]);rows=list(dm(BB[:,cols].T).rref()[1]);fp=list(dm(F[cs,fg.get(wt,[])]).rref()[1]);gr=len(dm(G[cs,gg.get(wt,[])]).rref()[1]);rD+=len(cols);rF+=len(fp);rG+=gr;basis=F[cs,fg.get(wt,[])][:,fp];rr=len(fp);ker=dm(BB).nullspace(divide_last=True).to_Matrix().T
  for j in range(ker.cols):
   candidate=ker[:,j];joined=basis.row_join(candidate);nr=len(dm(joined).rref()[1])
   if nr>rr:
    scale=s.ilcm(1,*[s.denom(a)for a in candidate]);candidate=scale*candidate;gcd=s.igcd(*[abs(int(a))for a in candidate if a]);candidate/=gcd
    v=s.zeros(150,1)
    for k,col in enumerate(cs):v[col]=candidate[k]
    normals.append(v);normalweights.append(wt-c);basis=basis.row_join(candidate);rr=nr
  ci=[cs[i]for i in cols];ri=[rs[i]for i in rows];implicit+=ci;selected+=ri;inverse=dm(D[ri,ci]).inv().to_Matrix()if ci else s.zeros(0)
  blocks.append(dict(column_weight=wt,relative_weight=wt-c,shape=[len(rs),len(cs)],D_rank=len(cols),G_rank=gr,F_rank=len(fp),implicit_columns=ci,implicit_rows=ri,implicit_inverse=[[str(a)for a in row]for row in inverse.tolist()],minor_determinant=str(dm(D[ri,ci]).det())if ci else'1'))
 NN=s.Matrix.hstack(*normals);assert D*NN==s.zeros(350,NN.cols)
 U=s.zeros(150,len(implicit))
 for j,col in enumerate(implicit):U[col,j]=1
 assert len(dm(F.row_join(U).row_join(NN)).rref()[1])==150
 return dict(canonical_point={str(i)+','+str(j):str(P[i,j])for i,j in PAIRS if P[i,j]!=0},weights=list(weights),canonical_tensor_weight=c,complete_D_shape=[350,150],D_rank=rD,actual_G_rank=rG,actual_full_source_F_rank=rF,full_kernel_dimension=150-rD,normal_count=len(normals),true_normal_relative_weights=normalweights,complete_blocks=blocks,true_normal_vectors=[{str(i):str(a)for i,a in enumerate(v)if a!=0}for v in normals],implicit_columns=implicit,implicit_rows=selected,all150_directsum_verified=True,full_DF_and_DN_zero=True,implicit_max_relative_weight=max(cw)-c,all_Jacobi_max_relative_weight=max(rw),seconds=time.monotonic()-start)
