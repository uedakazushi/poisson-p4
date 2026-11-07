"""Reusable completeQQ(b) graded ambient frame for actual Ga9 deformations."""
from pathlib import Path
import json,time,itertools
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,ROWS,jacobi
ROOT=Path(__file__).resolve().parents[1]
def build_frame():
 started=time.monotonic();prior_path=ROOT/'reports/CLASSIFICATION_GA9_CONSTANT_WEIGHT_BLOCKS_20261107.json';prior=json.loads(prior_path.read_text());b=s.Symbol('b');xx=s.symbols('h x y z w');E=s.Matrix(xx);weights=(0,4,1,3,2);c=-3;K=s.QQ.frac_field(b)
 normal_weights=(1,2,2,3);env={str(z):z for z in xx};env['b']=b;P=s.zeros(5)
 for ij,f in prior['canonical_family'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=env);P[j,i]=-P[i,j]
 def vector(Q):return s.Matrix([s.Poly(s.expand(Q[ij]),*xx).coeff_monomial(ex)for ij,ex in COLUMNS])
 def sparse(Q):return {ij:dict(s.Poly(s.expand(Q[ij]),*xx).terms())for ij in PAIRS if Q[ij]!=0}
 dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(K)
 ps=sparse(P);DD=[]
 for ij,ex in COLUMNS:
  variation={key:dict(val)for key,val in ps.items()};variation.setdefault(ij,{})[ex]=variation.get(ij,{}).get(ex,0)+1;DD.append(jacobi(variation).applyfunc(s.expand))
 D=s.Matrix.hstack(*DD);V=[];O=[]
 for i,j in itertools.product(range(5),repeat=2):
  A=s.zeros(5);A[i,j]=1;v=A*E;V.append(vector(E*v.T-v*E.T));O.append(vector((A*P+P*A.T-xx[j]*P.diff(xx[i])).applyfunc(s.expand)))
 F=s.Matrix.hstack(*V,*O,vector(P),vector(P.diff(b)))
 assert (D*F).applyfunc(s.expand)==s.zeros(350,F.cols)
 cw=[sum(a*w for a,w in zip(ex,weights))-weights[i]-weights[j]for(i,j),ex in COLUMNS];rw=[sum(a*w for a,w in zip(ex,weights))-sum(weights[i]for i in quad)-2*c for quad,ex in ROWS]
 Fgroups={}
 for j in range(F.cols):
  ws={cw[i]for i,a in enumerate(F[:,j])if a!=0};assert len(ws)<=1
  if ws:Fgroups.setdefault(next(iter(ws)),[]).append(j)
 F_rank=0;normals=[];nw=[];implicit=[];solvers={};allrows=[]
 for block in prior['complete_D_weight_blocks']:
  wt=block['column_weight'];cs=block['ambient_columns'];rs=block['ambient_rows'];BB=D[rs,cs];fc=Fgroups.get(wt,[]);FF=F[cs,fc];fp=list(dm(FF).rref()[1])if fc else [];F_rank+=len(fp)
  basis=FF[:,fp];rank_basis=len(fp);ker=dm(BB).nullspace(divide_last=True).to_Matrix().T
  for j in range(ker.cols):
   joined=basis.row_join(ker[:,j]);nr=len(dm(joined).rref()[1])
   if nr>rank_basis:
    v=s.zeros(150,1)
    for k,col in enumerate(cs):v[col]=ker[k,j]
    normals.append(v);nw.append(wt-c);basis=joined;rank_basis=nr
  assert rank_basis==len(cs)-block['rank']
  selected_cs=[cs[i]for i in block['columns']];selected_rs=[rs[i]for i in block['rows']];implicit.extend(selected_cs);allrows.extend(selected_rs)
  if selected_cs:
   assert all(rw[i]==wt-c for i in selected_rs)
   inverse=dm(D[selected_rs,selected_cs]).inv().to_Matrix();solvers[wt-c]=(selected_cs,selected_rs,[[K.from_sympy(z)for z in row]for row in inverse.tolist()])
 assert F_rank==48 and len(implicit)==98 and nw==list(normal_weights)
 NN=s.Matrix.hstack(*normals);assert dm(D).matmul(dm(NN)).is_zero_matrix
 # Directsum F48 + independent D-coordinate98 + genuine normal4 is all150.
 for block in prior['complete_D_weight_blocks']:
  cs=block['ambient_columns'];wt=block['column_weight'];ns=[j for j,v in enumerate(nw)if v==wt-c];us=[col for col in implicit if col in cs];coord=s.zeros(len(cs),len(us))
  for j,col in enumerate(us):coord[cs.index(col),j]=1
  assert len(dm(F[cs,Fgroups.get(wt,[])].row_join(coord).row_join(NN[cs,ns])).rref()[1])==len(cs)
 return dict(started=started,prior_path=prior_path,prior=prior,b=b,xx=xx,E=E,weights=weights,c=c,K=K,P=P,D=D,F=F,normals=normals,normal_weights=nw,implicit=implicit,solvers=solvers,selected_rows=allrows,column_weights=cw,row_relative_weights=rw)
