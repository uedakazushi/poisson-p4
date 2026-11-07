"""Complete graded frame at the EXISTING b=2 Ga9 constant-c anchor over QQ.

All ten actual source derivatives are retained.  The dummy rational-function
field returned to the shared native exporter contains only QQ constants.
"""
from pathlib import Path
import json,time,itertools
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,ROWS,jacobi
ROOT=Path(__file__).resolve().parents[1]
def build_frame():
 started=time.monotonic();prior_path=ROOT/'reports/CLASSIFICATION_GA9_CONSTANT_WEIGHT_BLOCKS_20261107.json';prior=json.loads(prior_path.read_text());b=s.Symbol('b');xx=s.symbols('h x y z w');E=s.Matrix(xx);weights=(0,4,1,3,2);c=-3;K=s.QQ.frac_field(b)
 source_path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL211_SHORT1_CONSTANT_C_20261107.json';source=json.loads(source_path.read_text())['nil_short_translation_constant_c'];pp=[s.Symbol(z)for z in source['parameters']];env={str(z):z for z in pp};env.update(dict(zip(['x','y','z','w'],xx[1:])));full=s.zeros(5)
 for(i,j),text in zip(itertools.combinations(range(1,5),2),source['universal_brackets']):
  f=s.sympify(text,locals=env);full[i,j]=s.expand(xx[0]**2*f.subs({a:a/xx[0]for a in xx[1:]},simultaneous=True));full[j,i]=-full[i,j]
 values={z:s.Integer(0)for z in pp};values.update({s.Symbol('r'):1,b:2});P=(8*full.subs(values)).applyfunc(s.expand)
 def vector(Q):return s.Matrix([s.Poly(s.expand(Q[ij]),*xx).coeff_monomial(ex)for ij,ex in COLUMNS])
 def sparse(Q):return {ij:dict(s.Poly(s.expand(Q[ij]),*xx).terms())for ij in PAIRS if Q[ij]!=0}
 dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ)
 ps=sparse(P);DD=[]
 for ij,ex in COLUMNS:
  variation={key:dict(val)for key,val in ps.items()};variation.setdefault(ij,{})[ex]=variation.get(ij,{}).get(ex,0)+1;DD.append(jacobi(variation))
 D=s.Matrix.hstack(*DD);V=[];O=[]
 for i,j in itertools.product(range(5),repeat=2):
  A=s.zeros(5);A[i,j]=1;v=A*E;V.append(vector(E*v.T-v*E.T));O.append(vector((A*P+P*A.T-xx[j]*P.diff(xx[i])).applyfunc(s.expand)))
 G=s.Matrix.hstack(*V,*O,vector(P));F=G.row_join(s.Matrix.hstack(*[vector(8*full.diff(z).subs(values).applyfunc(s.expand))for z in pp]));assert D*F==s.zeros(350,F.cols)
 cw=[sum(a*w for a,w in zip(ex,weights))-weights[i]-weights[j]for(i,j),ex in COLUMNS];rw=[sum(a*w for a,w in zip(ex,weights))-sum(weights[i]for i in quad)-2*c for quad,ex in ROWS]
 Fgroups={}
 for j in range(F.cols):
  ws={cw[i]for i,a in enumerate(F[:,j])if a!=0};assert len(ws)<=1
  if ws:Fgroups.setdefault(next(iter(ws)),[]).append(j)
 F_rank=0;normals=[];nw=[];implicit=[];solvers={};allrows=[];newblocks=[]
 for block in prior['complete_D_weight_blocks']:
  wt=block['column_weight'];cs=block['ambient_columns'];rs=block['ambient_rows'];BB=D[rs,cs];fp=list(dm(F[cs,Fgroups.get(wt,[])]).rref()[1]);F_rank+=len(fp);basis=F[cs,Fgroups.get(wt,[])][:,fp];rank_basis=len(fp);ker=dm(BB).nullspace(divide_last=True).to_Matrix().T
  for j in range(ker.cols):
   joined=basis.row_join(ker[:,j]);nr=len(dm(joined).rref()[1])
   if nr>rank_basis:
    local=ker[:,j];scale=s.ilcm(1,*[s.denom(a)for a in local]);local=scale*local;gcd=s.igcd(*[abs(int(a))for a in local if a]);local/=gcd
    v=s.zeros(150,1)
    for k,col in enumerate(cs):v[col]=local[k]
    normals.append(v);nw.append(wt-c);basis=basis.row_join(local);rank_basis=nr
  cols=list(dm(BB).rref()[1]);rows=list(dm(BB[:,cols].T).rref()[1]);assert rank_basis==len(cs)-len(cols)
  selected_cs=[cs[i]for i in cols];selected_rs=[rs[i]for i in rows];implicit.extend(selected_cs);allrows.extend(selected_rs)
  if selected_cs:
   inverse=dm(D[selected_rs,selected_cs]).inv().to_Matrix();solvers[wt-c]=(selected_cs,selected_rs,[[K.from_sympy(z)for z in row]for row in inverse.tolist()])
  newblocks.append(dict(column_weight=wt,ambient_columns=cs,ambient_rows=rs,columns=cols,rows=rows,rank=len(cols),shape=[len(rs),len(cs)]))
 assert F_rank==48 and len(implicit)==98 and nw==[1,2,2,3];NN=s.Matrix.hstack(*normals);assert D*NN==s.zeros(350,4)
 for block in newblocks:
  cs=block['ambient_columns'];wt=block['column_weight'];ns=[j for j,v in enumerate(nw)if v==wt-c];us=[col for col in implicit if col in cs];coord=s.zeros(len(cs),len(us))
  for j,col in enumerate(us):coord[cs.index(col),j]=1
  assert len(dm(F[cs,Fgroups.get(wt,[])].row_join(coord).row_join(NN[cs,ns])).rref()[1])==len(cs)
 prior=dict(prior);prior['complete_D_weight_blocks']=newblocks
 return dict(started=started,prior_path=prior_path,source_path=source_path,prior=prior,b=b,xx=xx,E=E,weights=weights,c=c,K=K,P=P,D=D,F=F,normals=normals,normal_weights=nw,implicit=implicit,solvers=solvers,selected_rows=allrows,column_weights=cw,row_relative_weights=rw)
