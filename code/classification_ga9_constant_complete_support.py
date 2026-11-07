#!/usr/bin/env python3
"""One complete positively graded IFT reduction of the Ga9 constant-c source.

The PRIOR weight theorem bounds implicit polynomials by weight10 and ALL
remaining equations by weight16. This is a finite exact construction, not
an unbounded normal-jet sequence. Failed resource runs make no support claim.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from sympy.polys.groebnertools import groebner
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,ROWS,TRIPLES,QUADS,jacobi
ROOT=Path(__file__).resolve().parents[1]

def main():
 started=time.monotonic();prior_path=ROOT/'reports/CLASSIFICATION_GA9_CONSTANT_WEIGHT_BLOCKS_20261107.json';prior=json.loads(prior_path.read_text());b=s.Symbol('b');xx=s.symbols('h x y z w');E=s.Matrix(xx);weights=(0,4,1,3,2);c=-3;K=s.QQ.frac_field(b)
 R,t0,t1,t2,t3=ring('t0,t1,t2,t3',K,order='grevlex');tt=(t0,t1,t2,t3);normal_weights=(1,2,2,3);env={str(z):z for z in xx};env['b']=b;P=s.zeros(5)
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
 colindex={col:i for i,col in enumerate(COLUMNS)};rowindex={row:i for i,row in enumerate(ROWS)}
 def add_column(A,index,co):
  if not co:return
  ij,ex=COLUMNS[index];f=A.setdefault(ij,{});f[ex]=f.get(ex,R.zero)+co
  if not f[ex]:del f[ex]
 def entry(A,i,j):return A.get((i,j),{})if i<j else {ex:-co for ex,co in A.get((j,i),{}).items()}if i>j else {}
 def schouten(A,B=None):
  triples={};terms=[(A,A)]if B is None else[(A,B),(B,A)]
  for i,j,k in TRIPLES:
   f={}
   for AA,BB in terms:
    for u,v,w in[(i,j,k),(j,k,i),(k,i,j)]:
     for ell in range(5):
      for ex,co in entry(AA,u,ell).items():
       for ee,cc in entry(BB,v,w).items():
        if ee[ell]:
         exponent=tuple(ex[q]+ee[q]-(q==ell)for q in range(5));f[exponent]=f.get(exponent,R.zero)+co*cc*ee[ell]
   triples[(i,j,k)]={ex:co for ex,co in f.items()if co}
  result=[R.zero for _ in range(350)]
  for quad in QUADS:
   for at,ell in enumerate(quad):
    for ex,co in triples[quad[:at]+quad[at+1:]].items():
     exponent=tuple(ex[q]+(q==ell)for q in range(5));i=rowindex[(quad,exponent)];result[i]+=(-1)**at*co
  return result
 dcols=[[(i,K.from_sympy(D[i,j]))for i in range(350)if D[i,j]!=0]for j in range(150)]
 def linear(A):
  result=[R.zero for _ in range(350)]
  for ij,terms in A.items():
   for ex,co in terms.items():
    for i,a in dcols[colindex[(ij,ex)]]:result[i]+=co*a
  return result
 def nonlinear_weight(d):
  result=[R.zero for _ in range(350)]
  for i in range(1,d//2+1):
   j=d-i
   if i not in W or j not in W:continue
   terms=schouten(W[i])if i==j else schouten(W[i],W[j])
   result=[a+bb for a,bb in zip(result,terms)]
  return result
 W={};implicit_maps={};residual=[];checkpoint=ROOT/'work/GA9_CONSTANT_COMPLETE_SUPPORT_CHECKPOINT.json'
 for degree in range(1,17):
  nonlin=nonlinear_weight(degree);Wd={}
  if degree<=10:
   for j,weight in enumerate(nw):
    if weight==degree:
     for i,a in enumerate(normals[j]):
      if a!=0:add_column(Wd,i,tt[j]*K.from_sympy(a))
   if degree in solvers:
    cs,rs,inverse=solvers[degree]
    for j,col in enumerate(cs):
     value=-sum((nonlin[row]*inverse[j][k]for k,row in enumerate(rs)),R.zero)
     if value:implicit_maps[col]=value;add_column(Wd,col,value)
   W[degree]=Wd
   ll=linear(Wd);current=[a+bb for a,bb in zip(nonlin,ll)]
  else:current=nonlin
  assert all(not current[i]for i in allrows)
  assert all(not f or rw[i]==degree for i,f in enumerate(current))
  for i,f in enumerate(current):
   if f:residual.append((i,f))
  state=dict(completed_weight=degree,implicit_map_count=len(implicit_maps),residual_equation_count=len(residual),seconds=time.monotonic()-started,complete=False)
  checkpoint.write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state),flush=True)
 # Every implicit polynomial is now exact: allselected98 equations vanish;
 # no ambientcoordinate above weight10 exists and all350 rows end at16.
 assert all(rw[i]<=16 for i in range(350))
 polys=[]
 for _,f in residual:
  monic=f.monic()
  if monic not in polys:polys.append(monic)
 gb=groebner(polys,R);leading=[p.LM for p in gb];pure=[]
 for i in range(4):pure.append(any(ex[i]>0 and sum(ex)==ex[i]for ex in leading))
 maximal_radical=all(pure)
 inputs=[prior_path,Path(__file__),ROOT/'code/classification_rank3_mixed_component_review.py']
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),
  coefficient_field='QQ(b)',whole_D_rank98=True,actual_F_rank48=True,normal_weights=nw,complete_weightstable_cut102=True,all150_directsum_verified=True,
  complete_normal_vectors=[{str(i):str(a)for i,a in enumerate(v)if a!=0}for v in normals],implicit_columns=implicit,implicit_rows=allrows,
  prior_bound_implicit_weight10=True,prior_bound_all_equations_weight16=True,complete_exact_implicit_polynomials={str(i):str(f.as_expr())for i,f in implicit_maps.items()},
  all_selected_98_equations_exactly_zero=True,complete_remaining_equations=[dict(row=i,polynomial=str(f.as_expr()))for i,f in residual],distinct_monic_generators=[str(f.as_expr())for f in polys],
  complete_grevlex_basis=[str(f.as_expr())for f in gb],leading_pure_power_for_each_normal=pure,maximal_radical_origin=maximal_radical,
  generic_whole_source_component23_proved=maximal_radical,local_smoothness_or_uniqueness_claimed=False,
  proof='The Cstar-stable genuine102-dimensional cut is complementary to actualF48. Its98 independentJacobi rows have a uniqueequivariant formalimplicit solution. Allnormalweights1,2,2,3 have one sign; eachimplicitcoordinate has relativeweight<=10, so this formal solution is a finite exact polynomial. AllJacobi rows have relativeweight<=16. The saved exact finite equations therefore give the complete formalcut, not a truncatedjet. A maximalradical would make the truecut0dimensional, and Krull gives localdimension<=48, matching actualF48.',seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_GA9_CONSTANT_COMPLETE_SUPPORT_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');checkpoint.write_text(json.dumps(dict(completed_weight=16,complete=True,maximal_radical_origin=maximal_radical,seconds=result['seconds']),indent=2)+'\n');print(json.dumps(dict(complete=True,maximal_radical_origin=maximal_radical,seconds=result['seconds'])))
if __name__=='__main__':main()
