#!/usr/bin/env python3
"""Complete QQ(b) ambient ranks of the weighted Ga9 c=1 canonical family.

The known (0,4,1,3,2) grading splits the whole350x150 differential,
avoiding a dense rational-function elimination. No normal jet is computed.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,ROWS,jacobi
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();source_path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL211_SHORT1_CONSTANT_C_20261107.json';source=json.loads(source_path.read_text())['nil_short_translation_constant_c']
 h,x,y,z,w=xx=s.symbols('h x y z w');b=s.Symbol('b');pp=[s.Symbol(v)for v in source['parameters']];env={str(v):v for v in pp};env.update(dict(zip(['x','y','z','w'],xx[1:])));full=s.zeros(5)
 for (i,j),f in zip(itertools.combinations(range(1,5),2),source['universal_brackets']):
  f=s.sympify(f,locals=env);full[i,j]=s.expand(h*h*f.subs({a:a/h for a in xx[1:]},simultaneous=True));full[j,i]=-full[i,j]
 values={v:s.Integer(0)for v in pp};values.update({s.Symbol('r'):1,b:b});P=(8*full.subs(values)).applyfunc(s.expand)
 def vector(Q):return s.Matrix([s.Poly(s.expand(Q[ij]),*xx).coeff_monomial(ex)for ij,ex in COLUMNS])
 def sparse(Q):return {ij:dict(s.Poly(s.expand(Q[ij]),*xx).terms())for ij in PAIRS if Q[ij]!=0}
 weights=(0,4,1,3,2);source_weight=-3
 column_weights=[sum(a*t for a,t in zip(ex,weights))-weights[i]-weights[j]for(i,j),ex in COLUMNS]
 row_weights=[sum(a*t for a,t in zip(ex,weights))-sum(weights[i]for i in quad)for quad,ex in ROWS]
 assert all(column_weights[i]==source_weight for i,a in enumerate(vector(P))if a!=0)
 ps=sparse(P);assert jacobi(ps).applyfunc(s.expand)==s.zeros(350,1)
 DD=[]
 for ij,ex in COLUMNS:
  variation={key:dict(val)for key,val in ps.items()};variation.setdefault(ij,{})[ex]=variation.get(ij,{}).get(ex,0)+1
  DD.append(jacobi(variation).applyfunc(s.expand))
 D=s.Matrix.hstack(*DD)
 assert all(row_weights[i]==column_weights[j]+source_weight for i in range(350)for j in range(150)if D[i,j]!=0)
 E=s.Matrix(xx);vertical=[];orbit=[]
 for i,j in itertools.product(range(5),repeat=2):
  A=s.zeros(5);A[i,j]=1;field=A*E;vertical.append(vector(E*field.T-field*E.T));orbit.append(vector((A*P+P*A.T-xx[j]*P.diff(xx[i])).applyfunc(s.expand)))
 V=s.Matrix.hstack(*vertical);G=V.row_join(s.Matrix.hstack(*orbit,vector(P)));directions=[8*full.diff(a).subs(values).applyfunc(s.expand)for a in pp];F=G.row_join(s.Matrix.hstack(*map(vector,directions)))
 assert (D*F).applyfunc(s.expand)==s.zeros(350,F.cols)
 field=s.QQ.frac_field(b);dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(field)
 def block_certificate(A):
  if not A.rows or not A.cols:return dict(rank=0,rows=[],columns=[],determinant='1')
  cols=list(dm(A).rref()[1]);rows=list(dm(A[:,cols].T).rref()[1]);det=dm(A[rows,cols]).det()if cols else 1
  assert det!=0
  return dict(rank=len(cols),rows=rows,columns=cols,determinant=str(s.factor(s.sympify(str(det)))))
 dblocks=[]
 for wt in sorted(set(column_weights)):
  cs=[i for i,v in enumerate(column_weights)if v==wt];rs=[i for i,v in enumerate(row_weights)if v==wt+source_weight];cert=block_certificate(D[rs,cs]);cert.update(column_weight=wt,ambient_rows=rs,ambient_columns=cs,shape=[len(rs),len(cs)]);dblocks.append(cert)
 def weighted_rank(A):
  groups={}
  for col in range(A.cols):
   ws={column_weights[i]for i,a in enumerate(A[:,col])if a!=0}
   assert len(ws)<=1,(col,ws)
   if ws:groups.setdefault(next(iter(ws)),[]).append(col)
  blocks=[]
  for wt,cs in sorted(groups.items()):
   rs=[i for i,v in enumerate(column_weights)if v==wt];cert=block_certificate(A[rs,cs]);cert.update(weight=wt,ambient_rows=rs,ambient_columns=cs,shape=[len(rs),len(cs)]);blocks.append(cert)
  return sum(a['rank']for a in blocks),blocks
 vr,vblocks=weighted_rank(V);gr,gblocks=weighted_rank(G);fr,fblocks=weighted_rank(F);dr=sum(a['rank']for a in dblocks)
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(source_path.relative_to(ROOT)):hashlib.sha256(source_path.read_bytes()).hexdigest(),str(Path(__file__).relative_to(ROOT)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'code/classification_rank3_mixed_component_review.py':hashlib.sha256((ROOT/'code/classification_rank3_mixed_component_review.py').read_bytes()).hexdigest()},producer_bytes=Path(__file__).read_text(),
  parameters='r=1,b free,alpha=p=q=A=D=G=H=k=0; k=1 point is the actual x-translation of this family.',weights=list(weights),canonical_tensor_weight=source_weight,canonical_family={str(i)+','+str(j):str(s.factor(P[i,j]))for i,j in PAIRS if P[i,j]!=0},
  coefficient_field='QQ(b)',complete_fullambient_D_shape=[350,150],whole350_Jacobi_zero=True,whole_DF_zero=True,ranks_QQ_b=dict(D=dr,vertical=vr,orbit_Euler_scalar=gr,all_ten_source=fr),full_projective_tangent_dimension=125-dr,actual_fullsource_projective_dimension_lower=fr-25,
  complete_D_weight_blocks=dblocks,complete_vertical_weight_blocks=vblocks,complete_orbit_weight_blocks=gblocks,complete_fullsource_weight_blocks=fblocks,
  maximal_D_block_columns=max(c['shape'][1]for c in dblocks),component_squeeze=dr+fr==150,scope='Complete generic fullambient ranks by an exact grading decomposition, not a chosen smaller submatrix. Source membership/componenthood requires the stated squeeze or actual knownfamily recognition; no higher normaljet is inferred.',seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_GA9_CONSTANT_WEIGHT_BLOCKS_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(ranks=result['ranks_QQ_b'],squeeze=result['component_squeeze'],max_block_cols=result['maximal_D_block_columns'],seconds=result['seconds'])))
if __name__=='__main__':main()
