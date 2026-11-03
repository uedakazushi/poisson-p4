#!/usr/bin/env python3
"""All-a Smith certificates for actual CA/CB/N source differentials at Qa.

These are differential IMAGE ranks, not Zariski tangent dimensions of targets.
"""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp
from ambient_CA_global_flag_model_20261103_v63 import X,raw as rawA,canonical,vector,action
from ambient_CB_global_flag_model_20261103_v63 import raw as rawB
ROOT=Path(__file__).resolve().parents[1]
h,x,y,z,w=X;a=s.Symbol('a');WEIGHT=[0,1,4,2,3];PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2));RW=[WEIGHT[i]+WEIGHT[j]-WEIGHT[m]-WEIGHT[n] for i,j in PAIRS for m,n in MONS]
def matrix(M):return [[str(v) for v in row] for row in M.tolist()]
def certify(M):
 R=s.QQ.poly_ring(a);B=DomainMatrix.from_Matrix(M).convert_to(R).to_dense();D,U,V=smith_normal_decomp(B);assert U*B*V==D;du,dv=U.det(),V.det();assert du and dv and not R.to_sympy(du).has(a) and not R.to_sympy(dv).has(a)
 d=D.to_Matrix();diag=[s.factor(d[i,i]) for i in range(min(d.shape))];assert all(d[i,j]==0 for i in range(d.rows) for j in range(d.cols) if i!=j)
 return dict(matrix=matrix(M),left=matrix(U.to_Matrix()),right=matrix(V.to_Matrix()),left_determinant=str(R.to_sympy(du)),right_determinant=str(R.to_sympy(dv)),diagonal=list(map(str,diag)),generic_rank=sum(v!=0 for v in diag))
def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic();H=4*h*x*w/3-h*z*z/3-x*x*z;Q=canonical(rawA(a*h*h,-s.Rational(2,9),H));assert rawA(a*h*h,-s.Rational(2,9),H)==Q==rawB(a*h*h,-2*h*h/9,H)
 cubics=[s.prod((h,x,z,w)[i]for i in m)for m in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in m)];assert len(cubics)==16
 A=[rawA(f,0,0)for f in [h*h,h*x,x*x]]+[rawA(0,1,0)]+[rawA(0,0,v)for v in cubics]
 B=[rawB(h*h,0,0),rawB(0,h*h,0),rawB(a*h*x,-2*h*x/9,0),rawB(a*x*x,-2*x*x/9,0)]+[rawB(0,0,v)for v in cubics]
 N=[rawA(f,0,0)for f in [h*h,h*x]]+[rawA(0,1,0)]+[rawA(0,0,v)for v in cubics]
 O=[action(s.eye(5)[:,i]*s.eye(5)[j,:],Q)for i,j in itertools.product(range(5),repeat=2)];families={};full={}
 for name,source in [('CA',A),('CB',B),('N',N)]:
  basis=[canonical(v)for v in source];M=s.Matrix.hstack(*[vector(v)for v in basis+O]);weights=[]
  for j in range(len(basis)):
   ws={RW[i]for i in range(150)if M[i,j]!=0};assert len(ws)==1;weights.append(ws.pop())
  weights += [5+WEIGHT[i]-WEIGHT[j]for i,j in itertools.product(range(5),repeat=2)]
  assert all(M[i,j]==0 for i in range(150)for j in range(M.cols)if RW[i]!=weights[j])
  blocks=[]
  for wt in sorted(set(weights)):
   rows=[i for i in range(150)if RW[i]==wt];cols=[j for j,v in enumerate(weights)if v==wt]
   if rows:blocks.append(dict(weight=wt,rows=rows,columns=cols,**certify(M[rows,cols])))
  nonunit=[s.factor(s.sympify(v))for b in blocks for v in b['diagonal']if s.sympify(v).has(a)];rank=sum(b['generic_rank']for b in blocks);ranks={str(v):sum(s.sympify(k).subs(a,v)!=0 for b in blocks for k in b['diagonal'])for v in [0,s.Rational(2,9),1]}
  families[name]=dict(source_fiber_columns=len(source),full_matrix_shape=list(M.shape),generic_image_rank=rank,nonunit_invariant_factors=list(map(str,nonunit)),special_image_ranks=ranks,all_a_ranks_certified_by_unimodular_blocks=True)
  full[name]=dict(full150_source_orbit_matrix=matrix(M),column_weights=weights,blocks=blocks)
  print(name,'image ranks',ranks,'nonunit',nonunit,flush=True)
 archive=out/'certificate.json.gz'
 with archive.open('wb')as f:
  with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0)as zz:zz.write((json.dumps(dict(physical_coordinate_order='(h,x,y,z,w)=(old0,old1,old4,old2,old3)',raw_Qa150=list(map(str,vector(Q))),cubic_monomials=list(map(str,cubics)),row_weights=RW,families=full),separators=(',',':'))+'\n').encode())
 inputs=[Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py']
 report=dict(id='Q_CA_CB_N_SOURCE_SMITH_20261103_v63',all_passed=True,field='QQ[a]',families=families,target_tangent_ranks_claimed=False,scope='Entire physical Qa line source differentials, including a0 as a separate rank2 point; no target normality follows from rank alone.',certificate_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),input_sha256={str(v.relative_to(ROOT)):hashlib.sha256(v.read_bytes()).hexdigest()for v in inputs},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('ALL-a CA/CB/N SOURCE SMITH PASS',report['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
