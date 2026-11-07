#!/usr/bin/env python3
"""Export the complete true transverse local Jacobi ideal, without a jet cutoff.

Usage: python this.py CASE [OUTPUT_DIRECTORY]. The existing reviewed complement
contains101 injective D-coordinate directions plus its one genuine normal.
A completed local standard basis is required for any support conclusion.
"""
import sys
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[1]/'code'))
import sympy as s,json,itertools,time,math,functools
from pathlib import Path
import classification_rank3_mixed_component_review as E
start=time.monotonic();xx=s.symbols('h x1:5');env=dict(zip(('h','x','y','z','w'),xx));envn=dict(zip(map(str,xx),xx));name=sys.argv[1]
source=json.load(open('reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json'));row=next(c for c in source['cases']if c['case']==name)
ns=json.load(open('reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_NORMALS_20261107.json'));normal=next(c for c in ns['cases']if c['case']==name)['normals'][0]
meta=next(c for c in json.load(open('reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_CUBIC_SUPPORT_20261107.json'))['cases']if c['case']==name)
def decode(data,en):
 A=s.zeros(5)
 for ij,f in data.items():
  i,j=map(int,ij.strip('()').split(','));A[i,j]=s.sympify(f,locals=en);A[j,i]=-A[i,j]
 return A
P=decode(row['point_brackets'],env)/s.Integer(row['integer_scale']);R=decode(normal,envn)
cs=[s.Rational(a)for a in R if a!=0 for a in s.Poly(a,*xx).coeffs()];g=s.igcd(*[int(a)for a in cs]);R=R/g
cols=meta['implicit_D_columns'];uv=s.symbols('v0:'+str(len(cols)));t=s.Symbol('t');variation=t*R
for j,c in enumerate(cols):
 (i,k),powers=E.COLUMNS[c];m=s.prod(a**p for a,p in zip(xx,powers));variation[i,k]+=uv[j]*m;variation[k,i]-=uv[j]*m
M=P+variation
sparse={ij:dict(s.Poly(s.expand(M[ij]),*xx).terms())for ij in E.PAIRS if M[ij]!=0}
J=E.jacobi(sparse).applyfunc(s.expand)
polys=[f for f in J if f!=0];folder=Path(sys.argv[2]) if len(sys.argv)>2 else Path('work/plane-local-support-20261107')/name;folder.mkdir(exist_ok=True)
vars=list(map(str,uv))+['t'];fmt=lambda f:str(f).replace('**','^')
body='ring R=0,('+','.join(vars)+'),ds;\nideal I='+','.join(map(fmt,polys))+';\n'
body+='print("INPUT_READY");\nideal G=std(I);\nprint("COMPLETE_STANDARD_BASIS");\nprint(dim(G));\nprint(vdim(G));\nprint(size(G));\n'
body+='write("'+str(folder/'basis.sing.txt')+'",G);\nwrite("'+str(folder/'leading.sing.txt')+'",lead(G));\nquit;\n'
(folder/'input.sing').write_text(body)
(folder/'meta.json').write_text(json.dumps(dict(case=name,columns=cols,normal_scale=str(g),point_scale=row['integer_scale'],normal=normal,polynomial_count=len(polys),seconds=time.monotonic()-start),indent=2))
print(json.dumps(dict(case=name,vars=len(vars),equations=len(polys),input_bytes=len(body),seconds=time.monotonic()-start,folder=str(folder))),flush=True)
