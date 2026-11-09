"""Independent first-order coordinate-change audit of the E8 jet models.

The generator uses affine component formulas and derivatives. This audit starts
from homogeneous rational coordinate functions and transforms P by their full
Jacobians, then extracts coefficients of formal normal displacements. It does
not import the generator. All comparisons are symbolic over QQ[parameters,t].
"""
from pathlib import Path
import argparse
import itertools
import json
import time
import sympy as s

def run(source_dir,jet_dir,out_path):
 reports=[]
 total=0
 for path in sorted(jet_dir.glob('*.json')):
  data=json.loads(path.read_text())
  if 'components' not in data:continue
  before=time.monotonic()
  model=json.loads((source_dir/path.name).read_text())['linear_model']
  x=s.symbols('x0:5'); t=s.Symbol('t'); u=s.symbols('u1:4')
  z=s.symbols('z0:'+str(model['linear_dimension']))
  local={str(a):a for a in x+z+(t,)}
  P=s.zeros(5)
  for ij,f in model['P_entries'].items():
   i,j=map(int,ij);P[i,j]=s.sympify(f,locals=local);P[j,i]=-P[i,j]
  checks=[]
  for k,component in enumerate(model['components']):
   recorded=data['components'][k]
   r=[s.sympify(v,locals=local) for v in component['parametrization']]
   pivot=next(i for i,ri in enumerate(r) if ri==1)
   tangent=next(i for i,ri in enumerate(r) if ri==t)
   normals=[i for i in range(5) if i not in (pivot,tangent)]
   point={x[pivot]:s.Integer(1),x[tangent]:t}
   point.update({x[i]:r[i]+u[j] for j,i in enumerate(normals)})
   # Homogeneous degree-zero curvilinear coordinates on the original cone.
   coord=[x[tangent]/x[pivot]]
   coord.extend(x[i]/x[pivot]-r[i].subs(t,x[tangent]/x[pivot]) for i in normals)
   jac=s.Matrix([[s.diff(f,a) for a in x] for f in coord])
   # Evaluation before multiplication keeps rational intermediate size small.
   J=jac.subs(point,simultaneous=True)
   PP=P.subs(point,simultaneous=True)
   first={}
   for i,j in itertools.combinations(range(4),2):
    transformed=s.expand(sum(J[i,a]*PP[a,b]*J[j,b] for a in range(5) for b in range(5)))
    poly=s.Poly(transformed,*u)
    c0=poly.coeff_monomial(1)
    checks.append(c0==0)
    cs=[s.expand(poly.coeff_monomial(v)) for v in u]
    first[(i,j)]=cs;first[(j,i)]=[-v for v in cs]
   H=s.Matrix(3,3,lambda i,j:first[(0,j+1)][i])
   B=s.Matrix(3,3,lambda i,j:first[[(2,3),(3,1),(1,2)][j]][i])
   Hr=s.Matrix([[s.sympify(v,locals=local) for v in row] for row in recorded['H']])
   Br=s.Matrix([[s.sympify(v,locals=local) for v in row] for row in recorded['B']])
   checks.extend(s.expand(v)==0 for v in H-Hr)
   checks.extend(s.expand(v)==0 for v in B-Br)
   checks.extend(s.expand(v)==0 for v in B-B.T)
   R=s.sympify(recorded['R'],locals=local)
   checks.append(s.expand(s.trace(H)+R)==0)
  assert all(checks), (path.name,[i for i,v in enumerate(checks) if not v])
  total+=len(checks)
  report={'name':data['name'],'components':len(model['components']),'checks':len(checks),'all_passed':True,'seconds':round(time.monotonic()-before,2)}
  reports.append(report);print(json.dumps(report),flush=True)
 out={'method':'Full homogeneous rational-coordinate Jacobian transformation, then symbolic normal Taylor coefficients; generator not imported',
      'models':reports,'total_components':sum(r['components'] for r in reports),'total_checks':total,'all_passed':all(r['all_passed'] for r in reports)}
 out_path.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'total_components':out['total_components'],'total_checks':total,'all_passed':out['all_passed']}),flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('source_directory',type=Path);ap.add_argument('jet_directory',type=Path);ap.add_argument('output',type=Path)
 a=ap.parse_args();run(a.source_directory,a.jet_directory,a.output)
