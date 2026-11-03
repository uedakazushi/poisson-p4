import sympy as s,json,time,hashlib
from pathlib import Path
root=Path('/workspace/poisson-p4');j=json.loads((root/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json').read_text());h,x,y,z,w,n0,n1,n2,r=s.symbols('h x y z w n0 n1 n2 r');X=[h,x,y,z,w];P=s.Matrix([[s.sympify(f)for f in row]for row in j['affine_P']]);raw=s.zeros(5)
for i in range(4):
 for k in range(4): raw[i+1,k+1]=s.expand(h*h*P[i,k].subs({x:x/h,y:y/h,z:z/h,w:w/h},simultaneous=True))
rows=[]
for e in [1,-1]:
 D=s.diag(1,r,r*r,r**3,e*r);D[3,0]=(r**3-e)/(6*n0);DI=D.inv();out=s.Matrix(X);back=dict(zip(X,DI*out));pushed=e/r**2*D*raw.subs(back,simultaneous=True)*D.T;target=raw.subs({n0:e*n0,n1:e*r*n1,n2:n2/r},simultaneous=True);diff=(pushed-target).applyfunc(s.cancel);assert diff==s.zeros(5);count=sum(len(s.Poly(diff[i,k],*X).terms())for i in range(5)for k in range(i+1,5));rows.append({'epsilon':e,'all_150_homogeneous_coefficients_zero':True,'z_translation':str(D[3,0]),'matrix':[[str(v)for v in row]for row in D.tolist()]})
print(json.dumps(rows,indent=2))
(root/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/transform_audit_homogeneous.json').write_text(json.dumps({'all_passed':True,'audit':'Independent homogeneous5x5 pushforward of frozen affine_P. Essential z-h shear retains constant yz bracket.','cases':rows},indent=2)+'\n')
