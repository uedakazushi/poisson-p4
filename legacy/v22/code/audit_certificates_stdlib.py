"""Independent exact verifier: Python standard library only.

Reconstruct conormal matrices directly from archived P, reconstruct necessary
E8 equations, and check every row-combination certificate without a CAS.
Rational coefficients are in Q[lam,lam^-1,(lam+1)^-1]; no numerical sampling.
"""
import ast,json,argparse,time
from fractions import Fraction as Q
from pathlib import Path
BASE = Path(__file__).resolve().parents[1]
from math import comb
from itertools import combinations
V=[f'z{i}' for i in range(60)]+['lam','t']+[f'x{i}' for i in range(5)]+[f'n{i}' for i in range(3)]
I={v:i for i,v in enumerate(V)};N=len(V);ZERO=(0,)*N;LI=I['lam'];TI=I['t'];NI=[I[f'n{i}'] for i in range(3)]
def addp(a,b):
 c=dict(a)
 for m,v in b.items():
  c[m]=c.get(m,Q(0))+v
  if not c[m]:del c[m]
 return c
def mulp(a,b):
 c={}
 for ma,va in a.items():
  for mb,vb in b.items():
   m=tuple(x+y for x,y in zip(ma,mb));c[m]=c.get(m,Q(0))+va*vb
 return {m:v for m,v in c.items() if v}
def shifts(a,l=0,b=0):
 c={}
 for m,v in a.items():
  for k in range(b+1):
   mm=list(m);mm[LI]+=l+k;mm=tuple(mm);c[mm]=c.get(mm,Q(0))+v*comb(b,k)
 return {m:v for m,v in c.items() if v}
def divplus(a):
 groups={}
 for m,c in a.items():
  p=m[LI];mm=list(m);mm[LI]=0;mm=tuple(mm);groups.setdefault(mm,{})[p]=c
 out={}
 for mm,p in groups.items():
  if sum(c*((-1)**k) for k,c in p.items()):return None
  p=dict(p);q={}
  for k in range(max(p),0,-1):
   c=p.get(k,Q(0))
   if c:q[k-1]=c;p[k-1]=p.get(k-1,Q(0))-c
  if p.get(0,Q(0)):return None
  for k,c in q.items():
   m=list(mm);m[LI]=k;out[tuple(m)]=c
 return out
def factoruni(a):
 if not a:raise ZeroDivisionError
 if any(any(e for j,e in enumerate(m) if j!=LI) for m in a):raise ValueError('non-lambda denominator')
 lo=min(m[LI] for m in a);p={}
 for m,c in a.items():
  mm=list(m);mm[LI]-=lo;p[tuple(mm)]=c
 b=0
 while len(p)>1 or next(iter(p))[LI]>0:
  p2=divplus(p)
  if p2 is None:raise ValueError('denominator outside lam and lam+1')
  p=p2;b+=1
 return p.get(ZERO,Q(0)),lo,b
class R:
 def __init__(self,num=None,den=(0,0)):
  p={m:Q(c) for m,c in (num or {}).items() if c};a,b=den
  if not p:self.p={};self.d=(0,0);return
  v=min(a,min(m[LI] for m in p));
  if v:
   p={tuple(e-v if j==LI else e for j,e in enumerate(m)):c for m,c in p.items()};a-=v
  while b:
   q=divplus(p)
   if q is None:break
   p=q;b-=1
  self.p=p;self.d=(a,b)
 @staticmethod
 def c(x):return R({ZERO:Q(x)})
 def __add__(self,o):
  if not isinstance(o,R):o=R.c(o)
  d=tuple(max(a,b) for a,b in zip(self.d,o.d));return R(addp(shifts(self.p,d[0]-self.d[0],d[1]-self.d[1]),shifts(o.p,d[0]-o.d[0],d[1]-o.d[1])),d)
 __radd__=__add__
 def __neg__(self):return R({m:-c for m,c in self.p.items()},self.d)
 def __sub__(self,o):return self+-o
 def __rsub__(self,o):return -self+o
 def __mul__(self,o):
  if not isinstance(o,R):o=R.c(o)
  return R(mulp(self.p,o.p),tuple(a+b for a,b in zip(self.d,o.d)))
 __rmul__=__mul__
 def __truediv__(self,o):
  if not isinstance(o,R):o=R.c(o)
  c,a,b=factoruni(o.p);num=shifts(self.p,*o.d);return R({m:v/c for m,v in num.items()},(self.d[0]+a,self.d[1]+b))
 def __pow__(self,k):
  if k<0:return R.c(1)/(self**(-k))
  r=R.c(1);v=self
  while k:
   if k%2:r=r*v
   v=v*v;k//=2
  return r
 def key(self):return self.d,tuple(sorted((m,c.numerator,c.denominator) for m,c in self.p.items()))
 def iszero(self):return not self.p

def var(name):
 m=list(ZERO);m[I[name]]=1;return R({tuple(m):Q(1)})
SY={v:var(v) for v in V}
CACHE={}
def parse(s,env=None):
 if s not in CACHE:CACHE[s]=ast.parse(s,mode='eval').body
 env=env or {}
 def ev(n):
  if isinstance(n,ast.Constant) and isinstance(n.value,int):return R.c(n.value)
  if isinstance(n,ast.Name):return env.get(n.id,SY[n.id])
  if isinstance(n,ast.UnaryOp):return -ev(n.operand) if isinstance(n.op,ast.USub) else ev(n.operand)
  if isinstance(n,ast.BinOp):
   a=ev(n.left)
   if isinstance(n.op,ast.Pow):
    k=ast.literal_eval(n.right);assert isinstance(k,int);return a**k
   b=ev(n.right)
   if isinstance(n.op,ast.Add):return a+b
   if isinstance(n.op,ast.Sub):return a-b
   if isinstance(n.op,ast.Mult):return a*b
   if isinstance(n.op,ast.Div):return a/b
  raise ValueError(ast.dump(n))
 return ev(CACHE[s])
def coeffs(p,ix):
 groups={}
 for m,c in p.p.items():
  k=m[ix];mm=list(m);mm[ix]=0;groups.setdefault(k,{})[tuple(mm)]=c
 return [R(v,p.d) for k,v in sorted(groups.items())]
def normal_linear(p,j):
 return R({tuple(0 if i in NI else e for i,e in enumerate(m)):c for m,c in p.p.items() if all(m[ix]==(1 if k==j else 0) for k,ix in enumerate(NI))},p.d)
def keys(es):return {e.key() for e in es if not e.iszero()}
def matmul(A,B):return [[sum((A[i][k]*B[k][j] for k in range(3)),R.c(0)) for j in range(3)] for i in range(3)]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--model',type=Path,default=(BASE/'previous_v21/data/tree_jacobi_models/22_node.json'));ap.add_argument('--data-directory',type=Path,default=BASE/'data/22');ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();root=args.data_directory
 orig=json.loads(args.model.read_text());heis=json.loads((root/'22_node_heisenberg.json').read_text());spec=json.loads((root/'22_node_spectral.json').read_text());checks=[]
 def ck(ok,label):
  if not ok:raise AssertionError(label)
  checks.append(label)
 allmins=[];alltrace=[];allcenter=[];t=SY['t'];n=[SY[f'n{i}'] for i in range(3)]
 for ci in range(2):
  xx=[R.c(1),t,t*t+n[0],n[1],n[2]] if ci==0 else [R.c(1),n[1],n[2],t,t*t+n[0]]
  env=dict(zip([f'x{i}' for i in range(5)],xx));P=[[R.c(0) for j in range(5)] for i in range(5)]
  for ij,e in orig['linear_model']['P_entries'].items():
   i,j=map(int,ij);P[i][j]=parse(e,env);P[j][i]=-P[i][j]
  p=[[P[i][j]-xx[i]*P[0][j]+xx[j]*P[0][i] for j in range(5)] for i in range(5)]
  if ci==0:hr=[p[1][2],p[1][3],p[1][4]];br=[p[3][4],p[4][2]-2*t*p[4][1],p[2][3]-2*t*p[1][3]]
  else:hr=[p[3][4],p[3][1],p[3][2]];br=[p[1][2],p[2][4]-2*t*p[2][3],p[4][1]-2*t*p[3][1]]
  H=[[normal_linear(e,j) for j in range(3)] for e in hr];B=[[normal_linear(e,j) for j in range(3)] for e in br]
  for i in range(3):
   for j in range(3):ck((H[i][j]-parse(spec['H'][ci][i][j])).iszero(),f'H{ci}{i}{j}');ck((B[i][j]-parse(heis['matrices'][ci][i][j])).iszero(),f'B{ci}{i}{j}')
  residue=t*SY[f'z{36+2*ci}']+t*t*SY[f'z{37+2*ci}']
  ck((sum((H[i][i] for i in range(3)),R.c(0))+residue).iszero(),f'trace{ci}')
  hess=P[0][1]*P[3][4]-P[0][3]*P[1][4]+P[0][4]*P[1][3] if ci==0 else P[0][1]*P[2][3]-P[0][2]*P[1][3]+P[0][3]*P[1][2]
  ck((normal_linear(hess,0)-residue*SY['z30' if ci==0 else 'z14']/5).iszero(),f'planeHessian{ci}')
  for ii in combinations(range(3),2):
   for jj in combinations(range(3),2):allmins+=coeffs(B[ii[0]][jj[0]]*B[ii[1]][jj[1]]-B[ii[0]][jj[1]]*B[ii[1]][jj[0]],TI)
  HH=matmul(H,H);alltrace+=coeffs(sum((HH[i][i] for i in range(3)),R.c(0))-Q(7,18)*residue*residue,TI)
  C=matmul(B,[[H[i][j]+(residue/2 if i==j else R.c(0)) for j in range(3)] for i in range(3)])
  for row in C:
   for e in row:allcenter+=coeffs(e,TI)
 ck(keys(allmins)==keys([parse(e) for e in heis['minors']]),'allHeisenbergCoefficients');ck(keys(alltrace)==keys([parse(e) for e in spec['trace_squared']]),'allSpectralTraceCoefficients');ck(keys(allcenter)==keys([parse(e) for e in spec['center']]),'allCenterCoefficients')
 relation_strings=['z30*z25-z35**2','z30*z21-z31*z35','z30*z26','z30*(3*z28+4*z35)-z31**2','z14*(z28+2*z35)-z25*(z21+z27)','z30*z14-z25*z31/2','z30*(z21+z27)-z31*z35']
 rel=[parse(e) for e in relation_strings]
 def monic_key(e):
  if e.iszero():return e.key()
  m=min(e.p);return (e/e.p[m]).key()
 available={monic_key(e) for e in allmins+[parse(e) for e in orig['jacobi_equations']] if not e.iszero()}
 for i,e in enumerate(rel):ck(monic_key(e) in available,f'handwrittenRelation{i}')
 ck((SY['z30']*rel[4]+SY['z25']*rel[6]-2*SY['z35']*rel[5]-SY['z30']*SY['z14']*SY['z28']).iszero(),'derive_r_zero')
 ck((rel[6]-rel[1]-SY['z30']*SY['z27']).iszero(),'derive_w_zero')
 cases=['22_node_NN_heis_both_certificate','22_node_SS_heis_both_SS1_certificate','22_node_NN_heis_both_slopes_certificate','22_node_NN_spectral_exception_certificate'];reports=[]
 for name in cases:
  d=json.loads((root/(name+'.json')).read_text());env={k:parse(v) for k,v in d['initial_substitution'].items()};source=orig['jacobi_equations']+heis['minors']
  if d['spectral_constraints_used']:source+=spec['trace_squared']+spec['center']
  expected=[parse(e,env) for e in source];start=len(checks);allowed=set()
  for ri,rr in enumerate(d['rounds']):
   inp=[parse(e) for e in rr['input_equations']];ck(keys(inp)==keys(expected),name+f':inputs{ri}')
   for k,c in enumerate(rr['certificates']):
    expr=R.c(0)
    for ix,v in c['combination'].items():
     factor=parse(v);allowed.update(a for a,b in zip(['lam','lam+1'],factor.d) if b);expr+=factor*inp[int(ix)]
    ck((expr-parse(c['target'])).iszero(),name+f':combination{ri}:{k}')
   new={k:parse(v) for k,v in rr['new_substitution'].items()}
   for k,l in enumerate(rr['linear_equations']):
    if l!='1':ck(parse(l,new).iszero(),name+f':substitution{ri}:{k}')
   expected=[parse(e,new) for e in rr['input_equations']]
  reports.append({'name':name,'checks':len(checks)-start,'denominator_factors':sorted(allowed),'passed':True});print(json.dumps(reports[-1]),flush=True)
 report={'all_passed':True,'checks':len(checks),'jet_and_handwritten_relation_checks':52,'certificates':reports,'implementation':'stdlib ast/Fraction sparse polynomial arithmetic, exact rational-function identities, no CAS, no sampling'}
 dest=args.output or root/'certificate_audit_stdlib.json';dest.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
