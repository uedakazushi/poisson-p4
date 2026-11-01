#!/usr/bin/env python3
"""T7 flag flat4 certificate: transport stage. Rebuild in supplied directory."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1]);WD.mkdir(parents=True,exist_ok=True)
import sys,json,time
from pathlib import Path
import sympy as s
sys.path.insert(0,str(ROOT/'code'))
from s61_model import S61,BASIS,X
st=time.time();m=S61();old=json.load(open(ROOT/'reports/S61_R0_MIXED_SURFACE_T7_FORMAL_20261101.json'))['universal'];H=s.Matrix(old['beta_kernel']).applyfunc(s.Rational);b=s.symbols('b0:11');a=s.Symbol('a');e=s.Symbol('e');one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==2];w=[0,1,-1,0,0];weights=[]
for k in m.certificate['raw_free_indices']:
 p,mon=BASIS[k];ex=s.Poly(mon,*X).monoms()[0];weights.append(sum(x*y for x,y in zip(ex,w))-sum(w[i]for i in p))
# H coordinates inherit a single torus weight.
bw=[]
for j in range(11):
 ww={weights[one[i]] for i in range(16)if H[i,j]};assert len(ww)==1;bw.append(ww.pop())
sub={b[j]:b[j]*a**(-bw[j])for j in range(11)}
q1=s.Matrix(json.load(open(ROOT/'reports/S61_FLAG_T7_BOUNDARY_20261101.json'))['projective_q']).applyfunc(s.sympify).subs(e,1)
oldfixed=next(i for i,v in enumerate(q1)if v);chart=[i for i in range(45)if i!=oldfixed]
coeff=[q1]+[s.zeros(45,1)for d in [1,2,3]]
for degree,key,col in [(1,'tangent',2),(2,'quadratic_correction',5),(3,'cubic_correction',9)]:
 mat=s.Matrix(old[key]).applyfunc(s.sympify);assert mat[-1,col]==0
 for j,i in enumerate(chart):coeff[degree][i]=mat[j,col]
# nu_old=a*nu and projective scale a^-3.
for degree in range(4):
 for i in range(45):
  value=s.cancel(coeff[degree][i].subs(sub,simultaneous=True)*a**(weights[two[i]]-3+degree));num,den=s.fraction(value);dp=s.Poly(den,a);assert len(dp.terms())==1
  (power,),dc=dp.terms()[0];value=0
  for (ex,),c in s.Poly(num,a).terms():
   k=ex-power;assert k%2==0,(degree,i,k);value+=c/dc*e**(-k//2)
  coeff[degree][i]=s.cancel(value)
print('transport complete',time.time()-st,flush=True)
# Change from old q24=3 to regular q25=-1 projective chart.
fixed=25;assert coeff[0][fixed]==-1
inverse=[s.S.One]
for degree in [1,2,3]:inverse.append(s.cancel(sum(coeff[k][fixed]*inverse[degree-k]for k in range(1,degree+1))))
new=[sum((coeff[k]*inverse[degree-k]for k in range(degree+1)),s.zeros(45,1)).applyfunc(s.cancel)for degree in range(4)]
assert new[0]==s.Matrix(json.load(open(ROOT/'reports/S61_FLAG_T7_BOUNDARY_20261101.json'))['projective_q']).applyfunc(s.sympify)
assert all(new[d][fixed]==0 for d in [1,2,3])
enc=lambda M:[[str(v)for v in row]for row in M.tolist()]
out=dict(variables=list(map(str,b)),beta_weights=bw,q=enc(new[0]),coefficients={str(d):enc(new[d])for d in [1,2,3]},fixed=fixed,seconds=time.time()-st)
(WD/'flag-T7-universal-transport.json').write_text(json.dumps(out,indent=2)+'\n');print('UNIVERSAL TRANSPORT PASS',time.time()-st,flush=True)
