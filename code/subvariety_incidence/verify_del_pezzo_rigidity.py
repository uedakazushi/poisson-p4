"""Exact certificate for generic extension rigidity of a smooth CI(2,2) surface.

The symbolic quadratic row relations are checked over Q(u,v), not sampled.
Only the two nonvanishing rank witnesses are specialized, to prove a nonempty
Zariski open set. No claim about every smooth surface/section is made.
Run from the bundle root: python code/subvariety_incidence/verify_del_pezzo_rigidity.py
Requires Python 3 and SymPy. Produces a small deterministic JSON certificate.
"""
import argparse, hashlib, itertools, json
from pathlib import Path
import sympy as s
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[2]/'certificates/del_pezzo_generic_rigidity.json')
args=parser.parse_args()
x=s.symbols('x0:5'); u,v=s.symbols('u v'); h=s.symbols('h0:5'); c=s.symbols('a0:10')+s.symbols('b0:10')
pairs=list(itertools.combinations(range(5),2)); triples=list(itertools.combinations(range(5),3)); quads=list(itertools.combinations(range(5),4)); xm=[x[i]*x[j]for i in range(5)for j in range(i,5)]; cm=[c[i]*c[j]for i in range(20)for j in range(i,20)];lam=[0,1,2,u,v]
def skew(vals):
 M=s.zeros(5)
 for (i,j),z in zip(pairs,vals):M[i,j]=z;M[j,i]=-z
 return M
A=skew(c[:10]);B=skew(c[10:]);q=[sum(z*z for z in x),sum(l*z*z for l,z in zip(lam,x))];P=s.zeros(5)
for i,j in pairs:
 for m in range(5):
  if m in (i,j):continue
  k,l=[z for z in range(5)if z not in(i,j,m)];word=(i,j,k,l,m);eps=(-1)**sum(word[a]>word[b]for a in range(5)for b in range(a+1,5))
  P[i,j]+=eps*h[m]*(lam[l]-lam[k])*x[k]*x[l]
 P[j,i]=-P[i,j]
def wedge(v,C,ij):
 i,j,k=ij;return v[i]*C[j,k]-v[j]*C[i,k]+v[k]*C[i,j]
def lin(C,ij):
 i,j,k=ij;return sum(C[i,l]*s.diff(P[j,k],x[l])+C[j,l]*s.diff(P[k,i],x[l])+C[k,l]*s.diff(P[i,j],x[l])for l in range(5))
rows=[]; labels=[]
for C in (A,B):
 w0=C*s.Matrix([s.diff(q[0],z)for z in x]);w1=C*s.Matrix([s.diff(q[1],z)for z in x])
 T={ijk:s.expand(lin(C,ijk)+wedge(w0,A,ijk)+wedge(w1,B,ijk))for ijk in triples}
 for ijkl in quads:
  pol=s.Poly(s.expand(sum((-1)**k*x[ijkl[k]]*T[tuple(i for i in ijkl if i!=ijkl[k])]for k in range(4))),*x)
  for mon in xm:rows.append(s.Poly(pol.coeff_monomial(mon),*c)); labels.append((0 if C==A else 1,ijkl,str(mon)))
assert all(p.total_degree() <= 2 and p.coeff_monomial(1) == 0 for p in rows)
L=s.Matrix([[p.coeff_monomial(z)for z in c]for p in rows]); Q=s.Matrix([[p.coeff_monomial(z)for z in cm]for p in rows]);Q0=Q.subs({u:3,v:4});L0=L.subs({u:3,v:4,**{z:1 for z in h}})
# Parameter monomials split by the (Z/2)^5 sign character.
groups={}
for j,mon in enumerate(cm):
 if not any(Q[:,j]):continue
 inds=[]
 for k,cc in enumerate(c):
  inds.extend(list(pairs[k%10])*s.degree(mon,cc))
 char=tuple(i for i in range(5)if inds.count(i)%2)
 groups.setdefault(char,[]).append(j)
for i in range(150):
    active=[char for char,cols in groups.items() if any(Q[i,j] for j in cols)]
    assert len(active)<=1, ('character block overlap',i,active)
N=[];block_records=[]; pivots=[]
for char,cols in sorted(groups.items()):
 ri=[i for i in range(Q.rows)if any(Q[i,j]for j in cols)]; sub=Q.extract(ri,cols); sub0=Q0.extract(ri,cols)
 cp=list(sub0.rref()[1]);rp=list(sub0.T.rref()[1]); assert len(cp)==len(rp)
 M=sub.extract(rp,cp);det=s.factor(M.det()); inv=M.inv(method='DM')
 for i in range(sub.rows):
  if i in rp:continue
  coeff=(sub.extract([i],cp)*inv).applyfunc(s.cancel)
  residual=sub[i,:]-coeff*sub.extract(rp,list(range(sub.cols)))
  assert all(s.cancel(z)==0 for z in residual),(char,i,residual)
  row=s.zeros(1,150);row[0,ri[i]]=1
  for k,r in enumerate(rp):row[0,ri[r]]=-coeff[k]
  N.append(row)
 block_records.append({'character':char,'rows':len(ri),'columns':len(cols),'rank':len(rp),'pivot_rows':[ri[z]for z in rp],'pivot_columns':[cols[z]for z in cp],'pivot_determinant':str(det)})
 print('block',char,len(ri),len(cols),len(rp),det,flush=True)
# any all-zero Q row provides another syzygy.
for i in range(150):
 if not any(Q[i,:]):
  row=s.zeros(1,150);row[i]=1;N.append(row)
N=s.Matrix.vstack(*N);print('syzygies',N.shape,flush=True)
N0=N.subs({u:3,v:4});R0=N0*L0;rp=list(R0.T.rref()[1]);cp=list(R0.rref()[1]);print('pure linear rank',len(rp),flush=True);assert len(rp)==18
M0=R0.extract(rp,cp);dep=cp;free=[i for i in range(20)if i not in dep];K0=-M0.inv()*R0.extract(rp,free); Y=s.zeros(20,2)
for k,i in enumerate(dep):Y[i,:]=K0[k,:]
for k,i in enumerate(free):Y[i,k]=1
z=s.symbols('z0:2'); submap=dict(zip(c,Y*s.Matrix(z)))
polys=[s.expand(p.as_expr().subs({u:3,v:4,**{a:1 for a in h}}).subs(submap,simultaneous=True))for p in rows]
M5=s.Matrix([[s.Poly(p,*z).coeff_monomial(m)for m in (z[0],z[1],z[0]**2,z[0]*z[1],z[1]**2)]for p in polys]);r5=list(M5.T.rref()[1]);assert len(r5)==5
print('reduced coefficient rank',M5.rank(),'minor',M5.extract(r5,range(5)).det(),flush=True)
assert N.shape == (90,150)
assert M0.det() == s.Rational(70875,88)
assert M5.extract(r5,range(5)).det() == s.Rational(4443337729552,7688671875)
record={'status':'PASS','normalization':'lambda=(0,1,2,u,v), Jacobian scaled by 1/4','sample':{'u':3,'v':4,'h':[1]*5},'quadratic_rows':150,'quadratic_parameter_columns':210,'block_certificates':block_records,'generic_quadratic_rank':sum(b['rank']for b in block_records),'linear_syzygy_count':N.rows,'pure_linear_witness':{'rank':18,'rows':rp,'columns':cp,'determinant':str(M0.det()),'free_columns':free},'post_elimination_witness':{'rank':5,'rows':r5,'determinant':str(M5.extract(r5,range(5)).det()),'matrix':[[str(a)for a in M5.row(i)]for i in r5]},'meaning':'These exact identities and nonzero minors certify a nonempty Zariski open set of pairs (smooth CI22 surface, nonzero Poisson section) with a unique reduced ambient Poisson extension. They do not certify every pair.'}
record['date']='2026-11-10'; record['version']='v1'; record['sympy']=s.__version__
record['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(record,indent=2)+'\n')
print('DONE',flush=True)
