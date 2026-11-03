import itertools,json,resource,signal,sys,time

from pathlib import Path
import sympy as s
ROOT=Path('/workspace/poisson-p4');sys.path.insert(0,str(ROOT/'code'))
from qb_source_series_generic_helper_20261103_v63 import make_qb_source
start=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
out=ROOT/'work/v62-review/CB-actualD-cubic-syzygy-generic-probe.json'
b=s.symbols('b');K=s.QQ.frac_field(b)
from sympy.polys.matrices import DomainMatrix
T=type(K.zero)
def Q(v=0):
 if isinstance(v,T):return v
 if isinstance(v,int):return K.convert(v)
 return K.from_sympy(s.sympify(v,locals={'b':b}))
D=json.loads((ROOT/'work/v62-review/CB-cubic-generic-setup.json').read_text())
signal.alarm(90)
assert D['all_passed'];R=D['base29_rows'];psel=D['base20_indices'];j0=D['f_row'];minimum=D['minimal8_linear_b_rows'];ri=D['rho_base_index'];cc=list(map(Q,D['constant_odd_row_normalization']))
I58=[[Q(c)for c in row]for row in D['inverse58']];I87=[[Q(c)for c in row]for row in D['inverse87']]
sparse=lambda A:[[(j,c)for j,c in enumerate(row)if c]for row in A]
S58,S87=map(sparse,[I58,I87]);J=[[[Q(c)for c in row]for row in A]for A in D['source_Jacobians150_by29']]
fact={n:[make_qb_source(k,parameter=b,order=n)for k in['CB1','CB2','CA']]for n in[2,3]}
for k,phi in enumerate(fact[3]):assert phi.jacobian()==J[k]
def mv(A,v):return[sum((c*v[j]for j,c in row),Q(0))for row in A]
def norm(v):return[q if i==j0 else q-cc[i]*v[j0]for i,q in enumerate(v)]
def constraints(mode,V):
 if mode==58:
  A,B=V;an,bn=norm(A),norm(B)
  return[A[i]-B[i]for i in R]+[B[R[i]]for i in psel]+[(an[i]-bn[i])/2 for i in minimum]+[B[R[ri]]]
 A,B,C=V;an,bn,cn=map(norm,V)
 return[A[i]-C[i]for i in R]+[B[i]-C[i]for i in R]+[C[R[i]]for i in psel]+[an[i]-cn[i]for i in minimum[1:]]+[cn[j0]-bn[j0],C[R[ri]]]
cache={}
def invert(mode,direction,n):
 key=(mode,tuple(direction),n)
 if key in cache:return cache[key]
 phis=fact[n][:2]if mode==58 else fact[n];inv=S58 if mode==58 else S87
 jets=[[Q(0)]*(n+1)for _ in range(mode)]
 first=mv(inv,direction)
 for i,q in enumerate(first):jets[i][1]=q
 for degree in range(2,n+1):
  # Unknown coefficient of degree d enters the actual source exactly via its
  # ORIGINAL Jacobian; update it linearly instead of reevaluating every jet.
  active=fact[degree][:2]if mode==58 else fact[degree]
  values=[phi.evaluate([row[:degree+1]for row in jets[29*k:29*(k+1)]])for k,phi in enumerate(active)]
  residual=constraints(mode,[[row[degree]for row in A]for A in values]);correction=mv(inv,[-q for q in residual])
  for i,q in enumerate(correction):jets[i][degree]+=q
  for k,A in enumerate(values):
   dq=mv(sparse(J[k]),correction[29*k:29*(k+1)])
   for i,q in enumerate(dq):A[i][degree]+=q
 for degree in range(1,n+1):assert constraints(mode,[[row[degree]for row in A]for A in values])==(direction if degree==1 else[Q(0)]*mode)
 cache[key]=(values,jets);return values,jets
slopes=[Q(0),Q(1),Q(-1),Q(2)]
interps={}
for degree in range(4):
 V=s.Matrix([[K.to_sympy(c)**j for j in range(degree+1)]for c in slopes[:degree+1]]).inv();interps[degree]=[[Q(str(q))for q in row]for row in V.tolist()]
def homogeneous(values,degree):
 coeff=mv(sparse(interps[degree]),values[:degree+1])
 for c,value in zip(slopes,values):assert sum((a*(c**j if j else Q(1)) for j,a in enumerate(coeff)),Q(0))==value
 return{(degree-j,j):a for j,a in enumerate(coeff)if a}
def add(*P):
 A={}
 for p in P:
  for k,q in p.items():A[k]=A.get(k,Q(0))+q
 return{k:q for k,q in A.items()if q}
def scale(P,c):return{k:q*c for k,q in P.items()if q*c}
def mul(P,V,n=3):
 A={}
 for(i,j),q in P.items():
  for(k,l),c in V.items():
   if i+j+k+l<=n:A[i+k,j+l]=A.get((i+k,j+l),Q(0))+q*c
 return{k:q for k,q in A.items()if q}
def truncate(P,n):return{k:q for k,q in P.items()if sum(k)<=n}
one={(0,0):Q(1)};x={(1,0):Q(1)};rho={(0,1):Q(1)}
def compose(P,X,n=3):
 A={}
 for(i,j),q in P.items():
  term=one
  for _ in range(i):term=mul(term,X,n)
  for _ in range(j):term=mul(term,rho,n)
  A=add(A,scale(term,q))
 return truncate(A,n)
def series(V,kind,degree):
 A,B=V[:2];C=V[2]if len(V)==3 else None
 if kind=='b':return[(q-r)/2 for q,r in zip(norm([v[degree]for v in A]),norm([v[degree]for v in B]))]
 if kind=='c':return[q-r for q,r in zip(norm([v[degree]for v in C]),norm([v[degree]for v in B]))]
 return[q-r for q,r in zip(norm([v[degree]for v in A]),norm([v[degree]for v in C]))]
def checkpoint(stage,**kw):
 out.write_text(json.dumps(dict(stage=stage,seconds=time.monotonic()-start,**kw),default=str)+'\n');print(stage,time.monotonic()-start,flush=True)
# Constant linear coefficient map b_j in conductor coordinates x0..x7.
beta=[[Q(0)]*8 for _ in range(150)]
for i in range(8):
 target=[Q(0)]*58;target[49+i]=1;u=mv(S58,target)
 a=mv(sparse(J[0]),u[:29]);b=mv(sparse(J[1]),u[29:]);linear=[(q-r)/2 for q,r in zip(norm(a),norm(b))]
 for j,q in enumerate(linear):beta[j][i]=q
assert all(beta[j][0]==(1 if j==j0 else 0)for j in range(150))
Bsamples=[]
for c in slopes:
 direction=[Q(0)]*58;direction[49]=1;direction[57]=c;V,_=invert(58,direction,3);Bsamples.append({d:series(V,'b',d)for d in range(1,4)});checkpoint('conductor_samples',count=len(Bsamples),samples=Bsamples)
bpoly=[add(*(homogeneous([B[d][j]for B in Bsamples],d)for d in range(1,4)))for j in range(150)]
assert bpoly[j0]==x and all(not bpoly[j]for j in minimum[1:])
j8=next(j for j in range(150)if j not in minimum and bpoly[j].get((0,2),Q(0)))
b8=add(bpoly[j8],scale(x,-beta[j8][0]));a8=b8[0,2];a83=b8.get((0,3),Q(0))
H8=[];H0=[]
for j in range(150):
 h0=bpoly[j].get((0,2),Q(0))/a8;h1=(bpoly[j].get((0,3),Q(0))-h0*a83)/a8
 H8.append({k:q for k,q in{(0,0):h0,(0,1):h1}.items()if q})
 numerator=add(bpoly[j],scale(mul(H8[-1],b8),Q(-1)))
 assert all(i>0 for i,k in numerator)
 H0.append({(i-1,k):q for(i,k),q in numerator.items()})
checkpoint('conductor_plane',quadratic_generator_row=j8)
mixed=[]
for i in range(1,8):
 e=[Q(0)]*58;e[49+i]=1;Vi,_=invert(58,e,2);bi=series(Vi,'b',2)
 ex=e.copy();ex[49]=1;Vx,_=invert(58,ex,2);bx=series(Vx,'b',2)
 er=e.copy();er[57]=1;Vr,_=invert(58,er,2);br=series(Vr,'b',2)
 mixed.append(([bx[j]-bi[j]-bpoly[j].get((2,0),Q(0))for j in range(150)],
               [br[j]-bi[j]-bpoly[j].get((0,2),Q(0))for j in range(150)]))
Hi=[]
for j in range(150):
 row=[]
 for i,(mx,mr)in enumerate(mixed,1):
  a=H8[j].get((0,0),Q(0));row.append({k:q for k,q in{(0,0):beta[j][i],(1,0):mx[j]-a*mx[j8],(0,1):mr[j]-a*mr[j8]}.items()if q})
 Hi.append(row)
checkpoint('mixed_second_derivatives',quadratic_generator_row=j8)
CA=[]
for c in slopes:
 direction=[Q(0)]*87;direction[85]=1;direction[86]=c;V,_=invert(87,direction,3);CA.append({kind:{d:series(V,kind,d)for d in range(1,4)}for kind in['b','c','z']})
bc=[add(*(homogeneous([A['b'][d][j]for A in CA],d)for d in range(1,4)))for j in range(150)]
cp=[add(*(homogeneous([A['c'][d+1][j]for A in CA],d)for d in range(3)))for j in range(150)]
zp=[add(*(homogeneous([A['z'][d][j]for A in CA],d)for d in range(1,4)))for j in range(150)]
assert cp[j0]==one and all(not p.get((0,0),Q(0))for j,p in enumerate(cp)if j!=j0)
assert all(not zp[j]for j in minimum[1:])
c8=add(cp[j8],*(scale(cp[j],-beta[j8][i])for i,j in enumerate(minimum)))
b8ca=add(bc[j8],*(scale(bc[j],-beta[j8][i])for i,j in enumerate(minimum)))
normal_cubic=[];all_d=[]
for j in range(150):
 h0=compose(H0[j],bc[j0],2);his=[compose(h,bc[j0],1)for h in Hi[j]]
 predicted=add(mul(h0,bc[j0]),mul(H8[j],b8ca),*(mul(h,bc[k])for h,k in zip(his,minimum[1:])))
 assert add(bc[j],scale(predicted,Q(-1)))=={},('b-syzygy',j)
 d=add(cp[j],scale(h0,Q(-1)),scale(mul(H8[j],c8,2),Q(-1)),*(scale(mul(h,cp[k],2),Q(-1))for h,k in zip(his,minimum[1:])))
 assert not any(sum(k)<2 for k in d),('linear-d',j,d)
 all_d.append(d);normal_cubic.append([d.get((2,0),Q(0)),d.get((1,1),Q(0)),d.get((0,2),Q(0))])
Z=s.Matrix([[K.to_sympy(p.get((2,0),Q(0))),K.to_sympy(p.get((1,1),Q(0))),K.to_sympy(p.get((0,2),Q(0)))]for p in zp]);rows=list(DomainMatrix.from_Matrix(Z.T).convert_to(K).rref()[1]);assert len(rows)==2
qb=Z[rows,:];ell=DomainMatrix.from_Matrix(qb).convert_to(K).nullspace().to_Matrix()[0,:].T;assert ell[0]==0 and ell[1]!=0
ef=[Q(q)for q in ell];values=[sum((q*c for q,c in zip(row,ef)),Q(0))/ef[1]for row in normal_cubic]
witness=next((j for j,q in enumerate(values)if q!=0),None)
def poly(p):return[[i,j,str(q)]for(i,j),q in sorted(p.items())]
result=dict(all_passed=True,seconds=time.monotonic()-start,all150_full_b_syzygy_identities_through3_checked=True,
 all150_zero_constant_and_linear_defects_checked=True,odd_nonbase_coordinate_count=121,
 quadratic_generator_row=j8,base_normalized_minimal9_rows=[*minimum,j8],actual_Z4_quadratic_rows=rows,actual_Z4_quadratic_matrix=[[str(q)for q in row]for row in Z.tolist()],
 all150_syzygy_defect_quadratics=[poly(p)for p in all_d],all150_socle_coefficients=list(map(str,values)),
 nonzero_cubic_socle_witness_row=witness,nonzero_cubic_socle_witness_value=str(values[witness])if witness is not None else None,
 conductor_plane_b150=[poly(p)for p in bpoly],CA_plane_b150=[poly(p)for p in bc],CA_plane_c150=[poly(p)for p in cp],
 division_h0=[poly(p)for p in H0],division_h8=[poly(p)for p in H8],division_h_linear=[[poly(p)for p in row]for row in Hi])
out.write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps({k:result[k]for k in['all_passed','quadratic_generator_row','actual_Z4_quadratic_rows','nonzero_cubic_socle_witness_row','nonzero_cubic_socle_witness_value','seconds']},indent=2),flush=True)
