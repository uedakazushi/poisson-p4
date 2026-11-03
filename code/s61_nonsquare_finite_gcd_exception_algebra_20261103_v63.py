from pathlib import Path
import json,hashlib,itertools,time,sys
P=32003;w=Path('work/nonsquare-v62/rank2-finite-model-v59');out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);start=time.monotonic()
def read(name):
 ts=iter((w/name).read_text().split());n,m,k=[int(next(ts))for _ in range(3)];assert k==0;out=[]
 for i in range(n):
  row=[]
  for j in range(m):
   q={}
   for _ in range(int(next(ts))):c=int(next(ts))%P;ex=(int(next(ts)),int(next(ts)));assert c;q[ex]=c
   row.append(q)
  out.append(row)
 return out
def add(q,r,scale=1):
 z=q.copy()
 for ex,c in r.items():
  c=(z.get(ex,0)+scale*c)%P
  if c:z[ex]=c
  elif ex in z:del z[ex]
 return z
def mul(q,r):
 z={}
 for (a,b),c in q.items():
  for (d,e),f in r.items():
   ex=(a+d,b+e);z[ex]=(z.get(ex,0)+c*f)%P
 return{k:v for k,v in z.items()if v}
def shift(q,ex,c):return{(a+ex[0],b+ex[1]):v*c%P for(a,b),v in q.items()if v*c%P}
order=lambda ex:(sum(ex),-ex[1],-ex[0])
G=read('finite-G-native2.txt')[0];lift=read('finite-G-lift-native2.txt');model=read('finite-model-polys-native2.txt')[0];wit=read('finite-U-bezout-native2.txt');assert len(model)==3 and len(lift)==2 and len(wit)==3
import sympy as sym
pp,ee=sym.symbols('P e');inp=json.loads((w/'primitive-input.json').read_text());qqN=sym.sympify(inp['N'],locals={'P':pp,'e':ee});qqM=sym.sympify(inp['primitive_minor'],locals={'P':pp,'e':ee});aa=-pp*(20*pp+3*ee+20)/180;qqU=aa*pp*(500*aa+90*ee**2+3*ee-24)*((50*aa+18*ee**2-9*ee-2)*pp+(60*ee-40)*aa)*(100*aa+18*ee**2-9*ee-2)*(5*pp-6*ee-4)
def qq_reduce(q):
 out={}
 for ex,c in sym.Poly(q,pp,ee,domain=sym.QQ).terms():
  den=int(sym.denom(c));assert den%P;v=int(sym.numer(c))*pow(den,-1,P)%P
  if v:out[ex]=v
 return out
assert model==[qq_reduce(q)for q in [qqN,qqM,qqU]]
LM=[max(q,key=order)for q in G];assert all(q[lm]==1 for q,lm in zip(G,LM))
def nf(q):
 q=q.copy();rem={}
 while q:
  ex=max(q,key=order);c=q[ex]
  for g,lm in zip(G,LM):
   if all(x>=y for x,y in zip(ex,lm)):
    q=add(q,shift(g,(ex[0]-lm[0],ex[1]-lm[1]),c),-1);break
  else:rem[ex]=c;del q[ex]
 return rem
for j,g in enumerate(G):assert add(add({},mul(model[0],lift[0][j])),mul(model[1],lift[1][j]))==g
assert all(not nf(q)for q in model[:2])
spairs=0
for i,j in itertools.combinations(range(len(G)),2):
 lcm=tuple(max(a,b)for a,b in zip(LM[i],LM[j]));pair=add(shift(G[i],tuple(a-b for a,b in zip(lcm,LM[i])),1),shift(G[j],tuple(a-b for a,b in zip(lcm,LM[j])),1),-1);assert not nf(pair);spairs+=1
for i,g in enumerate(G):
 for ex in g:
  for j,lm in enumerate(LM):
   if i!=j:assert not all(a>=b for a,b in zip(ex,lm))
unit={}
for q,c in zip(model,wit):unit=add(unit,mul(q,c[0]))
assert unit=={(0,0):1}
bp=min(a for a,b in LM if b==0);be=min(b for a,b in LM if a==0);basis=sorted([(a,b)for a in range(bp)for b in range(be)if not any(a>=c and b>=d for c,d in LM)],key=lambda ex:(sum(ex),*ex));assert len(basis)==28;ix={ex:i for i,ex in enumerate(basis)}
MP=[[0]*28 for _ in basis];ME=[[0]*28 for _ in basis]
for j,ex in enumerate(basis):
 for mat,offset in [(MP,(1,0)),(ME,(0,1))]:
  rem=nf({(ex[0]+offset[0],ex[1]+offset[1]):1});assert set(rem)<=set(ix)
  for key,c in rem.items():mat[ix[key]][j]=c
uinv=nf(wit[2][0]);assert not nf(add(mul(model[2],uinv),{(0,0):1},-1))
with(out/'multiplication-matrices.txt').open('w')as f:
 f.write('28 32003\n')
 for mat in [MP,ME]:
  for row in mat:f.write(' '.join(map(str,row))+'\n')
r=dict(all_exact_GF_checks_passed=True,all_literal_QQ_model_coefficients_matched=True,all_model_denominators_prime_units=True,prime=P,source_GB_generators=len(G),all_GB_source_lift_identities_verified=True,all_original_generators_in_GB_verified=True,all_S_pairs_verified=spairs,all_GB_monic_and_reduced=True,FULLU_exact_Bezout_identity_verified=True,basis=[list(v)for v in basis],algebra_dimension=28,U_inverse_reduced_terms=[[list(ex),c]for ex,c in sorted(uinv.items())],multiplication_matrix_sha256=hashlib.sha256((out/'multiplication-matrices.txt').read_bytes()).hexdigest(),seconds=time.monotonic()-start)
(out/'algebra-verification.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
