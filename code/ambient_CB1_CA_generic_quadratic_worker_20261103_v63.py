import gzip,itertools,json,signal,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path('/workspace/poisson-p4');start=time.monotonic();signal.alarm(120)
resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
out=ROOT/'work/v62-review/CB1-CA-two-normal-quadratic-generic-probe.json'
def rd(p):
 with gzip.open(ROOT/p,'rt') as f:return json.load(f)
h,x,y,z,w,t,b=s.symbols('h x y z w t b');X=s.Matrix([h,x,y,z,w])
K=s.QQ.frac_field(b)
piv=lambda A:list(DomainMatrix.from_Matrix(A).convert_to(K).rref()[1])
rank=lambda A:len(piv(A))
invmat=lambda A:DomainMatrix.from_Matrix(A).convert_to(K).inv().to_Matrix()
pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[u]*X[v]) for i,j in pairs for u,v in mons])
def canonical(P):
 d=s.Matrix([sum(s.diff(P[i,j],X[i]) for i in range(5)) for j in range(5)])
 return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
def rawCA(f,g,H):
 P=s.zeros(5);P[1,2]=f;P[3,4]=g*h*h;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z);return P-P.T
def rawCB(f,g,H):
 P=s.zeros(5);P[1,2]=f;P[3,4]=g;P[2,3]=-s.diff(f,x)*z/2-s.diff(H,w);P[2,4]=-s.diff(f,x)*w/2+s.diff(H,z);return P-P.T
def lie(i,j,P):
 E=s.zeros(5);E[i,j]=1;V=E*X
 return (E*P+P*E.T-sum((V[k]*P.diff(X[k]) for k in range(5)),s.zeros(5))).applyfunc(s.expand)
def push(G,U):return (G*U.subs(dict(zip(X,G.inv()*X)),simultaneous=True)*G.T).applyfunc(s.expand)
cubics=[s.prod((h,x,z,w)[i] for i in m) for m in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in m)]
fixed=[canonical(rawCB(f,0,0)) for f in [h*h,h*x,x*x]]+[canonical(rawCB(0,0,H)) for H in cubics]
with (ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json').open() as f:globalcb=json.load(f)
flags=[q for q in itertools.product(range(5),repeat=2) if list(q) not in globalcb['parabolic_GL_positions']]
source=rd('reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/certificate.json.gz')
ey=s.eye(5)[:,2];Db=s.Matrix([x*x+8*b*h*x/5,-2*b*x*x/5,0,h*x+3*b*x*z/5,h*h+z*z-7*b*x*w/5]);Pb=ey*Db.T-Db*ey.T;P=Pb
params=[];Js=[]
for datum in source['actual_two_CB_sources']:
 parse=lambda q:s.sympify(q,locals={'b':b,'h':h,'x':x,'y':y,'z':z,'w':w})
 G=s.Matrix([[parse(q) for q in row] for row in datum['G']]);F=parse(datum['f']);H=parse(datum['H']);Q=push(G,P)
 assert (Q-canonical(rawCB(F,0,H))).applyfunc(s.expand)==s.zeros(5)
 local=s.Matrix.hstack(*map(vec,[*fixed,canonical(rawCB(0,F,0)),*(lie(i,j,Q) for i,j in flags)]))
 J=s.Matrix.hstack(*(vec(push(G.inv(),U)) for U in [*fixed,canonical(rawCB(0,F,0)),*(lie(i,j,Q) for i,j in flags)]))
 assert rank(J)==29;Js.append(J);params.append((G,F,H,Q))

JCB=Js[0]
base=rd('reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz');C=s.Matrix(base['fixed_fiber_linear_basis150_by20']);caflags=base['flag_positions'];JCA=C.row_join(s.Matrix.hstack(*(vec(lie(i,j,P)) for i,j in caflags)))
A=JCB.row_join(-JCA);pc=piv(A);assert len(pc)==36
pr=piv(A[:,pc].T);inv=invmat(A[pr,pc])
O=s.Matrix.hstack(*(vec(lie(i,j,P)) for i,j in itertools.product(range(5),repeat=2)));S=O.row_join(vec(Pb.diff(b)));sc=piv(S);assert len(sc)==20;S=S[:,sc]
def left(J):
 rows=piv(J.T);return invmat(J[rows,:]),rows
lcb,rcb=left(JCB);lca,rca=left(JCA)
KS=(lcb*S[rcb,:]).col_join(lca*S[rca,:]);assert (A*KS).applyfunc(s.cancel)==s.zeros(150,20)
nn=DomainMatrix.from_Matrix(A).convert_to(K).nullspace().to_Matrix();ker=[nn[i,:].T for i in range(nn.rows)];assert len(ker)==22
for i,k in enumerate(ker):
 scale=next(q for q in reversed(k) if q!=0);ker[i]=k.applyfunc(lambda q:s.cancel(q/scale))
normals=[];span=S
for k in ker:
 v=JCB*k[:29,:]
 if span.row_join(v).subs(b,1).rank()>span.cols:normals.append(k);span=span.row_join(v)
 if len(normals)==2:break
assert len(normals)==2
F=KS.row_join(s.eye(58)[:,pc]).row_join(s.Matrix.hstack(*normals));detF=F.subs(b,1).det();assert detF!=0
result=dict(stage='frame',seconds=time.monotonic()-start,frame58_by58=[[str(q)for q in row]for row in F.tolist()],frame_determinant=str(detF),implicit36_target_rows=pr,implicit36_source_columns=pc,implicit36_inverse=[[str(q)for q in row]for row in inv.tolist()],normal58_by2=[[str(q)for q in row]for row in s.Matrix.hstack(*normals).tolist()])
out.write_text(json.dumps(result,indent=2)+'\n');print('frame',time.monotonic()-start,flush=True)

def source2(param,u,isCB):
 G,F0,H0,Q=param;positions=flags if isCB else caflags;Agroup=s.zeros(5)
 for q,(i,j) in zip(u[20:,:],positions):Agroup[i,j]+=q
 assert Agroup**5==s.zeros(5)
 E=s.eye(5)+t*Agroup+t*t*Agroup*Agroup/2;Ei=s.eye(5)-t*Agroup+t*t*Agroup*Agroup/2
 if isCB:
  dF=sum(q*m for q,m in zip(u[:3,:],[h*h,h*x,x*x]));dH=sum(q*m for q,m in zip(u[3:19,:],cubics));mu=t*u[19]
  U=canonical(rawCB(F0+t*dF,mu*(F0+t*dF),H0+t*dH))
 else:
  delta=s.zeros(5);vv=C*u[:20,:]
  for k,(i,j) in enumerate(pairs):
   delta[i,j]=sum(vv[15*k+m]*X[u0]*X[v0]for m,(u0,v0)in enumerate(mons));delta[j,i]=-delta[i,j]
  U=Q+t*delta
 pushed=(E*U.subs(dict(zip(X,Ei*X)),simultaneous=True)*E.T).applyfunc(s.expand)
 back=push(G.inv(),pushed)
 return [vec(back.applyfunc(lambda q:s.expand(q).coeff(t,k)))for k in range(3)]
def second(k):
 cb=source2(params[0],k[:29,:],True);ca=source2((s.eye(5),0,0,P),k[29:,:],False)
 assert all((U-V).applyfunc(s.cancel)==s.zeros(150,1)for U,V in[(cb[0],vec(P)),(ca[0],vec(P)),(cb[1],ca[1])])
 return cb[2]-ca[2]
q00,q11=map(second,normals);q01=second(normals[0]+normals[1])-q00-q11
raw=s.Matrix.hstack(q00,q01,q11);corrected=(raw-A[:,pc]*inv*raw[pr,:]).applyfunc(s.cancel);assert corrected[pr,:]==s.zeros(36,3)
rows=piv(corrected.T);rank=len(rows);print('rank',rank,'rows',rows,flush=True)
r,z0=s.symbols('r z');forms=[s.expand(corrected[i,0]*r*r+corrected[i,1]*r*z0+corrected[i,2]*z0*z0)for i in rows]
cubic=[]
for i in rows:
 aa,bb,cc=corrected[i,:];cubic.extend([[aa,bb,cc,0],[0,aa,bb,cc]])
M=s.Matrix(cubic);cp=piv(M.T);cubicrank=len(cp)
unit=s.factor(K.to_sympy(DomainMatrix.from_Matrix(M[cp,:]).convert_to(K).det()))if cubicrank==4 and len(cp)==4 else None
result.update(all_passed=True,stage='finished',seconds=time.monotonic()-start,quadratic_obstruction_rank=rank,obstruction_basis_rows=rows,obstruction_forms=list(map(str,forms)),degree3_multiple_rank=cubicrank,degree3_unit_minor=str(unit),raw150_by3=[[str(q)for q in row]for row in raw.tolist()],corrected150_by3=[[str(q)for q in row]for row in corrected.tolist()],normal_ratio_speeds=[str(k[19])for k in normals])
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in['quadratic_obstruction_rank','obstruction_forms','degree3_multiple_rank','degree3_unit_minor','normal_ratio_speeds','seconds']},indent=2),flush=True)
