import gzip,itertools,json,signal,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path('/workspace/poisson-p4');start=time.monotonic();signal.alarm(120)
resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
out=ROOT/'work/v62-review/CB-cubic-generic-setup.json'
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

pointsetup=json.loads((ROOT/'work/v62-review/CB-cubic-source-inversion-setup.json').read_text());R=pointsetup['base29_rows'];psel=pointsetup['base20_indices'];frow=pointsetup['f_row'];minimum=pointsetup['minimal8_linear_b_rows'];ri=pointsetup['rho_base_index']
def dm(A):return DomainMatrix.from_Matrix(A).convert_to(K)
def mm(A,B):return dm(A).matmul(dm(B)).to_Matrix()
J1,J2=Js;Jall=[J1,J2,JCA];V=[invmat(J[R,:])for J in Jall];Phi=[mm(J,v)for J,v in zip(Jall,V)]
B1,B2,BA=Phi;Nlin=(BA-B2).applyfunc(s.cancel);col=next(i for i,q in enumerate(Nlin[frow,:])if q!=0);cc=[s.cancel(Nlin[j,col]/Nlin[frow,col])for j in range(150)]
def transform(J):
 A=J.copy()
 for j in range(150):
  if j!=frow:A[j,:]=(J[j,:]-cc[j]*J[frow,:]).applyfunc(s.cancel)
 return A
TJB1,TJB2,TJA=map(transform,Jall);TB1,TB2,TBA=map(transform,Phi)
Psel=s.eye(29)[psel,:];xlin=(TB1-TB2)/2
Mx=Psel.col_join(xlin[minimum,:]).col_join(s.eye(29)[ri,:]);Iz=invmat(Mx)
B58=s.zeros(29,58);B58[:,29:]=s.eye(29);B58[20:28,:29]=-TB1[minimum,:]/2
Y58=mm(Iz,B58);E58=s.zeros(29,58);E58[:,:29]=s.eye(29)
I58=mm(V[0],Y58+E58).col_join(mm(V[1],Y58))
N58=J1[R,:].row_join(-J2[R,:]);N58=N58.col_join(s.zeros(20,29).row_join(J2[R,:][psel,:])).col_join((TJB1[minimum,:]/2).row_join(-TJB2[minimum,:]/2)).col_join(s.zeros(1,29).row_join(J2[R[ri],:]))
assert dm(N58).matmul(dm(I58))==dm(s.eye(58))
Mca=Psel.col_join((TB1-TBA)[minimum[1:],:]).col_join((TBA-TB2)[frow,:]).col_join(s.eye(29)[ri,:]);Ica=invmat(Mca)
B87=s.zeros(29,87);B87[:,58:]=s.eye(29);B87[20:27,:29]=-TB1[minimum[1:],:];B87[27,29:58]=TB2[frow,:]
Y87=mm(Ica,B87);E1=s.zeros(29,87);E1[:,:29]=s.eye(29);E2=s.zeros(29,87);E2[:,29:58]=s.eye(29)
I87=mm(V[0],Y87+E1).col_join(mm(V[1],Y87+E2)).col_join(mm(V[2],Y87))
N87=J1[R,:].row_join(s.zeros(29)).row_join(-JCA[R,:]);N87=N87.col_join(s.zeros(29).row_join(J2[R,:]).row_join(-JCA[R,:]))
N87=N87.col_join(s.zeros(20,58).row_join(JCA[R,:][psel,:])).col_join(TJB1[minimum[1:],:].row_join(s.zeros(7,29)).row_join(-TJA[minimum[1:],:])).col_join(s.zeros(1,29).row_join(-TJB2[frow,:]).row_join(TJA[frow,:])).col_join(s.zeros(1,58).row_join(JCA[R[ri],:]))
assert dm(N87).matmul(dm(I87))==dm(s.eye(87))
for key,A in [('inverse58',I58),('inverse87',I87)]:assert(A.subs(b,1)-s.Matrix(pointsetup[key])).applyfunc(s.cancel)==s.zeros(*A.shape)
mat=lambda A:[[str(q)for q in row]for row in A.tolist()]
result=dict(all_passed=True,field='QQ(b)',seconds=time.monotonic()-start,base29_rows=R,base20_indices=psel,f_row=frow,constant_odd_row_normalization=list(map(str,cc)),minimal8_linear_b_rows=minimum,rho_base_index=ri,normalizer_source_order=['CB1','CB2','CA'],Jacobian58=mat(N58),inverse58=mat(I58),Jacobian87=mat(N87),inverse87=mat(I87),source_Jacobians150_by29=[mat(J)for J in Jall])
out.write_text(json.dumps(result,separators=(',',':'))+'\n');print('GENERIC 58/87 SETUP PASS',time.monotonic()-start,flush=True)
