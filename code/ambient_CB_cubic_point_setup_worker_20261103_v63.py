import gzip,itertools,json,signal,resource,time
from pathlib import Path
import sympy as s
ROOT=Path('/workspace/poisson-p4');start=time.monotonic();signal.alarm(120)
resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
out=ROOT/'work/v62-review/CB-cubic-source-inversion-setup.json'
def rd(p):
 with gzip.open(ROOT/p,'rt') as f:return json.load(f)
h,x,y,z,w,t,b=s.symbols('h x y z w t b');X=s.Matrix([h,x,y,z,w])
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
ey=s.eye(5)[:,2];Db=s.Matrix([x*x+8*b*h*x/5,-2*b*x*x/5,0,h*x+3*b*x*z/5,h*h+z*z-7*b*x*w/5]);Pb=ey*Db.T-Db*ey.T;P=Pb.subs(b,1)
params=[];Js=[]
for datum in source['actual_two_CB_sources']:
 parse=lambda q:s.sympify(q,locals={'b':s.Integer(1),'h':h,'x':x,'y':y,'z':z,'w':w})
 G=s.Matrix([[parse(q) for q in row] for row in datum['G']]);F=parse(datum['f']);H=parse(datum['H']);Q=push(G,P)
 assert (Q-canonical(rawCB(F,0,H))).applyfunc(s.expand)==s.zeros(5)
 local=s.Matrix.hstack(*map(vec,[*fixed,canonical(rawCB(0,F,0)),*(lie(i,j,Q) for i,j in flags)]))
 J=s.Matrix.hstack(*(vec(push(G.inv(),U)) for U in [*fixed,canonical(rawCB(0,F,0)),*(lie(i,j,Q) for i,j in flags)]))
 assert J.rank()==29;Js.append(J);params.append((G,F,H,Q))

from sympy.polys.matrices import DomainMatrix
rk=lambda A:len(DomainMatrix.from_Matrix(A).convert_to(s.QQ).rref()[1])
inv=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ).inv().to_Matrix()
base=rd('reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz');C=s.Matrix(base['fixed_fiber_linear_basis150_by20']);caflags=base['flag_positions'];JCA=C.row_join(s.Matrix.hstack(*(vec(lie(i,j,P))for i,j in caflags)))
D=json.loads((ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json').read_text());Rrows=D['common28_linear_target_rows']+[D['common29th_target_row']]
J1,J2=Js;basemaps=[J*inv(J[Rrows,:])for J in[J1,J2,JCA]];B1,B2,BA=basemaps
O=s.Matrix.hstack(*(vec(lie(i,j,P))for i,j in itertools.product(range(5),repeat=2)));S=O.row_join(vec(Pb.diff(b).subs(b,1)));pc=list(S.rref()[1]);S=S[:,pc];SB=S[Rrows,:];psel=list(SB.T.rref()[1]);assert len(psel)==20
Psel=s.eye(29)[psel,:];Nlin=BA-B2;Zlin=B1-BA;assert rk(Nlin)==1 and rk(Zlin)==7
zrows=list(Zlin.T.rref()[1]);M=Psel.col_join(Zlin[zrows,:]);assert rk(M)==27
frow=next(i for i in range(150)if rk(M.col_join(Nlin[i,:]))==28);k=next(j for j,q in enumerate(Nlin[frow,:])if q);cc=[s.cancel(Nlin[i,k]/Nlin[frow,k])for i in range(150)];assert cc[frow]==1
def transform(J):
 A=J.copy()
 for i in range(150):
  if i!=frow:A[i,:]=J[i,:]-cc[i]*J[frow,:]
 return A
TJ1,TJ2,TJA=map(transform,[J1,J2,JCA]);Tb=(transform(B1)-transform(B2))/2
minimum=[frow];M=Psel.col_join(Tb[frow,:]);assert rk(M)==21
for i in range(150):
 if i==frow:continue
 if rk(M.col_join(Tb[i,:]))>M.rows:minimum.append(i);M=M.col_join(Tb[i,:])
 if len(minimum)==8:break
assert len(minimum)==8 and rk(M)==28
rhoindex=next(i for i in range(29)if rk(M.col_join(s.eye(29)[i,:]))==29)
zero=s.zeros(29);N58=(J1[Rrows,:].row_join(-J2[Rrows,:])).col_join(s.zeros(20,29).row_join(J2[Rrows,:][psel,:])).col_join((TJ1[minimum,:]/2).row_join(-TJ2[minimum,:]/2)).col_join(s.zeros(1,29).row_join(J2[Rrows[rhoindex],:]))
assert N58.shape==(58,58)and rk(N58)==58;I58=inv(N58);assert N58*I58==s.eye(58)
N87=(J1[Rrows,:].row_join(s.zeros(29)).row_join(-JCA[Rrows,:])).col_join(s.zeros(29).row_join(J2[Rrows,:]).row_join(-JCA[Rrows,:]))
N87=N87.col_join(s.zeros(20,58).row_join(JCA[Rrows,:][psel,:]))
N87=N87.col_join(TJ1[minimum[1:],:].row_join(s.zeros(7,29)).row_join(-TJA[minimum[1:],:]))
N87=N87.col_join(s.zeros(1,29).row_join(-TJ2[frow,:]).row_join(TJA[frow,:]))
N87=N87.col_join(s.zeros(1,58).row_join(JCA[Rrows[rhoindex],:]))
assert N87.shape==(87,87)and rk(N87)==87;I87=inv(N87);assert N87*I87==s.eye(87)
mat=lambda A:[[str(q)for q in row]for row in A.tolist()]
result=dict(all_passed=True,seconds=time.monotonic()-start,center150=list(map(str,vec(P))),base29_rows=Rrows,base20_indices=psel,f_row=frow,constant_odd_row_normalization=cc and list(map(str,cc)),minimal8_linear_b_rows=minimum,rho_base_index=rhoindex,normalizer_source_order=['CB1','CB2','CA'],Jacobian58=mat(N58),inverse58=mat(I58),Jacobian87=mat(N87),inverse87=mat(I87),source_Jacobians150_by29=[mat(J)for J in[J1,J2,JCA]])
out.write_text(json.dumps(result,separators=(',',':'))+'\n');print(json.dumps({k:result[k]for k in['f_row','minimal8_linear_b_rows','rho_base_index','seconds']},indent=2),flush=True)
