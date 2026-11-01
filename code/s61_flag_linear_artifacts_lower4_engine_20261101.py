import json,pickle,sys,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from sympy.polys.rings import ring
D=Path(sys.argv[1]);st=time.time();c=pickle.load(open(D/'flag-eight-flint-projection-input.pkl','rb'));qr=json.load(open(D/'flag-b1-universal-quadratic.json'));ms2=list(cwr(range(9),2));ms3=list(cwr(range(9),3));active=[4,5,7,8];rad=[0,1,2,3,6];pairs=list(cwr(active,2));actcols=[ms2.index(p)for p in pairs];rr=qr['pivot_rows']
if len(c['variables'])==7:
 for key in ['A2','W2','A3','T3','A4_surviving','T4_surviving']:
  for row in c[key]:
   for j,p in enumerate(row):row[j]=[[[0]+e,v]for e,v in p]
 c['variables']=['b2']+c['variables'];pickle.dump(c,open(D/'flag-eight-flint-projection-input.pkl','wb'))
v=list(map(s.Symbol,c['variables']));R,*_=ring(v,s.QQ);K=s.QQ.frac_field(*v)
def mat(key):return DM.from_dok({(i,j):K.field.new(R.from_dict({tuple(e):s.QQ(co)for e,co in p}))for i,row in enumerate(c[key])for j,p in enumerate(row)if p},(len(c[key]),len(c[key][0])),K)
A2,A3=mat('A2'),mat('A3');A=A2.extract(rr,actcols).to_Matrix();assert not any(x.free_symbols for x in A);qc=s.Matrix(qr['quadratic_rref']).applyfunc(s.sympify);mu=s.Matrix([-qc[i,44]for i in range(9)]+[1]);z=s.zeros(1,10);z[0,9]=1;pa=[actcols.index(i)for i in qr['pivot_columns']];Q=s.eye(10)-mu*z;inv=Q[:,pa]*A[:,pa].inv();assert A*mu==s.zeros(9,1)and A*inv==s.eye(9)and inv*A==Q
enc=lambda M:[[str(x)for x in row]for row in M.tolist()]
# All-active cubic lifts over QQ.
ac3=list(cwr(active,3));B=s.zeros(36,20)
for i in range(9):
 for j,a in enumerate(active):
  for k,m in enumerate(pairs):B[4*i+j,ac3.index(tuple(sorted(m+(a,))))]+=A[i,k]
sel=list(B.T.to_DM().rref()[1]);assert len(sel)==20;X=s.zeros(20,36);bi=B[sel,:].inv()
for j,col in enumerate(sel):X[:,col]=bi[:,j]
assert X*B==s.eye(20)
# Socle-radical classes, with no h inverse.
P={}
for(i,j),co in A3.to_dok().items():
 m=ms3[j];aa=tuple(a for a in m if a in active);bb=[a for a in m if a not in active]
 if len(aa)==2:
  key=(i,rad.index(bb[0]));P[key]=P.get(key,K.zero)+co*K.from_sympy(mu[pairs.index(aa)])
P=DM.from_dok(P,(135,5),K);C2=A2.extract(range(135),actcols)*inv.to_DM().convert_to(K);SH=(P-C2*P.extract(rr,range(5))).to_Matrix();h=s.Symbol('b8');pure=[i for i in range(135)if any(SH[i,:])and set().union(*(x.free_symbols for x in SH[i,:]))<={h}];S0=SH[pure,:].subs(h,0);rowpos=list(S0.T.to_DM().rref()[1]);pc=list(S0.to_DM().rref()[1]);assert len(rowpos)==len(pc)==3;sr=[pure[i]for i in rowpos];M=SH[sr,pc];det=s.factor(M.det());assert det.is_Rational and det;mi=M.inv().applyfunc(s.cancel);LR=(mi*SH[sr,:]).applyfunc(s.expand);assert(SH-SH[:,pc]*LR).applyfunc(s.expand)==s.zeros(135,5);L=s.zeros(3,135);
for j,col in enumerate(sr):L[:,col]=mi[:,j]
vn=s.symbols('v:9');ell=[str(sum(LR[i,j]*vn[rad[j]]for j in range(5)))for i in range(3)];free=[j for j in range(5)if j not in pc];subs={rad[j]:[rad[j],'1']for j in free}
for i,j in enumerate(pc):
 nonzero=[a for a in free if LR[i,a]];assert len(nonzero)<=1
 subs[rad[j]]=[rad[nonzero[0]],str(-LR[i,nonzero[0]])]if nonzero else[rad[free[-1]],'0']
q=dict(mu=enc(mu),socle=enc(z),basis=enc(A),right_inverse=enc(inv),row_contents=['1']*9,selected_raw_rows=rr,inverse_denominator='1',cubic_pivot_radical=pc,canonical_cubic_linear=ell,radical_substitutions=subs,free_radical=[rad[i]for i in free]);(D/'flag-global-low-quad.json').write_text(json.dumps(q)+'\n');(D/'flag-global-low-active-cubics.json').write_text(json.dumps(dict(lift_numerator=enc(X),denominator='1'))+'\n');(D/'flag-global-cubic-low-inverse.json').write_text(json.dumps(dict(left_inverse=enc(L),primitive_denominator='1'))+'\n');print('CONSTANT QUAD/CUBIC SPLIT PASS',D,ell,'seconds',time.time()-st,flush=True)
