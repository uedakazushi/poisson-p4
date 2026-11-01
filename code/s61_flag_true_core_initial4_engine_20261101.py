import json,pickle,time,sys
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from sympy.polys.rings import ring
st=time.time();wd=Path(sys.argv[1])
c=pickle.load(open(wd/'flag-true-core-fourth-flint-projection-input.pkl','rb'))
q=json.load(open('reports/S61_FLAG_TRUE_CORE_QUADRATIC_20261101.json'))
proof3=json.load(open('reports/S61_FLAG_TRUE_CORE_CUBIC_20261101.json'));p3=proof3['quadratic_product_space'];t3=proof3['cubic_space']
p4=json.load(open(wd/'flag-true-core-fourth-products-global.json'))
rawq=json.load(open(wd/'flag-b1-universal-quadratic.json'))
v=list(map(s.Symbol,c['variables']));R,*_=ring(v,s.QQ);K=s.QQ.poly_ring(*v)
ms={d:list(cwr(range(10),d))for d in [1,2,3,4]}
def ext(A,rr,cc):return A.extract(list(rr),list(cc))
def const(A):return s.Matrix(A).applyfunc(s.sympify).to_DM().convert_to(K)
def mat(key):
 return DM.from_dok({(i,j):K.ring.from_dict({tuple(e):s.QQ(co)for e,co in z})for i,row in enumerate(c[key])for j,z in enumerate(row)if z},(len(c[key]),len(c[key][0])),K)
def embed(A,columns,width):return DM.from_dok({(i,columns[j]):co for(i,j),co in A.to_dok().items()},(A.shape[0],width),K)
def projector(data,n):
 pc=data['constant_pivots'];fc=data['free'];pe=data['remaining_pivots'];fe=[j for j in range(len(fc))if j not in pe]
 RC,RE=const(data['constant_rref']),const(data['remaining_rref'])
 P1=DM.from_dok({(j,k):K.one for k,j in enumerate(fc)},(n,len(fc)),K)
 D=ext(RC,range(len(pc)),fc)
 P1=P1-DM.from_dok({(pc[i],j):co for(i,j),co in D.to_dok().items()},P1.shape,K)
 P2=DM.from_dok({(j,k):K.one for k,j in enumerate(fe)},(len(fc),len(fe)),K)
 D=ext(RE,range(len(pe)),fe)
 P2=P2-DM.from_dok({(pe[i],j):co for(i,j),co in D.to_dok().items()},P2.shape,K)
 return P1*P2,[fc[j]for j in fe]
P4,free4=projector(p4,715);P3,free3=projector(p3,220)
assert len(free4)==211 and len(free3)==106
active={4,5,7,8};assert ext(P4,[j for j,m in enumerate(ms[4])if sum(i in active for i in m)>=3],range(211)).is_zero_matrix
print('projection ready',time.time()-st,flush=True)
def multiples(A,d,project=False):
 ix={m:j for j,m in enumerate(ms[d+1])};out={}
 for(i,j),co in A.to_dok().items():
  for k in range(10):out[i*10+k,ix[tuple(sorted(ms[d][j]+(k,)))]]=co
 B=DM.from_dok(out,(A.shape[0]*10,len(ix)),K)
 return B*P4 if project else B
A2,A3=mat('A2'),mat('A3');A4=mat('A4_surviving')*ext(P4,c['A4_columns'],range(211))
traw4=mat('T4_surviving')*ext(P4,c['A4_columns'],range(211))
rows,cols=q['constant_minor']['rows'],q['constant_minor']['columns']
mi=ext(A2,rows,cols).to_field().inv().convert_to(K)
Q=const(q['polynomial_quadratic_basis']);F3=mi*ext(A3,rows,range(220));F4=mi*ext(A4,rows,range(211))
H3=A3-ext(A2,range(135),cols)*F3;H4=A4-ext(A2,range(135),cols)*F4
tc3=mat('T3')-ext(mat('W2'),[44],cols)*F3
tc4=traw4-ext(mat('W2'),[44],cols)*F4
G3=multiples(Q,2);G4=multiples(F3,3,True)
print('raw fourth/lower tails ready',time.time()-st,flush=True)
# Lift the constant product basis through the actual Q+F3+F4 relations.
rq=s.Matrix(rawq['quadratic_rref']).applyfunc(s.sympify)
C=rq[[i for i in range(15)if not any(x.free_symbols for x in rq[i,:])],:]
LC=const(C[:,cols]);assert LC*Q==const(C)
L=DM.from_dok({(i*10+k,j*10+k):co for(i,j),co in LC.to_dok().items()for k in range(10)},(130,150),K)
BC=L*G3;BCm=BC.to_Matrix();assert not any(x.free_symbols for x in BCm)
selected=list(BCm.T.to_DM().rref()[1]);pc=p3['constant_pivots'];fc=p3['free'];pe=p3['remaining_pivots']
bci=ext(BC,selected,pc).to_field().inv().convert_to(K)
RC=const(p3['constant_rref']);assert bci*ext(BC,selected,range(220))==RC
RT=bci*ext(L*G4,selected,range(211))
erows=[p3['remaining_rows'][i]for i in p3['minor']['rows']]
E3=ext(G3,range(150),fc)-ext(G3,range(150),pc)*ext(RC,range(100),fc)
ei=ext(E3,erows,pe).to_field().inv().convert_to(K)
RE=const(p3['remaining_rref']);assert ei*ext(E3,erows,range(120))==RE
ET=ei*(ext(G4,erows,range(211))-ext(G3,erows,pc)*RT)
B3=RC.vstack(embed(RE,fc,220));BT=RT.vstack(ET)
def decomp(A):
 first=ext(A,range(A.shape[0]),pc)
 reduced=ext(A,range(A.shape[0]),fc)-first*ext(RC,range(100),fc)
 second=ext(reduced,range(A.shape[0]),pe)
 return first.hstack(second),A*P3
GG,zero=decomp(G3);assert zero.is_zero_matrix and GG*B3==G3
GH,CH=decomp(H3);C3=const(t3['rref']);cp=t3['pivots'];sr=t3['selected_rows']
CT=ext(CH,range(135),cp);CMi=ext(CT,sr,range(5)).to_field().inv().convert_to(K)
TrueB=embed(C3,free3,220);assert GH*B3+CT*TrueB==H3
assert multiples(B3,3,True).is_zero_matrix and multiples(TrueB,3,True).is_zero_matrix
HC4=H4-GH*BT;TrueT=CMi*ext(HC4,sr,range(211))
extra=(G4-GG*BT).vstack(HC4-CT*TrueT)
GT,ST=decomp(tc3);ct=ext(ST,[0],cp);assert GT*B3+ct*TrueB==tc3
t4=tc4-GT*BT-ct*TrueT
print('ALL 285 TRUE FOURTH CANDIDATES',time.time()-st,flush=True)
enc=lambda A:[[str(x)for x in row]for row in A.to_Matrix().tolist()]
(wd/'flag-true-core-initial4-candidates.json').write_text(json.dumps(dict(variables=c['variables'],rows=enc(extra),corrected_t4=enc(t4),free_columns=free4,seconds=time.time()-st))+'\n')
M=extra.to_Matrix();M0=M.subs({x:0 for x in v});selected=list(M0.T.to_DM().rref()[1]);piv=list(M0[selected,:].to_DM().rref()[1]);print('zero candidate rank',len(selected),flush=True)
minor=ext(extra,selected,piv).to_field();det=s.factor(minor.domain.to_sympy(minor.det()));assert det.is_Rational and det;print('six constant minor verified',flush=True)
RR=minor.inv()*ext(extra,selected,range(211)).to_field();assert(extra.to_field()-ext(extra,range(285),piv).to_field()*RR).is_zero_matrix
tres=t4.to_field()-ext(t4,[0],piv).to_field()*RR
assert all(s.denom(x).is_Number for x in RR.to_Matrix());constant_t=[(j,str(x))for j,x in enumerate(tres.to_Matrix())if x and not x.free_symbols];assert constant_t
print('TRUE FOURTH RANK',len(piv),'corrected t nonzero',not tres.is_zero_matrix,'seconds',time.time()-st,flush=True)
(wd/'flag-true-core-initial4-result.json').write_text(json.dumps(dict(rank=len(piv),t_order=4,constant_t_witness=constant_t[0],selected_rows=selected,pivots=piv,minor_value=str(det),rref=enc(RR),corrected_t4_remainder=enc(tres),seconds=time.time()-st))+'\n')
