import json,pickle,time,sys
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix as DM
st=time.time();D=Path(sys.argv[1]);c=pickle.load(open(D/'flag-true-core-flint-projection-input.pkl','rb'));q=json.load(open('reports/S61_FLAG_TRUE_CORE_QUADRATIC_20261101.json'));p=json.load(open(D/'flag-true-core-cubic-products-global.json'));v=list(map(s.Symbol,c['variables']));R,*_=ring(v,s.QQ);K=s.QQ.frac_field(*v)
def mat(key):return DM.from_dok({(i,j):K.field.new(R.from_dict({tuple(e):s.QQ(co)for e,co in z}))for i,row in enumerate(c[key])for j,z in enumerate(row)if z},(len(c[key]),len(c[key][0])),K)
def ext(A,rr,cc):return A.extract(list(rr),list(cc))
def const(M):return s.Matrix(M).applyfunc(s.sympify).to_DM().convert_to(K)
A2,A3=mat('A2'),mat('A3');rows=q['constant_minor']['rows'];cols=q['constant_minor']['columns'];mi=ext(A2,rows,cols).inv();F3=mi*ext(A3,rows,range(220));H3=A3-ext(A2,range(135),cols)*F3;T3=mat('T3')-ext(mat('W2'),[44],cols)*F3;print('all true cubic/t2 tails',time.time()-st,flush=True)
rc=const(p['constant_rref']);re=const(p['remaining_rref']);pc=p['constant_pivots'];fc=p['free'];pe=p['remaining_pivots'];fe=[j for j in range(120)if j not in pe]
def project(A):
 B=ext(A,range(A.shape[0]),fc)-ext(A,range(A.shape[0]),pc)*ext(rc,range(100),fc)
 return ext(B,range(B.shape[0]),fe)-ext(B,range(B.shape[0]),pe)*ext(re,range(14),fe)
C=project(H3);t=project(T3);print('true cubic projected',time.time()-st,flush=True)
M=C.to_Matrix();sample={x:0 for x in v};sel=list(M.subs(sample).T.to_DM().rref()[1]);cp=list(M.subs(sample)[sel,:].to_DM().rref()[1]);minor=C.extract(sel,cp);det=K.to_sympy(minor.det());print('minor determinant',s.factor(det),flush=True);assert det.is_Rational and det
RR=minor.inv()*C.extract(sel,range(106));assert(C-C.extract(range(135),list(cp))*RR).is_zero_matrix;assert(t-t.extract([0],list(cp))*RR).is_zero_matrix
assert all(s.denom(x).is_Number for x in RR.to_Matrix())
print('TRUE CUBIC RANK',len(cp),'T3 IN IDEAL PASS','seconds',time.time()-st,flush=True)
enc=lambda M:[[str(x)for x in row]for row in M.to_Matrix().tolist()]
(D/'flag-true-core-cubic-global.json').write_text(json.dumps(dict(rank=len(cp),minor_value=str(det),selected_rows=sel,pivots=list(cp),rref=enc(RR),t3_after_t2_correction_in_ideal=True,seconds=time.time()-st))+'\n')
