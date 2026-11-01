#!/usr/bin/env python3
"""Exact ordinary-flag fourth-order core cover: initial4."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1])
import json,pickle,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix as DM
st=time.time();cache=pickle.load(open(WD/'flag-eight-flint-projected4.pkl','rb'));q=json.load(open(WD/'flag-global-low-quad.json'));act=json.load(open(WD/'flag-global-low-active-cubics.json'));tr=json.load(open(WD/'flag-global-cubic-low-inverse.json'));v=list(map(s.Symbol,cache['variables']));R,*_=ring(v,s.QQ);K=s.QQ.frac_field(*v);r=s.Symbol('b2');h=s.Symbol('b8');ms={d:list(cwr(range(9),d))for d in [2,3,4]};active=[4,5,7,8];rad=[0,1,2,3,6];pairs=list(cwr(active,2));ac3=list(cwr(active,3));rr=q['selected_raw_rows'];mu=s.Matrix(q['mu']).applyfunc(s.sympify);z=s.Matrix(q['socle']).applyfunc(s.sympify)
def ext(A,rr,cc):return A.extract(list(rr),list(cc))
def mat(key,project=False):
 out={};rows=cache[key]
 for i,row in enumerate(rows):
  for j,terms in enumerate(row):
   if not terms:continue
   p=R.from_dict({tuple(e):s.QQ(c)for e,c in terms})
   out[i,j]=K.field.new(p)
 return DM.from_dok(out,(len(rows),len(rows[0])),K)
def const(M):return s.Matrix(M).applyfunc(s.sympify).to_DM().convert_to(K)
A2,A3=mat('A2'),mat('A3');A4=mat('A4_projected',True);B2=DM.zeros((9,45),K);Braw=const(q['basis']);Bc={}
act2=[ms[2].index(p)for p in pairs]
for(i,j),c in Braw.to_dok().items():Bc[i,act2[j]]=c
B2=DM.from_dok(Bc,(9,45),K);Qinv=const(q['right_inverse']);contents=s.diag(*[1/s.sympify(c)for c in q['row_contents']]).to_DM().convert_to(K);F3=contents*ext(A3,rr,range(165));F4=contents*ext(A4,rr,range(213));C2=ext(A2,range(135),act2)*Qinv;assert C2*B2==A2
H3=A3-C2*F3;H4=A4-C2*F4;print('raw/canonical lower ready',time.time()-st,flush=True)
projection=[]
for p in cache['projection']:
 if p is None:projection.append(None);continue
 j,c=p;co=s.sympify(c);projection.append((j,K.from_sympy(co)))
def multiples(A,degree,project=False):
 ix={m:j for j,m in enumerate(ms[degree+1])};out={}
 for(i,j),c in A.to_dok().items():
  for k in range(9):
   dst=ix[tuple(sorted(ms[degree][j]+(k,)))];co=c
   if project:
    p=projection[dst]
    if p is None:continue
    dst,f=p;co*=f
   key=i*9+k,dst;out[key]=out.get(key,K.zero)+co
 return DM.from_dok(out,(A.shape[0]*9,213 if project else len(ms[degree+1])),K)
qi=Qinv.to_dok();muv=[K.from_sympy(x)for x in mu]
def decomp3(A):
 out={};soc={}
 for(i,j),c in A.to_dok().items():
  aa=tuple(a for a in ms[3][j]if a in active);bb=tuple(a for a in ms[3][j]if a not in active)
  if len(aa)==3:out[i,45+ac3.index(aa)]=c
  elif len(aa)==2:
   p=pairs.index(aa);k=rad.index(bb[0]);key=i,k;soc[key]=soc.get(key,K.zero)+c*muv[p]
   for a in range(9):
    key=i,a*5+k;co=qi.get((p,a),K.zero)
    if co:out[key]=out.get(key,K.zero)+c*co
 return DM.from_dok(out,(A.shape[0],65),K),DM.from_dok(soc,(A.shape[0],5),K)
G3=multiples(B2,2);G4=multiples(F3,3,True);radrows=[i*9+j for i in range(9)for j in rad];arows=[i*9+j for i in range(9)for j in active]
B3=ext(G3,radrows,range(165)).vstack(DM.from_dok({(i,ms[3].index(m)):K.one for i,m in enumerate(ac3)},(20,165),K))
GammaG,sg=decomp3(G3);assert sg.is_zero_matrix and GammaG*B3==G3
lift=const(act['lift_numerator']).scalarmul(K.one/K.from_sympy(s.sympify(act['denominator'])));assert lift*ext(G3,arows,[ms[3].index(m)for m in ac3])==DM.eye(20,K)
B3tail=ext(G4,radrows,range(213)).vstack(lift*ext(G4,arows,range(213)));print('global cubic product basis/tails',time.time()-st,flush=True)
GH,SH=decomp3(H3);CT=ext(SH,range(135),[1,2,3]);LC=const(tr['left_inverse']);assert LC*CT==DM.eye(3,K)
# Actual canonical socle*linear cubic basis, using the constant socle polynomial.
obs=[]
vn=s.symbols('v:9');zexpr=sum(z[0,k]*vn[a]*vn[b]for k,(a,b)in enumerate(pairs));ell=[vn[1]+(r-2)/h*vn[6],vn[2]+4*vn[6],vn[3]+2*vn[6]];dok={}
for i,l in enumerate(ell):
 p=s.Poly(s.expand(zexpr*l),*vn)
 for ex,c in p.terms():m=tuple(j for j,e in enumerate(ex)for _ in range(e));dok[i,ms[3].index(m)]=K.from_sympy(c)
TrueB=DM.from_dok(dok,(3,165),K);assert GH*B3+CT*TrueB==H3
assert multiples(B3,3,True).is_zero_matrix and multiples(TrueB,3,True).is_zero_matrix
HC4=H4-GH*B3tail;T3tail=LC*HC4;print('global true cubic basis/tails',time.time()-st,flush=True)
extra=(G4-GammaG*B3tail).vstack(HC4-CT*T3tail)
t3=mat('T3');Gt,St=decomp3(t3);ct=ext(St,[0],[1,2,3]);assert Gt*B3+ct*TrueB==t3
t4=mat('T4_projected',True)-Gt*B3tail-ct*T3tail
encp=lambda p:[[list(e),str(c)]for e,c in p.items()]
enc=lambda M:{str(i)+','+str(j):[encp(c.numer),encp(c.denom)]for(i,j),c in M.to_dok().items()if c}
den=set(q['row_contents']+[q['inverse_denominator'],act['denominator'],tr['primitive_denominator']]);den.update(str(s.denom(s.sympify(x)))for row in tr['left_inverse']for x in row if x!='0')
Path(WD/'flag-eight-initial4-field.json').write_text(json.dumps(dict(denominators=dict(global_basis=sorted(den)),projection_multiplier='b8**2',scaling_denominators=dict(ud=cache['ud'],td=cache['td']),variables=cache['variables'],shape=extra.shape,rows=enc(extra),corrected_t4=enc(t4),seconds=time.time()-st))+'\n');print('GLOBAL FOURTH CANDIDATES READY',time.time()-st,flush=True)
