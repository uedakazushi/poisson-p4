#!/usr/bin/env python3
"""Exact full quadratic obstruction and nonreduced tangent-cone certificate.

No network, cached matrices, or imports from previous producers. All 150
ambient linear columns and all 465 pairs of tangent directions are rebuilt.
The frame agrees with check_local_dimension28.py. Nonlinear normal directions
are eliminated only to order two; no claim about the full local ideal is made.
Run: python code/subvariety_incidence/check_isotropic_quadratic_cone.py --output FILE
"""
import argparse
import itertools as it
import json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

x=s.symbols('u0 u1 z w0 w1'); E=s.Matrix(x);u0,u1,z,w0,w1=x
pairs=list(it.combinations(range(5),2));triples=list(it.combinations(range(5),3));quads=list(it.combinations(range(5),4))
def mons(d):
 out=[]
 for t in it.combinations_with_replacement(range(5),d):
  ex=[0]*5
  for i in t:ex[i]+=1
  out.append(tuple(ex))
 return out
m2,m4=mons(2),mons(4)
rows=[(I,e) for I in quads for e in m4];ridx={a:i for i,a in enumerate(rows)}
def sign(t): return (-1)**sum(t[i]>t[j] for i in range(len(t)) for j in range(i+1,len(t)))
def bivector(q,xi):
 gam={(1,2):xi[0],(0,2):-xi[1],(0,1):xi[2]};B=s.zeros(5)
 for i,j in pairs:
  rest=tuple(k for k in range(5) if k not in(i,j))
  B[i,j]=s.expand(sign((i,j,*rest))*sum((-1)**r*s.diff(q,x[k])*gam.get(tuple(rest[t]for t in range(3)if t!=r),0)for r,k in enumerate(rest)))
  B[j,i]=-B[i,j]
 return B
q=u0**2+2*u1*w0+2*z*w1;xi=s.Matrix([u0,u1,-2*z]);P=bivector(q,xi);G=s.hessian(q,E)/2
Ms=[]
for i,j in it.product(range(5),repeat=2):
 A=s.zeros(5);A[i,j]=1;Ms.append(A)
def lie(M,B):
 v=M*E
 return s.Matrix([[s.expand(sum(v[k]*s.diff(B[i,j],x[k])-M[i,k]*B[k,j]-M[j,k]*B[i,k] for k in range(5)))for j in range(5)]for i in range(5)])
def gauge(M):
 v=M*E;return E*v.T-v*E.T
A1=s.zeros(5);A1[0,1]=1;A1[3,0]=-1
A2=s.zeros(5);A2[3,2]=1;A2[4,1]=-1
A3=s.diag(0,1,-1,-1,1)
def K(A):
 C=A*E; assert A.T*G+G*A==s.zeros(5)
 return (q*A*G.inv()/2+(E*C.T-C*E.T)/2).applyfunc(s.expand)
Mindices=[0,1,2,3,4,5,6,7,8,9,10,11,13,14,16,17,19,20,22,23]
Bs=[lie(Ms[i],P)for i in Mindices]
Bs += [bivector(f,xi)for f in [u0*w1,w0*w0,w0*w1,w1*w1]]
Bs += [bivector(q,s.Matrix([u1,0,0])),bivector(q,s.Matrix([0,u0,0])),bivector(q,s.Matrix([u0,0,-z]))]
Bs += [K(A) for A in [A1,A2,A3]]
def poly(f): return {e:F(int(c.p),int(c.q))for e,c in s.Poly(f,*x).terms()if c}
def sparse(B):return [[poly(B[i,j])for j in range(5)]for i in range(5)]
def diff(p,k):
 out={}
 for e,c in p.items():
  if e[k]:
   ee=list(e);ee[k]-=1;out[tuple(ee)]=c*e[k]
 return out
def addmul(out,p1,p2,scale=1):
 for e1,c1 in p1.items():
  for e2,c2 in p2.items():out[tuple(i+j for i,j in zip(e1,e2))]+=scale*c1*c2
def der(B):return [[[diff(B[i][j],k)for k in range(5)]for j in range(5)]for i in range(5)]
def cross(B,C,DB=None,DC=None):
 DB=DB or der(B);DC=DC or der(C);J={}
 for i,j,k in triples:
  v=defaultdict(F)
  for u,vv,w in [(i,j,k),(j,k,i),(k,i,j)]:
   for l in range(5):
    addmul(v,B[u][l],DC[vv][w][l]);addmul(v,C[u][l],DB[vv][w][l])
  J[(i,j,k)]={e:c for e,c in v.items()if c}
 out=defaultdict(F)
 for I in quads:
  for a,k in enumerate(I):
   for e,c in J[tuple(i for i in I if i!=k)].items():
    ee=list(e);ee[k]+=1;out[ridx[(I,tuple(ee))]]+=(-1)**a*c
 return {i:c for i,c in out.items()if c}
def Svec(v,n=350):return s.SparseMatrix(n,1,{(i,0):s.Rational(c.numerator,c.denominator)for i,c in v.items()})
def bcoords(B):
 return s.Matrix([s.Poly(B[i,j],*x).coeff_monomial(s.prod(x[k]**e[k]for k in range(5))) for i,j in pairs for e in m2])
def ratmatrix(M):return DomainMatrix.from_Matrix(M).convert_to(s.QQ)
def compute():
    p=sparse(P);dp=der(p)
    assert not cross(p,p,dp,dp)
    columns=[]
    for (i,j),e in it.product(pairs,m2):
     R=[[{}for _ in range(5)]for _ in range(5)];R[i][j]={e:F(1)};R[j][i]={e:F(-1)}
     columns.append(Svec(cross(p,R,dp)))
    M=s.SparseMatrix.hstack(*columns)
    
    RR,piv=ratmatrix(M.T).rref();RR=RR.to_Matrix()
    assert len(piv)==95
    
    nonp=[i for i in range(350)if i not in piv]
    residual=s.zeros(len(nonp),350)
    for k,i in enumerate(nonp):
     residual[k,i]=1
     for j,pv in enumerate(piv):residual[k,pv]=-RR[j,i]
    assert residual*M==s.zeros(255,150)
    bs=[sparse(B)for B in Bs];ds=[der(B)for B in bs]
    for b,d in zip(bs,ds):assert not cross(p,b,dp,d)
    
    C=s.Matrix.hstack(*[bcoords(B)for B in [gauge(A)for A in Ms]+[P]+Bs]);assert ratmatrix(C).rank()==55
    qcoeff=[];tmon=list(it.combinations_with_replacement(range(30),2))
    for l,(i,j)in enumerate(tmon):
     vv=cross(bs[i],bs[j],ds[i],ds[j]);w=residual*Svec(vv)
     if i==j:w=w/2
     qcoeff.append(w)
    
    Q=s.Matrix.hstack(*qcoeff); Qrr,Qpiv=ratmatrix(Q).rref();Qrr=Qrr.to_Matrix()[:len(Qpiv),:]
    
    t=s.symbols('t0:30');polys=[s.factor(sum(row[k]*t[i]*t[j]for k,(i,j)in enumerate(tmon)))for row in Qrr.tolist()]
    
    assert all(all(row[k]==0 for k,(i,j)in enumerate(tmon)if i<20)for row in Qrr.tolist())
    
    # The complete degree-two initial ideal in a minimal embedding.
    xx,yy,zz,ww,uu,vv,rr,aa,bb,cc=t[20:]
    expected=[xx*cc-zz*aa,yy*cc,cc*(zz+cc),ww*cc,
              uu*cc+rr*aa,vv*aa,vv*cc]
    assert all(s.expand(f-g)==0 for f,g in zip(polys,expected)) and len(polys)==7
    # Compact normalized witnesses. Each annihilates all 150 ambient columns.
    witness_specs=[(140,-s.Rational(1,12)),(52,s.Rational(1,9)),
                   (47,-s.Rational(10,9)),(145,-s.Rational(2,33)),
                   (110,s.Rational(1,2)),(95,-s.Rational(1,4)),
                   (88,s.Rational(1,2))]
    witnesses=[]
    for k,(row,factor) in enumerate(witness_specs):
        lam=factor*residual[row,:]
        assert lam*M==s.zeros(1,150)
        assert factor*Q[row,:]==Qrr[k,:]
        witnesses.append([{'four_indices':list(rows[j][0]),
                           'monomial_exponents':list(rows[j][1]),
                           'coefficient':str(lam[j])}
                          for j in range(350) if lam[j]])
    # Full radical decomposition and the nonzero square-zero homogeneous class.
    variables=(xx,yy,zz,ww,uu,vv,rr,aa,cc) # bb and the 20 orbit variables are free
    h=aa*(zz+cc)
    gb=s.groebner(expected,*variables,order='grlex')
    assert s.expand(gb.reduce(h)[1]-h)==0 and gb.reduce(h)[1]!=0
    assert gb.reduce(s.expand(h*h))[1]==0
    hrow=s.Matrix([[s.Poly(h,*t).coeff_monomial(t[i]*t[j])for i,j in tmon]])
    assert ratmatrix(Qrr.col_join(hrow)).rank()==8
    Ptop=[aa,cc];P26=[cc,zz,rr,vv]
    P24=[yy,ww,vv,zz+cc,xx+aa,uu*cc+rr*aa]
    def intersection(J,K):
        elim=s.Symbol('elim')
        g=s.groebner([elim*f for f in J]+[(1-elim)*f for f in K],
                     elim,*variables,order='lex')
        return [f.as_expr() for f in g.polys if not f.as_expr().has(elim)]
    radical=intersection(intersection(Ptop,P26),P24)
    radgb=s.groebner(radical,*variables,order='grlex')
    hgb=s.groebner(expected+[h],*variables,order='grlex')
    assert all(hgb.reduce(g)[1]==0 for g in radical)
    assert all(radgb.reduce(g)[1]==0 for g in expected+[h])
    # I:h = (yy,zz,ww,vv,rr,cc); computed using I intersect (h).
    colon=[]
    for f in intersection(expected,[h]):
        qout,rem=s.div(f,h,*variables)
        assert rem==0
        colon.append(qout)
    colongb=s.groebner(colon,*variables,order='grlex')
    assert list(colongb)==[yy,zz,ww,vv,rr,cc]
    # Mixed ac correction in the full ambient projective Jacobi equation.
    R=s.zeros(5)
    for (i,j),f in {
        (0,2):-(u0**2-u1*w0+2*w1*z)/2,
        (0,3):-3*w0*w1/s.Integer(2),
        (1,2):-3*u0*u1/s.Integer(4),
        (1,3):3*u0*w1/s.Integer(4),
        (2,3):3*u0*w0/s.Integer(4),
    }.items(): R[i,j]=f;R[j,i]=-f
    mixed=Svec(cross(bs[27],bs[29],ds[27],ds[29]))
    assert M*bcoords(R)+mixed==s.zeros(350,1)
    assert not cross(bs[27],bs[27],ds[27],ds[27])
    assert mixed!=s.zeros(350,1)
    
    # Independent 28-dimensional lower bound from the logarithmic source.
    qg=(u0**2+u1**2+z**2+w0**2+w1**2+u0*u1+u0*w0+u0*w1
        +2*u1*w0+u1*w1+z*w0+2*z*w1+w0*w1)
    fs=[u0,u1,z,qg]
    r1=s.Matrix([1,-1,0,0]);r2=s.Matrix([1,0,-1,0]);r3=s.Matrix([2,0,0,-1])
    residues=[{(i,j):v[i]*w[j]-v[j]*w[i]for i,j in it.combinations(range(4),2)}
              for v,w in [(r1,r2),(r1,r3),(r2,r3)]]
    cres={ij:residues[0][ij]+2*residues[1][ij]+5*residues[2][ij]for ij in residues[0]}
    def omega_source(f,c):
        out={ij:s.S.Zero for ij in pairs}
        for i,j in it.combinations(range(4),2):
            fac=s.prod(f[k]for k in range(4)if k not in (i,j))
            for h,l in pairs:
                out[h,l]+=c[i,j]*fac*(s.diff(f[i],x[h])*s.diff(f[j],x[l])
                                      -s.diff(f[i],x[l])*s.diff(f[j],x[h]))
        return {ij:s.expand(v)for ij,v in out.items()}
    def source_derivative(f,c,m,d):
        out={ij:s.S.Zero for ij in pairs}
        for i,j in it.combinations(range(4),2):
            oth=[k for k in range(4)if k not in (i,j)]
            fac=s.prod(f[k]for k in oth)
            if m in oth:
                tail=s.prod(f[k]for k in oth if k!=m)
                for h,l in pairs:
                    out[h,l]+=c[i,j]*tail*d*(s.diff(f[i],x[h])*s.diff(f[j],x[l])
                                            -s.diff(f[i],x[l])*s.diff(f[j],x[h]))
            if m in (i,j):
                for h,l in pairs:
                    dih=s.diff(d,x[h])if m==i else s.diff(f[i],x[h])
                    dil=s.diff(d,x[l])if m==i else s.diff(f[i],x[l])
                    djh=s.diff(d,x[h])if m==j else s.diff(f[j],x[h])
                    djl=s.diff(d,x[l])if m==j else s.diff(f[j],x[l])
                    out[h,l]+=c[i,j]*fac*(dih*djl-dil*djh)
        return {ij:s.expand(v)for ij,v in out.items()}
    variations=[]
    for j in range(3):
        variations.extend(source_derivative(fs,cres,j,x[i])for i in range(5)if i!=j)
    for i,j in it.combinations_with_replacement(range(5),2):
        if (i,j)!=(4,4):variations.append(source_derivative(fs,cres,3,x[i]*x[j]))
    variations.extend(omega_source(fs,cd)for cd in residues[1:])
    m3=mons(3)
    def form_vector(W):
        return s.Matrix([s.Poly(W[I],*x).coeff_monomial(s.prod(x[k]**e[k]for k in range(5)))
                         for I in pairs for e in m3])
    source_matrix=s.Matrix.hstack(*[form_vector(W)for W in variations])
    minor_rows=[2,3,4,6,7,8,9,10,11,16,19,20,21,25,26,27,28,29,30,38,39,42,43,45,46,52,53,55]
    source_det=source_matrix[minor_rows,:].det()
    assert source_det==-57177414000000
    assert ratmatrix(source_matrix).rank()==28
    assert ratmatrix(source_matrix.row_join(form_vector(omega_source(fs,cres)))).rank()==29
    # Check that the fixed point is itself in this logarithmic image.
    v1=s.Matrix([1,-1,0,0]);v2=s.Matrix([0,-4,-2,3])
    res0={(i,j):v1[i]*v2[j]-v1[j]*v2[i]for i,j in it.combinations(range(4),2)}
    omega0=omega_source([u0,u1,z,q],res0)
    gamma={(1,2):u0,(0,2):-u1,(0,1):-2*z}
    alpha={}
    for i,j,k in triples:
        alpha[i,j,k]=s.expand(s.diff(q,x[i])*gamma.get((j,k),0)
                             -s.diff(q,x[j])*gamma.get((i,k),0)
                             +s.diff(q,x[k])*gamma.get((i,j),0))
    for i,j in pairs:
        val=sum(sign((k,i,j))*x[k]*alpha[tuple(sorted((k,i,j)))]for k in range(5)if k not in(i,j))
        assert s.expand(val-omega0[i,j])==0
    # Numerical bookkeeping is exact; geometric conclusions are proved in the TeX.
    report={
        'date':'2026-11-10','version':'v12','status':'PASS','sympy':s.__version__,
        'coefficient_domain':'QQ','ambient_bivector_coefficients':150,
        'projective_Jacobi_coefficients':350,'linearization_rank':len(piv),
        'gauge_scale_and_tangent_rank':55,'projective_tangent_dimension':30,
        'tangent_basis':'20 orbit,4 q,3 xi,K1,K2,K3; indices as check_local_dimension28.py',
        'quadratic_pairs_checked':len(tmon),'complete_quadratic_relation_rank':len(Qpiv),
        'quadratic_generators':[str(f) for f in expected],
        'normalized_left_null_witnesses':witnesses,
        'radical_minimal_primes':[[str(f) for f in ideal]for ideal in (Ptop,P26,P24)],
        'quadratic_scheme_component_dimensions':[28,26,24],
        'nilradical_generator':str(h),'nilradical_square_zero':True,
        'nilradical_annihilator':[str(f)for f in colongb],
        'quadratic_rank_with_nilradical_class':8,
        'true_tangent_cone_nonreduced':True,
        'unique_top_tangent_plane':['t27','t29'],
        'local_dimension':28,'Hilbert_Samuel_multiplicity':1,
        'logarithmic_source_derivative_rank':28,'source_augmented_rank':29,
        'source_minor_rows':minor_rows,'source_minor_determinant':str(source_det),
        'mixed_ac_correction':{''.join(map(str,ij)):str(s.factor(R[ij]))for ij in pairs if R[ij]},
        'mixed_ac_corrected_modulo':'(a,c)^3 + (c^2)',
        'scope':'Complete degree-two initial ideal, not the full initial ideal. '
                'The explicit nonzero homogeneous degree-two square-zero class also '
                'survives in the true tangent cone. Ambient reducedness remains open. '
                'The local dimension 28 is independently rechecked via the logarithmic source.'
    }
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('work/isotropic_quadratic_cone.json'))
    parser.add_argument('--check',type=Path,help='Compare with an existing JSON certificate')
    args=parser.parse_args()
    if s.__version__ != '1.14.0':
        raise RuntimeError('Use the recorded SymPy 1.14.0 environment')
    report=compute()
    if args.check is not None:
        expected=json.loads(args.check.read_text(encoding='utf-8'))
        if report != expected:
            raise RuntimeError('Certificate mismatch')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('PASS: rank95; 7 complete quadratic relations; nonreduced true tangent cone; mixed ac cancellation')

if __name__=='__main__':
    main()
