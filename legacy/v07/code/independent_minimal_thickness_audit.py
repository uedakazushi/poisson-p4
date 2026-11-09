#!/usr/bin/env python3
"""Independent rational audit of the v7 minimal-nilpotent thickness-three claim.

No Poisson helper from any earlier investigation is imported.  Polynomial
Schouten brackets, Euler projection, and tangent matrices are assembled here.
All coefficients remain exact SymPy rationals.
"""
from v7_paths import ROOT, DATA, REFERENCE
from pathlib import Path
from itertools import combinations, combinations_with_replacement
from collections import defaultdict
import ast, hashlib, json, sys
import sympy as s
from sympy.polys.matrices import DomainMatrix

HERE=Path(__file__).resolve().parent

def rank(M):
    return len(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1])

class Calculus:
    def __init__(self,n):
        self.n=n; self.x=s.symbols('x0:'+str(n))
        self.pairs=list(combinations(range(n),2))
        self.triples=list(combinations(range(n),3))
        self.quads=list(combinations(range(n),4))
        self.m2=self.monomials(2); self.m3=self.monomials(3); self.m4=self.monomials(4)
        self.basis=[(I,m) for I in self.pairs for m in self.m2]
        self.bindex={z:i for i,z in enumerate(self.basis)}
        self.cubic_basis=[(I,m) for I in self.triples for m in self.m3]
        self.cindex={z:i for i,z in enumerate(self.cubic_basis)}
        self.projective_basis=[(I,m) for I in self.quads for m in self.m4]
        self.pindex={z:i for i,z in enumerate(self.projective_basis)}
    def monomials(self,d):
        return [tuple(v.count(i) for i in range(self.n)) for v in combinations_with_replacement(range(self.n),d)]
    def expr(self,f):return s.expand(sum(a*s.prod(x**k for x,k in zip(self.x,m)) for m,a in f.items()))
    def pol(self,f):return {m:s.Rational(a) for m,a in s.Poly(f,*self.x).terms() if a}
    def tensor(self,entries):return {I:self.pol(f) for I,f in entries.items() if f!=0}
    def as_matrix(self,P):
        M=s.zeros(self.n)
        for (i,j),f in P.items():M[i,j]=self.expr(f);M[j,i]=-M[i,j]
        return M
    def from_vector(self,v):
        P={}
        for (I,m),a in zip(self.basis,v):
            if a:P.setdefault(I,{})[m]=a
        return P
    def vector(self,P):
        return s.SparseMatrix(len(self.basis),1,{(self.bindex[I,m],0):a for I,f in P.items() for m,a in f.items() if a})
    def cubic_vector(self,C):
        return s.SparseMatrix(len(self.cubic_basis),1,{(self.cindex[I,m],0):a for I,f in C.items() for m,a in f.items() if a})
    def euler_vector(self,C):
        out=defaultdict(int)
        for I,f in C.items():
            for j in range(self.n):
                if j in I:continue
                sg=(-1)**sum(i<j for i in I);K=tuple(sorted((j,)+I))
                for m,a in f.items():
                    t=list(m);t[j]+=1;out[self.pindex[K,tuple(t)],0]+=sg*a
        return s.SparseMatrix(len(self.projective_basis),1,{k:v for k,v in out.items() if v})
    @staticmethod
    def entry(P,i,j):
        if i==j:return {}
        if i<j:return P.get((i,j),{})
        return {m:-a for m,a in P.get((j,i),{}).items()}
    @staticmethod
    def diff(f,j):
        out={}
        for m,a in f.items():
            if m[j]:
                t=list(m);t[j]-=1;out[tuple(t)]=a*m[j]
        return out
    @staticmethod
    def product(f,g):
        out=defaultdict(int)
        for m,a in f.items():
            for k,b in g.items():out[tuple(u+v for u,v in zip(m,k))]+=a*b
        return {m:a for m,a in out.items() if a}
    @staticmethod
    def add(out,f,mult=1):
        for m,a in f.items():out[m]+=mult*a
    def polar(self,P,R):
        out={}
        for I in self.triples:
            q=defaultdict(int)
            for i,j,k in [I,(I[1],I[2],I[0]),(I[2],I[0],I[1])]:
                for l in range(self.n):
                    self.add(q,self.product(self.entry(P,i,l),self.diff(self.entry(R,j,k),l)))
                    self.add(q,self.product(self.entry(R,i,l),self.diff(self.entry(P,j,k),l)))
            out[I]={m:a for m,a in q.items() if a}
        return out
    def jacobi(self,P):return {I:{m:a/2 for m,a in f.items()} for I,f in self.polar(P,P).items()}
    def derivative(self,P,projective=False):
        vf=self.euler_vector if projective else self.cubic_vector
        return s.Matrix.hstack(*[vf(self.polar(P,{I:{m:s.S.One}})) for I,m in self.basis])
    def lie(self,P,X):
        M=self.as_matrix(P);X=s.Matrix(X);A=X.jacobian(self.x)
        L=M.applyfunc(lambda q:sum(X[i]*s.diff(q,self.x[i]) for i in range(self.n)))-A*M-M*A.T
        return self.tensor({I:s.expand(L[I]) for I in self.pairs})
    def vertical(self,i,j):
        return self.tensor({(a,b):(self.x[a]*self.x[j] if b==i else 0)-(self.x[b]*self.x[j] if a==i else 0) for a,b in self.pairs})
    def lift(self,P,target):
        assert self.n+1==target.n
        return {(i+1,j+1):{(0,)+m:a for m,a in f.items()} for (i,j),f in P.items()}
    def encode(self,P):return {str(I):str(s.factor(self.expr(P.get(I,{})))) for I in self.pairs}

def constant_compatibility(C,J):
    # [J,R] is linear: use its natural16 coordinates, without degree-padding.
    cols=[]
    for I,m in C.basis:
        jj=C.polar(J,{I:{m:s.S.One}})
        cols.append(s.Matrix([jj.get(K,{}).get(tuple(int(i==j) for i in range(C.n)),0) for K in C.triples for j in range(C.n)]))
    M=s.Matrix.hstack(*cols)
    return M,s.Matrix.hstack(*M.nullspace())

def symplectic_orbit(C,J,Q):
    JM=C.as_matrix(J)
    H=[s.prod(x**a for x,a in zip(C.x,m)) for m in C.m2]
    return [C.lie(Q,-JM*s.Matrix([s.diff(h,x) for x in C.x])) for h in H]

def read_tensor(C,data):
    # The source labels x1..x4, while this affine calculus labels x0..x3.
    variables={f'x{i+1}':C.x[i] for i in range(4)}
    return C.tensor({tuple(ast.literal_eval(I)):s.sympify(f,locals=variables) for I,f in data.items()})

def main():
    source=REFERENCE/'minimal_thickness.json'
    record=json.loads(source.read_text());C=Calculus(4);P=Calculus(5)
    J=C.tensor({(0,1):1,(2,3):1});Q,R,T=[read_tensor(C,record[k]) for k in ('Q','R','T')]
    compatibility,B=constant_compatibility(C,J)
    assert rank(compatibility)==15 and rank(B)==45
    for A in [Q,R,T]:assert compatibility*C.vector(A)==s.zeros(16,1)
    assert not any(C.cubic_vector(C.jacobi(Q)))
    D=C.derivative(Q)*B
    assert not any(C.cubic_vector(C.polar(Q,R)))
    assert C.cubic_vector(C.polar(Q,T))+C.cubic_vector(C.jacobi(R))==s.zeros(80,1)
    O=symplectic_orbit(C,J,Q);actual=O+[Q]
    F=s.Matrix.hstack(*[C.vector(A) for A in actual]);K=F.row_join(C.vector(R))
    A=s.Matrix.hstack(D,*[C.cubic_vector(C.polar(R,Z)) for Z in actual+[R]])
    obstruction=C.cubic_vector(C.polar(R,T));A3=A.row_join(obstruction)
    qq={'compatibility_rank':rank(compatibility),'compatible_dimension':rank(B),
        'Jacobi_rank':rank(D),'Sp4_orbit_rank':rank(F[:,:10]),'actual_Sp4_times_scale_rank':rank(F),
        'actual_plus_R_rank':rank(K),'ambiguity_rank':rank(A),'third_rank':rank(A3)}
    expected=[34,9,10,11,34,35]
    assert [qq[k] for k in ['Jacobi_rank','Sp4_orbit_rank','actual_Sp4_times_scale_rank','actual_plus_R_rank','ambiguity_rank','third_rank']]==expected
    print('Q audit',qq,flush=True)
    LQ=C.lift(Q,P);LR=C.lift(R,P);LT=C.lift(T,P)
    for I in [(1,2),(3,4)]:LQ.setdefault(I,{})[(2,0,0,0,0)]=1
    DP=P.derivative(LQ,projective=True)
    VV=[P.vertical(i,j) for i in range(5) for j in range(5)]
    OP=[]
    for i in range(5):
        for j in range(5):
            X=s.zeros(5,1);X[i]=P.x[j];OP.append(P.lie(LQ,X))
    V=s.Matrix.hstack(*[P.vector(Z) for Z in VV]);FP=s.Matrix.hstack(V,*[P.vector(Z) for Z in OP+[LQ]])
    KP=FP.row_join(P.vector(LR))
    assert DP*KP==s.zeros(350,KP.cols)
    assert P.euler_vector(P.polar(LQ,LT))+P.euler_vector(P.jacobi(LR))==s.zeros(350,1)
    AP=s.Matrix.hstack(DP,*[P.euler_vector(P.polar(LR,Z)) for Z in VV+OP+[LQ,LR]])
    AP3=AP.row_join(P.euler_vector(P.polar(LR,LT)))
    pp={'Jacobi_rank':rank(DP),'vertical_rank':rank(V),'actual_vertical_orbit_scale_rank':rank(FP),
        'actual_plus_R_rank':rank(KP),'ambiguity_rank':rank(AP),'third_rank':rank(AP3)}
    assert list(pp.values())==[101,24,48,49,101,102]
    PM=P.as_matrix(LQ);pf=s.factor(PM[1,2]*PM[3,4]-PM[1,3]*PM[2,4]+PM[1,4]*PM[2,3])
    x0,x,y,z,w=P.x
    assert s.expand(pf-(x0**4+3*x0*x0*x*x+12*x**3*z))==0
    report={'independent_implementation':True,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'source':str(source),'rational_coefficients_preserved':True,
            'Q':C.encode(Q),'R':C.encode(R),'T':C.encode(T),'quadratic':qq,'P4':pp,
            'all_first_and_second_jet_identities_verified':True,'Pfaffian':str(pf),
            'formal_conclusion':{'quadratic_reduced_component_dimension':10,'full_projective_reduced_component_dimension':23,
                                 'generic_transverse_ring':'K[epsilon]/(epsilon^3)','generic_length':3},'all_passed':True}
    (DATA/'independent_minimal_thickness_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print('P4 audit',pp,flush=True);print('PASS',flush=True)

if __name__=='__main__':main()
