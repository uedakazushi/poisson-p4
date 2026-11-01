#!/usr/bin/env python3
"""Independent polynomial audit of the elliptic Cremona first-jet identity.

Derives every derivative/cofactor from the five quadratic forms; does not read
the proposed main certificate, saved adjugate, or M2 output.
"""
import argparse,hashlib,json
from itertools import combinations,permutations
from pathlib import Path
import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring
from elliptic_weight1_dominance import elliptic_family,projective_affine

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    R,lam,*y=ring('lam,y1,y2,y3,y4',QQ);zero=R.zero;one=R.one
    def matrix(n,m,func):return [[func(i,j) for j in range(m)] for i in range(n)]
    def mul(A,B):return matrix(len(A),len(B[0]),lambda i,j:sum((A[i][k]*B[k][j] for k in range(len(B))),zero))
    def trans(A):return list(map(list,zip(*A)))
    def determinant(A):
        n=len(A);out=zero
        for perm in permutations(range(n)):
            sign=(-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
            value=one*sign
            for i,j in enumerate(perm):value*=A[i][j]
            out+=value
        return out
    def adjugate(A):
        n=len(A)
        return matrix(n,n,lambda i,j:(-1)**(i+j)*determinant([[A[a][b] for b in range(n) if b!=i] for a in range(n) if a!=j]))
    xx,source_lam,B1=elliptic_family(1);_,_,B2=elliptic_family(2)
    zz=s.symbols('y1:5');ll=s.Symbol('lam')
    affine=projective_affine(B1,xx,zz)
    K=matrix(4,4,lambda i,j:R.from_expr(s.expand((source_lam**3*affine[i,j]).subs(source_lam,ll))))
    F=K[0][1]*K[2][3]-K[0][2]*K[1][3]+K[0][3]*K[1][2]
    N=[[zero,K[2][3],-K[1][3],K[1][2]],[-K[2][3],zero,K[0][3],-K[0][2]],
       [K[1][3],-K[0][3],zero,K[0][1]],[-K[1][2],K[0][2],-K[0][1],zero]]
    assert mul(K,N)==matrix(4,4,lambda i,j:-F if i==j else zero)
    # HomogeneousQ_i=lambda*x_i²+lambda²*x_(i-1)*x_(i+1)-x_(i-2)*x_(i+2).
    Q=s.Matrix([source_lam*xx[i]**2+source_lam**2*xx[(i-1)%5]*xx[(i+1)%5]
                 -xx[(i-2)%5]*xx[(i+2)%5] for i in range(5)])
    sub=dict(zip(xx,(s.Integer(1),*zz)))|{source_lam:ll}
    DQ=Q.jacobian(xx)
    jac=matrix(5,5,lambda i,j:R.from_expr(DQ[i,j].subs(sub,simultaneous=True)))
    adj=adjugate(jac);det=determinant(jac)
    assert det==2*F
    assert mul(jac,adj)==matrix(5,5,lambda i,j:det if i==j else zero)
    Ap=[[-y[i]]+[one if i==j else zero for j in range(4)] for i in range(4)]
    G=[]
    for a in range(4):
        H=matrix(5,4,lambda i,b:R.from_expr(s.diff(Q[i],xx[a+1],xx[b+1]).subs(sub,simultaneous=True)))
        G.append([[v/2 for v in row] for row in mul(mul(Ap,adj),H)])
    trace=[sum((G[a][i][i] for i in range(4)),zero) for a in range(4)]
    checked=0
    for a in range(4):
        gt=mul(trans(G[a]),N);ng=mul(N,G[a])
        for b in range(4):
            for c in range(4):
                intrinsic=-F*N[b][c].diff(y[a])+QQ(3,5)*F.diff(y[a])*N[b][c]
                intrinsic-=QQ(1,5)*N[b][a]*F.diff(y[c])
                intrinsic+=QQ(1,5)*F.diff(y[b])*N[c][a]
                connection=-gt[b][c]-ng[b][c]+QQ(2,5)*trace[a]*N[b][c]
                connection+=QQ(1,5)*N[b][a]*trace[c]-QQ(1,5)*trace[b]*N[c][a]
                assert connection-QQ(3,2)*intrinsic==zero
                checked+=1
    # Check the projectivePoisson multiplier without assuming an exact cone lift.
    # lambda³ clears every coefficient denominator in both bracket families.
    S,*sx=ring('lam,x0,x1,x2,x3,x4',QQ);sl=sx[0];x=sx[1:]
    toS=lambda e:S.from_expr(s.expand(e.subs(source_lam,s.Symbol('lam'))))
    Qh=[toS(v) for v in Q]
    JH=[[toS(v) for v in row] for row in DQ.tolist()]
    b1=[[toS(source_lam**3*v) for v in row] for row in B1.tolist()]
    replacement=dict(zip(xx,Q))
    b2=[[toS(source_lam**3*v.subs(replacement,simultaneous=True)) for v in row] for row in B2.tolist()]
    residual=[[sum((JH[i][a]*b1[a][b]*JH[j][b] for a in range(5) for b in range(5)),S.zero)+b2[i][j]
               for j in range(5)] for i in range(5)]
    for i,j,k in combinations(range(5),3):
        assert Qh[i]*residual[j][k]-Qh[j]*residual[i][k]+Qh[k]*residual[i][j]==S.zero
    report={'id':'ELLIPTIC_CREMONA_REVIEW_20261101','all_passed':True,'field':'QQ(lambda,p)',
        'derivatives_and_adjugate_derived_from_quadrics':True,
        'determinant_identity':'det D(lambdaQ)=2 Pf(lambda³ P1)',
        'projective_poisson_wedge_identities':10,'projective_poisson_multiplier':-1,
        'primitive_connection_polynomial_identities':checked,'primitive_connection_multiplier':'3/2',
        'primitive_firstjet_pullback_multiplier':'1/2',
        'normalized_parameter_transform':'lambda unchanged,p->Q(p),c->-4c,S correction from -2F2*Df*F1^-1',
        'scope':['Universal polynomial identities, not only point samples.',
                 'Birationality of Q is a separate required input for equality of parameter degrees.',
                 'With the independent Q inverse, rho1=rho2=3000 and d40=1, these identities imply d39=1.'],
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (Path(__file__),ROOT/'code/elliptic_weight1_dominance.py')}}
    with args.output.open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({'all_passed':True,'primitive_identities':checked,'projective_poisson_identities':10}))

if __name__=='__main__':main()
