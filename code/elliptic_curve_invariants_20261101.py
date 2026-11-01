#!/usr/bin/env python3
"""Fisher-Sadek invariants of the canonical first-jet Pfaffian curve.

Implements arXiv:1406.6953v1 section2/Theorem2.3 as exact tensor contractions.
The finite family samples support, but do not prove, a universal source-curve
identification. The input canonical Pfaffian model is verified separately.
"""
import argparse,hashlib,sys,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import sympy as s
from itertools import combinations,product
from elliptic_six_variable_fiber import normalized_jet

def invOmega(O,x):
    """c4,c6,Delta,j from an integer quadratic Omega, with FS14 normalization."""
    n=5
    W=[[[[int(s.diff(O[i,j],x[a],x[b])) for b in range(n)] for a in range(n)] for j in range(n)] for i in range(n)]
    Z=[[[[sum(W[i][r][u][a]*W[j][u][r][b]+W[i][r][u][b]*W[j][u][r][a]
               for r in range(n) for u in range(n)) for b in range(n)] for a in range(n)] for j in range(n)] for i in range(n)]
    c4=s.Rational(sum(Z[i][j][r][u]*Z[r][u][i][j] for i,j,r,u in product(range(n),repeat=4)),13440)
    T={}
    for i,j,k,a,b,c in product(range(n),repeat=6):
        T[i,j,k,a,b,c]=sum(Z[i][j][r][a]*W[r][k][b][c]+Z[i][j][r][b]*W[r][k][a][c]+Z[i][j][r][c]*W[r][k][a][b] for r in range(n))
    c6=-s.Rational(sum(v*T[r,u,t,i,j,k] for (i,j,k,r,u,t),v in T.items()),1036800)
    delta=(c4**3-c6**2)/1728
    return c4,c6,delta,s.cancel(c4**3/delta) if delta else None

def pfaff4(A):return A[0,1]*A[2,3]-A[0,2]*A[1,3]+A[0,3]*A[1,2]

def pfaffOmega(A,x):
    p=s.Matrix([(-1)**i*pfaff4(A.extract([j for j in range(5) if j!=i],[j for j in range(5) if j!=i])) for i in range(5)])
    assert (A*p).applyfunc(s.expand)==s.zeros(5,1)
    dp=[p.diff(t) for t in x];da=[A.diff(t) for t in x]
    O=s.zeros(5)
    for i,j in combinations(range(5),2):
        k,l,m=[a for a in range(5) if a not in (i,j)]
        perm=[i,j,k,l,m];sign=(-1)**sum(perm[a]>perm[b] for a in range(5) for b in range(a+1,5))
        value=sign*((dp[k].T*da[l]*dp[m])[0]+(dp[m].T*da[k]*dp[l])[0]+(dp[l].T*da[m]*dp[k])[0])
        O[i,j]=s.expand(value);O[j,i]=-O[i,j]
    return O,p

def canonical_phi(beta):
    """Phi_beta=[[star Z,a_beta],[-a_beta^T,0]], integral fifth row/column.

    Returns (coordinates,Phi,scale), where the fifth row and column were
    multiplied by scale. Invariants scale by scale^8,scale^12,scale^24.
    """
    x=s.symbols('x0:5');a,b,c,d,e=x
    pl=[a,b,c,d,e,-a]
    pairs=list(combinations(range(4),2))
    star=s.zeros(4)
    for (i,j),val in zip(pairs,[-a,-e,d,c,-b,a]):star[i,j]=val;star[j,i]=-val
    cov=s.Matrix([sum(beta[4*h+i]*pl[h] for h in range(6)) for i in range(4)])
    A=star.row_join(cov).col_join((-cov.T).row_join(s.zeros(1)))
    denominator=s.ilcm(*(s.denom(t) for t in beta))
    # Scale only the fifth row/column, preserving the represented curve.
    A[:4,4]=A[:4,4]*denominator;A[4,:4]=A[4,:4]*denominator
    return x,A,denominator

def curve_invariants_from_beta(beta):
    x,A,den=canonical_phi(beta);O,p=pfaffOmega(A,x)
    c4,c6,disc,j=invOmega(O,x)
    return c4/den**8,c6/den**12,disc/den**24,j

def hesse_matrix(aa,bb):
    """Fisher math/0610403v2 section4 u5(a,b), x1,...,x5 relabeledx0,...,x4."""
    x=s.symbols('x0:5');A=s.zeros(5)
    values=[aa*x[0],bb*x[1],-bb*x[2],-aa*x[3],aa*x[2],bb*x[3],-bb*x[4],aa*x[4],bb*x[0],aa*x[1]]
    for (i,j),v in zip(combinations(range(5),2),values):A[i,j]=v;A[j,i]=-v
    return x,A

def hesse_invariants():
    a,b=s.symbols('a b')
    D=a*b*(a**10-11*a**5*b**5-b**10)
    c4=s.expand(-s.det(s.hessian(D,(a,b)))/121)
    c6=s.expand((s.diff(D,a)*s.diff(c4,b)-s.diff(D,b)*s.diff(c4,a))/20)
    assert s.expand(c4**3-c6**2-1728*D**5)==0
    return a,b,D,c4,c6

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();a,b,D,H4,H6=hesse_invariants()
    hesse_checks=[]
    for aa,bb in [(2,-1),(1,2),(2,1),(3,-1),(1,3)]:
        x,A=hesse_matrix(aa,bb);O,_=pfaffOmega(A,x)
        c4,c6,disc,j=invOmega(O,x)
        assert c4==H4.subs({a:aa,b:bb}) and c6==H6.subs({a:aa,b:bb})
        assert disc==D.subs({a:aa,b:bb})**5
        hesse_checks.append({'a':aa,'b':bb,'c4':str(c4),'c6':str(c6),'discriminant':str(disc)})
    result=[]
    for kind in (1,2):
        for parameter in (2,3):
            for point in [(0,0,0,0),(1,2,3,4),(2,3,4,5)]:
                beta=normalized_jet(kind,s.Integer(parameter),list(map(s.Integer,point)))
                c4,c6,disc,j=curve_invariants_from_beta(beta)
                aa,bb=(parameter,-1) if kind==1 else (1,parameter)
                differential_scale=s.Rational(4 if kind==1 else 1,parameter**3)
                assert c4==H4.subs({a:aa,b:bb})*differential_scale**4
                assert c6==H6.subs({a:aa,b:bb})*differential_scale**6
                assert disc!=0
                result.append({'kind':kind,'lambda':parameter,'point':point,
                    'beta':[str(v) for v in beta],'c4':str(c4),'c6':str(c6),'discriminant':str(disc),'j':str(j)})
    # Raw tensor scaling tests the degree8/12 convention independently of c*pi.
    beta=s.Matrix([s.Rational(v) for v in result[1]['beta']])
    inv1=curve_invariants_from_beta(beta);inv3=curve_invariants_from_beta(3*beta)
    assert inv3[:3]==tuple(v*3**d for v,d in zip(inv1[:3],(8,12,24)))
    assert inv3[3]==inv1[3]
    report={'id':'ELLIPTIC_CURVE_INVARIANTS_20261101','all_passed':True,
        'field':'QQ','primary_formula':'https://arxiv.org/html/1406.6953v1#S2',
        'primary_formula_location':'cyclic Omega construction; Theorem2.3; Theorem2.1 smoothness criterion',
        'normalization_check_source':'https://arxiv.org/html/math/0610403v2',
        'normalization_check_location':'section3 discrete c4,c6,D and section4 Hesseu5',
        'hesse_normalization_checks':hesse_checks,'family_samples':result,
        'beta_scaling_degrees':[8,12,24],'sample_identifications_all_passed':True,
        'scope':['Exact invariants of the displayed canonical Pfaffian model.',
                 'The12 source samples agree with the corresponding source Hesse curve.',
                 'Finite samples alone do not prove source identification for everylambda,p.',
                 'No elliptic component W16 degree is claimed by this certificate.'],
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
             for p in (Path(__file__),ROOT/'code/elliptic_six_variable_fiber.py',ROOT/'code/elliptic_weight1_dominance.py')}}
    with args.output.open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({'all_passed':True,'source_family_samples':len(result),'hesse_model_checks':len(hesse_checks)}))

if __name__=='__main__':main()
