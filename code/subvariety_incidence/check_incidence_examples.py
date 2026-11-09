"""Exact polynomial checks used by poisson_subvariety_incidence_ja.tex.

Checks examples and finite rank witnesses, not the entire 147-component theorem.
Uses no old producers. Python 3 / SymPy only. Output is deterministic JSON.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json
from pathlib import Path
import sympy as s


def sign(word):
    return (-1)**sum(word[i]>word[j]for i in range(len(word))for j in range(i+1,len(word)))


def wedge_vectors(a,b):
    return a*b.T-b*a.T


def jacobi(P,x):
    return {ijk:s.expand(sum(P[ijk[0],l]*s.diff(P[ijk[1],ijk[2]],x[l])+
                            P[ijk[1],l]*s.diff(P[ijk[2],ijk[0]],x[l])+
                            P[ijk[2],l]*s.diff(P[ijk[0],ijk[1]],x[l]) for l in range(len(x))))
            for ijk in itertools.combinations(range(len(x)),3)}


def projective_pfaffian(P,x):
    ans=0
    for missing in range(5):
        i,j,k,l=[z for z in range(5)if z!=missing]
        ans+=(-1)**missing*x[missing]*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])
    return s.expand(ans)


def jacobian_three(F,G,H,x):
    P=s.zeros(5)
    for i,j in itertools.combinations(range(5),2):
        rest=[k for k in range(5)if k not in(i,j)]
        value=sign((i,j)+tuple(rest))*s.det(s.Matrix([[s.diff(f,x[k])for k in rest]for f in (F,G,H)]))
        P[i,j]=s.expand(value);P[j,i]=-P[i,j]
    return P


def centre_rank(P,x,degree):
    monomials=[s.prod(x[i]for i in inds)for inds in itertools.combinations_with_replacement(range(5),degree)]
    coeff=s.symbols(f'c0:{len(monomials)}'); f=sum(c*m for c,m in zip(coeff,monomials))
    eq=[]
    for out in P*s.Matrix([s.diff(f,z)for z in x]):
        eq.extend(s.Poly(s.expand(out),*x).coeffs())
    M=s.linear_eq_to_matrix(eq,coeff)[0]
    rows=list(M.T.rref()[1]);cols=list(M.rref()[1]); minor=M.extract(rows,cols).det() if rows else s.Integer(1)
    return {'dimension':len(monomials)-len(rows),'rank':len(rows),'nonzero_minor':str(minor),'row_indices':rows,'column_indices':cols}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'certificates/incidence_examples.json')
    args=p.parse_args(); checks={}
    def zero(name,exprs):
        values=list(exprs)
        for i,z in enumerate(values):
            if s.expand(z)!=0:raise AssertionError(f'{name}[{i}]: {s.expand(z)}')
        checks[name]=len(values)
    x=s.symbols('x0:5'); lam=(0,1,2,3,4)
    q0=sum(z*z for z in x);q1=sum(a*z*z for a,z in zip(lam,x));h=sum(x)
    P=jacobian_three(q0,q1,h,x)/4
    zero('CI22_Jacobi',jacobi(P,x).values());zero('CI22_projective_rank_two',[projective_pfaffian(P,x)])
    for i,f in enumerate((q0,q1,h)):zero(f'CI22_Casimir_{i}',P*s.Matrix([s.diff(f,z)for z in x]))
    centre={'CI22_degree1':centre_rank(P,x,1),'CI22_degree2':centre_rank(P,x,2)}
    assert centre['CI22_degree1']['dimension']==1 and centre['CI22_degree2']['dimension']==3
    F=sum(z**3 for z in x);P3=jacobian_three(F,x[0],x[1],x)
    zero('cubic_Jacobi',jacobi(P3,x).values());zero('cubic_projective_rank_two',[projective_pfaffian(P3,x)])
    centre['cubic_degree1']=centre_rank(P3,x,1);centre['cubic_degree3']=centre_rank(P3,x,3)
    assert centre['cubic_degree1']['dimension']==2 and centre['cubic_degree3']['dimension']==5
    u,a,b,c,v=x;q=2*u*v+a*a+b*b+c*c
    E=s.Matrix(x);U=[s.Matrix([0,u if k==0 else 0,u if k==1 else 0,u if k==2 else 0,-x[k+1]])for k in range(3)]
    beta=wedge_vectors(U[0],U[1]);T=s.zeros(5);T[4,3]=s.Rational(1,2);T[3,4]=-s.Rational(1,2)
    K=q*T-wedge_vectors(E,U[2])/2;t=s.Symbol('t');Pt=beta+t*K
    zero('quadric_family_strict_Jacobi',jacobi(Pt,x).values())
    zero('quadric_family_Casimir',Pt*s.Matrix([s.diff(q,z)for z in x]))
    zero('quadric_family_Pfaffian',[projective_pfaffian(Pt,x)+t*q*u**3/2])
    zero('quadric_family_same_restriction',list(K-q*T+wedge_vectors(E,U[2])/2))
    qaff=x[0]*x[4]-4*x[1]*x[3]+3*x[2]**2
    H=s.Matrix([2*x[0],x[1],0,-x[3],-2*x[4]])
    W=s.Matrix([0,x[0],2*x[1],3*x[2],4*x[3]])
    T0=s.zeros(5);T0[1,4]=-1;T0[4,1]=1;T0[2,3]=s.Rational(1,2);T0[3,2]=-s.Rational(1,2)
    K0=qaff*T0-wedge_vectors(E,W)/2;aff=wedge_vectors(H,W)
    Jt=jacobi(aff+t*K0,x)
    def one_wedge_two(vec,mat,ids):
        i,j,k=ids;return vec[i]*mat[j,k]-vec[j]*mat[i,k]+vec[k]*mat[i,j]
    zero('quadric_affine_full_Jacobi_pencil',[Jt[ids]+t*(t+1)*qaff*one_wedge_two(W,T0,ids)for ids in Jt])
    zero('quadric_affine_rank_four_factor',[projective_pfaffian(aff-K0,x)+s.Rational(5,2)*x[0]*qaff**2])
    aa,bb,cc=s.symbols('a b c');zz=(aa,bb,cc);dd=4*aa*cc-3*bb**2
    ss=s.symbols('s01 s02 s03 s04 s12 s13 s14 s23 s24 s34');pairs=list(itertools.combinations(range(5),2))
    Sm=s.zeros(5)
    for (i,j),value in zip(pairs,ss):Sm[i,j]=value;Sm[j,i]=-value
    Y=(Sm.T*s.Matrix([s.diff(qaff,vv)for vv in x])).subs(dict(zip(x,(1,aa,bb,cc,dd))),simultaneous=True)
    y=s.Matrix([Y[i+1]-zz[i]*Y[0]for i in range(3)])
    pb=s.zeros(3)
    for i,j,value in ((0,1,2*(bb-aa**2)),(0,2,3*(cc-aa*bb)),(1,2,6*(aa*cc-bb**2))):pb[i,j]=value;pb[j,i]=-value
    lie=[]
    for i,j in itertools.combinations(range(3),2):
        expr=s.expand(sum(y[k]*s.diff(pb[i,j],zz[k])-pb[k,j]*s.diff(y[i],zz[k])-pb[i,k]*s.diff(y[j],zz[k])for k in range(3)))
        lie.extend(s.Poly(expr,*zz).coeffs())
    stabilizer_matrix=s.linear_eq_to_matrix(lie,ss)[0]
    null=stabilizer_matrix.nullspace();assert len(null)==1 and stabilizer_matrix.rank()==9
    tcoeff=s.Matrix([T0[i,j]for i,j in pairs])
    assert stabilizer_matrix*tcoeff==s.zeros(stabilizer_matrix.rows,1)
    checks['quadric_affine_stabilizer_rank9']=1
    Y=[]
    for i in range(3):
        vec=s.zeros(5,1);vec[i]=x[i];vec[3]=-x[3];Y.append(vec)
    B=wedge_vectors(Y[0],Y[1]);E4=s.Matrix([*x[:4],0]);S=B+wedge_vectors(E4,Y[2])
    zero('hyperplane_suspension_Jacobi',jacobi(S,x).values())
    zero('hyperplane_suspension_Pfaffian',[projective_pfaffian(S,x)+4*s.prod(x)])
    out={'date':'2026-11-10','version':'v1','status':'PASS','sympy':s.__version__,
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'checks':checks,'check_count':sum(checks.values()),'low_degree_centres':centre,
         'quadric_family_divisor':'(-t/2) q u^3','hyperplane_divisor':'-4 x0 x1 x2 x3 x4',
         'scope':'Exact examples and rank witnesses. General cubic rigidity and hyperplane equivalence use written proofs; generic CI22 rigidity has its own symbolic certificate.'}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k]for k in ('status','check_count','low_degree_centres')},indent=2))

if __name__=='__main__':main()
