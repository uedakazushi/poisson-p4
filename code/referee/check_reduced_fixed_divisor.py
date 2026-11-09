#!/usr/bin/env python3
"""Independent whole-space checks for the reduced-fixed-divisor referee gate.

Uses cubic horizontal two-forms on C5, not the original150-cone adapter.
Exact QQ kernels and all pairwise Jacobi/Pfaffian identities are checked.
These three examples test the abstract proof; they do not prove the theorem.
"""
import argparse,hashlib,itertools,json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
X=s.symbols('a b c r t');PAIRS=list(itertools.combinations(range(5),2))
MON=[tuple(ij.count(k)for k in range(5))for ij in itertools.combinations_with_replacement(range(5),3)]
BASIS=list(itertools.product(PAIRS,MON))
def mon(e):return s.prod(x**n for x,n in zip(X,e))
def sign(seq):return (-1)**sum(seq[i]>seq[j]for i in range(len(seq))for j in range(i+1,len(seq)))
def kernel(F):
    f=s.Poly(F,*X);terms=f.terms();rows={}
    def add(key,col,value):
        row=rows.setdefault(key,{});row[col]=row.get(col,0)+value
    for col,((i,j),e)in enumerate(BASIS):
        for component,k,sgn in [(j,i,1),(i,j,-1)]:
            powers=list(e);powers[k]+=1;add(('horizontal',component,tuple(powers)),col,sgn)
        for k in set(range(5))-{i,j}:
            sg=sign((k,i,j));triple=tuple(sorted((k,i,j)))
            for power,c in terms:
                coefficient=sg*c*(e[k]-power[k])
                if coefficient:
                    exponent=tuple(e[l]+power[l]-int(l==k)for l in range(5))
                    add(('closed',triple,exponent),col,coefficient)
    rows=[{j:s.QQ.convert(v)for j,v in r.items()if v}for r in rows.values()];rows=[r for r in rows if r]
    dm=DomainMatrix({i:r for i,r in enumerate(rows)},(len(rows),len(BASIS)),s.QQ)
    null=dm.nullspace().to_Matrix();return null,len(rows)
def form(v):
    A=s.zeros(5)
    for value,((i,j),e)in zip(v,BASIS):A[i,j]+=value*mon(e)
    return A-A.T

def dual(A):
    B=A[1:,1:].subs(X[0],1);P=s.zeros(4)
    for (i,j),v in zip(itertools.combinations(range(4),2),[B[2,3],-B[1,3],B[1,2],B[0,3],-B[0,2],B[0,1]]):P[i,j]=v;P[j,i]=-v
    return P

def jac(P):
    x=X[1:]
    return [s.expand(sum(P[i,l]*s.diff(P[j,k],x[l])+P[j,l]*s.diff(P[k,i],x[l])+P[k,l]*s.diff(P[i,j],x[l])for l in range(4)))for i,j,k in itertools.combinations(range(4),3)]
def pf(P):return s.expand(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2])
def run(output):
    a,b,c,r,t=X
    fixtures=[('five_coordinate_hyperplanes',a*b*c*r*t,6,6),('smooth_quadric_times_elliptic_cubic_cone',(a*a+b*b+c*c+r*r+t*t)*(a**3+b**3+c**3),None,None),('irreducible_free_plane_quintic',a**5+a*a*b**3+b**4*c,10,4)]
    results=[]
    for name,F,dimension,qr in fixtures:
        null,nrows=kernel(F);n=null.rows
        if dimension is not None:assert n==dimension,(name,n)
        forms=[form(list(null.row(i)))for i in range(n)];poisson=[dual(A)for A in forms]
        x=X[1:];ff=s.expand(F.subs(a,1));point={v:1 for v in x};den=ff.subs(point);assert den!=0
        Q=s.zeros(n)
        for i in range(n):
            assert not any(jac(poisson[i]));q=s.cancel(pf(poisson[i]).subs(point)/den);assert s.expand(pf(poisson[i])-q*ff)==0;Q[i,i]=q
            for j in range(i):
                pair=poisson[i]+poisson[j];assert not any(jac(pair));qij=s.cancel(pf(pair).subs(point)/den);assert s.expand(pf(pair)-qij*ff)==0
                Q[i,j]=Q[j,i]=s.cancel((qij-Q[i,i]-Q[j,j])/2)
        rank=Q.rank()
        if qr is not None:assert rank==qr,(name,rank)
        row={'name':name,'F':str(F),'ambient_cubic_twoform_coefficients':len(BASIS),'linear_rows':nrows,'exact_kernel_dimension':n,'q_rank':rank,'all_kernel_self_and_pair_Jacobi':'PASS','all_kernel_Pfaffian_equals_scalar_F':'PASS','q_matrix':[[str(v)for v in row]for row in Q.tolist()],'kernel_basis':[{str((i,j)):str(s.expand(A[i,j]))for i,j in PAIRS if A[i,j]!=0}for A in forms]}
        results.append(row);print(name,'dimension',n,'q_rank',rank,flush=True)
    result={'scope':__doc__,'field':'QQ','examples':results,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);run(p.parse_args().output)
