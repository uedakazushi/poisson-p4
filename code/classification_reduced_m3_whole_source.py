#!/usr/bin/env python3
"""Whole-source shear reduction for remaining reduced nonconcurrent m3 symbols.

The saved complete K1 is used as a necessary space, not called Poisson.
All nine weight-two directions are retained before the original 350 equations.
No producer arithmetic or new first-symbol calculation is executed.
"""
from pathlib import Path
import hashlib,json,resource,time
from fractions import Fraction
import sympy as s
import classification_plane_m4_special_review as H

ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/plane_reduced_m3_whole_shear_source.json'
OUTPUT='reports/CLASSIFICATION_PLANE_REDUCED_M3_WHOLE_SHEAR_SOURCE_20261108.json'
INPUT='reports/CLASSIFICATION_PLANE_REDUCED_M3_NINE_FULLHIGHER_20261108.json'
x=s.symbols('a b c u v');a,b,c,u,v=x
mons=[s.prod(z**n for z,n in zip(x,e))for e in H.monomials(2)]
def scalar(z):return s.Rational(str(z))
def literal(P):
    return {ij:s.expand(sum(scalar(co)*s.prod(z**n for z,n in zip(x,e))for e,co in P.get(ij,{}).items()))for ij in H.PAIRS}
def from_vector(row):
    return {ij:sum(scalar(co)*m for ((kl,e),co,m)in zip(H.COLS,row,mons*10)if ij==kl)for ij in H.PAIRS}
def plus(*Ps):return {ij:s.expand(sum(P.get(ij,0)for P in Ps))for ij in H.PAIRS}
def jac(P):
    def at(i,j):return P.get((i,j),0)if i<j else -P.get((j,i),0)
    triples={ijk:s.expand(sum(at(i,l)*s.diff(at(j,k),x[l])for i,j,k in (ijk,(ijk[1],ijk[2],ijk[0]),(ijk[2],ijk[0],ijk[1]))for l in range(5)))for ijk in H.TRIPLES}
    return {quad:s.expand(sum((-1)**j*x[quad[j]]*triples[quad[:j]+quad[j+1:]]for j in range(4)))for quad in H.QUADS}
def jvector(J):
    pp={ij:s.Poly(poly,*x)for ij,poly in J.items()}
    return s.Matrix([pp[quad].coeff_monomial(s.prod(z**n for z,n in zip(x,e)))for quad,e in H.ROWS])
def actual(P):
    out=0
    for at in range(5):
        i,j,k,l=[t for t in range(5)if t!=at]
        out+=(-1)**at*x[at]*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])
    return s.expand(out)
def vec(P):
    pp={ij:s.Poly(poly,*x)for ij,poly in P.items()}
    return s.Matrix([pp[ij].coeff_monomial(s.prod(z**n for z,n in zip(x,e)))for ij,e in H.COLS])
def push(P,gs):
    inv={x[i]:x[i]-gs[2*i]*u-gs[2*i+1]*v for i in range(3)}
    M=s.eye(5)
    for i in range(3):M[i,3],M[i,4]=gs[2*i:2*i+2]
    T=s.zeros(5)
    for (i,j),pol in P.items():T[i,j]=pol.subs(inv,simultaneous=True);T[j,i]=-T[i,j]
    T=M*T*M.T
    return {ij:s.expand(T[ij])for ij in H.PAIRS}
def textvec(V):return [str(s.factor(z))for z in V]
def serialize(P):return {str(ij):str(s.factor(p))for ij,p in P.items()}

def main():
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(175,180));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    job=json.loads((ROOT/JOB).read_text())
    for item in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()==item['sha256'],item['path']
    assert not(ROOT/OUTPUT).exists()
    raw=json.loads((ROOT/INPUT).read_text());B1,B2=map(lambda seq:list(map(literal,seq)),H.sources())
    X=[H.a(1,i)for i in range(5)];G=[]
    for normal in (3,4):
        for i in range(3):
            Y=[{}for _ in range(5)];Y[i]=H.a(1,normal);G.append(literal(H.wedge(X,Y)))
    E=s.Matrix.hstack(*(vec(P)for P in B1+G));raw1=[i for i,(ij,e)in enumerate(H.COLS)if e[3]+e[4]-sum(k>=3 for k in ij)==1]
    Ei=E[raw1,:].inv();assert len(raw1)==36 and E[raw1,:]*Ei==s.eye(36)
    def project1(P):return Ei*vec(P)[raw1,:]
    out=[]
    for record in raw['cases']:
        if record['name']not in('M','SF_C3_UNUSED'):continue
        P0=from_vector(record['full150_P0']);K=s.Matrix(record['constrained_weight1_certificate']['complete_kernel']).applyfunc(scalar).T
        fields=[]
        for j in range(6):
            gs=[0]*6;gs[j]=1;full=project1(plus(push(P0,gs),{ij:-p for ij,p in P0.items()}));fields.append(full[:30,:])
        S=s.Matrix.hstack(*fields);assert S.rank()==6
        cols=list(fields);residual=[]
        for j in range(K.cols):
            if s.Matrix.hstack(*cols,K[:,j]).rank()>len(cols):cols.append(K[:,j]);residual.append(j)
        assert len(cols)==K.cols
        free=[j for j in range(30)if j not in record['constrained_weight1_certificate']['columns']]
        C=s.Matrix.hstack(*cols)[free,:];Ci=C.inv();assert K*s.Matrix.hstack(*cols)[free,:]==s.Matrix.hstack(*cols)
        assert K*C==s.Matrix.hstack(*cols) and C*Ci==s.eye(K.cols)
        pars=s.symbols('r0:'+str(len(residual)));p2=s.symbols('q0:9')
        P1={ij:s.expand(sum(pars[n]*sum(K[k,j]*B1[k][ij]for k in range(30))for n,j in enumerate(residual)))for ij in H.PAIRS}
        P=plus(P0,P1,{ij:sum(p2[k]*B2[k][ij]for k in range(9))for ij in H.PAIRS})
        JV=jvector(jac(P));rows,cc=record['weight2_certificate']['rows'],record['weight2_certificate']['columns'];rank=len(cc)
        A=s.Matrix(record['complete350x9']).applyfunc(scalar)
        base=JV.subs(dict.fromkeys(p2,0))
        assert JV[rows,:]==base[rows,:]+A[rows,:]*s.Matrix(p2)
        free2=[k for k in range(9)if k not in cc];sol={}
        chosen=A[rows,cc];assert chosen.det()==scalar(record['weight2_certificate']['minor_determinant'])
        values=chosen.inv()*(-base[rows,:]-A[rows,free2]*s.Matrix([p2[k]for k in free2])) if free2 else chosen.inv()*(-base[rows,:])
        sol.update({p2[k]:s.factor(z)for k,z in zip(cc,values)})
        final={ij:s.expand(pp.subs(sol))for ij,pp in P.items()};res=JV.subs(sol).applyfunc(s.factor)
        distinct=sorted(set(str(z)for z in res if z));D=actual(final)
        # Exact constant shear group and equivariance imply the inverse applies
        # globally; here verify its leading map with all six parameters at once.
        gs=s.symbols('s0:6');unit=project1(plus(push(P0,gs),{ij:-pp for ij,pp in P0.items()}))
        assert unit[:30,:]==S*s.Matrix(gs)
        assert all(not z for z in unit[30:,:]-sum((gs[j]*project1(plus(push(P0,[int(i==j)for i in range(6)]),{ij:-pp for ij,pp in P0.items()}))[30:,:]for j in range(6)),s.zeros(6,1)))
        case=dict(name=record['name'],K1_dimension=K.cols,shear_rank=6,residual_kernel_indices=residual,
            shear_30_by6=[textvec(row)for row in S.tolist()],kernel_free_coordinate_indices=free,
            global_first_source_matrix=[textvec(row)for row in C.tolist()],global_first_source_inverse=[textvec(row)for row in Ci.tolist()],global_constant_determinant=str(C.det()),
            normalized_P1=serialize(P1),all_nine_P2_variables=list(map(str,p2)),P2_pivot_variables=list(map(str,[p2[k]for k in cc])),
            P2_pivot_rows=rows,P2_pivot_columns=cc,P2_minor=str(chosen.det()),P2_solution={str(k):str(v)for k,v in sol.items()},P2_free_variables=list(map(str,[p2[k]for k in free2])),
            all350_after_P2_substitution=textvec(res),distinct_nonzero_remaining_Jacobi=distinct,
            normalized_complete_source=serialize(final),actual_divisor=str(D),actual_divisor_factorization=str(s.factor(D)),
            all_parameter_denominators_constant=all(not s.denom(z).free_symbols for z in list(sol.values())+list(Ci)))
        out.append(case);print(record['name'],'residual',residual,'P2',sol,'remainingJac',distinct,'D',s.factor(D),flush=True)
    result=dict(task_id=job['task_id'],whole_source_shear_reduction=True,cases=out,inputs=job['exact_inputs_and_hashes'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_bytes=Path(__file__).read_text(),
        no_new_first_symbol_or_jet=True,producer_not_imported=True,seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/OUTPUT).open('x')as f:json.dump(result,f,indent=2);f.write('\n')
if __name__=='__main__':main()
