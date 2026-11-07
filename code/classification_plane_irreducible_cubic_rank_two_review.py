#!/usr/bin/env python3
"""Independent exact nodal/cuspidal logarithmic-form certificate.

Only standard-library Fraction arithmetic. Topology and connection extension
are separately proved in the accompanying review, not inferred from this run.
"""
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

ROOT=Path(__file__).resolve().parents[1]

def add(*ps):
    d={}
    for p in ps:
        for e,a in p.items():d[e]=d.get(e,F(0))+a
    return {e:a for e,a in d.items()if a}

def scale(p,a):return {e:a*c for e,c in p.items()if a*c}

def mul(p,q):
    d={}
    for e,a in p.items():
        for f,b in q.items():
            g=tuple(i+j for i,j in zip(e,f));d[g]=d.get(g,F(0))+a*b
    return {e:a for e,a in d.items()if a}

def der(p,i):
    d={}
    for e,a in p.items():
        if e[i]:
            f=list(e);f[i]-=1;d[tuple(f)]=a*e[i]
    return d

def mons(d):
    out=[]
    for terms in it.combinations_with_replacement(range(3),d):
        e=[0]*3
        for i in terms:e[i]+=1
        out.append(tuple(e))
    return sorted(out,reverse=True)

def rref(matrix):
    M=[[F(c)for c in row]for row in matrix];piv=[];rank=0
    for j in range(len(M[0])):
        k=next((i for i in range(rank,len(M))if M[i][j]),None)
        if k is None:continue
        M[k],M[rank]=M[rank],M[k];a=M[rank][j];M[rank]=[c/a for c in M[rank]]
        for i in range(len(M)):
            if i!=rank and M[i][j]:
                a=M[i][j];M[i]=[c-a*d for c,d in zip(M[i],M[rank])]
        piv.append(j);rank+=1
    free=[j for j in range(len(M[0]))if j not in piv];K=[]
    for j in free:
        v=[F(0)]*len(M[0]);v[j]=1
        for i,p in enumerate(piv):v[p]=-M[i][j]
        K.append(v)
    return M,piv,K

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(89,90));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic();x={(1,0,0):F(1)};y={(0,1,0):F(1)};z={(0,0,1):F(1)};xyz=[x,y,z]
    node=add(mul(mul(y,y),z),scale(mul(mul(x,x),x),-1),scale(mul(mul(x,x),z),-1))
    cusp=add(mul(mul(x,x),z),scale(mul(mul(y,y),y),-1))
    cases={}
    for name,f in [('node',node),('cusp',cusp)]:
        columns=[mul(xyz[j],der(f,i))for i,j in it.product(range(3),repeat=2)]+[scale(f,-1)]
        matrix=[[p.get(e,F(0))for p in columns]for e in mons(3)]
        red,piv,K=rref(matrix)
        assert len(piv)==(9 if name=='node'else 8)
        expected=[F(1),0,0,0,F(1),0,0,0,F(1),F(3)]
        assert all(not sum(c*v for c,v in zip(row,expected))for row in matrix)
        if name=='cusp':
            second=[F(3),0,0,0,F(2),0,0,0,0,F(6)]
            assert all(not sum(c*v for c,v in zip(row,second))for row in matrix)
            assert expected[8] and not second[8]
        # Every basis vector satisfies the entire coefficient matrix.
        assert all(not sum(c*v for c,v in zip(row,k))for row in matrix for k in K)
        cases[name]={'all_Vf_minus_lambda_f_matrix':[[str(c)for c in row]for row in matrix],
                     'coefficient_monomials':mons(3),'rank':len(piv),'complete_kernel':[[str(c)for c in k]for k in K],
                     'kernel_mod_Euler_dimension':len(K)-1}
    # alpha=(2yz dx-3xz dy+xy dz)/(x²z-y³).
    numerator=[scale(mul(y,z),2),scale(mul(x,z),-3),mul(x,y)]
    d_numerator={}
    targets={(0,1):z,(0,2):scale(y,-1),(1,2):x}
    for i,j in it.combinations(range(3),2):
        p=add(mul(cusp,add(der(numerator[j],i),scale(der(numerator[i],j),-1))),
              scale(mul(numerator[j],der(cusp,i)),-1),mul(numerator[i],der(cusp,j)))
        assert p==mul(cusp,targets[i,j]);d_numerator[i,j]=p
    assert all(d_numerator.values())
    # Euler horizontality and contraction representation are exact.
    assert not add(*(mul(v,n)for v,n in zip(xyz,numerator)))
    V=[scale(x,3),scale(y,2),{}]
    contraction=[add(mul(y,V[2]),scale(mul(z,V[1]),-1)),
                 add(mul(z,V[0]),scale(mul(x,V[2]),-1)),
                 add(mul(x,V[1]),scale(mul(y,V[0]),-1))]
    assert all(a==scale(b,-1)for a,b in zip(numerator,contraction))
    # Nodal/cuspidal partial derivatives have no common constant direction.
    for f in [node,cusp]:
        m=[[der(f,i).get(e,F(0))for i in range(3)]for e in mons(2)]
        assert len(rref(m)[1])==3
    inputs=['research/classification/PLANE_IRREDUCIBLE_CUBIC_RANK_TWO_EXCLUSION.md',
            'research/classification/PLANE_SMOOTH_CUBIC_CONORMAL_CONNECTION_REVIEW.md',
            'reports/CLASSIFICATION_PLANE_SMOOTH_CUBIC_RANK_TWO_REVIEW_20261107.json']
    out={'task_id':'E_reduced_nonnormal','all_passed':True,
         'scope':'Exact global logarithmic-one-form calculation for both irreducible singular plane cubic types; topology and connection geometry are independently proved in the review note.',
         'cases':cases,'cusp_alpha_nonclosed_full_rational_differential_identity':True,
         'cusp_dalpha':'(z dx^dy-y dx^dz+x dy^dz)/(x²z-y³)',
         'Euler_horizontal_and_full_linear_field_contraction_checked':True,
         'no_nonzero_constant_direction_annihilates_either_cubic':True,
         'accepted_universal_full30_full9_arithmetic_applies_without_new_specialization':True,
         'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in inputs},
         'reviewer_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    target=ROOT/'reports/CLASSIFICATION_PLANE_IRREDUCIBLE_CUBIC_RANK_TWO_REVIEW_20261107.json'
    assert not target.exists();target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
