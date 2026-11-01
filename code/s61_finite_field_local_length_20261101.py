#!/usr/bin/env python3
"""Actual implicit jets and complete truncated Macaulay spaces over a prime field.

The caller must supply a rational point, all residual equations, and a tangent
slice. Stable consecutive lengths certify a finite special-fiber local ring.
Characteristic-zero consequences require the separate integral-model argument.
"""
import time
from itertools import combinations_with_replacement
import numpy as np
import sympy as s


def rref(matrix,prime):
    M=matrix.copy()%prime;row=0;pivots=[];ids=list(range(len(M)))
    for col in range(M.shape[1]):
        choices=np.flatnonzero(M[row:,col])
        if not len(choices):continue
        j=row+int(choices[0]);M[[row,j]]=M[[j,row]];ids[row],ids[j]=ids[j],ids[row]
        M[row]=M[row]*pow(int(M[row,col]),-1,prime)%prime
        ii=np.flatnonzero(M[:,col]);ii=ii[ii!=row]
        M[ii]=(M[ii]-M[ii,col,None]*M[row,None,:])%prime
        pivots.append(col);row+=1
        if row==len(M):break
    return M,pivots,ids


def implicit_jets(equations,J,T,chart_columns,maximum=6,prime=32003):
    def mod(value):
        value=s.Rational(value)
        assert int(value.q)%prime!=0
        return int(value.p)*pow(int(value.q),-1,prime)%prime
    def array(M):return np.array([[mod(v) for v in row] for row in M.tolist()],dtype=np.int64)
    fixed=next(i for i in range(45) if i not in chart_columns)
    j=array(J);t=array(T);n=T.cols
    _,columns,_=rref(j,prime);_,rows,_=rref(j[:,columns].T,prime)
    assert len(columns)==len(J.to_DM().rref()[1])
    rank=len(columns)
    inverse,_,_=rref(np.concatenate([j[np.ix_(rows,columns)],np.eye(rank,dtype=np.int64)],axis=1),prime)
    U=np.zeros((45,rank),dtype=np.int64);U[columns,:]=inverse[:,rank:];JU=j@U%prime
    assert np.array_equal(j[rows,:]@U%prime,np.eye(rank,dtype=np.int64))
    assert not np.any(j@t%prime)
    H=np.zeros((len(equations),45,45),dtype=np.int64)
    for row,F in enumerate(equations):
        for mon,coefficient in F.items():
            if fixed in mon:continue
            i,k=[chart_columns.index(v) for v in mon]+[44]*(2-len(mon));v=mod(coefficient)
            if i==k:H[row,i,k]+=v
            else:
                H[row,i,k]+=v*pow(2,-1,prime)%prime
                H[row,k,i]+=v*pow(2,-1,prime)%prime
    H%=prime
    monomials={d:list(combinations_with_replacement(range(n),d)) for d in range(1,maximum+1)}
    indices={d:{mon:i for i,mon in enumerate(monomials[d])} for d in monomials}
    W={1:t};jets={};homogenizer={}
    # Matrix contractions contain at most45 summands and fit in int64.
    assert 45*(prime-1)**2<2**63
    for degree in range(2,maximum+1):
        raw=np.zeros((len(equations),len(monomials[degree])),dtype=np.int64)
        for aa in range(1,degree//2+1):
            bb=degree-aa;factor=1 if aa==bb else 2
            dest=np.array([indices[degree][tuple(sorted(ma+mb))]
                           for ma in monomials[aa] for mb in monomials[bb]],dtype=np.int64)
            # bincount uses doubles; every integer sum is below2^53.
            assert len(dest)*(prime-1)<2**53
            for row in range(len(equations)):
                values=(factor*W[aa].T@(H[row]@W[bb]%prime)%prime).reshape(-1)
                raw[row]=(raw[row]+np.bincount(dest,weights=values,
                         minlength=len(monomials[degree])).astype(np.int64))%prime
        W[degree]=-U@raw[rows,:]%prime
        jets[degree]=(raw-JU@raw[rows,:])%prime
        homogenizer[degree]=W[degree][-1,:]
    return jets,homogenizer,dict(prime=prime,normal_variables=n,
        jacobian_rank=rank,jacobian_pivot_rows=rows,jacobian_pivot_columns=columns,
        all_rational_denominators_prime_units=True,implicit_pivot_prime_unit=True)


def truncated_lengths(jets,homogenizer,n,bounds=(5,6),prime=32003):
    started=time.time();maximum=max(bounds)
    mons={d:list(combinations_with_replacement(range(n),d)) for d in range(2,maximum+1)}
    index={d:{mon:i for i,mon in enumerate(mons[d])} for d in mons};profiles={}
    for bound in bounds:
        offsets={d:sum(len(mons[e]) for e in range(2,d)) for d in range(2,bound+1)}
        count=sum(len(mons[e]) for e in range(2,bound+1));rows=[];labels=[]
        for degree in range(bound-1):
            original=np.concatenate([jets[d] for d in range(2,bound-degree+1)],axis=1)
            _,basis,_=rref(original.T,prime)
            # These selected original rows span ALL original residual jets to
            # the required order. Multiplying a discarded row gives the same
            # linear combination after truncation, so no relation is omitted.
            for mon0 in combinations_with_replacement(range(n),degree):
                destinations={d:np.array([offsets[d+degree]+index[d+degree][tuple(sorted(mon0+mon))]
                    for mon in mons[d]],dtype=np.int64) for d in range(2,bound-degree+1)}
                for i in basis:
                    out=np.zeros(count,dtype=np.int64)
                    for d,loc in destinations.items():out[loc]=jets[d][i]
                    if np.any(out):rows.append(out);labels.append((i,mon0))
        M,pivots,ids=rref(np.array(rows),prime);rank=len(pivots)
        t=np.concatenate([homogenizer[d] for d in range(2,bound+1)])
        assert rank*(prime-1)**2<2**63
        reduced=(t-t[pivots]@M[:rank]%prime)%prime
        profiles[str(bound)]=dict(truncation_power=bound+1,rows=len(rows),columns=count,
            rank=rank,length=1+n+count-rank,homogenizer_zero=not bool(np.any(reduced)),
            pivot_columns=pivots,pivot_original_rows=[labels[i] for i in ids[:rank]])
        print('truncation',bound+1,'length',profiles[str(bound)]['length'],
              'seconds',time.time()-started,flush=True)
    return profiles
