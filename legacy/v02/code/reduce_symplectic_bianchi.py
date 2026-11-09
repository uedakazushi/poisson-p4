#!/usr/bin/env python3
"""Reduce the 61-variable ideal from154 to134 generators by exact Bianchi identities.

No Groebner basis or generic linear-system solve is used. All coefficients and
all20 generator-redundancy identities are verified over Q.
"""
from collections import defaultdict, Counter
from fractions import Fraction as F
from pathlib import Path
import itertools
import json
import sys
import time

HERE=Path(__file__).resolve().parent
SYM=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE.parent/'data'
import build_symplectic_slice as s
from build_universal_slice import (clean,add_to,mul,substitute,encode,rref,
                                   singular_polynomial)
p=s.p


def decode(q):
    return {tuple(m):F(n,d) for m,n,d in q}


def sparse_add(a,b,scalar=F(1)):
    for k,c in b.items():
        q=a.get(k,F(0))+scalar*c
        if q:a[k]=q
        else:a.pop(k,None)


class LinearCoordinates:
    def __init__(self,generators):
        self.basis={}
        for i,g in enumerate(generators):
            q=g.copy();rep={i:F(1)}
            while q:
                k=min(q,key=lambda m:(len(m),m))
                if k not in self.basis:
                    c=q[k]
                    q={m:v/c for m,v in q.items()}
                    rep={j:v/c for j,v in rep.items()}
                    self.basis[k]=(q,rep)
                    break
                b,r=self.basis[k];c=q[k]
                sparse_add(q,b,-c);sparse_add(rep,r,-c)
            assert q,'Generators must be linearly independent.'

    def coordinates(self,q):
        q=q.copy();out={}
        while q:
            k=min(q,key=lambda m:(len(m),m))
            assert k in self.basis,'Target outside the constant-linear span.'
            b,r=self.basis[k];c=q[k]
            sparse_add(q,b,-c);sparse_add(out,r,c)
        return out


def main():
    t0=time.time()
    cert=json.loads((SYM/'symplectic_slice_certificate.json').read_text())
    data=json.loads((SYM/'symplectic_slice_equations.json').read_text())
    weights=data['parameter_weights'];raw=cert['raw_slice_indices']
    free=cert['raw_free_indices']
    expr_raw={k:{(i,):F(1)} for i,k in enumerate(free)}
    expr_raw.update({q['raw_index']:decode(q['polynomial']) for q in data['eliminated_raw_variables']})
    expr={i:expr_raw[k] for i,k in enumerate(raw)}
    original=s.equations_from_raw(raw)
    allrows=[substitute(q,expr) for q in original]
    gens=[decode(q) for q in data['equations']]
    degree=lambda q:sum(weights[i] for i in next(iter(q)))
    inds4=[i for i,q in enumerate(gens) if degree(q)==4]
    inds5=[i for i,q in enumerate(gens) if degree(q)==5]
    E4=[gens[i] for i in inds4];E5=[gens[i] for i in inds5]
    assert len(E4)==70 and len(E5)==84
    co4=LinearCoordinates(E4);co5=LinearCoordinates(E5)
    triples=list(itertools.combinations(range(1,5),3))
    ys=lambda d:[a for a in p.homogeneous_monomials(d) if a[0]==0]
    Q4={};Q5={}
    for triple in triples:
        for beta in ys(3):
            alpha=(1,)+beta[1:]
            row=p.FOUR_VECTOR_INDEX[((0,)+triple,alpha)]
            Q4[triple,beta]=co4.coordinates(allrows[row])
        for beta in ys(4):
            row=p.FOUR_VECTOR_INDEX[((0,)+triple,beta)]
            Q5[triple,beta]=co5.coordinates(allrows[row])
    # Degree-one affine bivector coefficients after all linear normalizations.
    linear=[]
    for rawindex,q in expr_raw.items():
        pair,alpha=p.BIVECTOR_BASIS[rawindex]
        if pair[0]>0 and alpha[0]==1 and sum(alpha[1:])==1:
            beta=(0,)+alpha[1:]
            assert all(len(m)==1 and weights[m[0]]==1 for m in q)
            linear.append((pair,beta,{m[0]:c for m,c in q.items()}))
    gamma_list=ys(3);gamma_index={a:i for i,a in enumerate(gamma_list)}
    # d_(constant Pi) on the quartic Jacobi trivector.
    B=[[F(0) for _ in E5] for _ in gamma_list]
    constpairs={(1,2):1,(2,1):-1,(3,4):1,(4,3):-1}
    for (triple,beta),coords in Q5.items():
        omitted=next(k for k in range(1,5) if k not in triple)
        sign=p.permutation_sign((omitted,)+triple)
        for h in range(1,5):
            c=sign*constpairs.get((omitted,h),0)*beta[h]
            if c:
                gamma=list(beta);gamma[h]-=1;row=gamma_index[tuple(gamma)]
                for j,a in coords.items():B[row][j]+=c*a
    # [L,J4]: each entry is a coefficient of a linear parameter times E4_j.
    C=[defaultdict(F) for _ in gamma_list]
    for pair,alpha,lcoords in linear:
        for (triple,beta),qcoords in Q4.items():
            if not qcoords:continue
            for pos,h in enumerate(pair):
                sg=(-1)**(1-pos)*p.permutation_sign(pair[:pos]+pair[pos+1:]+triple)
                coefficient=sg*beta[h]
                if coefficient:
                    gamma=[a+b for a,b in zip(alpha,beta)];gamma[h]-=1
                    row=gamma_index[tuple(gamma)]
                    for e,a in qcoords.items():
                        for u,b in lcoords.items():C[row][e,u]+=coefficient*a*b
            for pos,h in enumerate(triple):
                sg=-(-1)**(2-pos)*p.permutation_sign(triple[:pos]+triple[pos+1:]+pair)
                coefficient=sg*alpha[h]
                if coefficient:
                    gamma=[a+b for a,b in zip(alpha,beta)];gamma[h]-=1
                    row=gamma_index[tuple(gamma)]
                    for e,a in qcoords.items():
                        for u,b in lcoords.items():C[row][e,u]+=coefficient*a*b
    C=[clean(q) for q in C]
    # Verify all20 Bianchi identities as full parameter polynomials.
    for row in range(20):
        q=defaultdict(F)
        for j,c in enumerate(B[row]):add_to(q,E5[j],c)
        for (j,u),c in C[row].items():add_to(q,mul({(u,):F(1)},E4[j]),c)
        assert not clean(q),('Bianchi failed',row,len(clean(q)))
    br,piv,transform=rref(B,transform=True)
    assert len(piv)==20
    transformed_C=[]
    for a in range(20):
        q=defaultdict(F)
        for b,c in enumerate(transform[a]):add_to(q,C[b],c)
        transformed_C.append(clean(q))
    keep=[i for i in range(84) if i not in piv]
    reduced=E4+[E5[i] for i in keep]
    assert len(reduced)==134
    # Directly verify each final eliminated-generator reconstruction identity.
    identities=[]
    for row,pivot in enumerate(piv):
        q=defaultdict(F)
        add_to(q,E5[pivot])
        for j in keep:add_to(q,E5[j],br[row][j])
        for (j,u),c in transformed_C[row].items():add_to(q,mul({(u,):F(1)},E4[j]),c)
        assert not clean(q)
        identities.append({
            'eliminated_old_generator':inds5[pivot],
            'constant_coefficients_on_retained_old_generators':
                [[inds5[j],(-br[row][j]).numerator,(-br[row][j]).denominator]
                 for j in keep if br[row][j]],
            'linear_coefficients_on_weight4_old_generators':
                [[inds4[j],u,(-c).numerator,(-c).denominator]
                 for (j,u),c in sorted(transformed_C[row].items())],
        })
    report={
        'old_generator_count':154,'new_generator_count':134,
        'new_weight_counts':{'4':70,'5':64},
        'quartic_trivector_differential_rank_exact':20,
        'Bianchi_identities_verified_over_Q':20,
        'generator_redundancy_identities_verified_over_Q':20,
        'retained_old_generator_indices':inds4+[inds5[j] for j in keep],
        'eliminated_old_generator_indices':[inds5[j] for j in piv],
        'generator_terms':sum(len(q) for q in reduced),
        'generator_max_ordinary_degree':max(len(m) for q in reduced for m in q),
        'runtime_seconds':round(time.time()-t0,3),
        'scope':'Same ideal over QQ; not a proof of minimality, radicality, or primary decomposition.'
    }
    (SYM/'symplectic_bianchi_certificate.json').write_text(json.dumps(report,indent=2)+'\n')
    (SYM/'symplectic_bianchi_identities.json').write_text(json.dumps({
        'parameter_count':61,'generator_redundancy_identities':identities,
        'Bianchi_constant_matrix':[[[c.numerator,c.denominator] for c in row] for row in B],
        'Bianchi_linear_matrix':[[[j,u,c.numerator,c.denominator] for (j,u),c in sorted(row.items())] for row in C]
    },separators=(',',':'))+'\n')
    (SYM/'symplectic_slice_reduced_equations.json').write_text(json.dumps({
        'parameter_count':61,'parameter_weights':weights,
        'equations':[encode(q) for q in reduced],
        'old_generator_indices':report['retained_old_generator_indices'],
        'reconstruction_file':'symplectic_slice_equations.json',
    },separators=(',',':'))+'\n')
    names=','.join('u'+str(i+1) for i in range(61))
    weights_s=','.join(str(w) for w in weights)
    text='// Universal rank4 transverse Poisson ideal:134 exact generators.\n'
    text+='// No primary decomposition or completeness assertion.\n'
    text+='ring R=0,('+names+'),wp('+weights_s+');\n'
    text+='ideal I=\n'+',\n'.join(singular_polynomial(q) for q in reduced)+';\n'
    text+='// Optional expensive command:\n// LIB "primdec.lib"; list minimalPrimes=minAssGTZ(I);\n'
    (SYM/'symplectic_slice_reduced.sing').write_text(text)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
