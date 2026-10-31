#!/usr/bin/env python3
"""Compute the 36 linear syzygies of the 70 weight-four leading quadrics.

Thirty-five are Schouten--Bianchi identities. The remaining one supplies a
universal weight-six, q-degree-at-most-two ideal consequence. At the saved
toric first jet it increases the quadratic span from 134 to 135. No fiber
degree or exhaustion assertion is made.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

import sympy as s

from elliptic_weight1_dominance import projective_affine
from s61_affine_fiber_macaulay import add_term, load_fiber
from s61_model import S61, X
from toric_algebra_degree import add_multiple, sparse, tracked_rows
from weight1_rank_one_fiber import row_reduce

ROOT=Path(__file__).resolve().parents[1]


def kernel(polynomials):
    basis={};relations=[]
    for j in range(45):
        for i,polynomial in enumerate(polynomials):
            row={tuple(sorted(mon+(j,))):c for mon,c in polynomial.items()}
            combination={j*70+i:Fraction(1)}
            while row:
                pivot=min(row)
                if pivot not in basis:
                    value=row[pivot]
                    basis[pivot]=({mon:c/value for mon,c in row.items()},
                                  {k:c/value for k,c in combination.items()})
                    break
                old,certificate=basis[pivot]
                value=-row[pivot]
                add_multiple(row,old,value);add_multiple(combination,certificate,value)
            if not row:relations.append(combination)
    assert len(basis)==3114 and len(relations)==36
    for relation in relations:
        check={}
        for index,c in relation.items():
            j,i=divmod(index,70)
            add_multiple(check,{tuple(sorted(mon+(j,))):v for mon,v in polynomials[i].items()},c)
        assert not check
    assert len(row_reduce([{(k,):c for k,c in relation.items()} for relation in relations]))==36
    return relations


def canonical_bianchi(model,variables,columns,leading):
    values=[s.S.Zero]*61
    for i,v in zip(columns,variables,strict=True):values[i]=v
    homogeneous=model.reconstruct(values)
    y=s.symbols('y0:4');p=projective_affine(homogeneous,X,y)
    j=s.zeros(4);j[0,1]=j[2,3]=1;j[1,0]=j[3,2]=-1
    quadratic=p-j
    assert all(s.Poly(quadratic[a,b],*y).is_homogeneous and
               s.Poly(quadratic[a,b],*y).total_degree()==2 for a,b in combinations(range(4),2))
    triples=list(combinations(range(4),3))
    mon3=[s.prod(y[i] for i in mon) for mon in combinations_with_replacement(range(4),3)]
    jacobi=[]
    for a,b,c in triples:
        expr=s.expand(sum(quadratic[a,h]*s.diff(quadratic[b,c],y[h])
                          +quadratic[b,h]*s.diff(quadratic[c,a],y[h])
                          +quadratic[c,h]*s.diff(quadratic[a,b],y[h]) for h in range(4)))
        jacobi.extend(s.Poly(expr,*y).coeff_monomial(mon) for mon in mon3)
    equations=[model.evaluate(eq,values) for eq in model.reduced['equations'][:70]]
    pivots=list(row_reduce(leading));assert len(pivots)==70
    a=s.Matrix([[poly.get(mon,0) for mon in pivots] for poly in leading])
    b=s.Matrix([[sparse(eq,variables).get(mon,0) for mon in pivots] for eq in jacobi])
    correspondence=b*a.inv()
    assert all(s.expand(expr)==0 for expr in correspondence*s.Matrix(equations)-s.Matrix(jacobi))
    formal=s.symbols('t0:80')
    components={triple:sum(formal[20*i+k]*mon for k,mon in enumerate(mon3))
                for i,triple in enumerate(triples)}
    def trivector(indices):
        if len(set(indices))!=3:return 0
        sign=(-1)**sum(indices[i]>indices[j] for i in range(3) for j in range(i+1,3))
        return sign*components[tuple(sorted(indices))]
    bracket=0
    for i in range(4):
        others=tuple(k for k in range(4) if k!=i)
        bracket+=(-1)**i*sum(quadratic[i,h]*s.diff(trivector(others),y[h]) for h in range(4))
    for i,k in combinations(range(4),2):
        others=tuple(a for a in range(4) if a not in (i,k))
        bracket+=(-1)**(i+k)*sum(s.diff(quadratic[i,k],y[h])*trivector((h,*others)) for h in range(4))
    mon4=list(combinations_with_replacement(range(4),4))
    bracket=s.Poly(s.expand(bracket),*y)
    operator=s.Matrix([[s.expand(bracket.coeff_monomial(s.prod(y[k] for k in mon))).coeff(t)
                        for t in formal] for mon in mon4])
    multipliers=operator*correspondence
    relations=[]
    for row in multipliers.tolist():
        relation={}
        for i,entry in enumerate(row):
            for powers,value in s.Poly(entry,*variables).terms():
                if value:
                    assert sum(powers)==1
                    relation[powers.index(1)*70+i]=Fraction(value)
        check={}
        for index,c in relation.items():
            k,i=divmod(index,70)
            add_multiple(check,{tuple(sorted(mon+(k,))):v for mon,v in leading[i].items()},c)
        assert not check
        relations.append(relation)
    assert len(row_reduce([{(k,):c for k,c in relation.items()} for relation in relations]))==35
    return relations,correspondence,mon4


def centered_consequence(relation,equations,center):
    result={}
    for index,c in relation.items():
        j,i=divmod(index,70)
        add_multiple(result,{tuple(sorted(mon+(j,))):v for mon,v in equations[i].items()},c)
        add_multiple(result,equations[i],c*center[j])
    assert all(len(mon)<=2 for mon in result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    model=S61();weights=model.reduced['parameter_weights']
    equation_weights=[{sum(weights[i] for i in mon) for mon,_,_ in eq}
                      for eq in model.reduced['equations']]
    assert equation_weights[:70]==[{4}]*70 and equation_weights[70:]==[{5}]*64
    columns=[i for i,w in enumerate(weights) if w==2];index={old:new for new,old in enumerate(columns)}
    leading=[]
    for equation in model.reduced['equations'][:70]:
        poly={}
        for mon,n,d in equation:
            if all(i in index for i in mon):add_term(poly,tuple(index[i] for i in mon),Fraction(n,d))
        assert all(len(mon)==2 for mon in poly)
        leading.append(poly)
    relations=kernel(leading)
    variables=s.symbols('q0:45')
    bianchi,correspondence,affine_mon4=canonical_bianchi(model,variables,columns,leading)
    bianchi_basis=row_reduce([{(k,):c for k,c in relation.items()} for relation in bianchi])
    candidates=[]
    for relation in relations:
        remainder={(k,):c for k,c in relation.items()}
        for pivot,row in sorted(bianchi_basis.items()):
            if pivot in remainder:add_multiple(remainder,row,-remainder[pivot])
        if remainder:candidates.append(relation)
    extra=min(candidates,key=len);assert len(extra)==125
    full,_,center,source,point_source=load_fiber()
    old_basis=tracked_rows(full)
    bianchi_redundancies=[]
    for relation in bianchi:
        consequence=centered_consequence(relation,full,center)
        remainder=dict(consequence);combination={}
        while remainder:
            pivot=min(remainder,key=lambda mon:(-len(mon),mon))
            assert pivot in old_basis
            row,certificate=old_basis[pivot];value=remainder[pivot]
            add_multiple(remainder,row,-value);add_multiple(combination,certificate,value)
        bianchi_redundancies.append(combination)
    universal={}
    for key,c in extra.items():
        j,i=divmod(key,70)
        for mon,n,d in model.reduced['equations'][i]:
            add_term(universal,tuple(sorted(mon+[columns[j]])),c*Fraction(n,d))
    assert len(universal)==9257
    assert all(sum(weights[i] for i in mon)==6 and sum(weights[i]==2 for i in mon)<=2 for mon in universal)
    new=centered_consequence(extra,full,center)
    combined=full+[new];pivots=list(row_reduce(combined))
    assert len(pivots)==135 and all(len(mon)==2 for mon in pivots)
    minor=s.Matrix([[poly.get(mon,0) for mon in pivots] for poly in combined])
    determinant=minor.to_DM().det();assert determinant
    # Exact requested flag test: u38=1, other weight-two variables zero.
    flag={}
    for mon,c in universal.items():
        q_part=[i for i in mon if weights[i]==2]
        if q_part==[37,37]:add_term(flag,tuple(i for i in mon if weights[i]==1),c)
    assert not flag
    encode_relation=lambda relation:[[k,str(c)] for k,c in sorted(relation.items())]
    report={
        'id':'S61_UNIVERSAL_WEIGHT6_20261031','all_passed':True,'field':'QQ',
        'weight_four_generators':70,'q_variables':45,
        'linear_syzygy_matrix_shape':[3150,16215],'linear_syzygy_matrix_rank':3114,
        'linear_syzygy_kernel_dimension':36,'canonical_Schouten_Bianchi_dimension':35,
        'multiplier_index_convention':'index=j*70+i represents q_j times weight4 generator i; all indices zero based',
        'q_to_u_indices_zero_based':columns,
        'all_36_kernel_vectors':[encode_relation(relation) for relation in relations],
        'canonical_35_Bianchi_vectors':[encode_relation(relation) for relation in bianchi],
        'affine_quartic_monomial_order':[list(mon) for mon in affine_mon4],
        'raw_Jacobi_to_stored70_sparse_rows':[[[i,str(value)] for i,value in enumerate(row) if value]
                                            for row in correspondence.tolist()],
        'selected_non_Bianchi_syzygy':encode_relation(extra),
        'universal_new_equation':[[list(mon),c.numerator,c.denominator] for mon,c in sorted(universal.items())],
        'universal_new_equation_weight':6,'universal_new_equation_q_degree':2,
        'universal_new_equation_monomial_count':len(universal),
        'sample_centered_new_equation':[[list(mon),str(c)] for mon,c in sorted(new.items())],
        'sample_Bianchi_constant_redundancies':[encode_relation(row) for row in bianchi_redundancies],
        'sample_augmented_quadratic_span_rank':135,
        'sample_augmented_constant_linear_consequence_rank':0,
        'sample_augmented_minor_monomials':[list(mon) for mon in pivots],
        'sample_augmented_minor_determinant':str(determinant),
        'degree_bound_for_augmented_minor_in_beta':66,
        'u38_flag_new_quadratic_leading_restriction':'0 identically in all 16 first-jet variables',
        'claim':'A non-Bianchi universal linear syzygy of the 70 weight4 quadratic leading parts gives a weight6 ideal consequence of q-degree at most two. On a nonempty first-jet open, adding it increases the quadratic span from 134 to 135 without changing the ideal.',
        'limitations':['No fiber degree, finiteness, component, or exhaustion claim.',
                       'The leading flag test is at the specified representative; no removal of the full infinity locus is asserted.',
                       'The new equation is already an ideal consequence, not a new Poisson condition.'],
        'proof_note':'research/notes/S61_UNIVERSAL_WEIGHT6_20261031.md',
        'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (Path(__file__),ROOT/'code/s61_model.py',source,point_source)},
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:
        json.dump(report,out,ensure_ascii=False,indent=2);out.write('\n')
    print(json.dumps({'all_passed':True,'kernel_dimension':36,'Bianchi_dimension':35,
                      'new_weight':6,'sample_quadratic_rank':135,'flag_restriction':0,'output':str(args.output)}))


if __name__=='__main__':main()
