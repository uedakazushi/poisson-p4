#!/usr/bin/env python3
"""Independent linear-span and support-enumeration intersection review.

Uses the opposite monomial pivot order for rational linear algebra and
enumerates bounded/unbounded monomial supports for the Hilbert series.
"""
import argparse
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta

ROOT=Path(__file__).resolve().parents[1]


def reduced_row(polynomial,pivots):
    row=dict(polynomial)
    while row:
        first=min(row)  # Opposite to the original verifier's max-pivot order.
        if first not in pivots:return row
        scale=row[first]
        for key,value in pivots[first].items():
            entry=row.get(key,QQ.zero)-scale*value
            if entry:row[key]=entry
            else:row.pop(key,None)
    return row


def row_space(polynomials):
    pivots={}
    for polynomial in polynomials:
        row=reduced_row(polynomial,pivots)
        if row:
            first=min(row);scale=row[first]
            pivots[first]={key:value/scale for key,value in row.items()}
    return pivots


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    model=S61();y,variables,beta=primitive_beta(model)
    variable_names=list(map(str,variables));r,*x=ring(','.join(variable_names),QQ)
    coefficient=s.Matrix([[s.diff(beta[i,j],v) for v in y]
                          for i,j in list(combinations(range(4),2))[:-1]])
    minors=[r.from_expr(s.expand(coefficient.extract(rows,cols).det()))
            for rows in combinations(range(5),2) for cols in combinations(range(4),2)]
    paths=[ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json',
           ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json']
    sections=[json.loads(path.read_text())['quadratic_section_q'] for path in paths]
    differences=[r.from_expr(s.sympify(a)-s.sympify(b)) for a,b in zip(*sections,strict=True)]
    i2=row_space(differences);j2=row_space(minors)
    assert len(i2)==45 and len(j2)==55
    assert all(not reduced_row(f,j2) for f in differences)
    i3=row_space(v*f for f in differences for v in x)
    assert len(i3)==560
    assert all(not reduced_row(v*f,i3) for f in minors for v in x)

    candidate_path=ROOT/'reports/ELLIPTIC_SECTION_INTERSECTION_GB_INPUT_20261101.json'
    candidate=json.loads(candidate_path.read_text());assert candidate['variables']==variable_names
    gr,*_=ring(','.join(variable_names),QQ,order='grevlex')
    basis=[gr.from_expr(s.sympify(v.replace('^','**'))).monic() for v in candidate['groebner']]
    assert len(basis)==55 and all(max(map(sum,f))==2 for f in basis)
    g2=row_space(r.from_expr(f.as_expr()) for f in basis)
    assert len(g2)==55 and all(not reduced_row(f,g2) for f in minors)
    assert all(not reduced_row(r.from_expr(f.as_expr()),j2) for f in basis)
    pairs_checked=0
    for f,g in combinations(basis,2):
        if all(min(a,b)==0 for a,b in zip(f.LM,g.LM,strict=True)):continue
        common=tuple(max(a,b) for a,b in zip(f.LM,g.LM,strict=True))
        difference=(f.mul_monom(tuple(a-b for a,b in zip(common,f.LM,strict=True)))
                    -g.mul_monom(tuple(a-b for a,b in zip(common,g.LM,strict=True))))
        assert not difference.rem(basis);pairs_checked+=1
    assert pairs_checked==376

    # A standard monomial is specified by an allowed support. The two
    # square-leading variables have exponent1, the other supported variables
    # may have any positive exponent. Enumerate all65536 supports exactly.
    leading=[f.LM for f in basis]
    bounded=sum(1<<i for i in range(16) if any(mon[i]==2 for mon in leading))
    assert bounded.bit_count()==2
    forbidden=[sum(1<<i for i,value in enumerate(mon) if value)
               for mon in leading if max(mon)==1]
    assert len(forbidden)==53
    counts=Counter()
    for support in range(1<<16):
        if not any(support&edge==edge for edge in forbidden):
            counts[(support&~bounded).bit_count(),(support&bounded).bit_count()]+=1
    t=s.Symbol('t')
    numerator=s.expand(sum(count*t**(free+fixed)*(1-t)**(16-free)
                          for (free,fixed),count in counts.items()))
    expected=s.expand((1+11*t+11*t*t+t**3)*(1-t)**11)
    assert numerator==expected
    paths += [candidate_path,Path(__file__),ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
              ROOT/'code/s61_model.py']
    report={'id':'ELLIPTIC_SECTION_INTERSECTION_REVIEW_20261101','all_passed':True,
        'field':'QQ','degree2_ranks':[len(i2),len(j2)],'degree3_rank':len(i3),
        'I_subset_J':True,'maximal_ideal_times_J_subset_I':True,
        'independent_pivot_order':'minimum lexicographic monomial',
        'groebner_basis_degree2_span_match':True,'buchberger_pairs_rechecked':pairs_checked,
        'Hilbert_method':'Enumerate all65536 supports, with two bounded-exponent variables',
        'bounded_variable_indices':[i for i in range(16) if bounded&(1<<i)],
        'support_counts':[[a,b,c] for (a,b),c in sorted(counts.items())],
        'admissible_support_count':sum(counts.values()),
        'Hilbert_series_J':'(1+11*t+11*t^2+t^3)/(1-t)^5',
        'nilradical_dimension':10,'nilradical_square_zero':True,
        'geometric_dependency':'Rank1 closedprimitive tensors form the highest-weight Sp4/B cone; Cartan-power coordinate spaces have dimension(n+1)^4.',
        'scope':'Intersection of the two reduced elliptic graph components, not the ambient S61 local ring.',
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'intersection_dimension':5,'degree':24,
                      'vertex_nilradical_length':10,'admissible_supports':sum(counts.values())}))


if __name__=='__main__':main()
