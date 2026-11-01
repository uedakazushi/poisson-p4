#!/usr/bin/env python3
"""Independent direct-Jacobi and character review of the q51 polynomial section.

Reads the proposed formulas as data; imports no section-derivation code.
The character check uses Weyl determinant identities, not weight alternants.
"""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import sympy as s
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from s61_model import S61,X

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json'


def evaluate(equation,values,r):
    total=r.zero
    for mon,n,d in equation:
        term=r(QQ(n,d))
        for i in mon:term*=values[i]
        total+=term
    return total


def independent_characters():
    x,y=s.symbols('x y')
    v=x+1/x+y+1/y
    wedge2=2+x*y+x/y+y/x+1/(x*y)
    # Closed primitive linear forms: V tensor wedge2(V) minus two copies of V.
    w=s.expand(v*wedge2-2*v)
    symw=s.expand((w*w+w.subs({x:x*x,y:y*y},simultaneous=True))/2)
    symv=s.expand((v*v+v.subs({x:x*x,y:y*y},simultaneous=True))/2)
    k=s.expand(symv*wedge2-v*v+1)

    def character(a,b):
        first=a+b+2;second=b+1
        numerator=((x**first-x**(-first))*(y**second-y**(-second))
                   -(y**first-y**(-first))*(x**second-x**(-second)))
        denominator=(x*x-x**-2)*(y-y**-1)-(y*y-y**-2)*(x-x**-1)
        return s.cancel(numerator/denominator)

    assert s.cancel(symw-2*character(2,0)-character(2,1)-character(2,2))==0
    assert s.cancel(k-character(2,0)-character(2,1))==0
    return {'method':'Exact Laurent-character identities using Weyl determinants',
            'Sym2W':'2*V(2,0)+V(2,1)+V(2,2)',
            'K45':'V(2,0)+V(2,1)','Hom_dimension':3}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    model=S61();u=s.symbols('u1:62')
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    source=json.loads(SOURCE.read_text())
    assert source['weight1_u_indices_zero_based']==one
    assert source['weight2_u_indices_zero_based']==two
    q=list(map(s.sympify,source['quadratic_section_q']));values=list(u)
    for i,v in zip(two,q,strict=True):values[i]=v
    assert all(s.Poly(v,*[u[i] for i in one]).total_degree()==2 for v in q if v)
    r,*gens=ring(','.join(map(str,[u[i] for i in one]+list(X))),QQ)
    x=gens[-5:];polyvalues=[r.from_expr(v) for v in values]
    assert all(not evaluate(eq,polyvalues,r) for eq in model.reduced['equations'])

    full=model.reconstruct(values)
    b={(i,j):r.from_expr(full[i,j]) for i in range(5) for j in range(5)}
    jacobi={}
    for i,j,k in combinations(range(5),3):
        jacobi[i,j,k]=sum((b[i,h]*b[j,k].diff(x[h])+b[j,h]*b[k,i].diff(x[h])
                           +b[k,h]*b[i,j].diff(x[h]) for h in range(5)),r.zero)
    for quad in combinations(range(5),4):
        residual=sum(((-1)**position*x[index]*jacobi[quad[:position]+quad[position+1:]]
                      for position,index in enumerate(quad)),r.zero)
        assert not residual

    known_path=ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'
    record=next(p for p in json.loads(known_path.read_text())['points'] if p['component_id']=='T39')
    point=list(map(s.Rational,record['S61_coordinates_u1_through_u61']))
    substitution={u[i]:point[i] for i in one}
    assert [v.subs(substitution) for v in q]==[point[i] for i in two]
    assert record['full_equation_jacobian_rank']==45
    q52_path=ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json'
    q52=s.Matrix(list(map(s.sympify,json.loads(q52_path.read_text())['quadratic_section_q'])))
    assert any(v.subs(substitution)!=point[i] for v,i in zip(q52,two,strict=True))

    # Independent coefficient-ansatz check at the single necessary test jet.
    # Universal identities above, and the previous q52 audit, supply sufficiency.
    covariants=[s.Matrix(list(map(s.sympify,v))) for v in source['covariants_in_q_coordinates']]
    sample=dict(zip([u[i] for i in one],range(1,17),strict=True))
    assert s.Matrix.hstack(*covariants).subs(sample).rank()==3
    a,beta,c=s.symbols('a b c')
    trial=(q52+a*covariants[0]+beta*covariants[1]+c*covariants[2]).subs(sample)
    trial_ring,*_=ring('a,b,c',QQ);trial_values=[None]*61
    for i in one:trial_values[i]=trial_ring(sample[u[i]])
    for i,v in zip(two,trial,strict=True):trial_values[i]=trial_ring.from_expr(s.expand(v))
    residuals=[evaluate(eq,trial_values,trial_ring).as_expr() for eq in model.reduced['equations']]
    groebner=s.groebner(residuals,a,beta,c,order='lex')
    expected=[a-20*c/3,beta-8*c,c*c+3*c/32]
    assert [v.as_expr() for v in groebner.polys]==expected
    character=independent_characters()

    paths=[Path(__file__),SOURCE,q52_path,known_path,ROOT/'code/s61_model.py']
    paths+=[model.source/name for name in ('symplectic_slice_certificate.json',
              'symplectic_slice_equations.json','symplectic_slice_reduced_equations.json')]
    report={'id':'Q51_UNIVERSAL_SECTION_REVIEW_20261101','all_passed':True,
        'field':'QQ; all16 first-jet variables symbolic',
        'independence':'No original section-derivation imports; direct homogeneous projective Jacobi; determinant-form character check.',
        'universal134':True,'direct_E_wedge_Jacobi_all5':True,
        'known_T39_all45_coordinates_match':True,'q52_distinct_at_known_T39':True,
        'smooth_witness_rank45_dependency':str(known_path.relative_to(ROOT)),
        'quadratic_monomials':sum(len(s.Poly(v,*[u[i] for i in one]).terms()) for v in q if v),
        'covariant_sample_rank':3,'ansatz_groebner':list(map(str,expected)),
        'character_review':character,
        'scope':['T39 reduced component closure is the polynomial A16 graph.',
                 'Exactly two graded Sp4-equivariant polynomial sections, using the natural covariance proof and the audited q52 section.',
                 'No classification of arbitrary rational sections, full fibers, or embedded structure.'],
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps({'all_passed':True,'direct_projective_Jacobi':True,'Hom_dimension':3,
                      'graded_equivariant_polynomial_sections':2}))


if __name__=='__main__':main()
