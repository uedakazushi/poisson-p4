#!/usr/bin/env python3
"""Exact first/second-order obstruction at the representative infinity flag.

Transverse means the homogenizing parameter itself is a uniformizer. Ramified
arcs t=s^e are not excluded. The second-order result is one rational jet only.
"""
from __future__ import annotations
import argparse,hashlib,json
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61
from toric_algebra_degree import sparse
from weight1_rank_one_fiber import row_reduce
ROOT=Path(__file__).resolve().parents[1]
EXTRA=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'


def specialize(model,beta,extra=False):
    weights=model.reduced['parameter_weights'];two=[i for i,w in enumerate(weights) if w==2];index={j:i for i,j in enumerate(two)}
    source=model.reduced['equations']+([json.loads(EXTRA.read_text())['universal_new_equation']] if extra else [])
    polynomials=[]
    for equation in source:
        out={}
        for mon,n,d in equation:
            coefficient=s.Rational(n,d);free=[]
            for i in mon:
                if weights[i]==1:coefficient*=beta[i]
                else:free.append(index[i])
            key=tuple(free);out[key]=out.get(key,0)+coefficient
        polynomials.append({k:s.expand(v) for k,v in out.items() if v})
    return two,polynomials


def chart(polynomials,two):
    k=two.index(37);columns=[j for j in range(45) if j!=k]
    assert all(not f.get((k,k),0) for f in polynomials)
    a=s.Matrix([[f.get(tuple(sorted((k,j))),0) for j in columns] for f in polynomials])
    b=s.Matrix([f.get((k,),0) for f in polynomials])
    return k,columns,a,b


def first_order(model):
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    beta={i:s.Symbol('u'+str(i+1)) for i in one};beta[3]=s.S.Zero
    two,polynomials=specialize(model,beta);_,_,a,b=chart(polynomials[:70],two)
    kernel=a.T.nullspace();assert a.rank()==15 and len(kernel)==55
    compatibility=[s.expand((v.T*b)[0]) for v in kernel]
    variables=[beta[i] for i in one if i!=3]
    basis=row_reduce([sparse(p,variables) for p in compatibility]);assert len(basis)==10
    equations=[s.expand(sum(s.Rational(c.numerator,c.denominator)*s.prod(variables[i] for i in mon)
                           for mon,c in f.items())) for f in basis.values()]
    u5=beta[4];delta=beta[13]-s.Rational(4,5)*beta[48]
    assert u5**2 in equations and s.expand(delta**2) in equations
    assert all(not s.expand(f.subs({u5:0,beta[13]:s.Rational(4,5)*beta[48]})) for f in equations)
    return {'coefficient_field':'QQ[u_i of weight1]/(u4)','matrix_rank':15,
            'left_kernel_dimension':55,'compatibility_span_rank':10,
            'chart_matrix':[[str(v) for v in row] for row in a.tolist()],
            'inhomogeneous_vector':list(map(str,b)),
            'left_kernel':[[str(v) for v in vector] for vector in kernel],
            'compatibility_polynomials':list(map(str,compatibility)),
            'compatibility_basis':list(map(str,equations)),
            'radical_generators':[str(u5),str(delta)],
            'radical_proof':'The compatibility ideal contains both squared generators, and every basis polynomial lies in their ideal.'}


def add(out,key,value):
    out[key]=out.get(key,F(0))+value
    if not out[key]:del out[key]


def multiply(a,b):
    out={}
    for m,c in a.items():
        for n,d in b.items():add(out,tuple(sorted(m+n)),c*d)
    return out


def second_order(model):
    beta={i:s.S.One for i,w in enumerate(model.reduced['parameter_weights']) if w==1}
    beta[3]=beta[4]=s.S.Zero;beta[13]=s.Rational(4,5)
    two,polynomials=specialize(model,beta,extra=True);k,columns,a,b=chart(polynomials,two)
    rr,pivots=a.row_join(-b).to_DM().rref();rr=rr.to_Matrix()
    assert len(pivots)==28 and all(j<44 for j in pivots)
    free=[j for j in range(44) if j not in pivots];assert len(free)==16
    solution=[{} for _ in range(45)]
    for i,j in enumerate(free):solution[columns[j]]={(i,):F(1)}
    for i,j in enumerate(pivots):
        solution[columns[j]]={():F(rr[i,-1]),**{(v,):F(-rr[i,q]) for v,q in enumerate(free)}}
        solution[columns[j]]={mon:c for mon,c in solution[columns[j]].items() if c}
    constants=s.Matrix([solution[j].get((),0) for j in columns])
    directions=s.Matrix([[solution[j].get((i,),0) for i in range(16)] for j in columns])
    assert a*constants+b==s.zeros(135,1) and a*directions==s.zeros(135,16) and directions.rank()==16
    # For q=qflag+t*q1+t²*q2, the order-t² term is A*q2+F(q1).
    expressions=[]
    for f in polynomials:
        out={}
        for mon,c in f.items():
            product={():F(c)}
            for j in mon:product=multiply(product,solution[j])
            for mm,cc in product.items():add(out,mm,cc)
        expressions.append(out)
    kernel=s.Matrix.hstack(*a.T.nullspace()).T
    monomials=[mon for degree in (2,1,0) for mon in combinations_with_replacement(range(16),degree)]
    coefficients=s.Matrix([[f.get(mon,0) for mon in monomials] for f in expressions])
    projected=(kernel.to_DM()*coefficients.to_DM()).to_Matrix()
    target=s.zeros(len(monomials),1);target[-1]=1
    reduced,pivs=projected.T.row_join(target).to_DM().rref();reduced=reduced.to_Matrix()
    assert all(j<projected.rows for j in pivs)
    combination=s.zeros(projected.rows,1)
    for i,j in enumerate(pivs):combination[j]=reduced[i,-1]
    assert combination.T*projected==target.T
    selected=[i for i,v in enumerate(combination) if v]
    assert selected==[6,8,11]
    return {'field':'QQ','beta':{str(i+1):str(v) for i,v in beta.items()},
            'flag_parameter':'u38=1; other weight2 directions zero',
            'first_order_rank':28,'first_order_solution_dimension':16,
            'q_to_u_indices_zero_based':two,'chart_columns':columns,
            'first_order_solution':[[[list(mon),str(c)] for mon,c in f.items()] for f in solution],
            'chart_matrix':[[str(v) for v in row] for row in a.tolist()],
            'inhomogeneous_vector':list(map(str,b)),
            'second_order_projected_span_rank':len(pivs),
            'selected_cokernel_rows':[[str(v) for v in kernel[i,:]] for i in selected],
            'selected_second_order_polynomials':[[[list(mon),str(projected[i,j])] for j,mon in enumerate(monomials) if projected[i,j]] for i in selected],
            'unit_combination':[str(combination[i]) for i in selected],
            'conclusion':'No solution modulo t^3 with q=qflag+t*q1+t²*q2 and chart qflag-coordinate=1 for this rational first jet.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    model=S61();first=first_order(model);second=second_order(model)
    report={'id':'S61_FLAG_TRANSVERSE_OBSTRUCTION_20261101','all_passed':True,
            'first_order':first,'second_order':second,
            'limitations':['Ramified arcs t=s^e for e>1 are not excluded.','The second-order unit certificate concerns one rational jet; generic nonexistence is not inferred from specialization.','No full-fiber dimension or exhaustion result.'],
            'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in (Path(__file__),ROOT/'code/s61_model.py',EXTRA,ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json')}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'first_order_compatibility_rank':10,'second_order_unit_certificate_terms':3,'output':str(args.output)}))

if __name__=='__main__':main()
