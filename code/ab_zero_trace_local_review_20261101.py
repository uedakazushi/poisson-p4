#!/usr/bin/env python3
"""Independent normal-cone, conormal-section and formal-lemma review for AB."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
from itertools import combinations
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def mat(data):return s.Matrix(data).applyfunc(s.sympify)
def run():
    source=ROOT/'reports/AB_ZERO_TRACE_LOCAL_20261101.json';data=json.loads(source.read_text());t=s.Symbol('t');T=mat(data['free_normal_transition']);D=mat(data['full_coordinate_differential']);assert s.factor(T.det())==144*t**5
    free=[0,1,2,3,4,7,10,11,12];Df=D[:,free];Di=(Df*T.inv()).applyfunc(s.cancel);at_inf=Di.applyfunc(lambda f:s.limit(f,t,s.oo));assert all(v.is_finite for v in at_inf)
    assert D.subs(t,0).rank()==13 and Df.subs(t,0).rank()==at_inf.rank()==9
    columns=[]
    for j in range(D.cols):
        for n in range(5):
            col=D[:,j].applyfunc(lambda f:s.Poly(f,t).coeff_monomial(t**n))
            if any(col):columns.append(col)
    span=s.Matrix.hstack(*columns);assert span.rank()==18
    vs=s.symbols('v0:36');phi=s.Matrix([sum(vs[4*j+n]*t**n for n in range(4)) for j in range(9)]).T;inf=(phi*T.inv()).applyfunc(s.expand);equations=[f.coeff(t,n) for f in inf for n in range(1,8)];C=s.linear_eq_to_matrix(equations,vs)[0];assert C.rank()==22
    # Derive the Rees equations directly from all nine minors, without reading
    # the main program's proposed lifted generating set.
    a,b,c,d,e,f,g,h,eps=s.symbols('A B C D E F G H eps');normals=(a,b,c,d,e,f,g,h);variables=(t,*normals,eps)
    N=s.Matrix([[eps*a,eps*b,t+eps*c],[eps*d,eps*e,eps*f],[eps*g,eps*h,-eps*(a+e)]]);minors=[s.expand(N.extract(rr,cc).det()) for rr in combinations(range(3),2) for cc in combinations(range(3),2)];divided=[]
    for p in minors:
        order=min(mon[0] for mon,coef in s.Poly(p,eps).terms());divided.append(s.cancel(p/eps**order))
    rees=s.groebner(divided,*variables,order='lex',domain=s.QQ);assert len(rees.polys)==10 and all(p.LM().exponents[-1]==0 for p in rees.polys)
    parse=lambda value:s.sympify(value,locals={str(v):v for v in variables})
    saved=s.groebner(list(map(parse,data['Rees_groebner_basis'])),*variables,order='lex',domain=s.QQ);assert list(rees)==list(saved)
    initial=[s.expand(p.as_expr().subs(eps,0)) for p in rees.polys];J=s.groebner(initial,t,*normals,order='lex',domain=s.QQ);expected=s.groebner(list(map(parse,data['normal_cone_core_ideal'])),t,*normals,order='lex',domain=s.QQ);assert list(J)==list(expected)
    assert all(J.reduce(t*q)[1]==0 for q in[d,e,g,h])
    assert all(s.expand(q.subs({d:0,e:0,g:0,h:0}))==0 for q in initial)
    vertical=[p.subs(t,0) for p in initial];V=s.groebner(vertical,*normals,order='grevlex',domain=s.QQ);leading=[p.LM(order='grevlex').exponents for p in V.polys];assert len(leading)<=16
    z=s.Symbol('z');num=0
    for mask in range(1<<len(leading)):
        selected=[leading[i] for i in range(len(leading)) if mask>>i&1];degree=sum(max(mon[j] for mon in selected) for j in range(8)) if selected else 0;num+=(-1)**len(selected)*z**degree
    assert s.expand(num-(1-z)**3*(1+3*z+z*z))==0
    free_series=(1+4*z)/(1-z)**10;torsion_series=(1+3*z+z*z)/(1-z)**10-1/(1-z)**9;assert s.cancel(free_series+torsion_series-(1+8*z+z*z)/(1-z)**10)==0
    paths=[Path(__file__),source,ROOT/'code/ab_zero_trace_local_20261101.py',ROOT/'research/notes/AB_ZERO_TRACE_LOCAL_20261101.md',ROOT/'research/notes/AB_INCIDENCE_20261101.md',ROOT/'research/notes/AB_ONE_VARIABLE_LOCAL_20261101.md']
    return {'id':'AB_ZERO_TRACE_LOCAL_REVIEW_20261101','all_passed':True,'field':'QQ','actual_coordinate_conormal_ranks':[13,9],'free_conormal_ranks':[9,9],'free_h0':14,'full_h0_and_coordinate_span':18,'Rees_basis_count':10,'eps_injective_from_leading_monomials':True,'torsion_generated_in_degree_one_and_killed_by_t':True,'vertical_grevlex_basis_count':len(leading),'tangent_hilbert_series':'(1+8z+z^2)/(1-z)^10','embedding_dimension':18,'multiplicity':10,'proof_review':'Actual fiber ideal, all-degree vanishing/generation, completion-generation lemma and graded-filtration identification accepted. No claim that completion equals tangent cone.','proof_note':'research/notes/AB_ZERO_TRACE_LOCAL_REVIEW_20261101.md','primary_sources':['https://stacks.math.columbia.edu/tag/03H0','https://stacks.math.columbia.edu/tag/02O7'],'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'embedding_dimension':18,'multiplicity':10,'output':str(args.output)}))
if __name__=='__main__':main()
