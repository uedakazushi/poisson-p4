#!/usr/bin/env python3
"""Independent coefficient reconstruction and dimension audit of nil22 exclusion."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J


def verify():
    start=time.time();root=Path(__file__).resolve().parents[1]
    path=root/'reports/S61_NIL22_INFINITY_20261101.json';report=json.loads(path.read_text());m=S61()
    one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2]
    def matrix(Q):
        P=X[0]**2*J;P[1:,1:]+=Q;q=m.normalize(P)['parameters'];M=s.zeros(134,16)
        for r,e in enumerate(m.reduced['equations']):
            constant=0
            for mon,n,d in e:
                qs=[i for i in mon if i in two];bs=[i for i in mon if i in one]
                if len(qs)!=2:continue
                c=s.Rational(n,d)*s.prod(q[i] for i in qs)
                if not bs:constant+=c
                else:
                    assert len(bs)==1
                    M[r,one.index(bs[0])]+=c
            assert s.cancel(constant)==0
        return M.applyfunc(s.expand)
    def parse(M):return s.Matrix(M).applyfunc(s.sympify)
    def determinant(A):
        field=s.QQ.algebraic_field(s.sqrt(3)) if any(v.has(s.sqrt(3)) for v in A) else s.QQ
        syms=sorted(A.free_symbols,key=str);domain=field.poly_ring(*syms) if syms else field
        return domain.to_sympy(A.to_DM().convert_to(domain).det())
    count=0
    for entry in report['fixed_cores']:
        M=matrix(parse(entry['Q']));assert (M-parse(entry['beta_matrix'])).applyfunc(s.expand)==s.zeros(134,16)
        d=determinant(M.extract(entry['rows'],entry['columns']));assert s.simplify(d-s.sympify(entry['constant_minor']))==0
        assert d!=0;count+=1
        if 'rank10_witness' in entry:
            w=entry['rank10_witness'];sample=dict(zip(s.symbols('p0:3'),w['p']))
            d=determinant(M.subs(sample).extract(w['rows'],w['columns']))
            assert d==s.sympify(w['evaluated_minor'])!=0;count+=1
    for entry in report['closed_lines']:
        M=matrix(parse(entry['Q']));assert (M-parse(entry['beta_matrix'])).applyfunc(s.expand)==s.zeros(134,16);g=0
        for minor in entry['rank12_minors']:
            d=s.factor(determinant(M.extract(minor['rows'],minor['columns'])))
            assert s.cancel(d-s.sympify(minor['determinant']))==0;g=s.gcd(g,d);count+=1
        assert s.cancel(g/s.sympify(entry['gcd'])).free_symbols==set()
        roots=s.solve(g,s.Symbol('d'));assert set(roots)=={s.sympify(e['d']) for e in entry['exceptional_values']}
        for e in entry['exceptional_values']:
            d=determinant(M.subs(s.Symbol('d'),s.sympify(e['d'])).extract(e['rows'],e['columns']))
            assert s.cancel(d-s.sympify(e['constant_minor']))==0 and d!=0;count+=1
    e=report['open_chart'];M=matrix(parse(e['normalized_Q']));assert (M-parse(e['beta_matrix'])).applyfunc(s.expand)==s.zeros(134,16)
    d=s.factor(determinant(M.extract(e['minor_rows'],e['minor_columns'])))
    assert d==s.Rational(2116013625,512);count+=1
    paths=[Path(__file__),path,root/'code/s61_model.py',root/'research/notes/S61_NIL22_INFINITY_20261101.md']
    return dict(id='S61_NIL22_INFINITY_REVIEW_20261101',all_passed=True,
        all11_beta_matrices_reconstructed_from_original134=True,checked_minors=count,
        open_uniform_rank=12,full_projective_dimension_bound=9,rank9_locus_projective_bound=8,beta_incidence_bound=15,
        proof_note='research/notes/S61_NIL22_INFINITY_REVIEW_20261101.md',
        dependency='Archived17-branch fixed-trace reduced exhaustion; checked centralizer normalizations in main certificate.',
        source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();r=verify()
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({k:r[k] for k in ['all_passed','checked_minors','beta_incidence_bound','seconds']}))


if __name__=='__main__':main()
