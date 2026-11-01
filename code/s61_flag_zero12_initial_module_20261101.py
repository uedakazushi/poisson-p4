"""Small exact module presentation of the ZERO12 degree-six initial subideal."""
from itertools import combinations_with_replacement

import sympy as s


def module_program(relations):
    polys = [{tuple(mon): s.Rational(c) for mon, c in row} for row in relations]
    assert len(polys) == 64 and all(0 not in mon for f in polys for mon in f)
    active = [4, 5, 7, 8]
    quadratic_mons = list(combinations_with_replacement(active, 2))
    quadratics = [f for f in polys if len(next(iter(f))) == 2]
    assert len(quadratics) == 8
    assert all(all(j in active for j in mon) for f in quadratics for mon in f)
    Q = s.Matrix([[f.get(mon, 0) for mon in quadratic_mons] for f in quadratics])
    reduced, pivots = Q.rref()
    nonpivots = [j for j in range(10) if j not in pivots]
    assert len(pivots) == 8 and [quadratic_mons[j] for j in nonpivots] == [(7, 7), (8, 8)]
    cubic_mons = list(combinations_with_replacement(active, 3))
    cubic_index = {mon: i for i, mon in enumerate(cubic_mons)}
    cubic = s.zeros(32, 20)
    for row in range(8):
        for variable, a in enumerate(active):
            for column, mon in enumerate(quadratic_mons):
                cubic[row * 4 + variable, cubic_index[tuple(sorted(mon + (a,)))]] += Q[row, column]
    assert cubic.rank() == 20
    base = [(), (4,), (5,), (7,), (8,), (7, 7), (8, 8)]
    reduction = {
        quadratic_mons[j]: ([(5 + k, -reduced[pivots.index(j), column])
                             for k, column in enumerate(nonpivots) if reduced[pivots.index(j), column]]
                            if j in pivots else [(5 + nonpivots.index(j), s.S.One)])
        for j in range(10)}

    def normal_form(polynomial, extra=()):
        out = {}
        for mon, coefficient in polynomial.items():
            mon = tuple(sorted(mon + extra))
            act = tuple(j for j in mon if j in active)
            if len(act) >= 3:
                continue
            rad = tuple(j for j in mon if j not in active)
            terms = reduction[act] if len(act) == 2 else [(base.index(act), s.S.One)]
            for index, value in terms:
                out[index, rad] = out.get((index, rad), 0) + coefficient * value
        return {key: value for key, value in out.items() if value}

    def close_generators(source):
        generators, degrees = [], []
        for f in source:
            degree = len(next(iter(f)))
            for extra in base:
                row = normal_form(f, extra)
                if row:
                    generators.append(row)
                    degrees.append(degree + len(extra))
        return generators, degrees

    generators, degrees = close_generators(polys)
    assert len(generators) == 278
    variable_names = {j: 'x_' + str(i) for i, j in enumerate([1, 2, 3, 6])}

    def expression(row, coordinate):
        return '+'.join('(' + str(c) + ')' + ''.join('*' + variable_names[j] for j in mon)
                        for (i, mon), c in row.items() if i == coordinate) or '0'

    def matrix(gens, degs):
        rows = ['{' + ','.join(expression(g, i) for g in gens) + '}' for i in range(7)]
        return ('map(F,S^{' + ','.join(str(-d) for d in degs) + '},matrix{'
                + ','.join(rows) + '})')

    # The expected ordinary-primary ideal is (v4,v7)+(v2,v3+2v6,v5,v8)^2.
    linear = [{(2,): s.S.One}, {(3,): s.S.One, (6,): s.Integer(2)},
              {(5,): s.S.One}, {(8,): s.S.One}]
    expected = [{(4,): s.S.One}, {(7,): s.S.One}]
    for i, f in enumerate(linear):
        for g in linear[i:]:
            out = {}
            for a, ca in f.items():
                for b, cb in g.items():
                    mon = tuple(sorted(a + b))
                    out[mon] = out.get(mon, 0) + ca * cb
            expected.append(out)
    expected_generators, expected_degrees = close_generators(expected)
    program = 'S=QQ[x_0..x_3];F=S^{0,-1,-1,-1,-1,-2,-2};\n'
    program += 'M=' + matrix(generators, degrees) + ';N=image M;A=coker M;\n'
    program += 'assert(dim A==2 and degree A==5);print("INITIAL MODULE DIM2 DEGREE5 PASS");stdio<<flush;\n'
    program += 'B=image(' + matrix(expected_generators, expected_degrees) + ');\n'
    program += 'assert(isSubset(N,B));assert(saturate(N,ideal(-3*x_0+x_2/2))==B);\n'
    program += 'BF=F/B;assert(dim BF==2 and degree BF==5 and degree(BF/(ideal(x_0,x_3)*BF))==5);\n'
    program += 'E=B/N;Ebar=E/(x_0*E);assert(dim E==1 and degree E==111);assert(dim Ebar==0 and degree Ebar==111);\n'
    program += 'assert(apply(toList(1..6),d->hilbertFunction(d,Ebar))=={2,15,28,38,18,10});\n'
    program += 'print("EMBEDDED FREE111 PASS");print("J_HILBERT",hilbertSeries(A,Reduce=>true));print("E_HILBERT",hilbertSeries(E,Reduce=>true));stdio<<flush;\n'
    return program, dict(active_variables=active, active_hilbert_function=[1, 4, 2],
                         active_cubic_rank=20, algebra_basis=[list(mon) for mon in base],
                         absent_variable=0, radical_coefficient_variables=[1, 2, 3, 6],
                         presentation_rows=7, presentation_columns=len(generators),
                         expected_primary_ideal='(v4,v7)+(v2,v3+2*v6,v5,v8)^2',
                         saturation_element='-3*v1+v3/2', embedded_free_base_variables=[0, 1],
                         embedded_free_generator_degrees={1: 2, 2: 15, 3: 28, 4: 38, 5: 18, 6: 10})
