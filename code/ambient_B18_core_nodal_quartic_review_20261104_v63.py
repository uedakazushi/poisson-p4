#!/usr/bin/env python3
"""Independent fixed QQ matrix/conic audit, without source-jet evaluation."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
ID = 'AMBIENT_B18_CORE_NODAL_QUARTIC_20261104_v63_REVIEW'


def main(destination):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2_000_000_000, 2_000_000_000))
    signal.alarm(58)
    output = Path(destination)
    output.mkdir(parents=True, exist_ok=False)
    inputs = [
        Path(__file__),
        ROOT / 'research/notes' / (ID + '.md'),
        ROOT / 'reports/AMBIENT_B18_CORE4_QUADRATIC_INITIAL_20261104_v63/report.json',
        ROOT / 'reports/AMBIENT_B18_CORE4_QUADRATIC_20261104_v63_REVIEW/report.json',
        ROOT / 'reports/AMBIENT_B18_CORE_NODAL_QUARTIC_20261104_v63/report.json',
    ]
    data = json.loads(inputs[2].read_text())
    assert data['all_passed']
    assert json.loads(inputs[3].read_text())['all_passed']
    x, y, a, b, t = s.symbols('x112 x117 tauA tauB t')
    variables = (x, y, a, b)
    names = {str(v): v for v in variables}
    f, g = [s.sympify(v, locals=names) for v in data['quadratic_expressions'][:2]]
    assert all(s.Poly(p, variables, domain=s.QQ).total_degree() == 2 for p in [f, g])
    A, B = [s.hessian(p, variables) / 2 for p in [f, g]]
    assert A.rank() == 4 and B.rank() == 3
    pencil = s.factor((A + t * B).det())
    assert s.expand(pencil - (25*t*t + 5*t + 3)/691200) == 0
    assert s.gcd(s.Poly(pencil, t), s.Poly(s.diff(pencil, t), t)).degree() == 0
    vertex = s.Matrix([s.Rational(7, 5), s.Rational(2, 5), s.Rational(7, 24), 1])
    assert B * vertex == s.zeros(4, 1)
    assert (vertex.T * A * vertex)[0] == 0
    assert (vertex.T * B * vertex)[0] == 0
    normal = vertex.T * A
    tangent = normal.nullspace()
    selected = []
    for column in tangent:
        before = s.Matrix.hstack(vertex, *selected).rank()
        if s.Matrix.hstack(vertex, *selected, column).rank() > before:
            selected.append(column)
    assert len(selected) == 2
    V = s.Matrix.hstack(*selected)
    node_matrix = V.T * B * V
    assert node_matrix.det() != 0

    # Projection from the sole possible singular point to the fixed plane b=0.
    w = s.Matrix([x, y, a, 0])
    Q = s.expand((w.T * B * w)[0])
    F = s.expand((w.T * A * w)[0])
    T = s.expand((vertex.T * A * w)[0])
    assert s.hessian(Q, (x, y, a)).rank() == 3
    direction = s.Matrix([8, 1, s.Rational(7, 6), 0])
    assert (direction.T * B * direction)[0] == 0
    assert (vertex.T * A * direction)[0] == s.Rational(1, 20)
    solved_a = s.solve(T, a)[0]
    q2 = s.Poly(s.expand(Q.subs(a, solved_a)), x, y, domain=s.QQ)
    f2 = s.Poly(s.expand(F.subs(a, solved_a)), x, y, domain=s.QQ)
    # No direction gives T=Q=F=0; check the finite chart and its omitted point.
    resultant = s.resultant(q2.as_expr().subs(x, 1), f2.as_expr().subs(x, 1), y)
    assert resultant != 0
    assert not (q2.eval({x: 0, y: 1}) == 0 and f2.eval({x: 0, y: 1}) == 0)
    node_discriminant = s.discriminant(q2.as_expr().subs(x, 1), y)
    assert node_discriminant != 0
    phi = 2*T*w - F*vertex
    assert s.expand((phi.T * A * phi)[0]) == 0
    assert s.expand((phi.T * B * phi)[0] - 4*T*T*Q) == 0

    proof = dict(
        source_variables=[str(v) for v in variables],
        matrix_A=[[str(v) for v in row] for row in A.tolist()],
        matrix_B=[[str(v) for v in row] for row in B.tolist()],
        ranks=[4, 3], determinant_pencil=str(pencil),
        unique_singular_point=[str(v) for v in vertex],
        node_tangent_basis=[[str(v) for v in column] for column in selected],
        node_tangent_determinant=str(node_matrix.det()),
        projection_conic=str(Q), projection_T=str(T),
        rational_conic_direction=[str(v) for v in direction],
        projection_map=[str(s.expand(v)) for v in phi],
        no_basepoint_binary_resultant=str(resultant),
        node_preimage_discriminant=str(node_discriminant),
        all_fixed_polynomial_identities_verified=True,
        conclusions=['geometrically irreducible nodal complete intersection (2,2)',
                     'QQ-rational normalization conic, hence P1 over QQ',
                     'first two actual core functions define an integral CM surface of multiplicity4'],
        excluded_conclusions=['full actual tangent cone', 'higher five core functions',
                              'normality of the actual surface', 'socle', 'D0 or E4 depth'],
    )
    evidence = output / 'independent-matrix-conic-proof.json'
    evidence.write_text(json.dumps(proof, indent=2) + '\n')
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    report = dict(id=ID, all_passed=True,
                  input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs},
                  evidence_sha256={evidence.name: sha(evidence)},
                  seconds=time.monotonic()-started,
                  max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  frozen=True)
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['id', 'all_passed', 'seconds', 'max_rss_kib']}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    main(parser.parse_args().output)
