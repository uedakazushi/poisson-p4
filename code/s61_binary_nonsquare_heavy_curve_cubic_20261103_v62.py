#!/usr/bin/env python3
"""Proper arithmetic certificate: every heavy-defect curve has generic cubic rank57."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import time

import numpy as np
import sympy as s
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]
PRIME = 32003
UNITS = "F=500*a+90*e^2+3*e-24;B=(50*a+18*e^2-9*e-2)*d+60*e-40;C=100*a+18*e^2-9*e-2;"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def matrix_text(rows):
    return "matrix{" + ",".join("{" + ",".join(str(c).replace("**", "^") for c in row) + "}" for row in rows) + "}"


def extract_matrix(log, end):
    start = log.index("matrix {{")
    return log[start:log.index("\n" + end, start)]


def run_m2(output, name, source, marker, timeout, reuse=None, old_name=None):
    path = output / (name + ".m2")
    path.write_text(source)
    if reuse is not None and old_name is not None:
        previous = reuse / (old_name + ".m2")
        assert previous.read_bytes() == path.read_bytes(), name
        log = (reuse / (old_name + ".log")).read_text()
        provenance = str(previous)
    else:
        cp = subprocess.run([shutil.which("M2"), "--script", str(path)],
                            cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        log = cp.stdout + cp.stderr
        (output / (name + ".log")).write_text(log)
        cp.check_returncode()
        provenance = "executed by this invocation"
    assert marker in log and "error:" not in log and "SIGSEGV" not in log, name
    (output / (name + ".log")).write_text(log)
    return log, provenance


def modular(coefficient):
    q = Fraction(str(coefficient))
    return q.numerator * pow(q.denominator, -1, PRIME) % PRIME


def show_mod(poly, variables):
    terms = []
    for mon, coefficient in poly.items():
        c = modular(coefficient)
        if c:
            terms.append(str(c) + "".join("*" + v + ("^" + str(n) if n > 1 else "")
                                         for v, n in zip(variables, mon) if n))
    return "+".join(terms) or "0"


class QuadraticField:
    def __init__(self, b, c):
        self.b, self.c = b % PRIME, c % PRIME
        assert pow((b*b-4*c) % PRIME, (PRIME-1)//2, PRIME) == PRIME-1

    def add(self, x, y):
        return ((x[0]+y[0]) % PRIME, (x[1]+y[1]) % PRIME)

    def mul(self, x, y):
        a, b = x
        c, d = y
        return ((a*c-self.c*b*d) % PRIME, (a*d+b*c-self.b*b*d) % PRIME)

    def pow(self, x, n):
        result = (1, 0)
        for _ in range(n):
            result = self.mul(result, x)
        return result

    def inv(self, x):
        a, b = x
        norm = (a*a-self.b*a*b+self.c*b*b) % PRIME
        assert norm
        inv = pow(norm, -1, PRIME)
        return ((a-self.b*b)*inv % PRIME, -b*inv % PRIME)

    def rank(self, values):
        A = np.array(values, dtype=np.int64)
        nr, nc, dim = A.shape
        assert dim == 2
        rank, pivots = 0, []
        for j in range(nc):
            ids = np.flatnonzero(np.any(A[rank:, j, :], axis=1))
            if not len(ids):
                continue
            i = int(ids[0])+rank
            A[[i, rank]] = A[[rank, i]]
            v = self.inv(tuple(map(int, A[rank, j])))
            x, y = A[rank, :, 0].copy(), A[rank, :, 1].copy()
            A[rank, :, 0] = (x*v[0]-self.c*y*v[1]) % PRIME
            A[rank, :, 1] = (x*v[1]+y*v[0]-self.b*y*v[1]) % PRIME
            for i in range(rank+1, nr):
                a, b = map(int, A[i, j])
                if not a and not b:
                    continue
                x, y = A[rank, :, 0], A[rank, :, 1]
                A[i, :, 0] = (A[i, :, 0]-a*x+self.c*b*y) % PRIME
                A[i, :, 1] = (A[i, :, 1]-a*y-b*x+self.b*b*y) % PRIME
            pivots.append(j)
            rank += 1
            if rank == nr:
                break
        return rank, pivots

    def det(self, values):
        A = [[tuple(map(int, v)) for v in row] for row in values]
        value = (1, 0)
        for j in range(len(A)):
            i = next(i for i in range(j, len(A)) if A[i][j] != (0, 0))
            if i != j:
                A[i], A[j] = A[j], A[i]
                value = ((-value[0]) % PRIME, (-value[1]) % PRIME)
            pivot = A[j][j]
            value = self.mul(value, pivot)
            inverse = self.inv(pivot)
            for i in range(j+1, len(A)):
                multiplier = self.mul(A[i][j], inverse)
                if multiplier == (0, 0):
                    continue
                for h in range(j, len(A)):
                    q = self.mul(multiplier, A[j][h])
                    A[i][h] = ((A[i][h][0]-q[0]) % PRIME, (A[i][h][1]-q[1]) % PRIME)
        return value


def run(output, reuse):
    began = time.time()
    output.mkdir(parents=True, exist_ok=False)
    assert shutil.which("M2")
    template_path = ROOT / "code/data/s61_binary_hyperplane_heavy_template_20261102.json"
    upstream_path = ROOT / "reports/S61_BINARY_HYPERPLANE_HEAVY_CUBIC_CODIMENSION_20261102.json"
    upstream = json.loads(upstream_path.read_text())
    assert upstream["all_passed"]
    for path, expected in upstream["input_sha256"].items():
        assert sha(ROOT / path) == expected
    template = json.loads(template_path.read_text())
    R, a, e, d = ring("a,e,d", s.QQ)
    known = [a, 500*a+90*e**2+3*e-24,
             (50*a+18*e**2-9*e-2)*d+60*e-40,
             100*a+18*e**2-9*e-2]
    original = [[R.from_dict({tuple(m): s.QQ(c) for m, c in p}) for p in row]
                for row in template["matrix"]]
    change = [[s.QQ(c) for c in row] for row in template["constant_row_change"]]
    assert s.Matrix(change).rank() == 15
    heavy = [[sum((q*original[k][j] for k, q in enumerate(row) if q), R.zero)
              for j in range(7)] for row in change]
    contents, retained_contents = [], []
    for row in heavy:
        content = R.zero
        for poly in row:
            content = content.gcd(poly)
        remainder = content
        orders = []
        for f in known:
            n = 0
            while not remainder.rem(f):
                remainder = remainder.exquo(f)
                n += 1
            orders.append(n)
        divisor = content.exquo(remainder)*remainder.LC
        row[:] = [poly.exquo(divisor) for poly in row]
        contents.append(orders)
        retained_contents.append(str(remainder/remainder.LC))
    module_source = ("R=QQ[a,e,d,z];\n" + UNITS + "\nM=" + matrix_text(heavy) + ";\n"
                     'E=transpose M | ((z*a*F*B*C-1)*id_(R^7));print "HEAVY_QQ_INPUT_READY"; <<flush; G=gb image E; print "HEAVY_QQ_GB_DONE"; <<flush; print(id_(R^7)%G==0);print(dim coker E);print toString gens G; print "HEAVY_QQ_DONE"; <<flush;\n')
    module_log, module_origin = run_m2(output, "module", module_source, "HEAVY_QQ_DONE", 720,
                                     reuse, "heavy-module-QQ-all")
    assert "\nfalse\n1\nmatrix" in module_log
    module_matrix = extract_matrix(module_log, "HEAVY_QQ_DONE")
    ann_source = ("R=QQ[a,e,d,z];\nG=" + module_matrix + ";\n"
                  'print(numRows G,numColumns G);<<flush;J=ann coker G;print "HEAVY_QQ_ANN_READY";<<flush;print(dim(R/J),degree(R/J));print toString gens gb J;print "HEAVY_QQ_ANN_DONE"; <<flush;\n')
    ann_log, ann_origin = run_m2(output, "ann", ann_source, "HEAVY_QQ_ANN_DONE", 240,
                               reuse, "heavy-ann-QQ")
    assert "\n(1, 16)\n" in ann_log
    ann_matrix = extract_matrix(ann_log, "HEAVY_QQ_ANN_DONE")
    Q, *avars = ring("a,e,d,z", s.QQ)
    ann = [Q.from_expr(s.sympify(text.replace("^", "**")))
           for text in ann_matrix[len("matrix {{"):-2].split(", ")]
    assert len(ann) == 25
    leading = []
    for poly in ann:
        for coefficient in poly.values():
            modular(coefficient)
        mon = max(poly, key=lambda m: (sum(m), tuple(-x for x in reversed(m))))
        value = modular(poly[mon])
        assert value
        leading.append(dict(monomial=mon, coefficient_modp=value, total_degree=sum(mon)))
    raw_generic = ",".join(show_mod(poly, ["aa", "ee", "dx1", "zx1"]) for poly in ann)
    raw_global = ",".join(show_mod(poly, ["a", "e", "d", "z"]) for poly in ann)
    proper_source = "\n".join([
        "T=(ZZ/32003)[aa];F=frac T;S=F[ee,dx1,zx1];J=ideal(" + raw_generic + ");",
        'G=gb(J,ChangeMatrix=>true);assert(dim(S/J)==0);assert(degree(S/J)==4);print "GENERIC_LENGTH4_PASS";',
        "CM62=(gens G)//(gens J);assert(gens J*CM62==gens G);hh62=1_T;",
        "scan(flatten apply((flatten entries gens G | flatten entries CM62),f->flatten entries last coefficients f),c->(hh62=lcm(hh62,denominator(leadCoefficient c))));",
        'assert(sub(hh62,{aa=>1})!=0);print("H_DENOM62="|toString hh62);print "DENOMINATOR_IDENTITY_PASS";',
        "R=(ZZ/32003)[a,e,d,z];Jall=ideal(" + raw_global + ");",
        "highestPart62=f->(q:=degree f;sum(select(terms f,t->degree t==q)));Lall=ideal apply(flatten entries gens Jall,highestPart62);",
        'assert(dim Jall==1);assert(degree Jall==16);assert(dim Lall==1);print "PROJECTIVE_INFINITY_DIM0_PASS";print toString gens G;',
    ]) + "\n"
    proper_log, _ = run_m2(output, "properness", proper_source, "PROJECTIVE_INFINITY_DIM0_PASS", 180)
    htext = proper_log.split("H_DENOM62=", 1)[1].splitlines()[0]
    aa = s.Symbol("aa")
    h = s.Poly(s.sympify(htext.replace("^", "**")), aa, modulus=PRIME)
    scalar, factors = s.factor_list(h)
    assert s.Poly(scalar*s.prod(f**n for f, n in factors), aa, modulus=PRIME) == h
    vertical = []
    vertical_source = []
    for i, (factor, multiplicity) in enumerate(factors):
        # SymPy's finite-field factorization is exact; record and recheck every
        # irreducible factor before constructing its residue field.
        assert factor.is_irreducible
        text = str(factor.as_expr()).replace("**", "^")
        vertical.append(dict(polynomial=str(factor.as_expr()), multiplicity=multiplicity,
                             degree=factor.degree(), irreducible=True))
        vertical_source += [
            "T=(ZZ/32003)[symbol aa];F=toField(T/ideal(" + text + "));S=F[symbol ee,symbol dx1,symbol zx1];aa=F_0;J=ideal(" + raw_generic + ");",
            'G=gb J;dv62=dim(S/J);assert(dv62<=0);print("ANN_VERTICAL_CASE_' + str(i) + '="|toString(dv62,degree(S/J)));<<flush;',
        ]
    vertical_source.append('print "ALL_ANN_VERTICAL_EXCEPTIONS_ZERO_DIMENSIONAL";<<flush;')
    vertical_log, _ = run_m2(output, "vertical", "\n".join(vertical_source)+"\n",
                             "ALL_ANN_VERTICAL_EXCEPTIONS_ZERO_DIMENSIONAL", 180)
    for i, record in enumerate(vertical):
        result = vertical_log.split("ANN_VERTICAL_CASE_"+str(i)+"=", 1)[1].splitlines()[0]
        record["fiber_dimension_and_length"] = result
    replay_path = ROOT / "reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/report.json"
    replay = json.loads(replay_path.read_text())
    assert replay["all_passed"]
    for name, expected in replay["generated_sha256"].items():
        assert sha(replay_path.parent/name) == expected
    with tarfile.open(replay_path.parent/"matrices.tar.gz") as tar:
        schur = tar.extractfile("schur").read()
    assert hashlib.sha256(schur).hexdigest() == replay["source"]["schur_matrix_sha256"]
    lines = schur.decode().splitlines()
    nr, nc, pivots = map(int, lines[0].split())
    assert (nr, nc, pivots) == (42, 70, 21)
    polynomials = []
    for line in lines[3:]:
        terms = line.split()
        count = int(terms[0])
        assert len(terms) == 1+4*count
        polynomials.append([(tuple(map(int, terms[2+4*k:5+4*k])), modular(terms[1+4*k])) for k in range(count)])
    assert len(polynomials) == nr*nc
    choice, points, e_minpolys = None, [], []
    for b, c, E, D in ((235, 1065, (6680, 24479), (223, 4274)),
                        (-2903, -7388, (7302, 20430), (26690, 6349))):
        field = QuadraticField(b, c)
        Fvalue = field.add((476, 0), field.add(field.mul((90, 0), field.pow(E, 2)), field.mul((3, 0), E)))
        Cvalue = field.add((98, 0), field.add(field.mul((18, 0), field.pow(E, 2)), field.mul((-9, 0), E)))
        Bvalue = field.add(field.mul(field.add((48, 0), field.add(field.mul((18, 0), field.pow(E, 2)), field.mul((-9, 0), E))), D), field.add(field.mul((60, 0), E), (-40, 0)))
        Z = field.inv(field.mul(field.mul(Fvalue, Bvalue), Cvalue))
        bases = ((1, 0), E, D, Z)
        maxpowers = [max(m[j] for poly in ann for m in poly) for j in range(4)]
        powers = [[field.pow(base, k) for k in range(n+1)] for base, n in zip(bases, maxpowers)]
        for poly in ann:
            result = (0, 0)
            for mon, coefficient in poly.items():
                value = (modular(coefficient), 0)
                for j, n in enumerate(mon):
                    value = field.mul(value, powers[j][n])
                result = field.add(result, value)
            assert result == (0, 0)
        trace = (2*E[0]-field.b*E[1]) % PRIME
        norm = (E[0]**2-field.b*E[0]*E[1]+field.c*E[1]**2) % PRIME
        assert E[1]
        e_minpolys.append([1, -trace % PRIME, norm])
        ep = [field.pow(E, k) for k in range(1+max(m[1] for poly in polynomials for m, _ in poly))]
        dp = [field.pow(D, k) for k in range(1+max(m[2] for poly in polynomials for m, _ in poly))]
        values = []
        for poly in polynomials:
            value = (0, 0)
            for mon, coefficient in poly:
                term = field.mul(ep[mon[1]], dp[mon[2]])
                value = field.add(value, (term[0]*coefficient % PRIME, term[1]*coefficient % PRIME))
            values.append(value)
        matrix = np.array(values, dtype=np.int64).reshape(nr, nc, 2)
        rank, columns = field.rank(matrix)
        assert rank == 36
        _, rows = field.rank(matrix[:, columns].transpose(1, 0, 2))
        if choice is None:
            choice = dict(rows=rows, columns=columns)
        determinant = field.det(matrix[np.ix_(choice["rows"], choice["columns"])])
        assert determinant != (0, 0)
        points.append(dict(quadratic=[b, c], e=E, d=D, z=Z,
                           ann_generators_vanish=25, e_minimal_polynomial=e_minpolys[-1],
                           full_cubic_residual_rank=rank, common_minor_determinant=determinant))
    assert e_minpolys[0] != e_minpolys[1]
    entry_degree = max(sum(mon) for i in choice["rows"] for j in choice["columns"]
                       for mon, _ in polynomials[i*nc+j])
    paths = [Path(__file__), template_path, upstream_path, replay_path]
    report = dict(
        all_passed=True, prime=PRIME, exact_QQ_heavy_cokernel_support_dimension=1,
        auxiliary_QQ_annihilator_degree=16, auxiliary_variables=["a", "e", "d", "z"],
        localization_relation="z*a*F*B*C-1", row_known_unit_content_orders=contents,
        unknown_row_contents_retained=retained_contents,
        integral_ann_generators=25, ann_grevlex_leading_terms_and_residues=leading,
        special_projective_infinity_dimension=0, special_generic_a_algebra_length=4,
        generic_GB_conversion_matrix_checked=True, denominator_polynomial=htext,
        denominator_degree=h.degree(), denominator_unit_at_a1=True,
        all_vertical_exception_fibers_zero_dimensional=True, vertical_factors=vertical,
        common_minor=choice, quadratic_Galois_pairs=points,
        four_special_points_are_distinct=True,
        common_entry_homogenization_degree=entry_degree,
        homogeneous_determinant_degree=36*entry_degree,
        conclusion="Over QQ, dim{rank(15x7 heavy block)<7 AND rank(full cubic Schur42x70)<36}<=0 on aFBC!=0.",
        proof_method="Fixed-degree homogenization in ordinary P4 and proper upper semicontinuity; no nonproper localization inference.",
        limitations=["No quartic surjectivity, transverse length65, or ambient t-nilpotence conclusion.",
                     "The QQ all-core quadratic-product rank41 assertion is separate."],
        earlier_same_session_calculations_reused={"module": module_origin, "ann": ann_origin},
        input_sha256={str(path.relative_to(ROOT)): sha(path) for path in paths},
        generated_sha256={path.name: sha(path) for path in output.iterdir() if path.suffix in (".m2", ".log")},
        seconds=time.time()-began,
    )
    (output/"report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({key: report[key] for key in ("all_passed", "conclusion", "seconds")}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reuse-dir", type=Path,
                        help="Reuse earlier same-session module/ann logs only when their generated M2 scripts match byte for byte.")
    args = parser.parse_args()
    run(args.output, args.reuse_dir)
