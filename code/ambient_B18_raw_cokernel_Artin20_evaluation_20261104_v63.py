#!/usr/bin/env python3
"""One bounded actual-source evaluation in the frozen 20-dimensional algebra.

Four fixed-origin-Jacobian corrections and final exact equations.  The source
constants have the original ['c',numerator,denominator] representation.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import resource
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = "AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63"
OUT = ROOT / "reports" / STEM
INPUTS = {
    "A20": "reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json",
    "source": "reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json",
    "joint105": "reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json",
    "preflight": "reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_SOURCE_PREFLIGHT_20261104_v63/actual-Artin20-source-and-cost-preflight.json",
}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def add(a, b, c=Q(1)):
    if not b or not c:
        return dict(a)
    out = dict(a)
    for i, x in b.items():
        v = out.get(i, Q(0)) + c * x
        if v:
            out[i] = v
        else:
            out.pop(i, None)
    return out


def scale(a, c):
    if not c:
        return {}
    return {i: c * x for i, x in a.items() if c * x}


def rec(a):
    return [[i, str(x)] for i, x in sorted(a.items())]


def read_rec(a):
    return {i: Q(x) for i, x in a if Q(x)}


def combine(columns, v):
    ans = {}
    for i, c in v.items():
        ans = add(ans, columns[i], c)
    return ans


def rank_basis(cols):
    piv = {}
    receipts = []
    for j, p in enumerate(cols):
        p = dict(p)
        steps = []
        while p:
            k = min(p)
            if k not in piv:
                c = p[k]
                piv[k] = scale(p, 1 / c)
                receipts.append({"source_column": j, "pivot": k, "divide": str(c), "subtract": steps})
                break
            c = p[k]
            steps.append([k, str(c)])
            p = add(p, piv[k], -c)
    return piv, receipts


class Algebra20:
    def __init__(self, d):
        self.data = d
        self.degrees = [sum(m) for m in d["quotient_basis_monomials"]]
        cols = {tuple(m): i for i, m in enumerate(d["monomial_order"])}
        red = [read_rec(r) for r in d["all495_monomial_reductions"]]
        self.table = {}
        for i, a in enumerate(d["quotient_basis_monomials"]):
            for j, b in enumerate(d["quotient_basis_monomials"]):
                m = tuple(x + y for x, y in zip(a, b))
                p = red[cols[m]] if sum(m) <= 4 else {}
                if p:
                    assert len(p) == 1
                    self.table[i, j] = next(iter(p.items()))
        assert len(self.table) == 80
        self.one = {0: Q(1)}

    def mul(self, a, b):
        if not a or not b:
            return {}
        if len(a) == 1 and 0 in a:
            return scale(b, a[0])
        if len(b) == 1 and 0 in b:
            return scale(a, b[0])
        out = {}
        for i, x in a.items():
            for j, y in b.items():
                t = self.table.get((i, j))
                if t:
                    k, c = t
                    z = out.get(k, Q(0)) + c * x * y
                    if z:
                        out[k] = z
                    else:
                        out.pop(k, None)
        return out

    def inv(self, a):
        c = a.get(0, Q(0))
        assert c
        n = dict(a)
        n.pop(0, None)
        n = scale(n, -1 / c)
        ans, term = dict(self.one), dict(self.one)
        for _ in range(4):
            term = self.mul(term, n)
            ans = add(ans, term)
        ans = scale(ans, 1 / c)
        assert self.mul(a, ans) == self.one
        return ans

    def order(self, a):
        return min((self.degrees[i] for i in a), default=5)


def dump_guard(p, d):
    s = json.dumps(d, separators=(",", ":")) + "\n"
    assert len(s.encode()) <= 32 << 20
    p.write_text(s)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw26-evidence", type=Path)
    parser.add_argument("--resume-saved-phase1", action="store_true",
                        help="complete action metadata from saved exact source values; never re-evaluate the source")
    args = parser.parse_args()
    start = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))
    OUT.mkdir(parents=True, exist_ok=True)
    inp = {k: json.loads((ROOT / p).read_text()) for k, p in INPUTS.items()}
    art = Algebra20(inp["A20"])
    if args.resume_saved_phase1:
        assert args.raw26_evidence
        phase_path = OUT / "actual-source105-Artin20-values.json"
        report = json.loads((OUT / "report.json").read_text())
        assert report["phase1_evidence"]["sha256"] == sha(phase_path)
        phase1 = json.loads(phase_path.read_text())
        assert phase1["all105_equations_exactly_zero"]
        finish_actions(args.raw26_evidence, art, phase1, report, start)
        report.update({"action_only_resume_after_relative_path_metadata_error": True,
                       "source_solve_rerun": False, "code_sha256": sha(Path(__file__)),
                       "elapsed_action_only_seconds": time.monotonic() - start,
                       "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
        dump_guard(OUT / "report.json", report)
        print(json.dumps(report), flush=True)
        return
    S, J, P = inp["source"], inp["joint105"], inp["preflight"]
    nodes = S["DAG_scalar_nodes"]
    ids = J["all105_unknown_original_node_order"]
    eqids = J["all105_equation_original_node_order"]
    order = P["needed_closure_order"]
    roots = P["requested_root_nodes"]
    uses0 = collections.Counter()
    for i in order:
        n = nodes[i]
        if n[0] in ["a", "m", "n", "i"]:
            uses0.update(n[1:])
    uses0.update(roots)
    inv_rows = [[(i, Q(x)) for i, x in enumerate(row) if Q(x)] for row in J["fixed105_Jacobian_inverse"]]
    assert sum(map(len, inv_rows)) == 866
    # Actual CA F9 variables: the first7 basis elements after1, then r twice.
    targets = [{i: Q(1)} for i in range(1, 8)] + [{8: Q(1)}, {8: Q(1)}]
    unknown = {i: {} for i in ids}

    def evaluate(u):
        val, uses = {}, uses0.copy()
        inverses = []
        for i in order:
            n, op = nodes[i], nodes[i][0]
            if op == "v":
                val[i] = u[i]
            elif op == "c":
                c = Q(n[1], n[2]) if len(n) == 3 else Q(n[1])
                val[i] = {0: c} if c else {}
            elif op == "a":
                val[i] = add(val[n[1]], val[n[2]])
            elif op == "m":
                val[i] = art.mul(val[n[1]], val[n[2]])
            elif op == "n":
                val[i] = scale(val[n[1]], -1)
            elif op == "i":
                val[i] = art.inv(val[n[1]])
                inverses.append([i, n[1], rec(val[n[1]]), rec(val[i])])
            else:
                raise ValueError(n)
            if op in ["a", "m", "n", "i"]:
                for j in n[1:]:
                    uses[j] -= 1
                    if uses[j] == 0:
                        del val[j]
        return {i: val[i] for i in roots}, inverses

    stages = []
    for stage in range(1, 5):
        values, inverses = evaluate(unknown)
        residual = [dict(values[i]) for i in eqids]
        for k, t in enumerate(targets):
            residual[48 + k] = add(residual[48 + k], t, -1)
        if stage == 1:
            assert all(not v.get(0, 0) for v in values.values() if v is not None)
            assert all(read_rec(q[2]).get(0) == 1 for q in inverses)
        assert min(map(art.order, residual)) >= stage
        correction = []
        for row in inv_rows:
            x = {}
            for i, c in row:
                x = add(x, residual[i], -c)
            correction.append(x)
        unknown = {i: add(unknown[i], x) for i, x in zip(ids, correction)}
        stages.append({"stage": stage, "residual105": [rec(x) for x in residual],
                       "correction105": [rec(x) for x in correction],
                       "unknown105_after": [rec(unknown[i]) for i in ids],
                       "minimum_residual_order": min(map(art.order, residual))})
        dump_guard(OUT / ("stage" + str(stage) + ".json"), stages[-1])
        print(json.dumps({"stage": stage, "elapsed_seconds": time.monotonic() - start}), flush=True)
    values, inverses = evaluate(unknown)
    final_eq = [dict(values[i]) for i in eqids]
    for k, t in enumerate(targets):
        final_eq[48 + k] = add(final_eq[48 + k], t, -1)
    assert all(not x for x in final_eq)
    assert all(not unknown[i].get(0, 0) for i in ids)
    coefficients = []
    for t in P["physical12_coefficient_roots"]:
        aA, aB, bA, bB = [values[t[k]] for k in ["aA_node", "aB_node", "bA_node", "bB_node"]]
        coefficients.append({"physical_row": t["physical_row"],
                             "roots": {k: t[k] for k in ["aA_node", "aB_node", "bA_node", "bB_node"]},
                             "aA": rec(aA), "aB": rec(aB), "delta": rec(add(aA, aB, -1)),
                             "bA": rec(bA), "bB": rec(bB)})
    nus = [values[i] for i in P["nu_roots_CA_CB"]]
    phase1 = {
        "A20_basis": inp["A20"]["quotient_basis_monomials"],
        "actual_CA_F9_targets": [rec(x) for x in targets],
        "joint105_unknown_order": ids, "joint105_equation_order": eqids,
        "all105_unknown_final_values": [rec(unknown[i]) for i in ids],
        "stages": stages, "all105_final_equation_residuals": [rec(x) for x in final_eq],
        "all105_equations_exactly_zero": True,
        "all_requested_original_source_root_values": [[i, rec(values[i])] for i in roots],
        "unit_inverse_exact_checks": inverses,
        "actual_physical_R0_rows": S["physical_zero_base_rows9"],
        "actual_physical_R0_readback": [rec(unknown[i]) for i in range(2, 11)],
        "nu_roots_CA_CB": P["nu_roots_CA_CB"], "nu_CA_CB_values": [rec(x) for x in nus],
        "all12_physical_coefficients": coefficients,
        "preflight_constant_representation_clarification": "source ['c',num,den] evaluated as Fraction(num,den); frozen preflight center-only reader omitted den; actual centers and inverses are rechecked here",
        "input_hashes": {k: {"path": p, "sha256": sha(ROOT / p)} for k, p in INPUTS.items()},
        "source_functions_are_actual_regular_normalizer_roots_not_old_r_g": True,
        "elapsed_seconds": time.monotonic() - start,
    }
    phase_path = OUT / "actual-source105-Artin20-values.json"
    dump_guard(phase_path, phase1)
    assert sum(len(x) for s in stages for x in s["residual105"] + s["correction105"] + s["unknown105_after"]) < 150000
    print(json.dumps({"phase1_saved": str(phase_path), "sha256": sha(phase_path), "elapsed_seconds": time.monotonic() - start}), flush=True)
    report = {"all_passed": True, "actual105_source_evaluation_passed": True,
              "all105_equations_exact": True, "four_fixed_J0_corrections": True,
              "phase1_evidence": {"path": str(phase_path.relative_to(ROOT)), "sha256": sha(phase_path), "bytes": phase_path.stat().st_size},
              "physical_actions_computed": False, "actual_C_special_dimension_computed": False,
              "generic_DVR_rank_or_K_torsion_computed": False, "input_hashes": phase1["input_hashes"],
              "code_sha256": sha(Path(__file__))}
    dump_guard(OUT / "report.json", report)
    if args.raw26_evidence:
        finish_actions(args.raw26_evidence, art, phase1, report, start)
    report["elapsed_seconds"] = time.monotonic() - start
    report["peak_rss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    dump_guard(OUT / "report.json", report)
    print(json.dumps(report), flush=True)


def finish_actions(path, art, phase1, report, start):
    path = path.resolve()
    R = json.loads(path.read_text())
    assert R["source_H20_sha256"] == phase1["input_hashes"]["A20"]["sha256"]
    assert R["H20_basis_monomials"] == phase1["A20_basis"]
    projections = {}
    for fam in ["CA", "CB"]:
        a = R["quotient_maps3_by20"][fam]
        projections[fam] = [{i: Q(a[i][j]) for i in range(3) if Q(a[i][j])} for j in range(20)]
    lift = [{i: Q(row[j]) for i, row in enumerate(R["chosen_common_lift20_by3"]) if Q(row[j])} for j in range(3)]
    assert R["chosen_lift_basis_indices"] == [0, 8, 16]
    assert all(combine(projections[f], lift[j]) == {j: Q(1)} for f in ["CA", "CB"] for j in range(3))
    coordinate_actions = [[read_rec(c) for c in A] for A in R["raw26_coordinate_actions8_column_sparse"]]

    def project(fam, p):
        return combine(projections[fam], p)

    def block(a, b, c):
        out = dict(a)
        out.update({i + 3: x for i, x in b.items()})
        out.update({i + 6: x for i, x in c.items()})
        return out

    nuA, nuB = [read_rec(x) for x in phase1["nu_CA_CB_values"]]
    actions, well_defined = [], []
    # Recheck all intrinsic row-annihilator relations with the actual functions.
    IA = [read_rec(dict(phase1["all_requested_original_source_root_values"])[i])
          for _, i in R["source_conductor_generators_rows_and_ids"]["CA"]]
    IB = [read_rec(dict(phase1["all_requested_original_source_root_values"])[i])
          for _, i in R["source_conductor_generators_rows_and_ids"]["CB"]]
    for p in IA:
        assert not project("CA", p)
    for p in IB:
        assert not project("CB", p)
    assert all(not art.mul(a, b) for a in IA for b in IB)
    for t in phase1["all12_physical_coefficients"]:
        delta, bA, bB = [read_rec(t[k]) for k in ["delta", "bA", "bB"]]
        assert not project("CA", bA) and not project("CB", bB)
        assert all(not art.mul(x, bB) for x in IA)
        assert all(not art.mul(x, bA) for x in IB)
        cols = []
        for x in lift:
            cols.append(block(project("CA", art.mul(delta, x)), {}, scale(art.mul(bB, x), -1)))
        for x in lift:
            cols.append(block({}, project("CB", art.mul(delta, x)), art.mul(bA, x)))
        coeffA = scale(art.mul(nuB, bB), -1)
        coeffB = art.mul(nuA, bA)
        for j in range(20):
            e = {j: Q(1)}
            cols.append(block(project("CA", art.mul(coeffA, e)), project("CB", art.mul(coeffB, e)), art.mul(delta, e)))
        assert len(cols) == 26
        actions.append(cols)
        well_defined.append({"physical_row": t["physical_row"], "bA_mod_IA": [], "bB_mod_IB": [],
                             "IA_times_bB_zero_in_A20": True, "IB_times_bA_zero_in_A20": True})
    for A in actions:
        for X in coordinate_actions:
            assert all(combine(A, X[j]) == combine(X, A[j]) for j in range(26))
    for i, A in enumerate(actions):
        for B in actions[i + 1:]:
            assert all(combine(A, B[j]) == combine(B, A[j]) for j in range(26))
    relation_columns = [c for A in actions for c in A]
    piv, receipts = rank_basis(relation_columns)
    free = [i for i in range(26) if i not in piv]
    free_pos = {i: j for j, i in enumerate(free)}
    reductions = {}
    for i in reversed(range(26)):
        if i in free_pos:
            reductions[i] = {free_pos[i]: Q(1)}
        else:
            q = {}
            for j, c in piv[i].items():
                if j != i:
                    q = add(q, reductions[j], -c)
            reductions[i] = q
    assert all(not combine([reductions[i] for i in range(26)], c) for c in relation_columns)
    quotient_actions = [[combine([reductions[i] for i in range(26)], X[j]) for j in free] for X in coordinate_actions]
    for i, A in enumerate(quotient_actions):
        for B in quotient_actions[i + 1:]:
            assert all(combine(A, B[j]) == combine(B, A[j]) for j in range(len(free)))
    evidence = {
        "source_phase1_file": report["phase1_evidence"]["path"], "source_phase1_sha256": report["phase1_evidence"]["sha256"],
        "raw26_file": str(path.relative_to(ROOT)), "raw26_sha256": sha(path),
        "raw26_basis": R["raw26_basis"], "sector_ranges": [[0, 3], [3, 6], [6, 26]],
        "physical_rows12": [x["physical_row"] for x in phase1["all12_physical_coefficients"]],
        "action_matrix_formula": "[[delta,0,-nuB*bB],[0,delta,nuA*bA],[-bB,bA,delta]]",
        "physical_actions12_column_sparse": [[rec(c) for c in A] for A in actions],
        "intrinsic_row_annihilator_checks": well_defined,
        "all64_actual_source_conductor_products_zero_in_A20": True,
        "all96_coordinate_intertwiners_exact": True, "all66_physical_action_commutators_exact": True,
        "all312_relation_columns_sparse": [rec(c) for c in relation_columns],
        "relation_rank": len(piv), "fixed_pivot_receipts": receipts,
        "echelon_rows": [[i, rec(c)] for i, c in sorted(piv.items())],
        "C_special_dimension": len(free), "C_basis_raw26_indices": free,
        "all26_projection_columns": [rec(reductions[i]) for i in range(26)],
        "C_coordinate_actions8_column_sparse": [[rec(c) for c in A] for A in quotient_actions],
        "right_exactness": "C/(ell)=C_raw/(ell) divided by all12 physical action images; no kernel is specialized",
        "no_K_torsion_or_generic_DVR_rank_claim": True,
    }
    p = OUT / "actual-physical12-actions-and-C-special.json"
    dump_guard(p, evidence)
    coeffs = sum(len(c) for A in evidence["physical_actions12_column_sparse"] for c in A)
    assert coeffs < 150000
    report.update({"physical_actions_computed": True, "actual_C_special_dimension_computed": True,
                   "actual_C_special_dimension": len(free), "physical_action_relation_rank": len(piv),
                   "physical_actions_nonzero_coefficients": coeffs,
                   "all_intrinsic_annihilators_and_intertwiners_checked": True,
                   "raw26_evidence": {"path": str(path.relative_to(ROOT)), "sha256": sha(path)},
                   "action_evidence": {"path": str(p.relative_to(ROOT)), "sha256": sha(p), "bytes": p.stat().st_size},
                   "elapsed_seconds": time.monotonic() - start})


if __name__ == "__main__":
    main()
