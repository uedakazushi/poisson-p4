#!/usr/bin/env python3
"""Actual conductor-product closed DVR fibre, using only frozen source <=3.

No source equation is solved here.  Two cached quadratic coordinate maps and
the cached CA-plane cubic restriction determine a legal unit-normalized IB
generating set modulo m^4; its IA products are therefore determined modulo m^5.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import resource
import time
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = "AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63"
OUT = ROOT / "reports" / STEM
FILES = {
    "full9_source2": "reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json",
    "joint_plane3": "reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json",
    "CB_plane3": "reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json",
    "original_source": "reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json",
    "length20_peer": "reports/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW/report.json",
}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def add(a, b, scale=Q(1)):
    c = dict(a)
    for m, x in b.items():
        y = c.get(m, Q(0)) + scale * x
        if y:
            c[m] = y
        else:
            c.pop(m, None)
    return c


def mul(a, b, cap):
    c = {}
    for m, x in a.items():
        for n, y in b.items():
            p = tuple(i + j for i, j in zip(m, n))
            if sum(p) <= cap:
                z = c.get(p, Q(0)) + x * y
                if z:
                    c[p] = z
                else:
                    c.pop(p, None)
    return c


def homogeneous(a, d):
    return {m: x for m, x in a.items() if sum(m) == d}


def scale(a, x):
    return {m: c * x for m, c in a.items() if c * x}


def var(n, i):
    return {tuple(int(j == i) for j in range(n)): Q(1)}


def compose(a, vs, cap):
    n = len(next(iter(vs[0])))
    one = {(0,) * n: Q(1)}
    c = {}
    powers = [[one] for _ in vs]
    for i, v in enumerate(vs):
        for _ in range(max((m[i] for m in a), default=0)):
            powers[i].append(mul(powers[i][-1], v, cap))
    for m, x in a.items():
        t = scale(one, x)
        for i, k in enumerate(m):
            t = mul(t, powers[i][k], cap)
        c = add(c, t)
    return c


def matinv(a):
    n = len(a)
    r = [[Q(v) for v in row] + [Q(i == j) for j in range(n)] for i, row in enumerate(a)]
    det = Q(1)
    for j in range(n):
        k = next(k for k in range(j, n) if r[k][j])
        if k != j:
            r[j], r[k] = r[k], r[j]
            det = -det
        p = r[j][j]
        det *= p
        r[j] = [v / p for v in r[j]]
        for k in range(n):
            if k != j and r[k][j]:
                c = r[k][j]
                r[k] = [v - c * w for v, w in zip(r[k], r[j])]
    return [row[n:] for row in r], det


def linear_apply(a, ps):
    out = []
    for row in a:
        p = {}
        for x, b in zip(row, ps):
            if x:
                p = add(p, b, x)
        out.append(p)
    return out


def exponent(n, indices):
    v = [0] * n
    for i in indices:
        v[i] += 1
    return tuple(v)


def poly_records(p):
    return [[list(m), str(c)] for m, c in sorted(p.items())]


def sparse_records(p):
    return [[i, str(c)] for i, c in sorted(p.items())]


def sparse_rank(rows):
    basis = {}
    for row in rows:
        row = dict(row)
        while row:
            p = min(row)
            if p not in basis:
                c = row[p]
                basis[p] = {i: x / c for i, x in row.items()}
                break
            row = add(row, basis[p], -row[p])
    return len(basis)


def main():
    t0 = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2 << 30, 2 << 30))
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = {k: json.loads((ROOT / v).read_text()) for k, v in FILES.items()}
    assert inputs["length20_peer"]["all_passed"]
    full, joint, cb = inputs["full9_source2"], inputs["joint_plane3"], inputs["CB_plane3"]
    for fam in ["CA", "CB"]:
        actual = inputs["original_source"]["families"][fam]["actual_zero_odd_I0_generators"]
        assert [p[1] for p in actual[:7]] == full["both57_source2"][fam]["F9_nodes"][:7]
        assert actual[7] == [63, full["both57_source2"][fam]["f63"]]
        assert full["both57_source2"][fam]["all57_equations_all55_coefficients_verified"]
    assert joint["all105_equations_zero_through_degree3"]
    assert cb["all57_equations_verified_through_degree3"]
    mons9 = [exponent(9, x) for x in full["source2_monomials"]]
    V9 = [var(9, i) for i in range(9)]
    xmaps = {}
    for fam in ["CA", "CB"]:
        f = full["both57_source2"][fam]
        assert f["unknown57"][-9:] == list(range(2, 11))
        jj = dict(f["all57_source2jets"])
        xmaps[fam] = [{m: Q(c) for m, c in zip(mons9, jj[i]) if Q(c)} for i in range(2, 11)]
    L = {f: [[p.get(exponent(9, [i]), Q(0)) for i in range(9)] for p in xmaps[f]] for f in xmaps}
    iLB, detLB = matinv(L["CB"])
    YB1 = linear_apply(iLB, [homogeneous(p, 1) for p in xmaps["CA"]])
    saved = [[Q(c) for c in row] for row in full["YB_to_YA_linear"]]
    assert [[p.get(exponent(9, [i]), Q(0)) for i in range(9)] for p in YB1] == saved
    YB2 = linear_apply(iLB, [add(homogeneous(a, 2), compose(homogeneous(b, 2), YB1, 2), -1)
                             for a, b in zip(xmaps["CA"], xmaps["CB"])])
    YB = [add(a, b) for a, b in zip(YB1, YB2)]
    assert [compose(p, YB, 2) for p in xmaps["CB"]] == xmaps["CA"]
    # The block-diagonal linear adapter is important for the cubic gB generator.
    assert all(saved[i][j] == 0 for i in range(7) for j in [7, 8])
    assert all(saved[i][j] == 0 for i in [7, 8] for j in range(7))
    A0 = [row[:7] for row in saved[:7]]
    iA0, detA0 = matinv(A0)
    A1 = [[{} for _ in range(7)] for _ in range(7)]
    plane2 = []
    for i, p in enumerate(YB2[:7]):
        b = {}
        for m, c in p.items():
            js = [j for j in range(7) if m[j]]
            if js:
                j = max(js)  # exact ordered-prefix coordinate Hadamard convention
                n = list(m)
                n[j] -= 1
                A1[i][j] = add(A1[i][j], {tuple(n): c})
            else:
                b[m] = c
        plane2.append(b)
        reconstructed = dict(b)
        for j in range(7):
            reconstructed = add(reconstructed, mul(A1[i][j], V9[j], 2))
        assert reconstructed == p
    plane = []
    plane_rows = []
    for record in joint["all8_CB_zero_restrictions"][:7]:
        plane_rows.append(record["CB_physical_zero_row"])
        p = {}
        for (r, s), c in zip(joint["monomials"], record["all_degree_le3_coefficients"]):
            if Q(c):
                p[(0,) * 7 + (r, s)] = Q(c)
        plane.append(p)
    assert plane_rows == [18, 19, 28, 29, 60, 67, 72]
    assert [homogeneous(p, 2) for p in plane] == plane2
    assert all(not homogeneous(p, 0) and not homogeneous(p, 1) for p in plane)
    beta2 = linear_apply(iA0, plane2)
    A1beta2 = [sum_poly([mul(A1[i][j], beta2[j], 3) for j in range(7)]) for i in range(7)]
    beta3 = linear_apply(iA0, [add(homogeneous(p, 3), t, -1) for p, t in zip(plane, A1beta2)])
    # Substitute ell=5(s-r)=0, retaining actual CA coordinates z0..z6,r.
    V8 = [var(8, i) for i in range(8)]
    restrict = V8[:7] + [V8[7], V8[7]]
    beta2_8 = [compose(p, restrict, 3) for p in beta2]
    beta3_8 = [compose(p, restrict, 3) for p in beta3]
    bnorm = [add(V8[i], add(beta2_8[i], beta3_8[i])) for i in range(7)]
    gA = {(r, s): Q(c) for (r, s), c in zip(joint["monomials"], joint["g_CA_all10_coefficients"]) if Q(c)}
    gB = {(r, s): Q(c) for (r, s), c in zip(cb["ordered_monomials"], cb["g_CB_all10_coefficients"]) if Q(c)}
    assert all(sum(m) == 3 for m in gA) and all(sum(m) == 3 for m in gB)
    gA8 = compose(gA, [V8[7], V8[7]], 4)
    gB8 = compose(gB, [compose(YB1[7], restrict, 1), compose(YB1[8], restrict, 1)], 4)
    assert len(gA8) == len(gB8) == 1 and next(iter(gA8.values())) and next(iter(gB8.values()))
    IA = V8[:7] + [gA8]
    IB = bnorm + [gB8]
    products = [mul(a, b, 4) for a in IA for b in IB]
    monomials = []
    for d in range(5):
        for comb in itertools.combinations_with_replacement(range(8), d):
            monomials.append(exponent(8, comb))
    monomials.sort(key=lambda m: (sum(m), -sum(m[:7]), tuple(-x for x in m)))
    col = {m: i for i, m in enumerate(monomials)}
    assert len(monomials) == 495
    rows, row_sources, aliases, row_keys = [], [], [], {}
    for k, p in enumerate(products):
        if not p:
            continue
        order = min(map(sum, p))
        for n in monomials:
            if sum(n) <= 4 - order:
                q = mul(p, {n: Q(1)}, 4)
                row = {col[m]: c for m, c in q.items()}
                lead = row[min(row)]
                key = tuple((i, c / lead) for i, c in sorted(row.items()))
                if key in row_keys:
                    rid, first_lead = row_keys[key]
                    aliases.append([k, list(n), rid, str(lead / first_lead)])
                else:
                    row_keys[key] = (len(rows), lead)
                    rows.append(row)
                    row_sources.append([k, list(n)])
    input_terms = sum(map(len, rows))
    assert len(rows) <= 2212 and input_terms <= 23000
    progress = {"source_preflight_passed": True, "rows": len(rows), "columns": 495,
                "input_nonzero_coefficients": input_terms, "elapsed_seconds": time.monotonic() - t0}
    (OUT / "progress.json").write_text(json.dumps(progress, indent=2) + "\n")
    print(json.dumps(progress), flush=True)
    basis, receipts = {}, []
    for rid, original in enumerate(rows):
        row, steps = dict(original), []
        while row:
            p = min(row)
            if p not in basis:
                c = row[p]
                basis[p] = {i: v / c for i, v in row.items()}
                receipts.append({"pivot": p, "source_row": rid, "subtract": [[q, str(v)] for q, v in steps],
                                 "divide_by": str(c)})
                break
            c = row[p]
            steps.append((p, c))
            row = add(row, basis[p], -c)
    free = [i for i in range(495) if i not in basis]
    assert len(free) == 20, f"truncated dimension {len(free)}; no inference beyond this point"
    qi = {c: i for i, c in enumerate(free)}
    reduced = {}
    for i in reversed(range(495)):
        if i in qi:
            reduced[i] = {qi[i]: Q(1)}
        else:
            p = {}
            for j, c in basis[i].items():
                if j != i:
                    p = add(p, reduced[j], -c)
            reduced[i] = p

    def reduce_poly(p):
        return sum_poly([scale(reduced[col[m]], c) for m, c in p.items()])

    assert all(not sum_poly([scale(reduced[j], c) for j, c in row.items()]) for row in rows)
    actions = []
    for v in V8:
        actions.append([reduce_poly(mul(v, {monomials[j]: Q(1)}, 4)) for j in free])

    def act(A, p):
        return sum_poly([scale(A[j], c) for j, c in p.items()])

    for i in range(8):
        for j in range(8):
            assert all(act(actions[i], actions[j][k]) == act(actions[j], actions[i][k]) for k in range(20))
    for p in products:
        for j in free:
            assert not reduce_poly(mul(p, {monomials[j]: Q(1)}, 4))
    mpowers = [sparse_rank([reduced[i] for i, m in enumerate(monomials) if sum(m) >= d]) for d in range(6)]
    hf = [mpowers[d] - mpowers[d + 1] for d in range(5)]
    assert hf == [1, 8, 8, 2, 1]
    evidence = {
        "variables": ["zA_" + str(i) for i in range(7)] + ["r_CA"],
        "ell_exact": "5*(s_CA-r_CA)",
        "input_hashes": {k: {"path": v, "sha256": sha(ROOT / v)} for k, v in FILES.items()},
        "actual_physical_rows9": inputs["original_source"]["physical_zero_base_rows9"],
        "source_unknown57_ids": {k: full["both57_source2"][k]["unknown57"] for k in ["CA", "CB"]},
        "source_F9_ids": {k: full["both57_source2"][k]["F9_nodes"] for k in ["CA", "CB"]},
        "source_conductor_generators_rows_and_ids": {k: inputs["original_source"]["families"][k]["actual_zero_odd_I0_generators"] for k in ["CA", "CB"]},
        "physical_coordinate_maps2": {k: [poly_records(p) for p in v] for k, v in xmaps.items()},
        "physical_CB_linear_inverse": [[str(c) for c in row] for row in iLB],
        "physical_CB_linear_determinant": str(detLB),
        "YB_in_YA_linear": [[str(c) for c in row] for row in saved],
        "YB_in_YA_quadratic": [poly_records(p) for p in YB2],
        "A0": [[str(c) for c in row] for row in A0],
        "A0_inverse": [[str(c) for c in row] for row in iA0], "A0_determinant": str(detA0),
        "A1_ordered_prefix": [[poly_records(p) for p in row] for row in A1],
        "b_plane3": [poly_records(p) for p in plane],
        "beta2_before_ell": [poly_records(p) for p in beta2],
        "beta3_before_ell": [poly_records(p) for p in beta3],
        "beta2_after_ell": [poly_records(p) for p in beta2_8],
        "beta3_after_ell": [poly_records(p) for p in beta3_8],
        "gA_plane3": poly_records(gA), "gB_own_plane3": poly_records(gB),
        "gA_after_ell": poly_records(gA8), "gB_after_ell": poly_records(gB8),
        "actual_generator_normalization": {
            "IA": "(zA_0,...,zA_6,gA(0,r,s)); exact coordinate-Hadamard replacement of gA",
            "IB": "(A(z,r,s)^-1 zB,gB(0,rB,sB)); each is an invertible actual ideal-generator change",
            "A": "zB=A*zA+b(r,s), exact ordered-prefix Hadamard; A(0)=A0 unit",
            "beta2": "A0^-1*b2", "beta3": "A0^-1*(b3-A1*A0^-1*b2)",
            "discarded_orders": "zB' unknown O(m4) times IA z lies in m5; gA/gB O(m4) times linear generator lies in m5; gA*gB lies in m6",
        },
        "IA_truncated_generators": [poly_records(p) for p in IA],
        "IB_truncated_generators": [poly_records(p) for p in IB],
        "all64_actual_pair_products_mod_m5": [poly_records(p) for p in products],
        "monomial_order": [list(m) for m in monomials],
        "relation_source_pair_and_multiplier": row_sources,
        "duplicate_source_pair_multiplier_row_and_scale": aliases,
        "all_relation_rows_sparse": [sparse_records(p) for p in rows],
        "pivot_receipts": receipts,
        "echelon_rows": [[i, sparse_records(p)] for i, p in sorted(basis.items())],
        "quotient_basis_monomials": [list(monomials[i]) for i in free],
        "all495_monomial_reductions": [sparse_records(reduced[i]) for i in range(495)],
        "multiplication8_column_sparse": [[sparse_records(p) for p in A] for A in actions],
        "m_power_dimensions": mpowers, "hilbert_function": hf,
        "all_coordinate_actions_commute": True,
        "all64_pair_products_kill_every_basis_vector": True,
        "all_relation_rows_reduce_to_zero": True,
    }
    path = OUT / "actual-special-fiber20-coefficients.json"
    data = json.dumps(evidence, separators=(",", ":")) + "\n"
    assert len(data.encode()) < 32 << 20
    path.write_text(data)
    report = {"all_passed": True, "source_bound_actual_special_fiber_length": 20,
              "truncated_quotient_dimension": 20, "relation_rank": len(basis),
              "relation_rows": len(rows), "input_nonzero_coefficients": input_terms,
              "hilbert_function": hf, "actual_m5_zero": True,
              "coordinate_multiplication_tables": 8,
              "proof": "actual H/(ell) surjects onto the computed quotient by m5; both lengths are 20 by the independently accepted DVR theorem, hence the surjection is an isomorphism",
              "no_generic_DVR_actions_or_K_torsion_computed": True,
              "evidence": {"path": str(path.relative_to(ROOT)), "sha256": sha(path), "bytes": path.stat().st_size},
              "input_hashes": evidence["input_hashes"],
              "code_sha256": sha(Path(__file__)),
              "note_sha256": sha(ROOT / "research" / "notes" / (STEM + ".md")),
              "elapsed_seconds": time.monotonic() - t0,
              "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (OUT / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)


def sum_poly(ps):
    out = {}
    for p in ps:
        out = add(out, p)
    return out


if __name__ == "__main__":
    main()
