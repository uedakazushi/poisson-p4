#!/usr/bin/env python3
"""Source/cost audit only: no Artin-valued source solve is executed."""
import collections
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = "AMBIENT_B18_RAW_COKERNEL_ARTIN20_SOURCE_PREFLIGHT_20261104_v63"
PATHS = {
    "A20": "reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json",
    "source": "reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json",
    "source2": "reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json",
    "joint105": "reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json",
    "free4": "reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json",
    "raw_coker": "reports/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63/DVR20-raw-cokernel26-and-comparison-preflight.json",
}


def digest(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    d = {k: json.loads((ROOT / p).read_text()) for k, p in PATHS.items()}
    A, S, F, J, N, P = (d[k] for k in ["A20", "source", "source2", "joint105", "free4", "raw_coker"])
    nodes = S["DAG_scalar_nodes"]
    mp = dict(N["original_node_to_pruned_node"])
    source_interfaces = []
    for fam in ["CA", "CB"]:
        eq48 = S["families"][fam]["all48_zero_IFT_quadratic_coefficient_equations"]
        cache = F["both57_source2"][fam]
        raw = next(x for x in N["normalizer_sources48_each"] if x["family"] == fam)
        assert eq48 == cache["equations57"][:48]
        assert [mp[x] for x in eq48] == raw["equations48"]
        assert [mp[x] for x in cache["unknown57"][:48]] == raw["unknown48"]
        assert all(Q(v[0]) == 0 for _, v in cache["all57_source2jets"])
        source_interfaces.append({"family": fam, "original48_unknowns": cache["unknown57"][:48],
                                  "original48_equations": eq48, "free4_pruned48_equations": raw["equations48"]})
    eq = J["all105_equation_original_node_order"]
    assert eq[:57] == F["both57_source2"]["CA"]["equations57"]
    assert eq[57:] == F["both57_source2"]["CB"]["equations57"][:48]
    assert eq[48:57] == F["both57_source2"]["CA"]["F9_nodes"]
    roots = eq[:]
    actions = P["all12_physical_difference_actions_on_raw_cokernel"]
    for t in actions:
        roots.extend(t[k] for k in ["aA_node", "aB_node", "bA_node", "bB_node"])
    nus = [S["families"][f]["actual_tracefree_nu"] for f in ["CA", "CB"]]
    roots.extend(nus)
    for fam in ["CA", "CB"]:
        roots.extend(x[1] for x in S["families"][fam]["actual_zero_odd_I0_generators"])
    closure, stack = set(), roots[:]
    while stack:
        i = stack.pop()
        if i in closure:
            continue
        closure.add(i)
        n = nodes[i]
        if n[0] in ["a", "m"]:
            stack.extend(n[1:])
        elif n[0] in ["n", "i"]:
            stack.append(n[1])
    order = sorted(closure)
    uses = collections.Counter()
    for i in order:
        n = nodes[i]
        if n[0] in ["a", "m", "n", "i"]:
            uses.update(n[1:])
    uses.update(roots)
    live, peak = set(), 0
    center, inverses = {}, []
    for i in order:
        n, op = nodes[i], nodes[i][0]
        if op == "v":
            center[i] = Q(0)
        elif op == "c":
            center[i] = Q(n[1])
        elif op == "a":
            center[i] = center[n[1]] + center[n[2]]
        elif op == "m":
            center[i] = center[n[1]] * center[n[2]]
        elif op == "n":
            center[i] = -center[n[1]]
        elif op == "i":
            assert center[n[1]] != 0
            center[i] = 1 / center[n[1]]
            inverses.append([i, n[1], str(center[n[1]])])
        else:
            raise ValueError(n)
        live.add(i)
        peak = max(peak, len(live))
        if op in ["a", "m", "n", "i"]:
            for a in n[1:]:
                uses[a] -= 1
                if not uses[a]:
                    live.remove(a)
    assert all(center[i] == 0 for i in eq)
    assert all(c == "1" for _, _, c in inverses)
    assert {i for i in closure if nodes[i][0] == "v"} == set(J["all105_unknown_original_node_order"])
    mc = {tuple(m): i for i, m in enumerate(A["monomial_order"])}
    red = [dict((i, Q(c)) for i, c in r) for r in A["all495_monomial_reductions"]]
    products = []
    for a in A["quotient_basis_monomials"]:
        for b in A["quotient_basis_monomials"]:
            m = tuple(x + y for x, y in zip(a, b))
            products.append(red[mc[m]] if sum(m) <= 4 else {})
    nonzero_products = sum(bool(p) for p in products)
    assert nonzero_products == 80 and max(map(len, products)) == 1
    kinds = collections.Counter(nodes[i][0] for i in order)
    inverse_nnz = sum(Q(c) != 0 for row in J["fixed105_Jacobian_inverse"] for c in row)
    assert inverse_nnz == 866
    result = {
        "source_interfaces": source_interfaces,
        "CA_F9_targets_in_A20": ["zA0", "zA1", "zA2", "zA3", "zA4", "zA5", "zA6", "rA", "rA"],
        "unknown105_order": J["all105_unknown_original_node_order"],
        "equation105_order": eq, "target_subtraction_positions": list(range(48, 57)),
        "fixed105_inverse_nonzeros": inverse_nnz,
        "needed_closure_order": order, "requested_root_nodes": sorted(set(roots)),
        "unit_inverse_center_checks": inverses,
        "original_node_to_free4_interface": [[i, mp[i]] for i in sorted(set(roots)) if i in mp],
        "nu_roots_CA_CB": nus, "physical12_coefficient_roots": actions,
        "cost": {"original_source_nodes": len(nodes), "needed_nodes": len(order), "operations": dict(kinds),
                 "peak_live_node_vectors": peak, "peak_live_coefficient_bound": 20 * (peak + 105),
                 "A20_nonzero_basis_products": nonzero_products,
                 "A20_max_basis_product_support": 1,
                 "four_pass_dense_product_coefficient_bound": 4 * kinds["m"] * nonzero_products,
                 "source_unknown_and_equation_stage_coefficient_bound": 4 * 2 * 105 * 20,
                 "serialized_nonzero_coefficient_ceiling": 150000,
                 "serialized_byte_ceiling": 32 << 20},
        "finite_algebra_solve_recipe": {
            "start": "all actual unknown offsets zero; retain CA F9 targets in A20",
            "correction": "u <- u - J0_inverse*(equations(u)-targets)",
            "passes": 4,
            "proof": "J(u)-J0 lies in m; residual in m^d becomes m^(d+1); m^5=0",
            "unit_inverse": "c^-1*sum(k=0..4,[-(x-c)/c]^k), exact in A20",
            "final_gates": ["all105 equations exact in A20", "all6 unit identities", "both actual conductor ideals replay", "all12 actions well-defined on the accepted raw26 quotients"]},
        "raw26_interface": {
            "owner": "connectivity", "components": ["A20/IA*tauA", "A20/IB*tauB", "A20*tauA*tauB"],
            "needed": ["piA20to3", "piB20to3", "chosen lifts3to20", "26basis and coordinate actions"],
            "physical_action_matrix": "[[delta,0,-nuB*bB],[0,delta,nuA*bA],[-bB,bA,delta]]",
            "right_exact_coker": "C/(ell)=raw26 / sum_j image(F_j on raw26)",
            "no_naive_kernel_basechange": True},
        "solve_executed": False, "actual_C_special_dimension_computed": False,
        "K_torsion_or_generic_DVR_rank_computed": False,
        "input_hashes": {k: {"path": p, "sha256": digest(ROOT / p)} for k, p in PATHS.items()},
    }
    out = ROOT / "reports" / STEM
    out.mkdir(parents=True, exist_ok=True)
    p = out / "actual-Artin20-source-and-cost-preflight.json"
    p.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    report = {"all_passed": True, "preflight_only": True, "cost": result["cost"],
              "source48_both_bound_to_actual_free4": True, "joint105_F9_adapter_bound": True,
              "all_unit_centers_one": True, "all_unknown_centers_zero": True,
              "no_finite_algebra_source_solve_executed": True,
              "input_hashes": result["input_hashes"],
              "evidence": {"path": str(p.relative_to(ROOT)), "sha256": digest(p), "bytes": p.stat().st_size},
              "code_sha256": digest(Path(__file__))}
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
