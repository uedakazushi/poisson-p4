#!/usr/bin/env python3
"""E0: inventory frozen image maps and export their exact graph presentations.

No factorization, rank extrapolation, elimination, or exhaustion assertion.
All localizations below live on the parameter source unless explicitly stated.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OLD = "code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json"
NEW = "code/data/known52_rank4_s61_catalogue_20261103_v63.json"
DIM = "reports/KNOWN51_FIRSTJET_DIMENSION_TABLE_20261103_v62.json"
MODEL = "legacy/v02/data/symplectic_slice_reduced_equations.json"
MODEL_CERT = "legacy/v02/data/symplectic_slice_certificate.json"


def digest(data):
    return hashlib.sha256(data).hexdigest()


class Inputs:
    def __init__(self):
        self.pins = {}

    def read(self, path, expected=None, compressed=False):
        raw = (ROOT / path).read_bytes()
        h = digest(raw)
        if expected is not None and h != expected:
            raise ValueError(f"Input digest mismatch: {path}")
        self.pins[path] = h
        decoded = gzip.decompress(raw) if compressed else raw
        return json.loads(decoded), digest(decoded)


def refs(node):
    if node[0] in {"q", "v"}:
        return []
    if node[0] == "inv":
        return [node[1]]
    if node[0] in {"+", "*"}:
        return node[1]
    raise ValueError(f"Unexpected node opcode: {node[0]}")


def load_catalogue(inputs):
    combined, _ = inputs.read(NEW)
    old, _ = inputs.read(combined["old51_catalogue"]["path"],
                         combined["old51_catalogue"]["sha256"])
    new, _ = inputs.read(combined["additional_record"]["path"],
                         combined["additional_record"]["sha256"])
    rows = old["records"] + [new]
    assert len(rows) == combined["component_count"] == 52
    assert len({x["component_id"] for x in rows}) == 52
    return rows


def examine(record, inputs, weight_one):
    archive, decoded = inputs.read(record["map"]["path"],
                                   record["map"]["sha256"], compressed=True)
    assert decoded == record["decoded_sha256"]
    assert archive["component_id"] == record["component_id"]
    names, nodes = archive["input_names"], archive["nodes"]
    assert names == record["input_names"]
    outputs = archive["output_u1_through_u61"]
    assert len(outputs) == 61
    source_support = []
    for i, node in enumerate(nodes):
        r = refs(node)
        assert all(type(j) is int and 0 <= j < i for j in r)
        if node[0] == "v":
            assert 0 <= node[1] < len(names)
            support = 1 << node[1]
        else:
            support = 0
            for j in r:
                support |= source_support[j]
        source_support.append(support)
    beta_outputs = [outputs[j] for j in weight_one]
    ancestors, todo = set(), list(beta_outputs)
    while todo:
        j = todo.pop()
        if j not in ancestors:
            ancestors.add(j)
            todo.extend(refs(nodes[j]))
    inverses = []
    for i, node in enumerate(nodes):
        if node[0] == "inv":
            j = node[1]
            inverses.append({
                "inverse_node": i, "argument_node": j,
                "argument_input_support": [name for k, name in enumerate(names)
                                           if source_support[j] & (1 << k)],
                "needed_for_firstjet_outputs": i in ancestors,
                "locus": "source_parameter_open",
                "base_polynomial_in_beta_identified": False,
                "source_complement": f"V(node_{j}) on the source before localization",
                "base_complement_disposition": "image/alternate-chart audit required; not V of an identified beta polynomial",
            })
    return archive, {
        "component_id": record["component_id"],
        "map": record["map"], "decoded_sha256": decoded,
        "coefficient_field": record.get("coefficient_field", "QQ"),
        "field_relation": record.get("field_relation"),
        "source_inputs": names,
        "declared_source_opens": record.get("opens", archive.get("normalization", {})),
        "source_inverse_nodes": inverses,
        "full_map_node_count": len(nodes),
        "firstjet_ancestor_node_count": len(ancestors),
        "firstjet_output_nodes": beta_outputs,
        "image_closed_ideal": {
            "construction": "contraction to C[beta] of the exact DAG graph ideal; see EXCEPTIONAL_BASE.md",
            "finite_graph_presentation_available": True,
            "expanded_elimination_generators_computed_here": False,
            "source_chart_boundary_covered_here": False,
            "whole_S61_fiber_above_image_covered_here": False,
        },
    }, sorted(ancestors)


def inventory(output):
    inputs = Inputs()
    model, _ = inputs.read(MODEL)
    inputs.read(MODEL_CERT)
    assert model["parameter_count"] == 61 and len(model["equations"]) == 134
    w1 = [i for i, w in enumerate(model["parameter_weights"]) if w == 1]
    assert [i + 1 for i in w1] == [4, 5, 12, 13, 14, 25, 26, 27, 34, 35, 36, 37, 48, 49, 50, 51]
    dimensions, _ = inputs.read(DIM)
    dims = {r["component_id"]: r for r in dimensions["records"]}
    records = []
    for row in load_catalogue(inputs):
        _, record, _ = examine(row, inputs, w1)
        cid = row["component_id"]
        if cid in dims:
            d = dims[cid]
            inputs.read(d["proof_record"], d["proof_sha256"])
            record.update({k: d[k] for k in ("s61_component_dimension", "firstjet_image_dimension",
                                            "generic_fiber_dimension", "general_firstjet_dominant")})
            record["dimension_evidence"] = {"path": d["proof_record"], "sha256": d["proof_sha256"]}
        else:
            p = "reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63/report.json"
            d, _ = inputs.read(p)
            assert d["all_passed"] and d["does_not_dominate_W16"]
            record.update(s61_component_dimension=16, firstjet_image_dimension=15,
                          generic_fiber_dimension=1, general_firstjet_dominant=False)
            record["dimension_evidence"] = {"path": p, "sha256": inputs.pins[p]}
        records.append(record)
    for p in ("reports/KNOWN51_FIRSTJET_IMAGE_DIMENSIONS_20261103_v62_REVIEW.json",
              "reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW/report.json",
              "reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63_REVIEW.json",
              "code/data/known51_coordinate_metadata_clarification_20261103_v62.json"):
        inputs.read(p)
    g28_peer_path = "reports/Q_NEW_FAMILY_DEGREE28_EXACT_EQUATION_20261104_v63_REVIEW/report.json"
    g28_peer, _ = inputs.read(g28_peer_path)
    assert g28_peer["all_passed"] and g28_peer["QQ_and_complex_defining_equation_up_to_nonzero_scalar"]
    assert g28_peer["actual_reduced_primitive_firstjet_B15_irreducible_equation_degree"] == 28
    g28_path = "reports/Q_NEW_FAMILY_DEGREE28_EXACT_EQUATION_20261104_v63/equation.json"
    g28_expected = g28_peer["input_sha256"][g28_path]
    g28, _ = inputs.read(g28_path, g28_expected)
    assert g28["name"] == "G28" and len(g28["terms"]) == 78
    records[-1]["image_closed_ideal"].update(
        known_accepted_generator="G28 in the primitive beta convention; transport to u via the audited linear adapter",
        known_accepted_generator_path=g28_path,
        known_accepted_generator_sha256=inputs.pins[g28_path],
        known_accepted_generator_degree=28,
        known_accepted_generator_review=g28_peer_path,
        stale_dense_map_record_does_not_override_completed_image_equation=True,
    )
    dominant = [r["component_id"] for r in records if r["general_firstjet_dominant"]]
    assert dominant == ["T01", "T39", "T40"]
    body = {
        "task_id": "E0", "version": "v64-classification-first", "date_JST": "2026-11-06",
        "scope": "source/dependency inventory; no new elimination/rank or exhaustion computation",
        "coefficient_ring": "C[beta_1,...,beta_16,q_1,...,q_45] with all 134 original equations",
        "base_coordinate_u_indices_1based": [j + 1 for j in w1],
        "component_count": 52, "known_base_dominant_components": dominant,
        "proper_known_image_count": 49,
        "image_dimension_histogram": dict(sorted(Counter(r["firstjet_image_dimension"] for r in records).items())),
        "base_only_generic_exhaustion_denominator_identified": False,
        "source_inversions_are_not_base_complements": True,
        "necessary_generic_five_base_denominator_condition": {
            "statement": "If full S61 over D(f(beta)) is covered by T01,T39,T40, then f belongs to every one of the49 proper known image ideals",
            "new25_consequence": "G28 divides f after the audited primitive-beta coordinate transport",
            "proof": "EXCEPTIONAL_BASE.md, denominator necessity lemma",
            "not_sufficient_for_exhaustion": True,
        },
        "E_exhaustion_proved": False,
        "records": records,
        "closed_support_queue": [
            {"id": "E.G-complement", "base_closed_set": "V(f_G(beta)) once G supplies base-only denominator",
             "status": "blocked: G has not produced f_G or generic exhaustion"},
            {"id": "E.known49", "base_closed_set": "49 individual closed image ideals defined by exact DAG contractions",
             "status": "image dimensions reported; no assertion their union exhausts exceptional bases"},
            {"id": "E.beta-zero", "base_closed_set": "V(beta_1,...,beta_16)",
             "status": "Q45_22 is repository-reported support exhaustion; complete22-to-known52 embedding crosswalk still to audit"},
            {"id": "E.discriminant", "base_closed_set": "V(Delta)",
             "status": "known cover degenerations do not prove full S61 support; new25 image is not Delta"},
            {"id": "E.fiber-localizations", "base_closed_set": "scheme-theoretic base images of retained unknown-coordinate chart complements",
             "status": "cannot substitute source denominator zeros for their projected base supports"},
        ],
        "input_sha256": inputs.pins,
        "code_sha256": digest(Path(__file__).read_bytes()),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as f:
        json.dump(body, f, ensure_ascii=False, indent=2)
        f.write("\n")
    print(json.dumps({"records": 52, "proper_images": 49, "dominant": dominant,
                      "sha256": digest(output.read_bytes()), "mathematical_exhaustion": False}))


def graph(component, output):
    """Export a finite ideal recipe, not an elimination result or coverage proof."""
    inputs = Inputs()
    model, _ = inputs.read(MODEL)
    w1 = [i for i, w in enumerate(model["parameter_weights"]) if w == 1]
    row = next(r for r in load_catalogue(inputs) if r["component_id"] == component)
    archive, record, ancestors = examine(row, inputs, w1)
    # DPM's archived input r is algebraic, so it MUST become a coefficient,
    # not a free source variable with a prose-only algebraicity condition.
    archived_names = record["source_inputs"]
    algebraic = {"r": "rho"} if component == "R4_DPM" else {}
    source_names = [name for name in archived_names if name not in algebraic]
    source_positions = {name: i for i, name in enumerate(source_names)}
    definitions = []
    for i in ancestors:
        node = archive["nodes"][i]
        if node[0] == "v":
            name = archived_names[node[1]]
            node = (["coeff", algebraic[name]] if name in algebraic
                    else ["v", source_positions[name]])
        definitions.append([i, node])
    field_specification = (
        {"field": "QQ(rho)", "generator": "rho", "minimal_polynomial": "rho^2-3",
         "complex_base_change": "choose one embedding QQ(rho)->C; rho is a fixed coefficient",
         "algebraic_source_variables": [], "conjugate_embeddings_unioned_here": False}
        if algebraic else {"field": "QQ", "complex_base_change": "QQ->C"}
    )
    body = dict(component_id=component, coefficient_field=record["coefficient_field"],
                field_relation=record["field_relation"],
                coefficient_field_specification=field_specification,
                archived_source_inputs=archived_names,
                algebraic_input_coefficient_specializations=algebraic,
                source_inputs=source_names,
                beta_coordinates=[f"u{i + 1}" for i in w1],
                beta_output_nodes=record["firstjet_output_nodes"],
                node_definitions=definitions,
                equations="q: t_i-c=0; coeff: t_i-c=0 over the specified coefficient field; v: t_i-x_j=0 indexed by source_inputs; +/*: t_i-op(children)=0; inv: t_i*t_j-1=0; beta_k-t_output_k=0",
                image_ideal="contract this graph ideal to C[beta]; no elimination has been executed",
                omitted_full_map_source_opens=record["declared_source_opens"],
                omission_justification="dense nonempty source opens do not change closed image on an irreducible source; needs recorded dense-map theorem and chosen algebraic-field embedding",
                input_sha256=inputs.pins)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as f:
        json.dump(body, f, ensure_ascii=False, indent=2)
        f.write("\n")


def beta_zero(output):
    """Restrict every original equation, retaining all row labels and zero rows."""
    inputs = Inputs()
    model, _ = inputs.read(MODEL)
    weights = model["parameter_weights"]
    qindices = [i for i, w in enumerate(weights) if w == 2]
    qmap = {old: new for new, old in enumerate(qindices)}
    restricted = []
    for i, equation in enumerate(model["equations"]):
        terms = {}
        for mon, numerator, denominator in equation:
            if any(weights[j] == 1 for j in mon):
                continue
            key = tuple(sorted(qmap[j] for j in mon))
            terms[key] = terms.get(key, Fraction(0)) + Fraction(numerator, denominator)
        terms = {mon: c for mon, c in terms.items() if c}
        assert all(len(mon) == 2 for mon in terms)
        restricted.append({"original_row_0based": i,
                           "equation": [[list(mon), c.numerator, c.denominator]
                                        for mon, c in sorted(terms.items())]})
    assert len(restricted) == 134 and len(qindices) == 45
    body = dict(task_id="E0", coefficient_ring="QQ[q0,...,q44]",
                parameter_q_u_indices_1based=[i + 1 for i in qindices],
                all134_original_rows_retained=True,
                nonzero_row_count=sum(bool(r["equation"]) for r in restricted),
                zero_row_count=sum(not r["equation"] for r in restricted),
                equations=restricted, input_sha256=inputs.pins,
                Q45_22_crosswalk_or_new_exhaustion_claim=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as f:
        json.dump(body, f, indent=2)
        f.write("\n")
    print(json.dumps({"rows": 134, "nonzero": body["nonzero_row_count"], "zero": body["zero_row_count"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["inventory", "image-graph", "beta-zero"])
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--component")
    args = parser.parse_args()
    if args.mode == "inventory":
        inventory(args.output)
    elif args.mode == "beta-zero":
        beta_zero(args.output)
    else:
        if not args.component:
            parser.error("image-graph requires --component")
        graph(args.component, args.output)
