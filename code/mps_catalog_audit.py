#!/usr/bin/env python3
"""Match the archived 40 records to the actual MPS Table 2 diagrams.

The optional --source-html re-extracts the factual diagram data from the arXiv
2010.08692v2 HTML (Table S5.T2). Otherwise the checked extraction is used.
All mathematical linear algebra is exact over QQ. SVG decimals only identify
vertices within each drawing; they are never used as mathematical coefficients.
This audits a literature crosswalk, not the exhaustion of all Poisson structures.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from fractions import Fraction
import hashlib
from itertools import combinations, permutations
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json"
EXTRACTION = ROOT / "code/mps_table2_diagrams.json"
SOURCE_URL = "https://arxiv.org/html/2010.08692v2#S5.T2"
PAIRS = list(combinations(range(5), 2))
PARAM_PAIRS = list(combinations(range(1, 5), 2))
PERMUTATIONS = list(permutations(range(5)))
NUMBER = r"-?\d+(?:\.\d+)?"


def path_vertices(path):
    return [tuple(Fraction(t) for t in pair) for pair in
            re.findall(r"[ML] (" + NUMBER + r") (" + NUMBER + r")", path)]


def extract_html(path):
    raw = path.read_bytes()
    html = raw.decode("utf-8")
    start = html.index('<table id="S5.T2.2"')
    table = ET.fromstring(html[start:html.index("</table>", start)+8])
    diagrams = []
    for row in table.iter("tr"):
        cells = [x for x in row if x.tag == "td"]
        if not cells:
            continue
        printed_dim = "".join(cells[-1].itertext()).strip()
        for col, svg in enumerate(row.iter("svg"), 1):
            centers = []
            for group in svg.iter("g"):
                if group.get("fill") == "#CCCCCC":
                    centers.extend(path_vertices(p.get("d"))[0]
                                   for p in group if p.tag == "path")
            vertices = sorted(set(centers))
            assert len(vertices) == 5, (svg.get("id"), vertices)
            labels = {point: i for i, point in enumerate(vertices)}
            clips = {x.get("id"): path_vertices(x.find("path").get("d"))
                     for x in svg.iter("clipPath")}
            edges = set()
            zero_edges = set()
            orders = Counter()
            for group in svg.iter("g"):
                # Blue solid line segments are smoothable edges; white dashed
                # overlays erase portions of a grey edge to indicate B_ij=0.
                if group.get("stroke") == "#0000FF" and group.get("stroke-width") == "1.2pt":
                    for p in group:
                        if p.tag == "path":
                            ends = path_vertices(p.get("d"))
                            assert len(ends) == 2
                            edges.add(tuple(sorted(labels[x] for x in ends)))
                if group.get("stroke-dasharray"):
                    assert group.get("stroke") == "#FFFFFF"
                    for p in group:
                        if p.tag == "path":
                            ends = path_vertices(p.get("d"))
                            assert len(ends) == 2
                            zero_edges.add(tuple(sorted(labels[x] for x in ends)))
                if group.get("clip-path"):
                    clip = group.get("clip-path")[5:-1]
                    triangle = clips[clip]
                    assert len(triangle) == 3
                    paths = list(group.iter("path"))
                    assert len(paths) == 1
                    vertex = path_vertices(paths[0].get("d"))[0]
                    assert vertex == triangle[1]
                    edge = tuple(sorted(labels[x] for x in (triangle[0], triangle[2])))
                    orders[(edge, labels[vertex])] += 1
            assert all(edge in edges for edge, vertex in orders)
            records = [{"edge": list(edge), "orders":
                        {str(k): orders[(edge, k)] for k in range(5) if k not in edge}}
                       for edge in sorted(edges)]
            assert all(sum(item["orders"].values()) == 2 for item in records)
            diagrams.append({
                "svg_id": svg.get("id"), "html_row_id": row.get("id"),
                "column_one_based": col,
                "printed_period_dimension": int(printed_dim) if printed_dim else None,
                "svg_vertex_coordinates": [[str(z) for z in x] for x in vertices],
                "edges": records, "forced_zero_entries": [list(e) for e in sorted(zero_edges)],
            })
    assert len(diagrams) == 40
    return {"source_url": SOURCE_URL, "source_html_sha256": hashlib.sha256(raw).hexdigest(),
            "extraction_method": "SVG blue segments, triangle clip paths and angle-fill multiplicities; white dashed zero-biresidue overlays",
            "vertex_labels": "Each SVG separately: ascending lexicographic pairs of its five grey vertex centres",
            "diagrams": diagrams}


def decorated_key(diagram, perm=tuple(range(5))):
    result = []
    for item in diagram["edges"]:
        edge = tuple(sorted(perm[i] for i in item["edge"]))
        orders = {perm[int(k)]: m for k, m in item["orders"].items()}
        result.append(edge + tuple(orders[k] for k in range(5) if k not in edge))
    return tuple(sorted(result))


def canonical_key(diagram):
    return min(decorated_key(diagram, p) for p in PERMUTATIONS)


def biresidue_matrix():
    variables = s.symbols("b12 b13 b14 b23 b24 b34")
    matrix = s.zeros(5)
    for (i, j), variable in zip(PARAM_PAIRS, variables, strict=True):
        matrix[i, j] = variable
        matrix[j, i] = -variable
    for i in range(1, 5):
        matrix[0, i] = sum(matrix[i, j] for j in range(1, 5))
        matrix[i, 0] = -matrix[0, i]
    assert matrix * s.ones(5, 1) == s.zeros(5, 1)
    return variables, matrix


def audit_record(record):
    variables, matrix = biresidue_matrix()
    equations = []
    weights = []
    for item in record["edges"]:
        i, j = item["edge"]
        weight = [-1 if k in (i, j) else item["orders"][str(k)] for k in range(5)]
        assert sum(weight) == 0
        weights.append(weight)
        equations.extend(matrix[j, int(k)] + matrix[int(k), i] - m*matrix[i, j]
                         for k, m in item["orders"].items())
    coefficient = s.Matrix([[s.expand(e).coeff(v) for v in variables] for e in equations])
    if not equations:
        coefficient = s.zeros(0, 6)
    kernel = s.Matrix.hstack(*coefficient.nullspace())
    assert kernel.cols == record["period_dimension"]
    archived_basis = s.Matrix(record["period_basis"]).applyfunc(s.Rational).T
    assert archived_basis.shape == kernel.shape
    assert coefficient * archived_basis == s.zeros(coefficient.rows, kernel.cols)
    assert archived_basis.rank() == kernel.cols
    substitutions = dict(zip(variables, (sum(kernel[i,j]*s.Symbol(f"t{j}") for j in range(kernel.cols)) for i in range(6)), strict=True))
    restricted = matrix.subs(substitutions)
    forced_zero = [list(pair) for pair in PAIRS if restricted[pair[0], pair[1]] == 0]
    assert forced_zero == record["forced_zero_entries"]
    sample = s.Matrix(record["biresidue_example"])
    assert sample + sample.T == s.zeros(5) and sample*s.ones(5,1) == s.zeros(5,1)
    assert sample.rank() == 4
    sampled_edges = []
    for i, j in PAIRS:
        if sample[i, j] == 0:
            continue
        orders = {str(k): (sample[j,k]+sample[k,i])/sample[i,j]
                  for k in range(5) if k not in (i, j)}
        if all(m.is_integer and m >= 0 for m in orders.values()):
            sampled_edges.append({"edge": [i,j], "orders": {k: int(m) for k,m in orders.items()}})
    assert sampled_edges == record["edges"]
    weight_rank = s.Matrix(weights).rank() if weights else 0
    assert weight_rank == record["smoothing_weight_rank"]
    period_dim, edge_count = kernel.cols, len(record["edges"])
    assert 19 + period_dim + edge_count == record["projective_component_dimension"]
    assert 20 + period_dim + edge_count == record["affine_component_dimension"]
    return {"period_dimension": period_dim, "edge_count": edge_count,
            "weight_rank": weight_rank, "linear_constraint_rank": coefficient.rank(),
            "period_rref": [[str(x) for x in row] for row in coefficient.rref()[0].tolist()],
            "s61_dimension": record["projective_component_dimension"] - 9}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-html", type=Path)
    parser.add_argument("--write-extraction", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "reports/MPS40_CROSSWALK_20261031.json")
    args = parser.parse_args()
    extracted = extract_html(args.source_html) if args.source_html else json.loads(EXTRACTION.read_text())
    if args.source_html and EXTRACTION.exists():
        assert extracted == json.loads(EXTRACTION.read_text()), "HTML extraction changed"
    if args.write_extraction:
        with args.write_extraction.open("x") as stream:
            json.dump(extracted, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    archived = json.loads(ARCHIVE.read_text())["diagrams"]
    primary = {canonical_key(d): d for d in extracted["diagrams"]}
    assert len(primary) == 40
    assert len({canonical_key(d) for d in archived}) == 40
    assert set(primary) == {canonical_key(d) for d in archived}
    csv_path = ARCHIVE.with_name("p4_toric_components_40.csv")
    with csv_path.open(newline="") as stream:
        csv_records = {r["id"]: r for r in csv.DictReader(stream)}
    output = []
    proposed_components = []
    for record in archived:
        exact = audit_record(record)
        drawing = primary[canonical_key(record)]
        perm = next(p for p in PERMUTATIONS if decorated_key(record, p) == decorated_key(drawing))
        mapped_zeros = sorted(sorted(perm[i] for i in e) for e in record["forced_zero_entries"])
        assert mapped_zeros == drawing["forced_zero_entries"]
        if drawing["printed_period_dimension"] is not None:
            assert drawing["printed_period_dimension"] == record["period_dimension"]
        csv_record = csv_records[record["id"]]
        for field in ("period_dimension", "edge_count", "smoothing_weight_rank", "projective_component_dimension", "affine_component_dimension"):
            assert int(csv_record[field]) == record[field]
        assert json.loads(csv_record["edge_orders"]) == record["edges"]
        assert json.loads(csv_record["biresidue_example"]) == record["biresidue_example"]
        output.append({
            "id": record["id"], "graph_type": record["graph_type"],
            "source_svg_id": drawing["svg_id"], "source_html_row_id": drawing["html_row_id"],
            "source_column_one_based": drawing["column_one_based"],
            "source_url": SOURCE_URL.replace("#S5.T2", "#" + drawing["svg_id"]),
            "archive_to_svg_vertex_permutation": list(perm),
            "edges": record["edges"], "forced_zero_entries": record["forced_zero_entries"],
            "projective_dimension": record["projective_component_dimension"],
            "named_family": record.get("named_family"),
            "biresidue_example": record["biresidue_example"],
            **exact,
            "generic_primitive_first_jet_dominance": "not_decided_by_this_audit" if exact["s61_dimension"] >= 16 else "excluded_by_dimension",
            "s61_coordinate_map_status": "not_constructed_by_this_audit",
        })
        if record["id"] in ("T39", "T40"):
            tilde_k = 2 if record["id"] == "T39" else 1
            first_row = [0] + [1]*tilde_k + [0]*(4-2*tilde_k) + [-1]*tilde_k
            cyclic = s.Matrix(5, 5, lambda i,j: first_row[(j-i) % 5])
            sample = s.Matrix(record["biresidue_example"])
            witnesses = [(p, scale) for p in PERMUTATIONS for scale in (1, -1)
                         if sample == scale*cyclic.extract(p, p)]
            assert witnesses
            p, scale = witnesses[0]
            k = (pow(1+tilde_k, -1, 5)-1) % 5
            assert record["named_family"] == f"q_5,{k}"
            output[-1]["feigin_odesskii_identification"] = {
                "source_url": "https://arxiv.org/html/2310.05284v1#S4.Thmtheorem3",
                "source_anchor": "Proposition 4.2 and Theorem 4.3",
                "cyclic_matrix_C_5_tilde_k": tilde_k,
                "permutation": list(p), "scale": scale, "q_5_k": k,
            }
        proposed_components.append({
            "id": record["id"],
            "name": record.get("named_family", "MPS smoothing diagram " + record["id"]),
            "counted_in_global_lower_bound": True, "count": 1,
            "category": "MPS_SNC_meeting_generic_rank_four_full_component",
            "generic_rank": 4,
            "projective_dimension": record["projective_component_dimension"],
            "s61_dimension": exact["s61_dimension"],
            "smoothing_diagram": {"graph_type": record["graph_type"],
                                  "edges": record["edges"],
                                  "forced_zero_entries": record["forced_zero_entries"]},
            "period_dimension": record["period_dimension"],
            "smoothing_parameter_count": len(record["edges"]),
            "toric_biresidue_example": record["biresidue_example"],
            "toric_biresidue_example_scope": "SNC boundary point, not a generic smoothed member",
            "meets_SNC_log_symplectic_locus": True,
            "claims": ["MPS40", "MPS40_CROSSWALK"],
            "proof_sources": [
                {"url": SOURCE_URL, "anchor": drawing["svg_id"]},
                {"path": "research/notes/MPS40_CROSSWALK_20261031.md"},
            ],
            "verification_records": ["reports/MPS40_CROSSWALK_20261031.json"],
            "s61_coordinate_map_status": "not_constructed_by_this_audit",
            "review_status": {"proof": "published_theorem", "primary_diagram_crosswalk": "exact_S5_match",
                              "full_independent_external_proof_review": False},
            "notes_ja": ["内部IDは文献の番号ではない。MPS40_BLOCKとの二重加算をしない。",
                         "個別図式・次元の台帳であり、全成分のS61定義イデアルを構築したわけではない。"],
        })
    candidates = [r["id"] for r in output if r["s61_dimension"] >= 16]
    assert candidates == ["T01", "T39", "T40"]
    report = {
        "all_passed": True, "field": "QQ", "id": "MPS40_CROSSWALK_20261031",
        "source_url": SOURCE_URL, "source_html_sha256": extracted["source_html_sha256"],
        "crosswalk_count": 40, "all_s5_classes_distinct": True,
        "primary_zero_biresidue_edges_match": True,
        "catalog_csv_and_json_agree": True,
        "projective_dimension_histogram": dict(sorted(Counter(r["projective_dimension"] for r in output).items())),
        "possible_dominant_mps_components_by_dimension": candidates,
        "mps40_records": output,
        "proposed_component_entries": proposed_components,
        "limitations": ["This checks the 40 published smoothing diagrams, not global exhaustion.",
                        "No complete S61 coordinate maps or component ideals are constructed.",
                        "Dimensions use the MPS unobstructed local model and the archived S61 dimension shift 9.",
                        "MPS v2 equation (5.11) has the opposite numerator to the zero-row-sum convention; see the note and Matviichuk 2310.05284v1 section 2."],
        "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (Path(__file__), ARCHIVE, csv_path, EXTRACTION)},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps({"all_passed": True, "matched_diagrams": 40, "candidates": candidates,
                      "output": str(args.output)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
