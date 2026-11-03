#!/usr/bin/env python3
"""Bind the new standalone paper to measured PDF and immutable proof receipts.

This changes neither the historical manuscript nor any claim/source register.
Run after the final LuaLaTeX pass and copying its PDF into papers/.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output):
    output.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'papers/positive_weight6_component_ja.tex'
    pdf = ROOT / 'papers/positive_weight6_component_ja.pdf'
    supplement = ROOT / 'papers/positive_weight6_initial21.m2'
    info = subprocess.check_output(['pdfinfo', str(pdf)], text=True)
    pages = int(re.search(r'^Pages:\s+(\d+)\s*$', info, re.M).group(1))
    before = json.loads((output / 'frozen_before.json').read_text())
    after = {path: sha(ROOT / path) for path in before}
    assert after == before, 'A historical frozen manuscript/index changed.'
    evidence = [
        'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW/report.json',
        'reports/Q_POSITIVE_WEIGHT6_SMALL_MODEL_20261103_v63_REVIEW/report.json',
        'reports/Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63/report.json',
        'reports/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json',
        'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json',
        'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63_REVIEW/report.json',
        'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/transform_audit_homogeneous.json',
        'reports/R4_POSITIVE_WEIGHT6_DISTINCTNESS_20261103_v63/report.json',
        'reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW/report.json',
        'reports/Q_NEW_FAMILY_GENERIC_JACOBIAN_RANK_20261103_v63_REVIEW.json',
        'reports/Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63_REVIEW.json',
        'reports/Q_NEW_FAMILY_GENERIC_QUADRATIC_RANK_20261103_v63_REVIEW.json',
        'reports/Q_NEW_FAMILY_TRANSVERSE_ARTIN_20261103_v63_REVIEW/report.json',
        'reports/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63_REVIEW/report.json',
        'reports/Q_NEW_FAMILY_LENGTH5_COMPLETION_20261103_v63.json',
        'reports/Q_POSITIVE_WEIGHT6_NONCM_20261103_v63.json',
        'reports/Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63/certificate.json.gz',
        'reports/Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63/certificate.json.gz',
    ]
    # Every receipt is inspected, not inferred merely from its presence.
    receipt_status = {}
    checked_input_hashes = {}
    for name in evidence:
        path = ROOT / name
        assert path.exists(), name
        if path.suffix == '.json':
            data = json.loads(path.read_text())
            assert data.get('all_passed') is True, name
            receipt_status[name] = True
            for dep, expected in data.get('input_sha256', {}).items():
                if isinstance(expected, str) and len(expected) == 64:
                    assert sha(ROOT / dep) == expected, (name, dep)
                    checked_input_hashes[dep] = expected
    manifest = dict(
        id='POSITIVE_WEIGHT6_COMPONENT_PAPER_20261103_v63',
        all_passed=True,
        source=str(source.relative_to(ROOT)),
        pdf=str(pdf.relative_to(ROOT)),
        supplement=str(supplement.relative_to(ROOT)),
        measured_pdf_pages=pages,
        manuscript_language='Japanese LuaLaTeX',
        manuscript_kind='New standalone computer-assisted research paper',
        established_scope=[
            'Explicit actual three-parameter Poisson family without N^2 constraint',
            'Actual projective irreducible component of dimension25',
            'Actual generic divisor3H+smooth rank5 quadric tangent toH',
            'Modular Jordan(3,1,1), generic orbit23 and connected stabilizerGa',
            'Shear-corrected two-parameter finite generic orbit cover',
            'Exact primitive firstjet image15 and generic S61 fiber dimension1',
            'Generic full-scheme nonreducedness, transverse embedding dimension2, exact length5',
            'At111: full projective completion C[[r25,u,v]]/(uv,u^3,v^3), CM type2',
            'AtP+: full Poisson local ring is not Cohen-Macaulay, from edim8/dim3/multiplicity<=5',
            'Distinct from previous57, hence global lower58, using accepted old catalogue',
        ],
        unresolved_scope=[
            'Global component exhaustion and any finite global upper bound',
            'Full finite equivalence relation on the two-parameter orbit cover',
            'Equation of the15-dimensional firstjet hypersurface and its fiber scheme',
            'Normality and singularities of the reduced component',
            'Complete branches atP+ and comprehensive literature priority',
        ],
        primary_claim_ids=[
            'R4_NEW25_POSITIVE_WEIGHT6_FULL_COMPONENT',
            'R4_NEW25_FIRSTJET_IMAGE_EXACT15_GENERIC_FIBER1',
            'Q_POSITIVE_RAY_TRUE_INITIAL_FOUR_LOCAL_UPPER25',
        ],
        source_sha256=sha(source), pdf_sha256=sha(pdf), supplement_sha256=sha(supplement),
        manifest_generator_sha256=sha(Path(__file__)),
        evidence_sha256={name: sha(ROOT / name) for name in evidence},
        receipt_all_passed=receipt_status,
        checked_receipt_input_sha256=checked_input_hashes,
        historical_frozen_sha256=after,
        historical_manuscripts_and_index_unchanged=True,
        auxiliary_initial21_check='Macaulay2 1.24.11: dim3, degree6, exact displayed Hilbert series',
        pdf_validation='LuaLaTeX completed; no undefined references, rerun request or overfull boxes; pages1 and10 visually inspected',
        source_index_mutated=False,
        commits_or_pushes_performed_by_this_agent=False,
    )
    (output / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    (output / 'pdfinfo.txt').write_text(info)
    request = dict(
        request='Register a new standalone v63 manuscript in RESEARCH_SOURCE_INDEX_v63.json; preserve historical SOURCE_INDEX.json',
        source=manifest['source'], pdf=manifest['pdf'], supplement=manifest['supplement'],
        language='Japanese', measured_pages=pages,
        source_sha256=manifest['source_sha256'], pdf_sha256=manifest['pdf_sha256'],
        primary_claim_ids=manifest['primary_claim_ids'],
        additional_accepted_review_receipts=[name for name in evidence if 'REVIEW' in name],
        frozen_historical_manuscripts_unchanged=True,
        scope=manifest['established_scope'], unresolved=manifest['unresolved_scope'],
    )
    (output / 'SOURCE_INDEX_REGISTRATION_REQUEST.json').write_text(json.dumps(request, ensure_ascii=False, indent=2) + '\n')
    print('NEW STANDALONE PAPER MANIFEST PASS; measured pages', pages)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output.resolve())
