#!/usr/bin/env python3
"""Freeze a measured new paper and its first-level immutable proof dependencies.

No CAS is run and no research ledger, source index, or old manuscript is edited.
Run after the final LuaLaTeX pass, from a checkout retaining that work build.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
ID = 'CA_CB_LOCAL_NORMALIZATION_PAPER_20261103_v63'
SOURCE = 'papers/ca_cb_local_normalization_ja.tex'
PDF = 'papers/ca_cb_local_normalization_ja.pdf'
SCOPE_NOTE = 'research/notes/CA_CB_LOCAL_NORMALIZATION_PAPER_20261103_v63.md'
RECEIPTS = [
    'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
    'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json',
    'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json',
    'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.json',
    'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json',
    'reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63_REVIEW.json',
    'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_CA_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_ORIGIN_BINARY_CUBIC_RECOGNITION_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_THREE_AXIS_CONDUCTOR_RIGIDITY_20261103_v63.json',
    'reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_ORIGIN_FORMAL_BAD_LOCI_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_ORIGIN_FORMAL99_M2_20261103_v63/report.json',
    'reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW.json',
    'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63_REVIEW.json',
    'reports/AMBIENT_CA_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW.json',
    'reports/AMBIENT_CA_SELF_PAIR_CB_WHOLE_FIBER_20261103_v63_REVIEW.json',
    'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_CA_SELF_PAIR_MEETING_DOUBLE_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW.json',
    'reports/AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json',
    'reports/AMBIENT_INTERSECTION_DIMENSION27_20261103_v62_REVIEW.json',
    'reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.json',
]
DIRECT_EVIDENCE = [
    'reports/AMBIENT_CB_ACTUAL_D_CUBIC_SYZYGY_20261103_v63_REVIEW/independent-source-curves.json.gz',
    'code/restore_CB_cubic_worker_inputs_20261103_v63.py',
    'research/notes/AMBIENT_CB_CUBIC_REPRODUCIBILITY_20261103_v63.md',
    'reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63/relations99.json.gz',
    'reports/AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63/certificate.json.gz',
]
ESTABLISHED = [
    'Individual reduced CA/CB full formal completions at the common rank2 point',
    'Conductor and quadratic-cover rigidity, G19 cover product, new Cohen base and intrinsic pushout',
    'Complete99 kernel, without a minimality assertion',
    'Explicit proper seminormalization A[wt,rt] and complete81 kernel',
    'Formal bad loci, normalization fibers and transverse invariants',
    'Two actual rank2 polynomial families Pa and Qb, with all-nonzero flag counts',
    'Full actual D types on nonempty parameter-and-group unit opens',
    'Pa nil module generic length1 and Qb length4, with distinct primary quotient lengths',
    'Two distinct projective19 global embedded associated-component closures of actual D',
]
EXCLUDED = [
    'No new reduced full-Poisson component; lower58 is unchanged',
    'No full relative CA/CB intersection type at the common rank2 base point',
    'No global exhaustion or finite component upper bound',
    'No convergent holomorphic equivalence or literature-priority assertion',
    'No all-nonzero-parameter extension of unit-open actual D normal forms',
]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def main(build, output):
    output.mkdir(parents=True, exist_ok=True)
    frozen = json.loads((build / 'frozen_before.json').read_text())
    after = {name: sha(ROOT / name) for name in frozen}
    assert after == frozen, 'An old frozen manuscript/PDF/historical index changed.'
    source = ROOT / SOURCE
    pdf = ROOT / PDF
    log_path = build / 'ca_cb_local_normalization_ja.log'
    log = log_path.read_text(errors='replace')
    failures = {
        'overfull_boxes': bool(re.search(r'Overfull \\[hv]box', log)),
        'undefined_reference_or_citation': bool(re.search(
            r'LaTeX Warning: (?:Reference|Citation).*undefined|There were undefined (?:references|citations)', log)),
        'changed_labels_or_rerun_warning': bool(re.search(
            r'Label\(s\) may have changed|Rerun to get|Package rerunfilecheck Warning', log)),
        'latex_error': bool(re.search(r'^!|Emergency stop|Fatal error', log, re.M)),
        'hyperref_warning': 'Package hyperref Warning' in log,
    }
    assert not any(failures.values()), failures
    assert 'Output written on ca_cb_local_normalization_ja.pdf (17 pages' in log
    # Copy only the final PDF belonging to this new manuscript.
    shutil.copyfile(build / 'ca_cb_local_normalization_ja.pdf', pdf)
    info = subprocess.check_output(['pdfinfo', str(pdf)], text=True)
    pages = int(re.search(r'^Pages:\s+(\d+)\s*$', info, re.M).group(1))
    assert pages == 17 and 'Form:            none' in info and 'JavaScript:      no' in info
    assert pdf.stat().st_size == 428342
    subprocess.run(['pdftotext', '-layout', str(pdf), str(output / 'paper-text.txt')], check=True)
    text = (output / 'paper-text.txt').read_text()
    assert text.count('\f') == pages
    checked = {}
    statuses = {}
    issues = []
    for name in RECEIPTS:
        receipt = ROOT / name
        data = json.loads(receipt.read_text())
        statuses[name] = data.get('all_passed') is True
        if not statuses[name]:
            issues.append([name, 'all_passed is not true'])
        for key in ('input_sha256', 'evidence_sha256'):
            for dep, expected in data.get(key, {}).items():
                if not isinstance(expected, str) or not re.fullmatch(r'[0-9a-f]{64}', expected):
                    continue
                candidates = [ROOT / dep, receipt.parent / dep]
                matches = [p for p in candidates if p.is_file() and sha(p) == expected]
                if not matches:
                    issues.append([name, key, dep, 'missing or hash mismatch'])
                else:
                    checked[relative(matches[0])] = expected
    assert not issues, json.dumps(issues, ensure_ascii=False)
    for name in DIRECT_EVIDENCE:
        assert (ROOT / name).is_file(), name
    for name in ('frozen_before.json',):
        shutil.copyfile(build / name, output / name)
    shutil.copyfile(log_path, output / 'lualatex-final.log')
    shutil.copyfile(build / 'lualatex-pass4.stdout', output / 'lualatex-final.stdout')
    for name in ('page1.png', 'page16.png'):
        shutil.copyfile(build / name, output / name)
    (output / 'pdfinfo.txt').write_text(info)
    lualatex = subprocess.check_output(['lualatex', '--version'], text=True).splitlines()[0]
    font_warnings = [line.strip() for line in log.splitlines() if
                     'LaTeX Font Warning:' in line or
                     'Package luatexja-preset Warning:' in line or
                     'Package luatexja-fontspec Warning:' in line]
    manifest = {
        'id': ID, 'all_passed': True,
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'repository_head_when_manifest_created': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'source': SOURCE, 'pdf': PDF, 'scope_note': SCOPE_NOTE,
        'source_sha256': sha(source), 'pdf_sha256': sha(pdf),
        'scope_note_sha256': sha(ROOT / SCOPE_NOTE),
        'manifest_generator_sha256': sha(Path(__file__)),
        'language': 'Japanese LuaLaTeX', 'new_standalone_paper': True,
        'measured_pdf_pages': pages, 'pdf_bytes': pdf.stat().st_size,
        'lualatex_version': lualatex, 'pdf_has_forms': False, 'pdf_has_javascript': False,
        'established_scope': ESTABLISHED, 'excluded_scope': EXCLUDED,
        'receipt_all_passed': statuses,
        'receipt_sha256': {name: sha(ROOT / name) for name in RECEIPTS},
        'checked_first_level_dependency_sha256': checked,
        'direct_evidence_sha256': {name: sha(ROOT / name) for name in DIRECT_EVIDENCE},
        'historical_frozen_sha256': after, 'historical_frozen_file_count': len(after),
        'historical_frozen_files_unchanged': True,
        'saved_M2_validation': {
            'version': '1.24.11', 'backend': 'default GB and full source elimination',
            'entire99_kernel_checked': True, 'rerun_for_this_paper': False,
            'not_a_geometric_application_proof': True,
            'relation_minimality_asserted': False,
        },
        'heavy_CAS_rerun_for_paper': False,
        'build_problem_flags': failures,
        'retained_font_warnings': font_warnings,
        'font_information_scope': 'Japanese italic shape INFO and small math-font substitution; not unresolved references',
        'visually_checked_pages': [1, 16],
        'build_artifact_sha256': {name: sha(output / name) for name in (
            'lualatex-final.log', 'lualatex-final.stdout', 'pdfinfo.txt',
            'paper-text.txt', 'page1.png', 'page16.png', 'frozen_before.json')},
        'source_index_mutated_by_this_agent': False,
        'ledgers_mutated_by_this_agent': False,
        'external_publication_performed': False,
        'commit_or_push_performed_by_this_agent': False,
    }
    (output / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    registration = {
        'request': 'Register a new standalone manuscript in RESEARCH_SOURCE_INDEX_v63.json; preserve historical SOURCE_INDEX.json',
        'source': SOURCE, 'pdf': PDF, 'language': 'Japanese',
        'measured_pages': pages, 'source_sha256': manifest['source_sha256'],
        'pdf_sha256': manifest['pdf_sha256'],
        'claim_ids_to_be_selected_by_integration_writer': True,
        'accepted_receipts': RECEIPTS, 'scope': ESTABLISHED, 'excluded_scope': EXCLUDED,
        'old_nine_frozen_files_unchanged': True,
    }
    (output / 'SOURCE_INDEX_REGISTRATION_REQUEST.json').write_text(
        json.dumps(registration, ensure_ascii=False, indent=2) + '\n')
    report = f'''# CA/CB 局所正規化・交叉原稿の保存

新規原稿: `{SOURCE}` / `{PDF}`。最終 PDF は **{pages} 頁**、{pdf.stat().st_size} bytes。
LuaLaTeX: `{lualatex}`。A4、form/JavaScript なし。

- source SHA256: `{manifest['source_sha256']}`
- PDF SHA256: `{manifest['pdf_sha256']}`
- 受理済み有限入力・独立監査: {len(RECEIPTS)} report、全て `all_passed=true`。
- 第一階層の入力・証拠 SHA256: {len(checked)} ファイル一致。
- 旧原稿4組と歴史的 index、計 {len(after)} ファイルは作成前と同一。
- 全参照・引用解決、overfull なし、再実行要求なし。1・16頁を画像確認。

個別 CA/CB の全形式型、完全な99関係、seminormalizationと81関係、形式bad loci、
二つの実自己旗族、非空 unit open 上の全 actual D と相異なる射影19次元embedded閉包を収録。
二つの D の nil 長1/4と embedded primary 商長2/6を分けた。
共通基点の相対的な全 D、global分類、収束同値、優先性は主張しない。下界58は不変。

Macaulay2 1.24.11 default GB による既存の全99 kernel検査を保存証明書から参照し、
原稿作成で重いCASは再実行していない。独立三次 raw source jets と fresh-checkout用
復元 helper も manifest に直接結び付けた。font警告は小さな数式字体の代替等であり、
最終 log に保持した。既存論文、台帳、indexはこの担当で編集せず、登録依頼だけを保存した。
外部投稿は行っていない。

再検査入口: `.venv/bin/python code/ca_cb_local_paper_manifest_20261103_v63.py`。
これは指定済み work 内の完成 PDF/build を保存し、入力を検査する。新しい checkoutでの
組版は、既設環境を activate して source を LuaLaTeX で work 出力へビルドする。
数学的証明とCASの再実行手順は原稿付録および各不変証明書を参照する。
'''
    (output / 'REPORT.md').write_text(report)
    print(json.dumps({'all_passed': True, 'pages': pages,
                      'receipt_count': len(RECEIPTS), 'checked_dependency_count': len(checked),
                      'source_sha256': manifest['source_sha256'],
                      'pdf_sha256': manifest['pdf_sha256'],
                      'output': relative(output)}, ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path,
                        default=ROOT / 'work/ca-cb-local-paper-20261103')
    parser.add_argument('--output', type=Path, default=ROOT / 'reports' / ID)
    args = parser.parse_args()
    main(args.build.resolve(), args.output.resolve())
