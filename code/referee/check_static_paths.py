#!/usr/bin/env python3
"""Resolve the finite path candidates reported by the dependency scanner.

Checks source AST sites and real target bytes. This is a bounded review of
reported candidates, not an exhaustive analysis of arbitrary generated code.
No path is waived merely because its name contains RUNTIME or CHECKPOINT.
"""
import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Source:
    def __init__(self, path):
        self.path = path
        self.text = (ROOT/path).read_text()
        self.tree = ast.parse(self.text)
        self.parents = {child: node for node in ast.walk(self.tree)
                        for child in ast.iter_child_nodes(node)}
        self.strings = [n for n in ast.walk(self.tree)
                        if isinstance(n, ast.Constant) and isinstance(n.value, str)]
        self.scope_cache = {}

    def sites(self, path):
        return [n for n in self.strings
                if n.value == path or path.endswith('/'+n.value)]

    def enclosing(self, node, types):
        while not isinstance(node, types) and node in self.parents:
            node = self.parents[node]
        return node

    def assignment_methods(self, node):
        assignment = self.enclosing(node, ast.Assign)
        assert isinstance(assignment, ast.Assign), (self.path, node.lineno)
        names = {n.id for n in assignment.targets if isinstance(n, ast.Name)}
        assert names, (self.path, assignment.lineno)
        scope = self.enclosing(assignment, (ast.FunctionDef, ast.Module))
        if scope not in self.scope_cache:
            calls = [n for n in ast.walk(scope) if isinstance(n, ast.Call)
                     and isinstance(n.func, ast.Attribute)
                     and isinstance(n.func.value, ast.Name)]
            self.scope_cache[scope] = calls
        methods = {n.func.attr for n in self.scope_cache[scope] if n.func.value.id in names}
        return assignment, names, methods


def resolve(path, source, original):
    sites = source.sites(path)
    assert sites, (source.path, path)
    lines = sorted({n.lineno for n in sites})
    target = None
    if '_RUNTIME_' in path or '_CHECKPOINT_' in path:
        assert original == 'code/classification_plane_rank_zero_fullrank_weights.py'
        for node in sites:
            _, _, methods = source.assignment_methods(node)
            assert 'write_text' in methods and methods <= {'write_text', 'exists', 'unlink'}, (path, methods)
        kind, reason = 'producer_output', 'All bound-variable file calls in the enclosing scope are write_text, exists or unlink; there is no read call.'
    elif path.startswith('audits/'):
        assert original == 'code/normal_actual_eight_support_audit_20261103_v62.py'
        assert "base=ROOT/'legacy/v22'" in source.text
        target = 'legacy/v22/'+path
        assert all(isinstance(source.enclosing(n, ast.Call), ast.Call)
                   and source.enclosing(n, ast.Call).func.attr == 'read_text' for n in sites)
        kind, reason = 'resolved_input', 'The runtime base is ROOT/legacy/v22, not the repository root.'
    elif path in ['code/audit_tree_jacobi_sparse.py', 'code/build_p4_tree_linear.py']:
        assert original == 'legacy/v21/code/run_checks.py'
        assert 'bundle=Path(__file__).resolve().parents[1]' in source.text
        target = 'legacy/v21/'+path
        kind, reason = 'resolved_input', 'The subprocess script is relative to the legacy/v21 bundle root.'
    elif path == 'code/poisson_p4_components_check.py':
        assert original == 'code/classification_plane_quintic_saito.py'
        assert all(isinstance(source.parents[n], ast.keyword) and source.parents[n].arg == 'archive_member'
                   for n in sites)
        kind, reason = 'archive_locator', 'Only archive_member metadata; the imported helper uses the separate full legacy/v03 path.'
    elif path.startswith('legacy/v02/'):
        assert original in ['legacy/v02/code/build_symplectic_slice.py', 'legacy/v02/code/build_universal_slice.py']
        for node in sites:
            call = source.enclosing(node, ast.Call)
            assert isinstance(call.func, ast.Attribute)
            assert call.func.attr == 'write_text' or (call.func.attr == 'open'
                and call.args and isinstance(call.args[0], ast.Constant) and call.args[0].value == 'w')
        kind, reason = 'producer_output', 'The path expression is written with write_text or open(w); ROOT is an output directory selected by CLI.'
    elif path in ['legacy/v22/22_node_heisenberg.json', 'legacy/v22/22_node_spectral.json']:
        assert original == 'legacy/v22/code/audit_certificates_stdlib.py'
        assert "default=BASE/'data/22'" in source.text and 'root=args.data_directory' in source.text
        target = 'legacy/v22/data/22/'+Path(path).name
        kind, reason = 'resolved_input', 'The runtime data-directory default is legacy/v22/data/22.'
    elif path == 'legacy/v22/certificate_audit_stdlib.json':
        assert original == 'legacy/v22/code/audit_certificates_stdlib.py'
        for node in sites:
            _, _, methods = source.assignment_methods(node)
            assert methods == {'write_text'}, methods
        kind, reason = 'producer_output', 'Only the destination write; --output is required by this CLI.'
    elif path == 'reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SQUAREFREE_20261107.json':
        assert original == 'code/classification_plane_rank_zero_fullrank_weights.py'
        for node in sites:
            _, _, methods = source.assignment_methods(node)
            assert 'write_text' in methods and methods <= {'write_text', 'exists'}, methods
        kind, reason = 'producer_output', 'Destination of the bounded squarefree batch; checked for nonexistence, then written.'
    elif path == 'work/old-plane-v62/weight6-review-fresh/report.json':
        assert original == 'code/q_weight6_two_axes_review_20261103_v63.py'
        kind, reason = 'generated_prerequisite', 'A real read input, generated by the documented weight6 producer before the review command; not an external or unrecoverable asset.'
    elif path.startswith('work/ga9-'):
        if Path(path).name in ['failure.txt', 'resource.json']:
            for node in sites:
                call = source.enclosing(node, ast.Call)
                # One duplicate is inside an exception message, not an I/O call.
                assert (isinstance(call.func, ast.Attribute) and call.func.attr == 'write_text') or (
                    isinstance(call.func, ast.Name) and call.func.id == 'str')
            kind, reason = 'producer_output', 'Diagnostic text or timeout/failure report written by this runner.'
        else:
            for node in sites:
                _, names, _ = source.assignment_methods(node)
                assert 'output_path' in names
            assert 'subprocess.run' in source.text
            assert 'str(output_path)' in source.text and ('write(' in source.text or 'write(' in source.text.replace(' ', ''))
            kind, reason = 'generated_native_output', 'The runner emits a native input, executes it, then reads its generated text output; no preexisting coefficient file is consumed at this path.'
    else:
        raise ValueError(('Unreviewed path', source.path, path))
    result = {'disposition': kind, 'reason': reason, 'source_lines': lines}
    if target:
        assert (ROOT/target).is_file(), target
        result['required_file'] = {'path': target, 'sha256': sha(ROOT/target)}
    if kind == 'generated_prerequisite':
        producer = 'code/q_global_weight6_two_axes_20261103_v63.py'
        result['producer'] = {'path': producer, 'sha256': sha(ROOT/producer),
            'command': 'python '+producer+' --output work/old-plane-v62/weight6-review-fresh'}
    return result


def main(audit_path, output):
    audit = json.loads(audit_path.read_text())
    pairs = sorted({(r['origin'], r['path']) for r in audit['excluded_nonpremises']
                    if r.get('reason', '').startswith('Unresolved syntactic path candidate')})
    assert pairs
    frozen = json.loads((ROOT/'research/reproduction/frozen_producer_index.json').read_text())['records']
    originals = {r['materialized_path']: r['original_path'] for r in frozen}
    sources, rows = {}, []
    for origin, path in pairs:
        if origin not in sources:
            sources[origin] = Source(origin)
        resolution = resolve(path, sources[origin], originals.get(origin, origin))
        rows.append({'origin': origin, 'path': path, **resolution})
    result = {'scope': __doc__, 'reported_paths': len({p for _,p in pairs}),
              'origin_path_pairs': len(pairs), 'unreviewed_pairs': 0,
              'source_sha256': {p: sha(ROOT/p) for p in sources}, 'resolutions': rows,
              'disposition_counts_by_unique_path': dict(Counter(
                  next(r['disposition'] for r in rows if r['path'] == p) for p in sorted({p for _,p in pairs}))),
              'script_sha256': sha(Path(__file__))}
    # Two historical failed-execution records also occur in JSON metadata.
    # Resolve these by exact origin/digest and the replacement certificate;
    # do not waive all missing files whose names happen to contain RUNTIME.
    manifest_path = 'research/classification/inputs/plane_m4_two_line_casimir_manifest_20261107.json'
    job_path = 'research/classification/jobs/plane_m4_two_line_casimir.json'
    manifest = json.loads((ROOT/manifest_path).read_text())
    job = json.loads((ROOT/job_path).read_text())
    snapshot = manifest['producer_snapshot']
    assert sha(ROOT/snapshot['path']) == snapshot['sha256']
    failed = []
    for key in ['pre_execution_implementation_correction', 'serialization_correction']:
        record = job[key]
        path, digest = record['runtime_receipt'], record['runtime_receipt_sha256']
        archived = next(r for r in manifest['stable_files'] if r['path'] == path)
        assert archived['sha256'] == digest
        assert 'No source matrix or mathematical computation changed' in record['reason'] or 'Failed runtime producer' in record['reason']
        failed.append({'path': path, 'sha256': digest, 'origins': [manifest_path, job_path],
            'reason': record['reason'], 'replacement_certificate': {'path': snapshot['path'], 'sha256': snapshot['sha256']}})
    result['diagnostic_metadata_exclusions'] = failed
    result['source_sha256'].update({p: sha(ROOT/p) for p in [manifest_path, job_path]})
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['reported_paths', 'origin_path_pairs', 'unreviewed_pairs',
                                           'disposition_counts_by_unique_path']}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    main(args.audit, args.output)
