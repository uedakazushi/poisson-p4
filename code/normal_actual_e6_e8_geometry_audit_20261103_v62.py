#!/usr/bin/env python3
"""Independent integer and incidence-tree coverage for normal actual divisors.

The geometry is proved in the accompanying note.  This bounded enumeration
checks the finite combinatorial steps, hashes audited archived inputs, and reads
saved exact certificate passes without rerunning unchanged older calculations.
"""
import argparse
import hashlib
import itertools as it
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]


def canonical(matrix,degrees):
    r,n=len(degrees),len(matrix[0])
    orders=[p for p in it.permutations(range(r)) if tuple(degrees[i] for i in p)==degrees]
    return min(tuple(matrix[i][j] for i in p for j in q)
               for p in orders for q in it.permutations(range(n)))


def enumerate_trees(degrees):
    r=len(degrees);answer=[]
    for n in range(1,r):
        found={}
        for bits in it.product([0,1],repeat=r*n):
            if sum(bits)!=r+n-1:continue
            matrix=[bits[i*n:(i+1)*n] for i in range(r)]
            if any(not 1<=sum(row)<=2 for row in matrix):continue
            if any(sum(matrix[i][j] for i in range(r))<2 for j in range(n)):continue
            adjacency=[[] for _ in range(r+n)]
            for i in range(r):
                for j in range(n):
                    if matrix[i][j]:adjacency[i].append(r+j);adjacency[r+j].append(i)
            visited={0};pending=[0]
            while pending:
                for vertex in adjacency[pending.pop()]:
                    if vertex not in visited:visited.add(vertex);pending.append(vertex)
            if len(visited)!=r+n:continue
            key=canonical(matrix,degrees)
            found.setdefault(key,matrix)
        answer.extend(dict(component_degrees=degrees,singular_points=n,
                           incidence=[list(row) for row in matrix],
                           component_mark_counts=[sum(row) for row in matrix],
                           singular_branch_counts=[sum(matrix[i][j] for i in range(r)) for j in range(n)])
                      for _,matrix in sorted(found.items()))
    return answer


def verify(output):
    start=time.monotonic();inputs={};archive_count=0
    def record(path):
        nonlocal archive_count
        p=ROOT/path;value=hashlib.sha256(p.read_bytes()).hexdigest();inputs[path]=value
        if path.startswith('legacy/'):
            version=Path(path).parts[1];local='/'.join(Path(path).parts[2:])
            manifest=ROOT/'legacy'/version/'SHA256SUMS.txt'
            hashes=dict((line.split(maxsplit=1)[1].lstrip('*'),line.split(maxsplit=1)[0])
                        for line in manifest.read_text().splitlines() if line.strip())
            assert hashes[local]==value,path
            archive_count+=1
        return p
    integer_solutions=[(a,b,c) for a in range(6) for b in range(5) for c in range(5)
                       if 8*a+9*b+10*c==40]
    assert integer_solutions==[(0,0,4),(5,0,0)]
    trees={''.join(map(str,degrees)):enumerate_trees(degrees)
           for degrees in [(3,1),(2,2),(2,1,1),(1,1,1,1)]}
    assert [len(trees[x]) for x in ['31','22','211','1111']]==[1,1,3,3]
    for item in it.chain.from_iterable(trees.values()):
        assert sum(item['singular_branch_counts'])-item['singular_points']==len(item['component_degrees'])-1
        item['complete_linear_system_dimension']=sum(d+1 for d in item['component_degrees'])-(len(item['component_degrees'])-1)
        assert item['complete_linear_system_dimension']==5
    for version in [20,21]:
        for p in sorted((ROOT/f'legacy/v{version}/sections').glob('*.tex')):record(str(p.relative_to(ROOT)))
        for base in ['proofs','code']:
            for p in sorted((ROOT/f'legacy/v{version}/{base}').glob('**/*.py')):record(str(p.relative_to(ROOT)))
    # The previous independent scope audit identifies the exact line theorem input.
    record('reports/NORMAL_ACTUAL_INACTIVE_LINE_SCOPE_AUDIT_20261103_v62.json')
    for path in ['legacy/v14/sections/normal_quintics.tex','legacy/v15/sections/normal_geometry.tex',
                 'legacy/v16/sections/smooth_curves.tex',
                 'legacy/v20/p4_poisson_completeness_v20.tex','legacy/v21/p4_poisson_completeness_v21.tex']:
        record(path)
    saved=[]
    for version in [20,21]:
        for p in sorted((ROOT/f'legacy/v{version}/data/reproduced').glob('**/*.json')):
            data=json.loads(record(str(p.relative_to(ROOT))).read_text())
            if 'all_passed' in data:
                assert data['all_passed'],p
                saved.append(str(p.relative_to(ROOT)))
    primary=[]
    for name,url,role in [
        ('pym-1507.05668.pdf','https://arxiv.org/pdf/1507.05668',
         'Definition4.2,Remark4.4,Theorem4.5 pp17-18 and fullJacob length statement p21; Theorem1.3'),
        ('ams-secant.pdf','https://www.cambridge.org/core/services/aop-cambridge-core/content/view/ACDCB8F4BE19E7640B3DEB608F66564F/S0013091513000497a.pdf/div-class-title-a-degree-formula-for-secant-varieties-of-curves-div.pdf',
         'Proposition2.5 p309; Theorem3.5 p312, used only for smooth nondegenerate components'),
        ('pym-survey.html','https://pmc.ncbi.nlm.nih.gov/articles/PMC5818673/',
         'Theorem5.12 explicitly requires modularresidue nowherezero; not invoked from mere genericactivity')]:
        cached=ROOT/'work/cloud-d0-20261103'/name
        primary.append(dict(url=url,sha256=hashlib.sha256(cached.read_bytes()).hexdigest(),
                            bytes=cached.stat().st_size,read_for=role,cache_path=str(cached.relative_to(ROOT))))
    record(str(Path(__file__).relative_to(ROOT)))
    report=dict(id='NORMAL_ACTUAL_E6_E8_GEOMETRY_AUDIT_20261103_v62',all_passed=True,
                integer_degree_solutions=integer_solutions,eight_incidence_trees=trees,
                incidence_tree_counts={key:len(value) for key,value in trees.items()},
                archived_checksum_matches=archive_count,saved_exact_reports_read=saved,
                unchanged_old_certificates_rerun=False,
                primary_sources_opened_and_read=primary,
                mathematical_scope='Normal actual divisor: inactive exclusion through E6completion and E8reduction to the eightP4trees; v22finaleight exclusions separate.',
                clarified_dependencies=['Pymactivepoint means normalD+nonzeroZ in dimension4',
                                        'actual local product gives FULL Jacobian lengths8,9,10',
                                        'smooth secant formula used after singular directrix is excluded',
                                        'derivation lift in characteristic0 plus conductor preservation gives zeros at all singular preimages',
                                        'Hilbert-function bounds do not assume ACM',
                                        'same Hilbert-function argument replaces the extra Castelnuovo/RR input for singular integral space quartics',
                                        'genus0 and delta equality imply equal-constant ordinary multibranch points',
                                        'markedP1components at most2 marks and tree Picard gluing imply one projective orbit per incidence type'],
                unresolved_geometry_gap=None,
                limits=['The report does not freshly execute or re-prove the archived coefficient certificates.',
                        'The main note supplies the human proof audit; the enumeration only proves its finite combinatorial steps.',
                        'The v16 independent exclusion of irreducible rational normal quartic support is retained as an explicit earlier input.',
                        'No assertion about all nonnormal divisors or total Poisson component count.'],
                input_sha256=inputs,seconds=time.monotonic()-start)
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    r=verify(p.parse_args().output);print('PASS',r['id'],r['incidence_tree_counts'])
