#!/usr/bin/env python3
"""Certify an enlarged known-five section and its finite-normal-module bounds.

Uses the original134 homogeneous equations, an integral extension of the frozen
P25 normal plane, and exact F32003 Groebner arithmetic. Does not claim full
S61 affine exhaustion. Reuses the existing Cloud M2, with its supported
`groebnerBasis(...,Strategy=>\"F4\")` entrypoint.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np
import sympy as sp

from s61_affine_finite_normal_module_20261102 import slice_model
from s61_affine_normal_projection_20261102 import ROOT, P, load_fiber, residue, pivot_minor


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(normal_count):
    section, _, axes, inputs = slice_model()
    known_path = ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json'
    known = json.loads(known_path.read_text())
    equations, ids, center, _, _ = load_fiber()
    n = 4+normal_count
    plane_variables = sp.symbols('v0:3')
    coordinates = np.zeros((46,n),dtype=np.int64)
    coordinates[0,0] = 1
    coordinates[1:,1:4] = [[residue(sp.sympify(e).coeff(x)) for x in plane_variables]
                           for e in known['affine_plane_substitution']]
    normals = np.column_stack([np.array(section['fixed_normal_plane_integer_basis'],dtype=np.int64),
                               np.eye(42,dtype=np.int64)[:,axes[:normal_count-22]]])
    assert len(pivot_minor(normals)[0]) == normal_count
    coordinates[4:,4:] = normals%P
    quadrics = np.zeros((134,46,46),dtype=np.int64)
    for i,equation in enumerate(equations):
        for mon,c in equation.items():
            if len(mon)==1:
                quadrics[i,0,mon[0]+1] = residue(c)
            else:
                assert len(mon)==2
                quadrics[i,mon[0]+1,mon[1]+1] = residue(c)
    restricted = coordinates.T@quadrics@coordinates%P
    monomials = list(itertools.combinations_with_replacement(range(n),2))
    rows = np.column_stack([restricted[:,i,j]+(restricted[:,j,i] if i!=j else 0)
                            for i,j in monomials])%P
    t = sp.Symbol('t')
    conjugate = int(known['toric_quadratic_conjugate_witness']['t_mod_prime'])
    points = []
    for root in [0,conjugate,-conjugate]:
        points.append([1]+[residue(sp.sympify(known['toric_S61_coordinates'][ids[i]]).subs(t,root)-center[i])
                           for i in range(3)])
    for name in ['q51_S61_coordinates','q52_S61_coordinates']:
        points.append([1]+[residue(sp.sympify(known[name][ids[i]])-center[i]) for i in range(3)])
    assert len(set(map(tuple,points)))==5
    evaluation=np.array([[point[i]*point[j]%P for i,j in itertools.combinations_with_replacement(range(4),2)]
                         for point in points],dtype=np.int64)
    assert len(pivot_minor(evaluation)[0])==5
    for point in points:
        full = point+[0]*normal_count
        assert not np.any(rows@np.array([full[i]*full[j]%P for i,j in monomials])%P)
    text = 'R=(ZZ/32003)['+','.join(f'x{i}' for i in range(n))+'];\n'
    text += 'I=ideal('+','.join('+'.join(str(int(c))+f'*x{i}*x{j}'
                                       for (i,j),c in zip(monomials,row) if c) or '0' for row in rows)+');\n'
    return text,normals,points,inputs+[known_path],n


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--normal-count',type=int,default=24,choices=range(23,39))
    parser.add_argument('--timeout',type=int,default=120)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    started = time.monotonic()
    text,normals,points,inputs,n = build(args.normal_count)
    script,log = args.output/'section.m2',args.output/'section.log'
    leading = args.output/'leading-exponents.txt'
    text += 'print "SECTION START";stdio << flush;\n'
    text += 'G=groebnerBasis(I,Strategy=>"F4");H=monomialIdeal leadTerm G;\n'
    text += 'print("HF="|toString apply(0..7,d->hilbertFunction(d,R/H)));stdio << flush;\n'
    text += 'assert(dim(R/H)==1 and degree(R/H)==5);\n'
    text += 'print("DIM="|toString dim(R/H));print("DEG="|toString degree(R/H));\n'
    text += 'out=openOut '+json.dumps(str(leading.resolve()))+';\n'
    text += 'scan(flatten entries leadTerm G,f->out << toString first exponents f << endl);close out;\n'
    text += 'print "ENLARGED SECTION PASS";stdio << flush;exit 0;\n'
    script.write_text(text)
    executable = shutil.which('M2')
    assert executable
    with log.open('w') as stream:
        result = subprocess.run([executable,'--script',str(script.resolve())],stdout=stream,stderr=subprocess.STDOUT,
                                check=True,timeout=args.timeout)
    content = log.read_text()
    assert 'ENLARGED SECTION PASS' in content
    match = re.search(r'HF=\(([^\n]+)\)',content)
    assert match
    hf = list(map(int,match.group(1).split(',')))
    equality = min(d for d in range(2,len(hf)) if hf[d]==5)
    assert all(h==5 for h in hf[equality:])
    retained = 42-args.normal_count
    known_hf = [1,4]+[5]*6
    bounds = [retained*known_hf[d-1]+hf[d]-known_hf[d] for d in range(1,equality+1)]
    (args.output/'normal-basis.json').write_text(json.dumps(normals.tolist(),indent=2)+'\n')
    inputs += [Path(__file__),ROOT/'code/s61_affine_finite_normal_module_20261102.py',
               ROOT/'code/s61_affine_normal_projection_20261102.py',ROOT/'code/s61_affine_fiber_macaulay.py']
    report = dict(id='S61_AFFINE_ENLARGED_SECTION_20261103_v62',all_passed=True,prime=P,
                  projective_section_dimension=n-1,normal_direction_count=args.normal_count,
                  retained_normal_count=retained,coefficient_equations=134,
                  hilbert_function_mod_32003=hf,first_full_five_ideal_degree=equality,
                  known_five_mod_32003=points,known_five_degree_two_evaluation_rank=5,section_reduced_length=5,
                  generic_QQ_section_reduced_length=5,
                  finite_normal_module_generator_bounds_by_degree=bounds,
                  finite_normal_module_generator_bound=sum(bounds),
                  finite_normal_module_generation_degree=equality,
                  limitations=['Finite-field degree equality lifts to QQ/general beta as explained in proof note.',
                               'Full S61 affine exhaustion and zero residual module are not proved.',
                               'The finite-module bounds are upper bounds; the action matrix is not constructed here.'],
                  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(inputs)},
                  evidence_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()},
                  seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',n-1,retained,hf,bounds,sum(bounds),report['seconds'])


if __name__=='__main__':
    main()
