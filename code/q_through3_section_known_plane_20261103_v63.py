#!/usr/bin/env python3
"""Exact12-linear-section intersection with the known13-plane cubic divisor."""
import argparse
from functools import reduce
import hashlib
import json
from math import gcd
from pathlib import Path
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--section',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    source=args.section/'section.json';section=json.loads(source.read_text())
    plane_path=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/coordinate-plane.json'
    plane=json.loads(plane_path.read_text());assert plane['all_passed']
    cubic_path=ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63/plane-certificate.json'
    cubic=json.loads(cubic_path.read_text());assert cubic['all_passed']
    names=section['original_variables'];active=plane['active_variables'];assert len(active)==13
    indices=[names.index(name) for name in active];T=section['integer_matrix_T']
    M=s.Matrix(12,13,lambda i,j:(int(indices[j]==i) if indices[j]<12 else -T[i][indices[j]-12]))
    # Signed maximal minors give an explicit integer kernel without numerical
    # or finite-field nullspace reconstruction.
    minors=[]
    for j in range(13):
        block=M[:,[k for k in range(13) if k!=j]]
        minors.append((-1)**j*int(DomainMatrix.from_Matrix(block).det()))
    assert any(minors) and M*s.Matrix(minors)==s.zeros(12,1)
    content=reduce(gcd,minors);v=[value//content for value in minors]
    assert M*s.Matrix(v)==s.zeros(12,1)
    unit=next(i for i,value in enumerate(minors) if value%P)
    values=dict(zip(active,v))
    factor=12*values['z14']+values['z17']+6*values['z3']
    f=values['z8']**2*factor;assert f and f%P
    report=dict(id='Q_THROUGH3_SECTION_KNOWN_PLANE_20261103_V63',all_passed=True,
                coefficient_ring='Z, localized at32003 when a unit is used',
                active_plane_variables=active,restricted_linear_matrix=[list(map(int,row)) for row in M.tolist()],
                signed_maximal_minors=minors,primitive_integer_kernel_vector=v,
                chosen_unit_minor_column=unit,chosen_signed_minor=minors[unit],
                restricted_linear_rank_QQ_and_GF32003=12,
                z8_on_kernel=values['z8'],ell_on_kernel=factor,
                cubic_f_on_kernel=f,cubic_f_mod32003=f%P,
                all_required_linear_and_cubic_coefficients_are32003_units=True,
                intersection_ring_over_V='V[t]/(t^3)',
                known13plane_intersection_is_Artinian_length3=True,
                known13plane_projective_section_empty_over_QQ_and_GF32003=True,
                full23variable_section_Artinian_not_inferred=True,
                input_sha256={str(path.resolve().relative_to(ROOT)):sha(path) for path in
                              [Path(__file__),source,plane_path,cubic_path]},seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS known13plane section V[t]/t^3; fmodp',f%P,'seconds',report['seconds'])


if __name__=='__main__':main()
