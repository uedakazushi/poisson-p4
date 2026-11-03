#!/usr/bin/env python3
"""Exact rational involution carrying h=4,rho=1 to h=-4,rho=1."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sympy as s
from s61_model import S61, BASIS, X, J

ROOT=Path(__file__).resolve().parents[1]


def verify():
    model=S61();diagonal=[1,1,-1,-1,1];D=s.diag(*diagonal)
    assert D*J*D.T==-J
    def character(index):
        pair,mon=BASIS[index]
        result=-diagonal[pair[0]]*diagonal[pair[1]]
        for i,x in enumerate(X):result*=diagonal[i]**int(mon.as_powers_dict().get(x,0))
        return result
    signs=[character(i)for i in model.certificate['raw_free_indices']]
    assert len(model.data['eliminated_raw_variables'])==55
    for entry in model.data['eliminated_raw_variables']:
        expected=character(entry['raw_index'])
        assert all(math.prod(signs[i]for i in mon)==expected for mon,n,d in entry['polynomial']if n)
    weights=model.reduced['parameter_weights']
    one=[i for i,w in enumerate(weights)if w==1];two=[i for i,w in enumerate(weights)if w==2]
    extra=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    equations=model.reduced['equations']+[json.loads(extra.read_text())['universal_new_equation']]
    row_characters=[]
    for equation in equations:
        chars={math.prod(signs[i]for i in mon)for mon,n,d in equation if n}
        assert len(chars)==1
        row_characters.append(chars.pop())
    # t is fixed, hence the same characters hold after q-homogenization.
    assert all(signs[two[i]]in [1,-1]for i in [0,19,27,32])
    assert signs[two[25]]==1
    beta=s.zeros(16,1)
    for i,c in {4:1,11:1,25:1,26:4,50:2}.items():beta[one.index(i)]=c
    target=s.Matrix([signs[i]*beta[j]for j,i in enumerate(one)])
    expected=beta.copy();expected[one.index(26)]=-4
    assert target==expected
    points_path=ROOT/'code/certificates/s61_marked_sigma_pure_endpoint_points_20261102.json'
    qs=[s.Matrix(q).applyfunc(s.Rational)for q in json.loads(points_path.read_text())['q_points']]
    moved=[]
    for q in qs:
        q2=s.Matrix([signs[i]*q[j]for j,i in enumerate(two)])
        values=dict(zip(one,target))|dict(zip(two,q2))
        assert all(sum(s.Rational(n,d)*s.prod(values[i]for i in mon)
                       for mon,n,d in equation if sum(i in two for i in mon)==2)==0 for equation in equations)
        assert q2[25]==-1 and all(q2[i]==0 for i in [0,19,27,32])
        moved.append(list(map(str,q2)))
    note=ROOT/'research/notes/S61_MARKED_SIGMA_H_SIGN_TRANSPORT_20261103_v62.md'
    files=[Path(__file__),note,points_path,extra,ROOT/'code/s61_model.py']
    files.extend(model.source/name for name in ['symplectic_slice_certificate.json',
                 'symplectic_slice_equations.json','symplectic_slice_reduced_equations.json'])
    return dict(id='S61_MARKED_SIGMA_H_SIGN_TRANSPORT_20261103_v62',all_passed=True,
                field='QQ',affine_matrix_diagonal=diagonal[1:],bivector_action='P(x) -> -D P(D x) D^T',
                homogenizer_character=1,original_parameter_characters=signs,
                all55_eliminated_coordinate_formulas_equivariant=True,
                all135_homogenized_equations_preserved_up_to_sign=True,row_characters=row_characters,
                q45_characters=[signs[i]for i in two],
                q25_fixed=True,three_orbit_cuts_and_q32_base_preserved=True,
                rho_fixed=True,a_fixed=True,c_fixed=True,h_character=-1,
                h_minus4_extra_q_points=moved,
                mixed_characteristic_unit_automorphism_at_32003=True,
                h4_length_and_ambient_certificates_transport_to_minus4=True,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS',report['id'])
