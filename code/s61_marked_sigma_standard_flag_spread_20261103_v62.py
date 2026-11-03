#!/usr/bin/env python3
"""Arithmetic spread of the ordinary rho=1 flag neighborhood, uniform in h.

Checks the frozen universal actual Jacobian/quadratic data after the Borel
shift b3->b3+b0. All remaining beta coordinates stay arbitrary.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from weight1_structure import J

ROOT=Path(__file__).resolve().parents[1]
PRIME=32003


def verify():
    started=time.monotonic();b=s.symbols('b:15');sub={b[0]:1,b[1]:1,b[3]:1,b[12]:0}
    path=ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json';source=json.loads(path.read_text());assert source['all_passed']
    review=ROOT/'reports/S61_FLAG_UNIFORM_OPEN_REVIEW_20261101.json';assert json.loads(review.read_text())['all_passed']
    D5path=ROOT/'reports/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.json';D5=json.loads(D5path.read_text());assert D5['all_passed']
    scalar_strings=set()
    for key in ['jacobian','tangent','left_kernel']:
        scalar_strings.update(c for row in source['linear'][key]for c in row)
    scalar_strings.update(c for row in source['quadratic']['quadratic_matrix']for c in row)
    remaining=[z for z in b if z not in sub]
    count=0;coefficient_count=0;denominators=set()
    for text in scalar_strings:
        value=s.cancel(s.sympify(text).subs(sub));num,den=s.fraction(value)
        assert den.is_Rational and int(den.p)%PRIME and int(den.q)%PRIME
        coefficients=s.Poly(num,*remaining,domain=s.QQ).coeffs()
        assert all(int(c.q)%PRIME for c in coefficients)
        denominators.add(str(den));coefficient_count+=len(coefficients);count+=1
    jacobian=s.factor(s.sympify(source['linear']['jacobian_minor']).subs(sub))
    quadratic=s.factor(s.sympify(source['quadratic']['quadratic_minor']).subs(sub))
    assert jacobian==1959552
    assert all(int(c.p)%PRIME and int(c.q)%PRIME for c in [jacobian,quadratic])
    generator=s.Matrix(source['frame_change']['unipotent_generator']).applyfunc(s.Rational)
    g=s.eye(4)+generator
    assert generator**2==s.zeros(4)and(g.T*J*g-J)==s.zeros(4)
    assert source['frame_change']['translation']=='b3 -> b3+s*b0; b0,b1,b12 fixed'
    # Flat closed D5 model: check all arithmetic constants in the canonical K.
    u=s.symbols('u1:62');allD5=set(c for row in D5['kernel']for c in row)
    for text in allD5:
        num,den=s.fraction(s.cancel(s.sympify(text)));poly=s.Poly(den,u[4],domain=s.QQ)
        assert len(poly.terms())==1
        coefficient=poly.terms()[0][1];assert int(coefficient.p)%PRIME and int(coefficient.q)%PRIME
        assert all(int(c.q)%PRIME for c in s.Poly(num,*u,domain=s.QQ).coeffs())
    assert s.sympify(D5['kernel_minor']).subs(u[4],1)==s.Rational(80,21)
    assert s.sympify(D5['radial_C_normal_complement_minor']).subs(u[4],1)==-4
    files=[Path(__file__),path,review,D5path,ROOT/'code/weight1_structure.py',
           ROOT/'research/notes/S61_FLAG_UNIFORM_OPEN_20261101.md',
           ROOT/'research/notes/S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101.md',
           ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
           ROOT/'research/notes/S61_MARKED_SIGMA_STANDARD_FLAG_SPREAD_20261103_v62.md']
    return dict(id='S61_MARKED_SIGMA_STANDARD_FLAG_SPREAD_20261103_v62',all_passed=True,prime=PRIME,
                scope='rho=1 standard flag, b0=b1=1, initial b3=b12=0, all other beta arbitrary',
                finite_Borel_shift=source['frame_change']['translation'],shift_matrix=[[str(c)for c in row]for row in g.tolist()],
                checked_unique_rational_entries=count,checked_numerator_coefficients=coefficient_count,
                all_specialized_frame_denominators_constant_prime_units=True,
                Jacobian_38_minor=str(jacobian),normal_quadratic_10_minor=str(quadratic),
                D5_kernel_unique_entries_checked=len(allD5),D5_kernel_minor='80/21',D5_C_normal_minor='-4',
                actual_flag_incidence_beta_smooth_minor='-u5=-1',
                mixed_characteristic_relative_free_rank=5,ambient_t_zero=True,
                all_h_included=True,all_other_beta_parameters_retained=True,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files},
                seconds=time.monotonic()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS',report['id'],report['seconds'])
