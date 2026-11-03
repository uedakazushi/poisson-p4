#!/usr/bin/env python3
"""Identify the exact resonance curve with the known q51 component directly."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from s61_model import X
from s61_rational_family_adapter_20261103_v62 import normalized_field,verify_field_map

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();start=time.monotonic()
    curve=ROOT/'reports/Q_RESONANCE_INTEGRABLE_DIRECTION_20261103_v62.json'
    d=json.loads(curve.read_text());assert d['all_passed'];e=s.Symbol('e');P=s.zeros(5)
    entries={(1,4):s.Rational(2,9)*X[0]**2,(2,3):-s.Rational(2,9)*X[0]**2,
             (2,4):s.Rational(4,3)*X[0]*X[1],(3,4):s.Rational(2,3)*X[0]*X[2]+X[1]**2}
    for key,text in d['integrable_R'].items():
        ij=tuple(map(int,key.split(',')));entries[ij]=entries.get(ij,0)+e*s.sympify(text)
    entries[0,1]-=75*e**2*X[4]**2
    for (i,j),v in entries.items():P[i,j]=v;P[j,i]=-v
    model,K,values,normal,metadata=normalized_field(P,[e])
    actual=verify_field_map(model,K,values,normal)
    graph=ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json'
    D=json.loads(graph.read_text());assert D['all_passed'];u=s.symbols('u1:62')
    one=D['weight1_u_indices_zero_based'];two=D['weight2_u_indices_zero_based'];beta=[u[i]for i in one]
    for j,text in zip(two,D['quadratic_section_q']):
        expected=K.zero
        for mon,c in s.Poly(s.sympify(text),*beta,domain=s.QQ).terms():
            term=K.convert(c)
            for i,ex in zip(one,mon):
                if ex:term*=values[i]**ex
            expected+=term
        assert values[j]==expected,j
    files=[Path(__file__),curve,graph,ROOT/'reports/Q51_UNIVERSAL_SECTION_REVIEW_20261101.json',
           ROOT/'code/s61_rational_family_adapter_20261103_v62.py',ROOT/'code/s61_model.py',
           ROOT/'code/data/mps_s61_coordinate_catalogue_20261103_v62.json']
    report=dict(id='Q_RESONANCE_Q51_INCIDENCE_20261103_V62',all_passed=True,
                coefficient_field='QQ(e)',normalizer_metadata=metadata,
                full61_normalized_parameters=list(map(str,values)),actual_model_checks=actual,
                all45_q51_graph_identities=True,known_component='T39, q_(5,1)',
                entire_resonance_curve_in_component=True,resonant_Q_in_component=True,
                closure_argument='The generic curve point is PGL-equivalent to an actual point of the proved Q51 graph. The known component is closed and PGL-invariant; its preimage under the polynomial resonance curve is closed and contains the generic point, hence the whole A1, including e=0.',
                normal_actual_divisor_classification_used=False,new_component_claimed=False,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files},
                seconds=time.monotonic()-start)
    with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS direct q51 incidence, all45 graph equations and original134+150',report['seconds'])


if __name__=='__main__':main()
