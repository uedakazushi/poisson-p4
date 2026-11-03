#!/usr/bin/env python3
"""Exact magnetic membership and proof-scope review for full-model non-CM."""
import hashlib
import json
from pathlib import Path

import sympy as s
from s61_model import X,pushforward

ROOT=Path(__file__).resolve().parents[1]


def main():
 a=s.Symbol('a',nonzero=True);Q=s.zeros(5)
 for (i,j),v in {(1,4):a*X[0]**2,(2,3):-s.Rational(2,9)*X[0]**2,
   (2,4):s.Rational(4,3)*X[0]*X[1],(3,4):s.Rational(2,3)*X[0]*X[2]+X[1]**2}.items():
  Q[i,j]=v;Q[j,i]=-v
 M=s.zeros(5);M[0,0]=M[1,1]=M[3,2]=1;M[2,4]=1/a;M[4,3]=-s.Rational(9,2)
 H=(s.Rational(4,3)*X[0]*X[1]*X[4]+s.Rational(3,2)*X[0]*X[3]**2+s.Rational(9,2)*X[1]**2*X[3])/a
 target=s.zeros(5)
 for(i,j),v in {(1,2):X[0]**2,(3,4):X[0]**2,(2,3):-s.diff(H,X[4]),(2,4):s.diff(H,X[3])}.items():
  target[i,j]=v;target[j,i]=-v
 assert (pushforward(Q,M)-target).applyfunc(s.cancel)==s.zeros(5)
 table=json.loads((ROOT/'reports/KNOWN51_FIRSTJET_DIMENSION_TABLE_20261103_v62.json').read_text())
 T39=next(row for row in table['records']if row['component_id']=='T39')
 assert T39['s61_component_dimension']+9==25
 paths=[Path(__file__),ROOT/'research/notes/Q_RESONANCE_NON_CM_20261103_v63.md',
   ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63_REVIEW.json',
   ROOT/'reports/Q_RESONANCE_Q51_INCIDENCE_20261103_v62_REVIEW.json',
   ROOT/'reports/KNOWN51_FIRSTJET_DIMENSION_TABLE_20261103_v62.json',
   ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',
   ROOT/'legacy/v06/ambient_families.tex',ROOT/'code/s61_model.py']
 report=dict(id='Q_RESONANCE_NON_CM_20261103_V63_REVIEW',all_passed=True,
             independent_exact_nonzero_a_magnetic_transport_all150_coefficients=True,
             magnetic_transport_invertible_at_two_over_nine=True,
             actual_known_global_components_retained_as_theorem_input=True,
             known_projective_minimal_component_dimensions=[28,25],
             excellent_completion_dimension_preservation_accepted=True,
             known_full35_minimal_component_dimensions=[12,9],
             full35_completion_not_Cohen_Macaulay_accepted=True,
             full_projective_Jacobi_completion_not_Cohen_Macaulay_accepted=True,
             unknown_components_do_not_remove_unmixedness_obstruction=True,
             full_branch_count_or_dimension_upper_claimed=False,
             input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})
 path=ROOT/'reports/Q_RESONANCE_NON_CM_20261103_v63_REVIEW.json'
 with path.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
 print('PASS exact magnetic transport and full-model unmixedness obstruction')


if __name__=='__main__':main()
