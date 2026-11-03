#!/usr/bin/env python3
"""Freeze the exact missing premise in the proposed integral Nakayama shortcut."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert not args.output.exists()
    stems=['S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102','S61_FIXED_PAIR_SINGLE_COLON_20261102',
           'S61_AFFINE_FIXED_PAIR_LOCALIZATION_20261102','S61_AFFINE_HOMOGENIZER_FINITE_ACTION_20261102']
    checked={};reports={}
    for stem in stems:
        path=ROOT/'reports'/(stem+'.json');report=json.loads(path.read_text());assert report['all_passed']
        checked[str(path.relative_to(ROOT))]=sha(path);reports[stem]=report
        for name,value in report['input_sha256'].items():
            source=ROOT/name;assert sha(source)==value,name;checked[name]=value
    finite=reports[stems[2]];action=reports[stems[3]]
    assert finite['retained_base_variables']==20 and not finite['all_twenty_parameters_empty_proved']
    assert action['homogenizer_candidate_size']==191 and not action['all_directions_excluded_claimed']
    assert reports[stems[1]]['colon_exponent']==2 and not reports[stems[1]]['fixed_c_affine_algebra_length_determined']
    p=32003;y,t=sp.symbols('Y t');affine=p*y-1;homogeneous=t-p*y
    assert affine.subs(y,sp.Rational(1,p))==0
    assert sp.Poly(affine,y,modulus=p).as_expr()==-1
    assert homogeneous.subs({t:p,y:1})==0
    assert sp.Poly(homogeneous,t,y,modulus=p).as_expr()==t
    # In the special fiber t^2=t*t belongs to(t), while1 does not. The
    # relative quotient V[t,Y]/(t-pY)=V[Y] is a domain, so relative colon=t-pY.
    assert sp.expand(t*t-t**2)==0
    checked[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    result=dict(id='S61_AFFINE_INTEGRAL_FINITE_PREMISE_AUDIT_20261103_v62',all_passed=True,
                scope='source/report hash audit and explicit obstruction; preserves the valid generic geometric projection theorem',
                checked_file_count=len(checked),old_fixed_pair_model_finite_over='R_beta,c[Y1,...,Y20]',
                finite_over_parameter_local_ring_at_frozen_point_certified=False,
                generic_center_free_open_contains_frozen_point_certified=False,
                relative_integral_saturation_or_complete_infinity_nilpotence_cover_certified=False,
                counterexample=dict(prime=p,homogeneous_ideal='(t-pY)',monic_unknown_pivot=1,
                                    finite_over_V_Y_generators=1,affine_quotient='V[Y]/(pY-1)',
                                    special_affine_quotient_zero=True,generic_affine_quotient_nonzero=True,
                                    generic_point_Y=str(sp.Rational(1,p)),relative_colon='(t-pY):t^2=(t-pY)',
                                    specialized_colon='(t):t^2=(1)',projective_model_finite_over_V=True,
                                    projective_t_open_not_finite_over_V=True),
                sufficient_extra_premise='Actual affine quotient finite over the p-integral parameter local ring; or a valid relative-saturation/center-free certificate implying that finiteness.',
                generic_rank_or_unit_claimed=False,input_sha256=checked)
    args.output.write_text(json.dumps(result,indent=2)+'\n');print('PASS',len(checked),'checked receipts; missing finite premise recorded')


if __name__=='__main__':main()
