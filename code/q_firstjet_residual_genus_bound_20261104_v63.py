#!/usr/bin/env python3
"""Proof-gated coarse normalization-genus bound; no CAS or points."""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_RESIDUAL_GENUS_BOUND_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic()
    primary=ROOT/'reports/Q_FIRSTJET_RESIDUAL_RECIPROCAL_QUOTIENT_20261104_v63'
    report=json.loads((primary/'report.json').read_text())
    cert=json.loads((primary/'certificate.json').read_text())
    assert report['all_passed']
    assert report['actual_plane_total_degree_r_m']==4
    assert cert['double_cover']=='h^2=r^2+4*f*u*v'
    table=[]
    for d in range(1,5):
        plane=(d-1)*(d-2)//2
        poles=3*d
        table.append({'plane_component_degree':d,'plane_genus_upper':plane,
                      'pole_divisor_degree_upper':poles,'ramification_degree_upper':2*poles,
                      'connected_quadratic_genus_upper':2*plane-1+poles})
    assert max(row['connected_quadratic_genus_upper'] for row in table)==17
    paths=[Path(__file__),primary/'report.json',primary/'certificate.json',
           ROOT/'research/notes/Q_FIRSTJET_RESIDUAL_RECIPROCAL_QUOTIENT_20261104_v63.md']
    data={'id':STEM,'all_passed':True,
          'proof_only_plane_genus_and_char0_Riemann_Hurwitz_not_CAS_genus':True,
          'scope':'each m-dominant geometric component after ordered eigenbase and existing relative-sign quotient',
          'coarse_normalization_genus_upper_bound':17,'degree_table':table,
          'exact_genus_smoothness_irreducibility_selected_degree8_or_Weyl_descent_claim':False,
          'no_new_beta_or_coefficient_point_CAS_GB_or_invariant':True,
          'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
          'seconds':time.monotonic()-start}
    out.mkdir(parents=True,exist_ok=False)
    (out/'report.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'genus_upper_bound':17,'seconds':data['seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
