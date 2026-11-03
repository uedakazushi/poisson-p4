#!/usr/bin/env python3
"""Preserve frozen firstjet evidence; clarify raw vs adjusted coordinate ranks."""
import argparse,gzip,hashlib,json
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def matrix(rows):return s.Matrix([[s.Rational(v)for v in row]for row in rows])
def main(output):
    folder=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63'
    path=folder/'certificate.json.gz';data=json.loads(gzip.decompress(path.read_bytes()));results={}
    for name,record in data['CA_CB_full_firstjet_data'].items():
        A=matrix(record['actual']);E=matrix(record['tangent']);normal=matrix(record['normal']);L=matrix(record['L']);dnu=matrix(record['discriminant'])
        change=s.eye(29);change[28,:28]=L
        assert A*change==E.row_join(normal)
        ker=s.Matrix.hstack(*dnu.nullspace());assert ker.shape==(29,28)
        assert A.rank()==7 and A[:,:28].rank()==7 and E.rank()==6 and(A*ker).rank()==6
        results[name]=dict(actual_all29_rank=7,raw_Y28_columns_at_fixed_rawrow61_rank=7,
            adjusted_Y28_at_fixed_eta_equals_row61_minus_LY_rank=6,
            intrinsic_restriction_to_kernel_of_discriminant_differential_rank=6)
    N=data['N_full_firstjet_data'];A=matrix(N['actual']);E=matrix(N['tangent']);dnu=matrix(N['discriminant']);ker=s.Matrix.hstack(*dnu.nullspace())
    assert A.rank()==7 and E.rank()==6 and(A*ker).rank()==6
    note=ROOT/'research/notes/AMBIENT_B18_ODD_FIRSTJET_COORDINATE_SCOPE_20261103_v63.md'
    result=dict(id='AMBIENT_B18_ODD_FIRSTJET_COORDINATE_SCOPE_20261103_v63',all_passed=True,CA_CB=results,
        N=dict(full_actual28_rank=7,adjusted27_rank=6,discriminant_tangent_restriction_rank=6),
        frozen_original_evidence_unchanged=True,original_rank6_label_means_adjusted_tangent_not_fixed_raw29th_coordinate=True,
        N_projection_pinch_and_discriminant_independence_unaffected=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),note,path,folder/'report.json']})
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print('ODD FIRSTJET COORDINATES: RAW7 / ADJUSTED6 / INTRINSIC6 PASS')
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
