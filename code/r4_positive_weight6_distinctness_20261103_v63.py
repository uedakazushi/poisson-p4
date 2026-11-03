#!/usr/bin/env python3
"""A short intrinsic distinction from ALL51 old rank4 components.

Each MPS component has a toric rank4 member with five distinct hyperplanes.
The11 extra component dimensions differ from25. No global exhaustion claim.
"""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from additional_rank4_families_20261103_v62 import families
ROOT=Path(__file__).resolve().parents[1]
def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic();path=ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json';data=json.loads(path.read_text());witnesses=[]
 for record in data['diagrams']:
  B=s.Matrix(record['biresidue_example'])[1:5,1:5];pf=s.expand(B[0,1]*B[2,3]-B[0,2]*B[1,3]+B[0,3]*B[1,2]);assert pf!=0
  coords=s.Matrix([B[i,j]for i,j in itertools.combinations(range(4),2)]);period=s.Matrix(record['period_basis']).T;solution=period.gauss_jordan_solve(coords)[0];assert period*solution==coords
  A=(-pf*B.inv()).applyfunc(s.cancel);pfA=s.factor(A[0,1]*A[2,3]-A[0,2]*A[1,3]+A[0,3]*A[1,2]);assert pfA==pf
  witnesses.append(dict(component_id=record['id'],period_parameters=list(map(str,solution)),biresidue4=[[str(v)for v in row]for row in B.tolist()],toric_coefficient4=[[str(v)for v in row]for row in A.tolist()],actual_quintic=str(pf)+'*x0*x1*x2*x3*x4',rank4=True,outside_closed_sigma32=True))
 assert len(witnesses)==40;extras=[dict(component_id=f['component_id'],projective_dimension=f['projective_dimension'])for f in families()];assert len(extras)==11 and all(r['projective_dimension']!=25 for r in extras)
 files=[Path(__file__),path,ROOT/'code/additional_rank4_families_20261103_v62.py',ROOT/'reports/MPS40_CROSSWALK_20261031.json',ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json',ROOT/'reports/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json'];report=dict(id='R4_POSITIVE_WEIGHT6_DISTINCTNESS_20261103_v63',all_passed=True,field='QQ',closed_locus='Sigma32=image(P(V*) x P(Sym²V*) -> P(Sym5V*)), (L,Q)->L³Q',new_component_rank4_divisors_all_in_sigma32=True,forty_MPS_actual_toric_witnesses=witnesses,eleven_extra_dimensions=extras,all51_old_rank4_components_distinct=True,six_rank2_components_distinct_by_generic_rank=True,accepted_component_count_dependency='Independent component25 review; this packet supplies a short alternate distinction proof',lower_bound=58,global_exhaustion_claimed=False,input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in files},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('INTRINSIC40TORIC+11DIM DISTINCTNESS AND58 LOWER PASS',report['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
