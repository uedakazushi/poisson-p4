#!/usr/bin/env python3
"""A closed, free-rank4 T7 surface quotient through a principal flag-boundary open."""
import argparse,hashlib,json,subprocess,sys,tempfile,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();start=time.time()
 engines=[f'code/s61_flag_T7_flat4_{stage}_engine_20261101.py'for stage in ['transport','identities','lattice','geometry']]
 with tempfile.TemporaryDirectory(prefix='flag-T7-flat4-')as directory:
  wd=Path(directory)
  for engine in engines:subprocess.run([sys.executable,str(ROOT/engine),str(wd)],cwd=ROOT,check=True)
  identities=json.loads((wd/'flag-T7-universal-transport-identities.json').read_text());geometry=json.loads((wd/'flag-T7-universal-lattice-geometry.json').read_text());lattice=json.loads((wd/'flag-T7-universal-lattice-ring.json').read_text())
  assert identities['all_passed']and geometry['all_passed'];digest=hashlib.sha256((wd/'flag-T7-universal-lattice-ring.json').read_bytes()).hexdigest()
 paths=engines+[str(Path(__file__).relative_to(ROOT)),'reports/S61_FLAG_T7_BOUNDARY_20261101.json','reports/S61_R0_MIXED_SURFACE_T7_FORMAL_20261101.json','code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py','reports/S61_UNIVERSAL_WEIGHT6_20261031.json','research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
 result=dict(id='S61_FLAG_T7_FLAT4_20261101',all_passed=True,identities=identities,geometry=geometry,selected_q_coordinates=lattice['selected_q_rows'],lattice_payload_sha256=digest,lattice_payload_rebuilt_not_used_as_input=True,flat_rank=4,reduced_base_parameters=['epsilon','u'],endpoint_nil_cotangent_dimension=2,homogenizer_zero_in_quotient=True,closed_formal_quotient=True,all45_coordinate_and_all6_product_coefficients_regular=True,scope='The stated principal open in H11, with formal epsilon and inverse-beta-compensated shear parameter u.',limitations=['This is a closed flat4 surface quotient of the actual completed equations, not an identification of the ambient local ring.','No ambient t-nilpotence or saturation removal is asserted.','The complement of the explicit principal open, and the other flag strata, remain separate.'],proof_note='research/notes/S61_FLAG_T7_FLAT4_20261101.md',seconds=time.time()-start,input_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in paths})
 with a.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('CLOSED T7 FLAG FLAT4 CERTIFICATE PASS',result['seconds'],flush=True)
if __name__=='__main__':main()
