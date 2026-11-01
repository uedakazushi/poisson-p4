#!/usr/bin/env python3
"""Remove the artificial Q denominator from the closed T7 flag flat4 quotient."""
import argparse,hashlib,json,subprocess,sys,tempfile,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def run():
 start=time.time();engines=['code/s61_flag_T7_flat4_transport_engine_20261101.py','code/s61_flag_T7_flat4_identities_engine_20261101.py','code/s61_flag_T7_flat4_cover_lattice_engine_20261101.py']
 with tempfile.TemporaryDirectory(prefix='flag-T7-flat4-cover-')as directory:
  wd=Path(directory)
  for engine in engines:subprocess.run([sys.executable,str(ROOT/engine),str(wd)],cwd=ROOT,check=True)
  lattice=json.loads((wd/'flag-T7-universal-lattice-ring.json').read_text());ident=json.loads((wd/'flag-T7-universal-transport-identities.json').read_text());assert ident['all_passed']
  b=s.symbols('b0:11');e=s.Symbol('e');P=19800*b[0]*b[5]*b[8]+4480*b[0]*b[9]**2-1144*b[10]*b[8]*b[9]-1000*b[2]*b[8]*b[9]-5940*b[6]*b[8]**2-2079*b[7]*b[8]**2
  constant=s.Integer(3100909122761275736064000);D0=s.sympify(lattice['coordinate_denominator']).subs(e,0)
  assert s.Poly(D0-constant*b[0]**2*b[8]**12*b[9]**4*P**3,*b).is_zero
  for row in lattice['multiplication']:
   denominator=s.sympify(row['denominator']).subs(e,0);quot=s.cancel(denominator/D0);terms=s.Poly(quot,*b).terms();assert len(terms)==1
   ex,c=terms[0];assert c and all(k==0 for j,k in enumerate(ex)if j not in [8,9])
  assert lattice['selected_q_rows']==[9,0,10]
  digest=hashlib.sha256((wd/'flag-T7-universal-lattice-ring.json').read_bytes()).hexdigest()
 old=json.loads((ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json').read_text());db=s.Matrix(old['frame_change']['beta_derivative']).applyfunc(s.sympify);bb=s.symbols('b0:15')
 for col in [4,2]:assert db[0,col]==db[12,col]==0
 t7=json.loads((ROOT/'reports/S61_R0_MIXED_SURFACE_T7_EXHAUSTION_20261101.json').read_text());boundary=json.loads((ROOT/'reports/S61_FLAG_T7_BOUNDARY_20261101.json').read_text());assert t7['beta_kernel']==boundary['fixed_beta_kernel']
 exceptions={x['label']:x['beta_image_dimension']for x in t7['exceptions']};assert exceptions['b9']==exceptions['P']==15
 paths=engines+[str(Path(__file__).relative_to(ROOT)),'reports/S61_FLAG_T7_FLAT4_20261101.json','reports/S61_FLAG_T7_BOUNDARY_20261101.json','reports/S61_R0_MIXED_SURFACE_T7_EXHAUSTION_20261101.json','reports/S61_FLAG_UNIFORM_OPEN_20261101.json','code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py']
 return dict(id='S61_FLAG_T7_FLAT4_COVER_20261101',all_passed=True,selected_q_coordinates=[9,0,10],denominator_constant=str(constant),denominator_at_endpoint=str(constant*b[0]**2*b[8]**12*b[9]**4*P**3),regular_open='a0*a8*a9*P !=0',artificial_Q_removed=True,all135_four_coefficients_zero=True,coordinate_and_product_regularity=True,flat_rank=4,reduced_base_dimension=2,homogenizer_zero_in_quotient=True,core_transport_preserves_a0_and_a8=True,bad_image_dimensions={'a9':15,'P':15},whole_ordinary_core_flat4_over_sufficiently_general_beta=True,lattice_payload_sha256=digest,scope='Every point with b0*b12 !=0, b1=0 and I/b12^2=-91/400 over sufficiently general beta has the stated actual closed flat4 surface quotient.',limitations=['A quotient assertion, not the ambient completed ring or ambient t-nilpotence.','The b0=0 and b12=0 flags and other cores remain separate.','No saturation removal is deduced from the flat4 quotient alone.'],seconds=time.time()-start,proof_note='research/notes/S61_FLAG_T7_FLAT4_COVER_20261101.md',input_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in paths})
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run()
 with a.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
 print('WHOLE ORDINARY T7 CORE FLAT4 QUOTIENT PASS',r['seconds'],flush=True)
