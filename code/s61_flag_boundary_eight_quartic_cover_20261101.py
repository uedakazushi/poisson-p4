#!/usr/bin/env python3
"""Rebuild the genuine fourth initial core minor cover on the eight-variable flag chart."""
import argparse,hashlib,json,subprocess,sys,tempfile,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
 p.add_argument('--flint-include',type=Path,default=ROOT/'work/flint-headers/flint/3.5.0/include');p.add_argument('--mpfr-include',type=Path,default=ROOT/'work/flint-headers/mpfr/4.2.2/include');p.add_argument('--flint-library',type=Path,default=ROOT/'work/m2-runtime/lib/libflint.24.0.dylib')
 args=p.parse_args();start=time.time();early=[f'code/s61_flag_boundary_eight_{stage}_engine_20261101.py'for stage in ['frame','inverse','quadratic']];late=[f'code/s61_flag_boundary_eight_{stage}_engine_20261101.py'for stage in ['raw4','project4','initial4','cover4']];native=['code/parametric_implicit_jets_20261101.cpp','code/univariate_polynomial_minor_gcd_20261101.cpp','code/univariate_polynomial_gcd_factor_20261101.cpp']
 with tempfile.TemporaryDirectory(prefix='flag-eight-quartic-')as directory:
  wd=Path(directory);binaries=[wd/'raw-jets',wd/'polynomial-minors',wd/'polynomial-factors']
  for source,binary in zip(native,binaries):subprocess.run(['c++','-O3','-std=c++17','-I'+str(args.flint_include),'-I'+str(args.mpfr_include),'-I/opt/homebrew/include',str(ROOT/source),str(args.flint_library),'-L/opt/homebrew/lib','-lgmp','-o',str(binary)],check=True)
  for engine in early:subprocess.run([sys.executable,str(ROOT/engine),str(wd)],cwd=ROOT,check=True)
  for i,engine in enumerate(late):
   command=[sys.executable,str(ROOT/engine),str(wd)]
   if i==0:command.append(str(binaries[0]))
   if i==3:command.extend(map(str,binaries[1:]))
   subprocess.run(command,cwd=ROOT,check=True)
  result=json.loads((wd/'flag-eight-quartic-cover-result.json').read_text())
 paths=early+late+native+[str(Path(__file__).relative_to(ROOT)),'code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','reports/S61_FLAG_UNIFORM_OPEN_20261101.json','reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json','reports/S61_FLAG_BOUNDARY_EIGHT_QUARTIC_20261101.json']
 result.update(seconds=time.time()-start,raw135_rebuilt=True,consistent_raw_four_jets=True,quadratic_native_sympy_crosscheck=True,scope='Eight-variable ordinary flag chart over its function field and the explicitly certified core principal open.',limitations=['No full initial ideal equality or ambient t-nilpotence.','The zero sets of every recorded core factor must be covered separately; divided row contents are included.','The b8=0 boundary and four linear core exceptions are not covered by this certificate.'],proof_note='research/notes/S61_FLAG_BOUNDARY_EIGHT_QUARTIC_COVER_20261101.md',input_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in paths})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('EIGHT VARIABLE QUARTIC CORE COVER PASS',result['seconds'],flush=True)
if __name__=='__main__':main()
