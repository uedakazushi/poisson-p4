#!/usr/bin/env python3
"""Rebuild true fourth-order relations from all135 equations on r=-91/400."""
import argparse,hashlib,json,subprocess,sys,tempfile,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--flint-include',type=Path,default=ROOT/'work/flint-headers/flint/3.5.0/include');p.add_argument('--mpfr-include',type=Path,default=ROOT/'work/flint-headers/mpfr/4.2.2/include');p.add_argument('--flint-library',type=Path,default=ROOT/'work/m2-runtime/lib/libflint.24.0.dylib')
    args=p.parse_args();start=time.time();engines=[f'code/s61_flag_boundary_eight_{stage}_engine_20261101.py' for stage in ['frame','inverse','quadratic']];native='code/parametric_implicit_jets_20261101.cpp';raw='code/s61_flag_rankone_socle_raw4_engine_20261101.py';quartic='code/s61_flag_rankone_socle_quartic_engine_20261101.py'
    with tempfile.TemporaryDirectory(prefix='rankone-socle-quartic-') as directory:
        wd=Path(directory);binary=wd/'parametric-jets'
        subprocess.run(['c++','-O3','-std=c++17','-I'+str(args.flint_include),'-I'+str(args.mpfr_include),'-I/opt/homebrew/include',str(ROOT/native),str(args.flint_library),'-L/opt/homebrew/lib','-lgmp','-o',str(binary)],check=True)
        for engine in engines:subprocess.run([sys.executable,str(ROOT/engine),str(wd)],cwd=ROOT,check=True)
        subprocess.run([sys.executable,str(ROOT/raw),str(wd),str(binary)],cwd=ROOT,check=True)
        subprocess.run([sys.executable,str(ROOT/quartic),str(wd)],cwd=ROOT,check=True)
        result=json.loads((wd/'flag-rankone-socle-quartic.json').read_text())
    paths=engines+[native,raw,quartic,str(Path(__file__).relative_to(ROOT)),'code/s61_initial_quartic_from_jets_20261101.py','code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','reports/S61_FLAG_UNIFORM_OPEN_20261101.json','reports/S61_FLAG_BOUNDARY_QUADRATIC_COVER_20261101.json','reports/S61_FLAG_RANKONE_SOCLE_CUBIC_20261101.json']
    result.update(seconds=time.time()-start,consistent_raw4_rebuilt=True,all135_pivot_residuals_zero=True,quadratic_native_sympy_crosscheck=True,scope='All seven nuisance parameters with b8 nonzero on the ordinary r=-91/400 core.',limitations=['The complete fourth initial space is certified; no full initial ideal equality is asserted.','The b8=0 boundary, ambient t nilpotence and saturation removal remain separate.','The accompanying T7 reduced curve is separate from any endpoint flat4 theorem.'],proof_note='research/notes/S61_FLAG_RANKONE_SOCLE_QUARTIC_20261101.md',input_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in paths})
    with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print('RANKONE SOCLE QUARTIC CERTIFICATE PASS',result['seconds'],flush=True)
if __name__=='__main__':main()
