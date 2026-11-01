#!/usr/bin/env python3
"""Rebuild the ordinary flag four-jet cover with b8=0 and b2 different from 2."""
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--flint-include', type=Path, default=ROOT/'work/flint-headers/flint/3.5.0/include')
    parser.add_argument('--mpfr-include', type=Path, default=ROOT/'work/flint-headers/mpfr/4.2.2/include')
    parser.add_argument('--flint-library', type=Path, default=ROOT/'work/m2-runtime/lib/libflint.24.0.dylib')
    args = parser.parse_args()
    start = time.time()
    old = {stage: f'code/s61_flag_boundary_eight_{stage}_engine_20261101.py' for stage in ['frame','inverse','quadratic','raw4']}
    new = {stage: f'code/s61_flag_boundary_global_{stage}_engine_20261101.py' for stage in ['quadratic','active3','project4','cubic','initial4','rank4','cover4']}
    new.update({stage:f'code/s61_flag_boundary_zero8_{stage}_engine_20261101.py' for stage in ['project4','cubic','initial4','cover4']})
    special='code/s61_flag_boundary_zero8_specialize_engine_20261101.py'
    native = ['code/parametric_implicit_jets_20261101.cpp','code/rational_coefficient_linear_solution_20261101.cpp','code/flag_core_polynomial_rank_20261101.cpp','code/univariate_polynomial_minor_gcd_20261101.cpp','code/univariate_polynomial_gcd_factor_20261101.cpp']
    with tempfile.TemporaryDirectory(prefix='flag-global-fourth-proof-') as directory:
        wd = Path(directory)
        binaries = [wd/f'native-{i}' for i in range(len(native))]
        for source, binary in zip(native, binaries):
            subprocess.run(['c++','-O3','-std=c++17','-I'+str(args.flint_include),'-I'+str(args.mpfr_include),'-I/opt/homebrew/include',str(ROOT/source),str(args.flint_library),'-L/opt/homebrew/lib','-lgmp','-o',str(binary)],check=True)
        stages = [(old[x],[]) for x in ['frame','inverse','quadratic']]
        stages += [(new['quadratic'],[]),(new['active3'],[binaries[1]]),(special,[]),(old['raw4'],[binaries[0]])]
        stages += [(new[x],[]) for x in ['project4','cubic','initial4']]
        stages += [(new['rank4'],[binaries[2]]),(new['cover4'],binaries[3:])]
        for engine, extra in stages:
            print('STAGE', engine, time.time()-start, flush=True)
            subprocess.run([sys.executable,str(ROOT/engine),str(wd),*map(str,extra)],cwd=ROOT,check=True)
        result = json.loads((wd/'flag-eight-quartic-cover-result.json').read_text())
        rank = json.loads((wd/'flag-global-field-rank-result.json').read_text())
        lower = {name: json.loads((wd/file).read_text()) for name,file in [('quadratic','flag-global-low-quad.json'),('active_cubic','flag-global-low-active-cubics.json'),('true_cubic','flag-global-cubic-low-inverse.json')]}
    paths = [special]+list(old.values())+list(new.values())+native+[str(Path(__file__).relative_to(ROOT)),'code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','reports/S61_FLAG_UNIFORM_OPEN_20261101.json','reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json','reports/S61_FLAG_BOUNDARY_QUADRATIC_COVER_20261101.json']
    result.update(seconds=time.time()-start,raw135_rebuilt=True,consistent_raw_four_jets=True,quadratic_native_sympy_crosscheck=True,field_rank_certificate=rank,lower_basis_certificates=lower,scope='Every value of the seven parameters with b8=0 and (b2-2)*D(b2)*(400*b2+91) nonzero.',limitations=['No full initial ideal equality, ambient t-nilpotence, or whole flag removal.','The b2=2 locus, four linear D zeros and rank-one socle core remain separate.'],proof_note='research/notes/S61_FLAG_BOUNDARY_ZERO8_QUARTIC_COVER_20261101.md',input_sha256={file:hashlib.sha256((ROOT/file).read_bytes()).hexdigest() for file in paths})
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print('ZERO8 ORDINARY FLAG FOURTH COVER PASS',result['seconds'],flush=True)
if __name__ == '__main__':
    main()
