#!/usr/bin/env python3
"""True fourth jets on the three ordinary-flag linear frame artifact cores."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.time()
    frame='code/s61_flag_boundary_eight_frame_engine_20261101.py'
    raw='code/s61_flag_boundary_eight_raw4_engine_20261101.py'
    rank='code/s61_flag_boundary_global_rank4_engine_20261101.py'
    stages={name:f'code/s61_flag_linear_artifacts_{name}_engine_20261101.py' for name in ['frame','quadratic','lower4','project4','initial4','cover4','zero8']}
    native=['code/parametric_implicit_jets_20261101.cpp','code/flag_core_polynomial_rank_20261101.cpp','code/univariate_polynomial_minor_gcd_20261101.cpp','code/univariate_polynomial_gcd_factor_20261101.cpp']
    cases=[]
    with tempfile.TemporaryDirectory(prefix='flag-linear-fourth-proof-')as temp:
        wd=Path(temp);bins=[wd/f'native-{i}'for i in range(4)]
        for source,binary in zip(native,bins):
            subprocess.run(['c++','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),'-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(ROOT/source),str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(binary)],check=True)
        def run(engine,directory,extra=()):
            print('STAGE',engine,directory.name,time.time()-started,flush=True)
            subprocess.run([sys.executable,str(ROOT/engine),str(directory),*map(str,extra)],cwd=ROOT,check=True)
        run(frame,wd);run(stages['frame'],wd)
        for k,core in enumerate(['-6','-12/5','-119/20']):
            case=wd/f'core-{k}';run(stages['quadratic'],case);run(raw,case,[bins[0]])
            for name in ['lower4','project4','initial4']:run(stages[name],case)
            run(rank,case,[bins[1]]);run(stages['cover4'],case,bins[2:]);run(stages['zero8'],case)
            saved={name:json.loads((case/file).read_text()) for name,file in [('field_rank','flag-global-field-rank-result.json'),('core_cover','flag-eight-quartic-cover-result.json'),('zero8','flag-linear-zero8-result.json'),('lower_basis','flag-global-low-quad.json')]}
            cases.append(dict(core=core,**saved))
    paths=[frame,raw,rank,*stages.values(),*native,str(Path(__file__).relative_to(ROOT)),'code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','reports/S61_FLAG_UNIFORM_OPEN_20261101.json','reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json']
    result=dict(all_passed=True,cases=cases,hilbert_function=[1,9,36,97,201],t_order=4,raw135_rebuilt=True,all_nuisance_values_when_b8_nonzero=True,zero8_nonempty_open=True,all_sufficiently_general_beta_points=True,scope='Three cores -6,-12/5,-119/20 in the normalized b0=b12=1 flag chart; true jets only through order4.',limitations=['The true linear rank-drop core -4/25 is separate.','No full initial equality, ambient t-nilpotence or flag saturation removal.'],proof_note='research/notes/S61_FLAG_LINEAR_ARTIFACTS_QUARTIC_20261101.md',seconds=time.time()-started,input_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest()for path in paths})
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('THREE LINEAR FLAG CORE FOURTH COVERS PASS',result['seconds'],flush=True)
if __name__=='__main__':main()
