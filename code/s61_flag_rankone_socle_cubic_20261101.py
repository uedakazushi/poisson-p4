#!/usr/bin/env python3
"""Rebuild genuine cubic initial classes on the core r=-91/400 flag stratum."""
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
    parser.add_argument('--output',type=Path,default=ROOT/'reports/S61_FLAG_RANKONE_SOCLE_CUBIC_20261101.json')
    args=parser.parse_args();start=time.time()
    engines=[f'code/s61_flag_boundary_eight_{stage}_engine_20261101.py' for stage in ['frame','inverse','quadratic']]
    engines.append('code/s61_flag_rankone_socle_cubic_engine_20261101.py')
    with tempfile.TemporaryDirectory(prefix='flag-rankone-socle-') as directory:
        wd=Path(directory)
        for engine in engines:subprocess.run([sys.executable,str(ROOT/engine),str(wd)],cwd=ROOT,check=True)
        result=json.loads((wd/'flag-rankone-socle-cubic.json').read_text())
    assert result['quadratic_rank']==9 and result['quadratic_linear_multiple_rank']==64
    assert result['cubic_rank']==2 and result['rref_symbols']==['b13','b8']
    inputs=engines+['code/s61_flag_rankone_socle_cubic_20261101.py','code/s61_model.py',
        'code/s61_flag_transverse_obstruction_20261101.py','code/primitive_jet_elliptic_curve_20261101.py',
        'code/weight1_structure.py','reports/S61_FLAG_UNIFORM_OPEN_20261101.json',
        'reports/S61_FLAG_BOUNDARY_QUADRATIC_COVER_20261101.json']
    result.update(all_passed=True,seconds=time.time()-start,
        hilbert_function_through_degree_three=[1,9,36,99],
        scope='Function field of all seven nuisance parameters on r=-91/400; no all-point cubic rank assertion.',
        limitations=['All active cubics were retained before the 64-dimensional quadratic-product quotient.',
            'The two new cubic classes are exact over the function field, with all135 row relations checked.',
            'Fourth and higher initial relations, the full formal ring and saturation are unresolved.'],
        input_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('RANK-ONE SOCLE TRUE CUBIC PASS',result['seconds'],flush=True)

if __name__=='__main__':main()
