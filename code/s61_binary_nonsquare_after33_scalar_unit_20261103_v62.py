#!/usr/bin/env python3
"""Exact known-unit scalar row after the33-pivot primary cubic reduction."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import tempfile
import time
ROOT=Path(__file__).resolve().parents[1];STEM='S61_BINARY_NONSQUARE_AFTER33_SCALAR_UNIT_20261103_v62';DIR=ROOT/'reports'/STEM

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def polynomials(it,n):
    out=[]
    for i in range(n):
        q={}
        for j in range(int(next(it))):
            c=Fraction(next(it));ex=tuple(int(next(it))for _ in range(3));assert c and ex not in q;q[ex]=c
        out.append(q)
    return out

def verify():
    t=time.monotonic();info=json.loads((DIR/'source.json').read_text());prior=ROOT/'reports/S61_BINARY_NONSQUARE_SECOND_SCALAR_PIVOT_20261103_v62/report.json';P=json.loads(prior.read_text());assert P['all_passed'];source=ROOT/info['source'];assert info['source_sha256']==P['generated_output']['sha256'] and sha(source)==info['source_sha256']
    rb=gzip.decompress((DIR/'source_row.txt.gz').read_bytes());wb=gzip.decompress((DIR/'cofactor.txt.gz').read_bytes());it=iter(rb.decode().split());assert[int(next(it))for _ in range(3)]==[1,58,0];row=polynomials(it,58);assert next(it,None)is None
    with source.open()as f:
        assert list(map(int,f.readline().split()))==[24,58,1];f.readline();f.readline()
        for i in range(8):
            raw=[f.readline()for _ in range(58)]
            if i==7:
                parsed=[]
                for s in raw:parsed.extend(polynomials(iter(s.split()),1))
                assert parsed==row
    it=iter(wb.decode().split());assert int(next(it))==58;L=polynomials(it,58);assert next(it,None)is None
    nt=sum(len(q)for q in L);degree=max(sum(ex)for q in L for ex in q);bits=max(max(abs(c.numerator).bit_length(),c.denominator.bit_length())for q in L for c in q.values());assert(nt,degree,bits)==(975,7,3299)
    cpp=ROOT/'code/polynomial_exact_scalar_unit_check_20261103_v62.cpp';units='a\n500*a+90*e^2+3*e-24\n(50*a+18*e^2-9*e-2)*d+(60*e-40)*a\n100*a+18*e^2-9*e-2\n5*d-6*e-4\nd\n'
    with tempfile.TemporaryDirectory(prefix='nonsquare-scalar-')as tmp:
        w=Path(tmp);rp=w/'row';wp=w/'witness';up=w/'units';rp.write_bytes(rb);wp.write_bytes(wb);up.write_text(units);binary=w/'check'
        flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-O2','-std=c++17',str(cpp),'-o',str(binary),*flags],check=True,capture_output=True,text=True)
        result=subprocess.run([str(binary),str(rp),str(wp),str(up)],check=True,capture_output=True,text=True);assert'EXACT_FULL_QQ_SCALAR_UNIT_IDENTITY_PASS'in result.stdout and'EXPONENTS 1 0 0 0 1 2'in result.stdout
    files=[Path(__file__),cpp,ROOT/'code/polynomial_bounded_degree_unit_lift_20261103_v62.cpp',prior,*[DIR/s for s in ['source.json','source_row.txt.gz','cofactor.txt.gz','D3-unit.log','D4-unit.log','independent-unit.log']],ROOT/'research/notes'/f'{STEM}.md']
    return dict(id=STEM,all_passed=True,field='QQ',source=info['source'],source_sha256=info['source_sha256'],source_shape=[24,58],target_product_row=7,
        coefficient_identity='source_row7*cofactor=a*P^2*(5P-6e-4)',known_unit_exponents_in_a_F_B0_C_J_P=[1,0,0,0,1,2],known_unit_scalar='1',
        cofactor_terms=nt,active_components=sum(bool(q)for q in L),cofactor_degree=degree,max_coefficient_bits=bits,
        discovery_degree_bound_on_primitive_entry_cofactors=4,exact_QQ_Dixon_square_size=1008,
        independent_exact_original_QQ_row_identity_verified=True,
        next_unit_row_available=True,executed_native_unit_pivots=33,native34th_pivot_executed=False,
        expected_append_pivot_shape_if_no_unit_coefficient=[23,58],replacement_not_yet_certified=True,
        all_point_product41_or_full57_or_quartic_coverage_claimed=False,
        input_sha256={str(q.relative_to(ROOT)):sha(q)for q in files},seconds=time.monotonic()-t)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=verify()
    with args.output.open('x')as f:json.dump(out,f,indent=2);f.write('\n')
    print('PASS',STEM,out['seconds'])
