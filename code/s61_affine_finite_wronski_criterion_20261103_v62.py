#!/usr/bin/env python3
"""A finite87-value and one algebraic-degree87 UNIT criterion for surfaces."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    previous=ROOT/'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62'
    report=json.loads((previous/'report.json').read_text());assert report['all_passed']
    for name,value in report['input_sha256'].items():assert sha(ROOT/name)==value
    for name,value in report['evidence_sha256'].items():assert sha(previous/name)==value
    n=45;bound=2*n-4;count=bound+1;assert(bound,count)==(86,87)
    values=list(range(count));vandermonde=math.prod(j-i for i in values for j in values if i<j)
    assert vandermonde==math.prod(math.factorial(j)for j in range(count)) and vandermonde>0
    factorization={str(p):int(sum(sum(j//(p**k)for k in range(1,10)if p**k<=j)for j in range(count)))for p in sp.primerange(2,count)}
    assert math.prod(int(p)**e for p,e in factorization.items())==vandermonde
    z=sp.Symbol('z');rational_test=sp.Poly(sp.prod(z-v for v in values),z)
    assert rational_test.degree()==87 and sp.gcd(rational_test,rational_test.diff()).degree()==0
    eisenstein=sp.Poly(z**87-2,z);coefs=eisenstein.all_coeffs()
    assert coefs[0]==1 and all(c%2==0 for c in coefs[1:]) and coefs[-1]%4!=0
    assert sp.expand(z*sp.diff(eisenstein.as_expr(),z)-87*eisenstein.as_expr())==174
    source=(previous/'wronski-exact.sing').read_text();raw=source.split('ideal I=',1)[1].split(';\n',1)[0].replace('zz1','zz')
    setup='A=QQ[ss1,ss2,ss3,ss4,ss5,cc1,cc2];K=frac A;'
    full=setup+'R=K[zz,'+','.join(f'x{i}'for i in range(45))+'];\nI=ideal(zz^87-2,'+raw+');\n'
    full+='print "ORDINARY87_EXACT_MODEL_START";stdio << flush;G=gb(I,ChangeMatrix=>true);print "ORDINARY87_EXACT_GB_DONE";stdio << flush;\n'
    full+='if 1_R%G==0 then (C=matrix{{1_R}}//(gens I);assert((gens I)*C==matrix{{1_R}});print "ORDINARY87_EXACT_UNIT_COFACTOR_PASS";stdio << flush;) else (print("NONUNIT_DIM="|toString dim(R/I));stdio << flush;);exit 0;\n'
    (args.output/'algebraic87-polynomial-exact.m2').write_text(full)
    toy=setup+'R=K[zz,x,y];I=ideal(zz^87-2,x*x-4/5*ss1*x,y-cc1);G=gb(I,ChangeMatrix=>true);C=(gens G)//(gens I);assert((gens I)*C==gens G);assert(numColumns basis(R/I)==174);print "ORDINARY87_FUNCTION_FIELD_TRACE_PASS";stdio << flush;exit 0;\n'
    probe=args.output/'coefficient-probe.m2';log=args.output/'coefficient-probe.log';probe.write_text(toy)
    with log.open('w')as stream:subprocess.run([shutil.which('M2'),'--script',str(probe.resolve())],stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=30)
    assert 'ORDINARY87_FUNCTION_FIELD_TRACE_PASS' in log.read_text()
    historical=ROOT/'work/affine-infinity-lift-v62/finite87-incomplete-nested/coefficient-probe.log'
    shutil.copyfile(historical,args.output/'nested-field-unsupported.log')
    assert 'expected coefficient ring of the form ZZ/n, ZZ, QQ, or GF' in historical.read_text()
    data=dict(wronskian_degree_bound=86,rational_test_values=values,
              rational_test_polynomial_coefficients=list(map(str,rational_test.all_coeffs())),
              vandermonde_determinant_prime_factorization=factorization,
              vandermonde_determinant_bit_length=vandermonde.bit_length(),
              single_algebraic_test_minimal_polynomial='aa^87-2',eisenstein_prime=2,
              separability_identity='aa*fprime-87*f=174',coefficient_field='QQ(aa)(ss1,...,ss5,cc1,cc2)',
              ordinary_polynomial_presentation='QQ(ss1,...,ss5,cc1,cc2)[zz,x0,...,x44]/(zz^87-2)',
              nested_field_GB_status='unsupported by CloudM2; retained failure log',
              ordinary_presentation_GB_status='exact174-dimensional toy quotient and change matrix passed')
    (args.output/'finite-criterion.json').write_text(json.dumps(data,indent=2)+'\n')
    inputs=[Path(__file__),previous/'report.json',previous/'wronski-exact.sing',
            ROOT/'research/notes/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62.md']
    answer=dict(id='S61_AFFINE_FINITE_WRONSKI_CRITERION_20261103_v62',all_passed=True,field='QQ',
                transcendental_parameter_count=7,algebraic_constant_degree=87,original_equations=134,linear_projection_equations=2,
                all87_rational_unit_tests_equivalent_to_no_surfaces=True,
                single_algebraic87_unit_test_equivalent_to_no_surfaces=True,
                irreducibility_Eisenstein_checked=True,separability_checked=True,finite_evaluation_injectivity_checked=True,
                exact_ordinary_polynomial_number_field_trace_probe=True,nested_number_field_GB_unsupported=True,
                actual_full_unit_result_obtained=False,surface_exclusion_claimed=False,five_point_exhaustion_claimed=False,
                limitations=['The finite criterion is proved; the full87-degree algebraic model was not solved by this producer.',
                             'This is a UNIT criterion. The degree of a nonunit Wronski fiber is not identified with projective surface degree.'],
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(answer,indent=2)+'\n');print('PASS',count,'values; algebraic degree87;7 transcendental parameters',answer['seconds'])


if __name__=='__main__':main()
