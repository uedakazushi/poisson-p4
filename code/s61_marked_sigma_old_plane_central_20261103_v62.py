#!/usr/bin/env python3
"""Inspect full rho=0 marked-square central algebras at h^2+2, not just heavy quotients."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
from collections import defaultdict
from fractions import Fraction

from s61_marked_sigma_pencil_rank_20261103 import read_rows, ROOT, PRIME, WEIGHTS


def read_exact(path):
    """Reduce h^2=-2 over QQ before constructing the full normal ideal."""
    records=[]
    with path.open() as stream:
        assert tuple(map(int,stream.readline().split()))==(14,136)
        for row in range(136):
            weight,count=map(int,stream.readline().split())
            terms=defaultdict(lambda:defaultdict(Fraction))
            for _ in range(count):
                value,*powers=stream.readline().split()
                h,*exponents=map(int,powers)
                if row==135:
                    continue
                terms[tuple(exponents)][h%2]+=Fraction(value)*(-2)**(h//2)
            records.append((row,weight,{ex:{j:c for j,c in cs.items()if c}
                                        for ex,cs in terms.items()if any(cs.values())}))
        assert not stream.read().strip()
    return records


def program(records, root, include_rho=False,prime=PRIME,flat_only=False):
    count = 14 if include_rho else 13
    def expression(row):
        terms = []
        for exponents, coefficients in row.items():
            if not include_rho and exponents[-1]:
                continue
            if not prime:
                value='('+ '+'.join(f'({c.numerator}/{c.denominator})'+('*a' if h else '')
                                   for h,c in coefficients.items())+')'
            else:
                value = sum(c * pow(root, h, PRIME) for h, c in coefficients.items()) % PRIME
            if value:
                terms.append(str(value) + ''.join(f'*v{i}^{power}'
                             for i,power in enumerate(exponents[:count]) if power))
        return '+'.join(terms) or '0'
    setup=('A=QQ[h]; K=toField(A/ideal(h^2+2)); a=promote(h,K);\n' if not prime else '')
    base='K' if not prime else '(ZZ/32003)'
    source = (setup+'R='+base+'[' + ','.join(f'v{i}' for i in range(count))
              + ',Degrees=>{' + ','.join(map(str,WEIGHTS[:count])) + '}];\n'
              + 'I=ideal(' + ','.join(expression(row) for i,w,row in records if i<135) + ');\n')
    generators = [f'v{i}' for i in range(9)] + [f'v9-{"a" if not prime else root}*v10','v11']
    if flat_only:
        assert not include_rho
        source += ('J=ideal('+','.join(generators)+'); assert(isSubset(I,J));\n'
                   'print "FLAT PLANE QUOTIENT VERIFIED"; exit 0;\n')
        return source
    source += ('print "INPUT READY"; flush stdio; G=gb I; print "GB READY";\n'
              + 'print("DIM="|toString dim(R/I));\n'
              + 'print("HS="|toString hilbertSeries(R/I,Reduce=>true));\n'
              + 'print("GB_GENERATORS="|toString numColumns gens G);\n')
    if not include_rho:
        source += 'J=ideal(' + ','.join(generators) + '); assert(isSubset(I,J));\n'
        source += ('powers=apply(flatten entries gens J,f->select(1..12,n->f^n%G==0));\n'
                   'print("RADICAL_POWERS="|toString powers);\n'
                   'print("RADICAL_PROVED="|toString all(powers,v->#v>0));\n')
    source += ('print("T_POWER="|toString select(1..12,n->v11^n%G==0));\n'
               'print "CENTRAL INSPECTION COMPLETE"; exit 0;\n')
    return source


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--root',type=int,default=3083)
    parser.add_argument('--prime',type=int,choices=[0,PRIME],default=PRIME)
    parser.add_argument('--include-rho',action='store_true')
    parser.add_argument('--flat-quotient-only',action='store_true')
    parser.add_argument('--timeout',type=int,default=180)
    args=parser.parse_args()
    assert (args.root**2+2)%PRIME==0
    args.output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    expected=json.loads((ROOT/'reports/S61_MARKED_SIGMA_PENCIL_MODEL_20261102.json').read_text())['row_output_sha256']
    assert hashlib.sha256(args.rows.read_bytes()).hexdigest()==expected
    records=read_rows(args.rows) if args.prime else read_exact(args.rows)
    script=args.output/'central.m2';script.write_text(program(records,args.root,args.include_rho,args.prime,args.flat_quotient_only))
    executable=shutil.which('M2');assert executable
    log=args.output/'central.log'
    returncode=None;timed_out=False
    with log.open('w') as stream:
        try:
            result=subprocess.run([executable,'--script',str(script.resolve())],
                                  stdout=stream,stderr=subprocess.STDOUT,timeout=args.timeout)
            returncode=result.returncode
        except subprocess.TimeoutExpired:
            timed_out=True
    report=dict(id='S61_MARKED_SIGMA_OLD_PLANE_CENTRAL_20261103_V62',
                execution_completed=returncode==0,returncode=returncode,timed_out=timed_out,
                field='GF(32003)' if args.prime else 'QQ[h]/(h^2+2)',root=args.root,include_rho=args.include_rho,
                central_support_proved=(returncode==0 and 'RADICAL_PROVED=true' in log.read_text()),
                flat_quotient_verified=(returncode==0 and 'FLAT PLANE QUOTIENT VERIFIED' in log.read_text()),
                full_ambient_removal_proved=False,
                global_rows_sha256=hashlib.sha256(args.rows.read_bytes()).hexdigest(),
                input_sha256={str(Path(__file__).relative_to(ROOT)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                              'code/s61_marked_sigma_pencil_rank_20261103.py':hashlib.sha256((ROOT/'code/s61_marked_sigma_pencil_rank_20261103.py').read_bytes()).hexdigest()},
                seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
    if returncode or timed_out:
        raise SystemExit(1)


if __name__=='__main__':
    main()
