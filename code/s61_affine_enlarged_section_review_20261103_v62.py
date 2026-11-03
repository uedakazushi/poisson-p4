#!/usr/bin/env python3
"""Independently count the saved leading ideal and verify all report hashes.

This does not rerun Groebner arithmetic or independently establish its basis
membership. It checks the finite-degree rank used by the section proof.
"""
import argparse
from collections import Counter
import hashlib
from itertools import combinations_with_replacement
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=json.loads(args.report.read_text());directory=args.report.parent
    assert report['all_passed']
    for name,value in report['input_sha256'].items():
        assert sha(ROOT/name)==value,name
    for name,value in report['evidence_sha256'].items():
        assert sha(directory/name)==value,name
    n=report['projective_section_dimension']+1
    source=directory/'leading-exponents.txt'
    leading=[]
    for line in source.read_text().splitlines():
        exponents=json.loads(line.replace('{','[').replace('}',']'))
        assert len(exponents)==n and all(x>=0 for x in exponents)
        leading.append(tuple(i for i,e in enumerate(exponents) for _ in range(e)))
    degree=report['first_full_five_ideal_degree'];ranks={};hf=[1]
    for d in range(1,degree+1):
        eliminated=set()
        for monomial in leading:
            if len(monomial)>d:
                continue
            for shift in combinations_with_replacement(range(n),d-len(monomial)):
                eliminated.add(tuple(sorted(monomial+shift)))
        ranks[d]=len(eliminated)
        hf.append(math.comb(n+d-1,d)-len(eliminated))
    assert hf==report['hilbert_function_mod_32003'][:degree+1]
    assert hf[-1]==5
    answer=dict(all_passed=True,scope='saved leading monomial count and report hash verification',
                leading_generator_degree_counts=dict(Counter(map(len,leading))),
                degree_space_ranks=ranks,hilbert_function_through_equality=hf,
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in [Path(__file__),args.report,source]})
    with args.output.open('x') as stream:
        json.dump(answer,stream,indent=2);stream.write('\n')
    print('PASS',hf,ranks)


if __name__=='__main__':
    main()
