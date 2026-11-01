#!/usr/bin/env python3
"""Independent QQ Groebner check of the one quotient that proves all base-adic cutoffs."""
import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.time()
    source='reports/S61_FLAG_ZERO12_RATIONAL_INITIAL_20261101.json'
    r=json.loads((ROOT/source).read_text())
    v=s.symbols('v:9');base={0,1,3};free=[v[i] for i in range(9) if i not in base]
    polys=[sum(s.Rational(c)*s.prod(v[i] for i in mon) for mon,c in row if not(base&set(mon))) for row in r['initial_relations']]
    G=s.groebner([p for p in polys if p],*free,domain=s.QQ,order='grevlex')
    lm=[p.LM(order='grevlex').exponents for p in G.polys];hf=[]
    for d in range(7):
        n=0
        for inds in itertools.combinations_with_replacement(range(6),d):
            ex=tuple(inds.count(i) for i in range(6))
            n+=not any(all(a>=b for a,b in zip(ex,lead)) for lead in lm)
        hf.append(n)
    assert hf==[1,6,13,15,10,6,0]
    assert all(G.reduce(p)[1]==0 for p in polys)
    paths=[str(Path(__file__).relative_to(ROOT)),source,
           'reports/S61_FLAG_ZERO12_BASE_ADIC_BOUNDS_20261102.json',
           'research/notes/FINITE_BASE_ADIC_CUTOFF_20261102.md',
           'research/notes/S61_FLAG_ZERO12_BASE_ADIC_BOUNDS_20261102.md']
    report=dict(all_passed=True,accepted=True,seconds=time.time()-started,
                independent_implementation='SymPy QQ grevlex after setting the three base variables to zero; no primary M2 invocation.',
                reduced_variables=list(map(str,free)),Groebner_basis=list(map(str,G.polys)),
                leading_monomials=lm,Hilbert_function_through_six=hf,
                graded_generation_bound=5,all_r_argument='Graded induction then actual local Nakayama; see review note.',
                limitations=['Only r=1 is independently recomputed; this suffices for the stated all-r bound.',
                             'Saved true initial relations are dependencies, not independently rederived from raw135.',
                             'No actual quotient length equality or ambient t nilpotence is inferred.'],
                input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print('ZERO12 BASE ADIC REVIEW PASS',report['seconds'])

if __name__=='__main__':main()
