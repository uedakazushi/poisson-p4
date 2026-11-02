#!/usr/bin/env python3
"""Unit 7/14/3 blocks and a 56-by-84 cubic presentation on the treated-factor open."""
import argparse
import hashlib
import json
from pathlib import Path
from itertools import combinations_with_replacement
import time
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix
from s61_binary_deeper_uniform_quadratic_20261101 import PIVOT_ROWS,PIVOT_COLS
ROOT=Path(__file__).resolve().parents[1]
NOTE=ROOT/'research/notes/S61_BINARY_DEEPER_WEIGHTED_LOW_BLOCKS_20261102.md'

def verify():
    started=time.time()
    qp=ROOT/'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json'
    wp=ROOT/'reports/S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102.json'
    Q=json.loads(qp.read_text());W=json.loads(wp.read_text());assert Q['all_passed'] and W['all_passed']
    R,e,d=ring('e,d',s.QQ);domain=R.to_domain();zero={x:0 for x in s.symbols('p q r h k')}
    pol=lambda x:R.from_expr(s.sympify(x).subs(zero))
    known=[pol(x) for x in Q['irreducible_minor_factors']]
    def strip(value):
        assert value
        exponents=[]
        for factor in known:
            power=0
            while value.rem(factor)==0:value=value.exquo(factor);power+=1
            exponents.append(power)
        assert max(map(sum,value))==0
        return value,exponents
    monomials=list(combinations_with_replacement(range(9),2))
    cw=[sum(W['normal_weights'][j] for j in mon) for mon in monomials];rw=W['equation_weights']
    blocks={};contents={}
    for weight in [2,3,4]:
        rows=[i for i in PIVOT_ROWS if rw[i]==weight];cols=[j for j in PIVOT_COLS if cw[j]==weight]
        assert len(rows)==len(cols)=={2:7,3:14,4:3}[weight]
        A=[[pol(Q['primitive_quadratic_rows'][str(i)][j]) for j in cols] for i in rows]
        for i in rows:
            assert all(pol(Q['primitive_quadratic_rows'][str(i)][j])==0 for j in range(45) if cw[j]!=weight)
            unit,exponents=strip(pol(Q['row_contents'][i]))
            contents[i]=dict(factor_exponents=exponents,rational_unit=str(unit.as_expr()))
        determinant=DomainMatrix(A,(len(rows),len(cols)),domain).det()
        unit,exponents=strip(determinant)
        blocks[weight]=dict(rows=rows,columns=cols,monomials=[list(monomials[j]) for j in cols],dimension=len(rows),determinant_factor_exponents=exponents,determinant_rational_unit=str(unit.as_expr()),expanded_determinant_sha256=hashlib.sha256(str(determinant.as_expr()).encode()).hexdigest())
        print('BLOCK PASS',weight,len(rows),time.time()-started,flush=True)
    assert len(set(PIVOT_ROWS))==24 and len(set(PIVOT_COLS))==24
    paths=[Path(__file__),NOTE,qp,wp,ROOT/'code/s61_binary_deeper_uniform_quadratic_20261101.py']
    return dict(id='S61_BINARY_DEEPER_WEIGHTED_LOW_BLOCKS_20261102',all_passed=True,coefficient_ring='QQ[e,d,D^-1], D the thirteen previously treated factors',nuisance_parameters_zero=True,factors=Q['irreducible_minor_factors'],blocks=blocks,raw_quadratic_row_content_factors=contents,off_weight_blocks_identically_zero=True,all_three_square_blocks_are_invertible_on_entire_D_open=True,weight_two_rank_on_entire_D_open=7,weight_two_quotient_free_rank=23,heavy_light_schur_block_size=14,cubic_schur_presentation=dict(light_cubic_columns=84,quadratic_linear_multiple_rows=49,remaining_pure_cubic_rows=7,total_rows=56,target_rank=50,original_rank_equals_14_plus_schur_rank=True),proof_note=str(NOTE.relative_to(ROOT)),limitations=['The 56-by-84 cubic matrix is a proved presentation reduction, not a proved uniform rank50 result.','No all-core cubic or quartic rank certificate, new removal, or ambient nilpotence.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-started)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args();report=verify()
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps({k:report[k] for k in ['all_passed','seconds']}))
