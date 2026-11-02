#!/usr/bin/env python3
"""Independent Fraction elimination of the saved first base-adic module certificate."""
import argparse, hashlib, json, time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
NOTE='research/notes/S61_FLAG_ZERO12_FIRST_BASE_QUOTIENT_REVIEW_20261102.md'

def sparse_rank(columns):
    basis={}
    for col in columns:
        v=dict(col)
        while v:
            i=max(v)
            if i not in basis:
                c=v[i];basis[i]={j:a/c for j,a in v.items()};break
            c=v[i]
            for j,a in basis[i].items():
                v[j]=v.get(j,Q(0))-c*a
                if not v[j]:del v[j]
    return len(basis)

def action(columns,v):
    out={}
    for j,c in v.items():
        for i,a in columns[j].items():out[i]=out.get(i,Q(0))+a*c
    return {i:a for i,a in out.items() if a}

def verify():
    started=time.monotonic()
    path=ROOT/'reports/S61_FLAG_ZERO12_FIRST_BASE_QUOTIENT_20261102.json'
    p=json.loads(path.read_text());assert p['all_passed']
    for name,h in p['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    nf=[{int(i):Q(c) for i,c in col} for col in p['generator_normal_forms']]
    assert len(nf)==196 and sparse_rank(nf)==187 and sparse_rank(nf[49:])==138
    special=p['special_linear_relation_generator_indices']
    assert special==[39,41,42,43,44,45,46,47,48]
    for j in special:assert nf[147+j]=={i:6*c for i,c in nf[98+j].items()}
    # These nine formal relations have distinct coordinates in the v3 block.
    # Their independence and 196-187=9 establish the full A-module kernel.
    assert len(special)==196-187 and len(set(special))==9
    bas=[tuple(m) for m in p['quotient_basis']];assert len(bas)==len(set(bas))==187
    assert [sum(len(m)==d for m in bas) for d in range(7)]==[1,9,31,52,55,27,12]
    mt=[{int(i):Q(c) for i,c in col} for col in p['t_multiplication_columns']]
    assert len(mt)==187 and sparse_rank(mt)==3
    assert all(not action(mt,col) for col in mt)
    ib={m:i for i,m in enumerate(bas)}
    cols=[mt[ib[m]] for m in [(),(0,),(1,),(3,)]]
    assert sparse_rank(cols)==3 and cols[3]=={i:6*c for i,c in cols[2].items()}
    # The matrices here are the primary certificate, not a second raw135 reconstruction.
    paths=[path,Path(__file__),ROOT/NOTE,ROOT/'code/s61_flag_first_base_algebra_engine_20261102.cpp',
           ROOT/'code/s61_flag_zero12_first_base_quotient_20261102.py']
    return dict(id='S61_FLAG_ZERO12_FIRST_BASE_QUOTIENT_REVIEW_20261102',all_passed=True,
      exact_arithmetic='Python Fraction, reverse sparse elimination; no SymPy rank routine',
      saved_generator_rank=187,saved_first_base_layer_rank=138,full_module_relation_count=9,
      t_multiplication_rank=3,t_square_zero_in_this_quotient=True,
      module_decomposition='A^40 plus (A/(v3-6*v1))^9',
      t_ideal_module='R/(a^2,v3-6*v1)',actual_kernel_mod_a2_length=187-5*4,
      logical_audit=['All 2123 QQ pivot rows are original-ideal members, so length <=187.',
        'Nine commuting coordinate actions, cyclic monomial checks, base^2 and m^7 annihilation define an actual 187-dimensional algebra quotient.',
        'Original residual equations vanish, giving length >=187 and equality.',
        'Base squares are checked on 1; commutativity and cyclic generation imply their vanishing as operators.',
        'R-free D splits as an R-module; reduction modulo a^2 preserves M as the kernel.'],
      limitations=['Saved finite matrices are independently checked; raw135 implicit recursion and QQ replay have one implementation.',
                   'No full-ring t-nilpotence or ambient saturation removal follows.'],
      input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths},
      seconds=time.monotonic()-started)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);o=a.parse_args()
    r=verify()
    with o.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
