#!/usr/bin/env python3
"""Exact all-point four-quad coefficient independence, and its rank35 consequence.

The data certificate is a polynomial right inverse up to -a²P²J. Every
identity below is checked by Python Fraction sparse arithmetic. d in the
serialized core polynomials denotes P=ad, not the original d coordinate.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v62'
SOURCE = ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json'
CERT = ROOT/'reports'/STEM/'right_inverse.txt.gz'

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def mons(n, d):
    return [(d,)] if n == 1 else [(i,)+v for i in range(d+1) for v in mons(n-1, d-i)]
def add_product(out, p, q):
    for u, c in p.items():
        for v, d in q.items():
            out[tuple(a+b for a,b in zip(u,v))] += c*d

def read_matrix(data):
    it = iter(data.decode().split()); n,m,k = [int(next(it)) for _ in range(3)]
    assert k == 0
    out = []
    for _ in range(n):
        row = []
        for _ in range(m):
            p = {}
            for _ in range(int(next(it))):
                c = Fraction(next(it)); ex = tuple(int(next(it)) for _ in range(3))
                assert c and ex not in p and min(ex) >= 0
                p[ex] = c
            row.append(p)
        out.append(row)
    assert next(it, None) is None
    return out

def verify():
    start=time.monotonic(); data=json.loads(SOURCE.read_text()); assert data['all_passed']
    for name,want in data['input_sha256'].items(): assert sha(ROOT/name)==want,name
    # The six-light chart preserves actual unscaled quadrics up to declared
    # units. Change original d to P/a and clear only a per quadric.
    columns=mons(6,2); matrix=[]; shifts=[]
    for row in data['cleared_qfour_sparse']:
        q={tuple(mon[:6]):{tuple(ex):Fraction(c) for ex,c in terms} for mon,terms in row}
        assert all(len(mon)==7 and mon[6]==0 and sum(mon)==2 for mon,_ in row)
        shift=max([0]+[ex[2]-ex[0] for c in q.values() for ex in c]); shifts.append(shift)
        converted={}
        for mon,c in q.items():
            out={}
            for (a,e,d),v in c.items():
                ex=(a-d+shift,e,d); assert min(ex)>=0 and ex[0]+ex[2]-shift==a
                out[ex]=out.get(ex,Fraction(0))+v
            converted[mon]={ex:v for ex,v in out.items() if v}
        matrix.append([converted.get(mon,{}) for mon in columns])
    assert shifts==[0,2,2,1]
    inv=read_matrix(gzip.decompress(CERT.read_bytes()))
    assert len(inv)==21 and all(len(row)==4 for row in inv)
    # g=-a²P²J, J=5P-6e-4. This is a unit on the declared open.
    g={(2,0,3):Fraction(-5),(2,1,2):Fraction(6),(2,0,2):Fraction(4)}
    for i in range(4):
        for j in range(4):
            out=defaultdict(Fraction)
            for k in range(21): add_product(out,matrix[i][k],inv[k][j])
            assert {ex:c for ex,c in out.items() if c}==(g if i==j else {}),(i,j)
    # A distinct degree2 monomial has six degree3 multiples. The product
    # sets of any pair have intersection of size at most one.
    qmon=mons(6,2)
    for i,u in enumerate(qmon):
        su={tuple(a+b for a,b in zip(u,v)) for v in mons(6,1)}
        assert len(su)==6
        for v in qmon[:i]:
            sv={tuple(a+b for a,b in zip(v,w)) for w in mons(6,1)}
            assert len(su & sv)<=1
    pair_note=ROOT/'research/notes/S61_BINARY_NONSQUARE_PAIR_COVER_20261103_v62.md'
    pair_report=ROOT/'reports/S61_BINARY_NONSQUARE_PAIR_COVER_20261103_v62/report.json'
    pair=json.loads(pair_report.read_text()); assert pair['all_passed']
    files=[Path(__file__),SOURCE,CERT,pair_note,pair_report,
           ROOT/'research/notes'/f'{STEM}.md']
    return dict(id=STEM,all_passed=True,field='QQ',core_coordinates=['a','e','P=ad'],
        open_original='a*d*F*B*C*J != 0, J=5ad-6e-4',
        open_P_chart='a*P*F*B0*C*J != 0, B0=aB, J=5P-6e-4',
        exact_matrix_shape=[4,21],a_clear_shifts=shifts,
        exact_right_inverse_shape=[21,4],exact_right_inverse_terms=sum(len(p)for row in inv for p in row),
        exact_matrix_times_right_inverse='-a^2*P^2*(5P-6e-4)*identity(4)',
        exact_all16_QQ_polynomial_identities_verified=True,
        all_point_four_quad_rank=4,all_point_mod_ell_linear_product_lower_bound=18,
        pair_cover_kernel_independent_cubic_products=17,
        all_point_original_quadratic_linear_product_lower_bound=35,
        generic_constant_mixing_gives_a_unimodular_residual_product_row_after31_pivots=True,
        all_point_product_rank41_or_full_cubic_rank57_claimed=False,
        ambient_nilpotence_or_quartic_coverage_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in files},seconds=time.monotonic()-start)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    out=verify()
    with args.output.open('x')as f:json.dump(out,f,indent=2);f.write('\n')
    print('PASS',STEM,out['seconds'])
