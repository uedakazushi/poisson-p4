#!/usr/bin/env python3
"""G: common-eigenvector obstruction for the affine incidence rank-three branch.

Only a nonzero minor is transported from GF(32003) to characteristic zero.
No zero/rank-deficient modular result is called a generic upper bound.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import resource
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
P = 32003
CROSS = 'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/cross-input.txt.gz'


def mod(c):
    c = Fraction(c)
    assert c.denominator % P
    return c.numerator * pow(c.denominator, -1, P) % P


def rref(a):
    a = a.copy() % P
    ids = list(range(a.shape[0]))
    pivots, rank = [], 0
    for j in range(a.shape[1]):
        k = next((k for k in range(rank, a.shape[0]) if a[k,j]), None)
        if k is None:
            continue
        a[[rank,k]] = a[[k,rank]]
        ids[rank], ids[k] = ids[k], ids[rank]
        a[rank] = a[rank] * pow(int(a[rank,j]), -1, P) % P
        for k in range(a.shape[0]):
            if k != rank and a[k,j]:
                a[k] = (a[k] - a[k,j]*a[rank]) % P
        pivots.append(j)
        rank += 1
        if rank == a.shape[0]:
            break
    return a, pivots, ids[:rank]


def nullspace(a):
    rr, pivots, _ = rref(a)
    free = [j for j in range(a.shape[1]) if j not in pivots]
    v = np.zeros((a.shape[1],len(free)), dtype=np.int64)
    for k,j in enumerate(free):
        v[j,k] = 1
        for i,pivot in enumerate(pivots):
            v[pivot,k] = -rr[i,j] % P
    assert not (a @ v % P).any()
    return v


def determinant(a):
    a = a.copy() % P
    answer = 1
    for j in range(a.shape[0]):
        k = next((k for k in range(j, a.shape[0]) if a[k,j]), None)
        if k is None:
            return 0
        if k != j:
            a[[j,k]] = a[[k,j]]
            answer = -answer % P
        value = int(a[j,j])
        answer = answer*value % P
        a[j] = a[j]*pow(value,-1,P) % P
        for k in range(j+1,a.shape[0]):
            a[k] = (a[k]-a[k,j]*a[j]) % P
    return answer


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(30,31))
    resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,)*2)
    start=time.monotonic()
    lines=gzip.decompress((ROOT/CROSS).read_bytes()).decode().splitlines()
    assert lines[0]=='129 126 945' and len(lines)==259
    C=np.array([[mod(c) for c in line.split()] for line in lines[1:130]],dtype=np.int64)
    D=np.array([[mod(c) for c in line.split()[:42]] for line in lines[130:]],dtype=np.int64)
    _,_,ids=rref(C)
    assert len(ids)==126
    selected=C[ids]
    determinant_C=determinant(selected)
    assert determinant_C
    aug, pivots, _=rref(np.concatenate((selected,np.eye(126,dtype=np.int64)),axis=1))
    assert pivots==list(range(126))
    inverse=aug[:,126:]
    assert np.array_equal(selected @ inverse % P,np.eye(126,dtype=np.int64))
    L=np.zeros((126,129),dtype=np.int64)
    L[:,ids]=inverse
    assert np.array_equal(L @ C % P,np.eye(126,dtype=np.int64))
    T=(-L @ D) % P
    operators=[T[42*h:42*(h+1)] for h in range(3)]
    H=(D+C @ T) % P
    pairs=[(0,1),(0,2),(1,2)]
    commutators=[(operators[i] @ operators[j]-operators[j] @ operators[i]) % P for i,j in pairs]
    comm=np.concatenate(commutators,axis=0)
    _,_,selected_comm_rows=rref(comm)
    obstruction=np.concatenate((comm,H),axis=0)
    V=nullspace(obstruction)
    dims=[V.shape[1]]
    while V.shape[1]:
        annihilator=nullspace(V.T).T
        constraints=np.concatenate([annihilator @ op @ V % P for op in operators],axis=0)
        smaller=V @ nullspace(constraints) % P
        if smaller.shape[1]==V.shape[1]:
            break
        V=smaller
        dims.append(V.shape[1])
    minor=None
    if len(selected_comm_rows)==42:
        minor=determinant(comm[selected_comm_rows])
        assert minor
    result={
        'task_id':'G','prime':P,
        'exact_definition':'C=X129x126 and D=D1_129x42 over QQ; S=selected rows; L=C_S^(-1)pi_S; T_h=-(LD)_h',
        'source_cross_path':CROSS,
        'source_cross_sha256':hashlib.sha256((ROOT/CROSS).read_bytes()).hexdigest(),
        'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'all_cross_and_D1_denominators_p_units':True,
        'selected_cross_rows':ids,'selected_cross_determinant_mod_p':determinant_C,
        'left_inverse_mod_p':L.tolist(),'operators_T3_mod_p':[x.tolist() for x in operators],
        'pure_linear_remainder_H_mod_p':H.tolist(),
        'commutator_pairs':pairs,'commutator_rank_lower_bound':len(selected_comm_rows),
        'selected_commutator_rows':selected_comm_rows,'selected42_minor_mod_p':minor,
        'common_invariant_subspace_dimensions_mod_p':dims,
        'final_invariant_subspace_mod_p':V.tolist(),
        'rank3_branch_excluded_over_QQ_by_direct42_minor':minor is not None,
        'generic_transport_requires_same_rational_X_D1_specialization':True,
        'all_generic_residual_excluded':False,
        'seconds':time.monotonic()-start,
        'scope':'A nonzero42minor proves no common eigenvector over algebraic closure of QQ for this source and any generic rational model specializing to it. Rank-deficient modular output alone proves no characteristic-zero upper bound.',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:
        json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ['selected_cross_determinant_mod_p','commutator_rank_lower_bound','selected42_minor_mod_p','common_invariant_subspace_dimensions_mod_p','rank3_branch_excluded_over_QQ_by_direct42_minor','seconds']},indent=2))


if __name__=='__main__':
    main()
