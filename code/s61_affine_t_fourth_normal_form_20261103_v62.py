#!/usr/bin/env python3
"""Expand the exact t^4 row on all1031 symbols and all17 parameters.

The optional dense scratch inverse cache accelerates the calculation only.
The output original134 cofactor identity is checked coefficient by coefficient
and can be verified independently of that cache. No t-nilpotence is asserted.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time

import numpy as np

from s61_affine_seventeen_base_lift_20261103_v62 import build,verify_coefficients
from s61_affine_t_action_circuit_20261103_v62 import ordered_monomials,read_native
from s61_affine_normal_projection_20261102 import ROOT,P


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrices(path):
    with path.open('rb')as stream:
        header=np.fromfile(stream,dtype='<u4',count=8).tolist()
        assert header==[0x53363154,29,17,134,435,4495,3795,700]
        result=[]
        for shape in [(134,134),(134,301),(3795,3795),(3795,700)]:
            result.append(np.fromfile(stream,dtype='<u2',count=int(np.prod(shape))).reshape(shape).astype(np.int64))
        assert stream.read()==b''
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--circuit',type=Path,default=ROOT/'reports/S61_AFFINE_T_ACTION_CIRCUIT_20261103_v62/native.txt')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    high,cross,base,pairs,cubics,quartics,bpairs,v,inputs,_,_=build()
    summary,piv2,rows3,piv3,_,_,_,_,_=read_native(args.circuit)
    assert summary['GENERATORS']==1031
    inv2,R2,inv3,R3=matrices(args.cache)
    free2=[i for i in range(435) if i not in set(piv2)]
    gcubics=ordered_monomials(29,3);ci={m:i for i,m in enumerate(cubics)}
    free3=[i for i in range(4495) if i not in set(piv3)]
    piv3lex=[ci[gcubics[i]] for i in piv3];free3lex=[ci[gcubics[i]] for i in free3]
    lift=ROOT/'reports/S61_AFFINE_SEVENTEEN_BASE_LIFT_20261103_v62/native.txt'
    H=np.zeros_like(high);pi={m:i for i,m in enumerate(pairs)}
    for line in lift.read_text().splitlines()[4:]:
        r,i,j,c=map(int,line.split());assert H[r,pi[i,j]]==0;H[r,pi[i,j]]=c
    A3,A2=verify_coefficients(high,cross,base,pairs,cubics,quartics,H,0)
    # The first quartic reduction gives -A3*Y-A2*YY in the full algebra.
    C3=(-A3[piv3lex].T)%P
    NF3=(-A3[free3lex].T+C3@R3)%P
    original3=C3@inv3%P
    L=np.zeros((17,134,29),dtype=np.int64)
    for k,row in enumerate(rows3):L[:,row//29,row%29]=(L[:,row//29,row%29]+original3[:,k])%P
    raw_lc=np.einsum('brk,ric->bcki',L,cross,optimize=True)%P
    lc_unknown=np.array([raw_lc[:,:,k,i]+(raw_lc[:,:,i,k] if k!=i else 0) for k,i in pairs]).transpose(1,2,0)%P
    LC=np.array([lc_unknown[b,c]+(lc_unknown[c,b] if b!=c else 0) for b,c in bpairs])%P
    poly2=(-A2.T-LC)%P
    C2=poly2[:,piv2]
    NF2=(poly2[:,free2]+C2@R2)%P
    D=C2@inv2%P
    LB=np.einsum('brk,rq->qbk',L,base,optimize=True)%P
    DC=np.einsum('qr,rkb->qbk',D,cross,optimize=True)%P
    btriples=list(itertools.combinations_with_replacement(range(17),3));bquartics=list(itertools.combinations_with_replacement(range(17),4))
    bt={m:i for i,m in enumerate(btriples)};bq={m:i for i,m in enumerate(bquartics)}
    triple_product=np.array([[bt[tuple(sorted((*pair,b)))] for b in range(17)] for pair in bpairs])
    fourth_product=np.array([[bq[tuple(sorted((*a,*b)))] for b in bpairs] for a in bpairs])
    NF1=np.zeros((len(btriples),29),dtype=np.int64)
    np.add.at(NF1,triple_product.ravel(),-(LB+DC).reshape(-1,29));NF1%=P
    DB=D@base%P
    NF0=np.zeros(len(bquartics),dtype=np.int64)
    np.add.at(NF0,fourth_product.ravel(),-DB.ravel());NF0%=P
    # Check all formal parameter/unknown coefficients against original134.
    # This does not infer the identity from numerical parameter samples.
    products23=np.array([[ci[tuple(sorted((*pair,k)))] for k in range(29)] for pair in pairs])
    LH=np.einsum('brk,rq->bkq',L,high,optimize=True)%P
    check3=A3.T.copy()
    for b in range(17):np.add.at(check3[b],products23.T.ravel(),LH[b].ravel())
    expected3=np.zeros_like(check3);expected3[:,free3lex]=-NF3
    assert np.array_equal(check3%P,expected3%P)
    check2=(A2.T+LC+D@high)%P
    expected2=np.zeros_like(check2);expected2[:,free2]=-NF2
    assert np.array_equal(check2%P,expected2%P)
    check1=np.zeros_like(NF1);np.add.at(check1,triple_product.ravel(),(LB+DC).reshape(-1,29))
    assert np.array_equal(check1%P,(-NF1)%P)
    check0=np.zeros_like(NF0);np.add.at(check0,fourth_product.ravel(),DB.ravel())
    assert np.array_equal(check0%P,(-NF0)%P)
    # Independent direct original-row accumulation for the mixed L*cross
    # tensor, which is the most index-sensitive part of the expansion.
    direct=np.zeros_like(raw_lc)
    for r in range(134):direct+=L[:,r,None,:,None]*cross[r].T[None,:,None,:]
    assert np.array_equal(direct%P,raw_lc)
    # Evaluate this global identity at the known-five fiber Y=v.
    v=v%P
    y2=np.array([int(v[i])*int(v[j])%P for i,j in bpairs],dtype=np.int64)
    y3=np.array([int(np.prod(v[list(m)],dtype=np.int64))%P for m in btriples],dtype=np.int64)
    y4=np.array([int(np.prod(v[list(m)],dtype=np.int64))%P for m in bquartics],dtype=np.int64)
    coeff3=v@NF3%P;coeff2=y2@NF2%P;coeff1=y3@NF1%P;coeff0=int(y4@NF0)%P
    known_path=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    known=json.loads(known_path.read_text())
    for point in known['known_five_mod_32003']:
        Z=np.array(point+[0]*25,dtype=np.int64)
        eval3=np.array([int(np.prod(Z[list(gcubics[i])],dtype=np.int64))%P for i in free3],dtype=np.int64)
        eval2=np.array([int(Z[pairs[i][0]])*int(Z[pairs[i][1]])%P for i in free2],dtype=np.int64)
        assert (int(coeff3@eval3)+int(coeff2@eval2)+int(coeff1@Z)+coeff0)%P==1
    identity=dict(convention='t^4=normal_form+sum_i(H_i(Z)+sum_b Y_b L_bi(Z)+sum_bc Y_bY_c D_bci)*F_i(Z,Y)',
                  basis='1,29 singletons,301 free lex pairs,700 free grevlex cubics; pivots in circuit certificate',
                  base_order=[f'Y{i}' for i in range(1,18)],
                  cofactor_quadratic=[[int(r),*pairs[k],int(H[r,k])] for r,k in zip(*np.nonzero(H))],
                  cofactor_base_linear_unknown_linear=[[int(b),int(r),int(k),int(L[b,r,k])] for b,r,k in zip(*np.nonzero(L))],
                  cofactor_base_quadratic_constant=[[*bpairs[q],int(r),int(D[q,r])] for q,r in zip(*np.nonzero(D))],
                  normal_base_linear_unknown_cubic=[[int(b),*gcubics[free3[j]],int(NF3[b,j])] for b,j in zip(*np.nonzero(NF3))],
                  normal_base_quadratic_unknown_quadratic=[[*bpairs[q],*pairs[free2[j]],int(NF2[q,j])] for q,j in zip(*np.nonzero(NF2))],
                  normal_base_cubic_unknown_linear=[[*btriples[q],int(k),int(NF1[q,k])] for q,k in zip(*np.nonzero(NF1))],
                  normal_base_quartic_constant=[[*bquartics[q],int(NF0[q])] for q in np.flatnonzero(NF0)])
    with gzip.open(args.output/'identity.json.gz','wt')as stream:json.dump(identity,stream,separators=(',',':'))
    counts=[int(np.count_nonzero(a)) for a in [NF3,NF2,NF1,NF0]]
    inputs += [Path(__file__),ROOT/'code/s61_affine_t_action_circuit_20261103_v62.py',
               ROOT/'code/s61_affine_seventeen_base_lift_20261103_v62.py',args.circuit,lift,known_path]
    report=dict(id='S61_AFFINE_T_FOURTH_NORMAL_FORM_20261103_v62',all_passed=True,prime=P,
                all17_parameter_identity_verified=True,original_quadrics=134,symbol_count=1031,
                nonzero_normal_terms_by_parameter_degree_1_2_3_4=counts,normal_term_count=sum(counts),
                nonzero_original_cofactor_terms=dict(quadratic=int(np.count_nonzero(H)),base_linear_unknown_linear=int(np.count_nonzero(L)),base_quadratic_constant=int(np.count_nonzero(D))),
                verified_coefficient_space_dimensions_by_parameter_degree_0_1_2_3_4=[35960,17*4495,153*435,969*29,4845],
                dense_mixed_tensor_independent_accumulation=True,known_five_evaluation=True,
                scratch_cache_sha256=sha(args.cache),scratch_cache_committed=False,
                t_nilpotence_global_claimed=False,full_affine_exhaustion_claimed=False,
                limitations=['Exact frozen-beta F32003 identity. GenericQQ analogous construction uses the same nonzero finite inverses.',
                             'The nonzero normal form is an actual all-parameter reduction; it does not remove the uncomputed relation kernel.'],
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(inputs)},
                evidence_sha256={p.name:sha(p) for p in args.output.iterdir() if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',counts,sum(counts),report['seconds'])


if __name__=='__main__':main()
