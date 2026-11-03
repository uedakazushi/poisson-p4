#!/usr/bin/env python3
"""Certify the global degree<=1031 t-determinant annihilator and its constant.

Uses the exact global lift phi*T=t*phi; no detK localization is involved.
The central block matrices are reconstructed from the saved exact reduction.
No billion-degree resultant or symbolic dense1031 matrix need be expanded.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import numpy as np

from s61_affine_t_action_circuit_20261103_v62 import ordered_monomials,read_native
from s61_affine_t_fourth_normal_form_20261103_v62 import matrices
from s61_affine_normal_projection_20261102 import ROOT,P,pivot_minor


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    action=ROOT/'reports/S61_AFFINE_T_ACTION_CIRCUIT_20261103_v62'
    previous=json.loads((action/'report.json').read_text());assert previous['all_passed']
    for name,value in previous['input_sha256'].items():assert sha(ROOT/name)==value
    for name,value in previous['evidence_sha256'].items():assert sha(action/name)==value
    summary,piv2,_,piv3,_,_,_,_,_=read_native(action/'native.txt')
    assert summary['GENERATORS']==1031
    _,R2,_,R3=matrices(args.cache)
    pairs=list(itertools.combinations_with_replacement(range(29),2));pi={m:i for i,m in enumerate(pairs)}
    cubes=ordered_monomials(29,3);ci={m:i for i,m in enumerate(cubes)}
    free2=[i for i in range(435)if i not in set(piv2)];free3=[i for i in range(4495)if i not in set(piv3)]
    pr2={j:i for i,j in enumerate(piv2)};fr2={j:i for i,j in enumerate(free2)}
    pr3={j:i for i,j in enumerate(piv3)};fr3={j:i for i,j in enumerate(free3)}
    L=np.zeros((301,29),dtype=np.int64);C=np.zeros((700,301),dtype=np.int64)
    for i in range(29):
        k=pi[0,i]
        if k in pr2:L[:,i]=R2[pr2[k]]
        else:L[fr2[k],i]=1
    for j,k in enumerate(free2):
        mon=tuple(sorted((0,*pairs[k])));index=ci[mon]
        if index in pr3:C[:,j]=R3[pr3[index]]
        else:C[fr3[index],j]=1
    CL=C@L%P
    rL=len(pivot_minor(L)[0]);rC=len(pivot_minor(C)[0]);rCL=len(pivot_minor(CL)[0])
    rank1=1+rL+rC;rank2=int(np.any(L[:,0]))+rCL;rank3=int(np.any(CL[:,0]))
    ranks=[1031,rank1,rank2,rank3,0]
    assert rank3==1
    # Degree blocks are0,1,2,3. Multiplication by t atY=0 has only
    # (1,0),(2,1),(3,2) blocks; hence its fourth power is identically zero.
    np.savez_compressed(args.output/'central-action-blocks.npz',linear_to_quadratic=L.astype(np.uint16),
                        quadratic_to_cubic=C.astype(np.uint16),linear_to_cubic=CL.astype(np.uint16))
    definition=dict(operator='Exact corrected T_t(Y) on B^1031, phi*T_t=t*phi',
                    global_annihilator='G(Y)=det(T_t(Y)-I_1031)',
                    global_cofactor_recipe=['u=adj(T_t-I)*e0','t*poly_E(u)-poly_E(T_t*u)=sum_i K_i(u)*F_i',
                                           'G=(t-1)*poly_E(u)-sum_i K_i(u)*F_i'],
                    generic_map_recipe='Fix any full quartic minor nonzero at the saved point; use the exact p-integral correction in the spread note.',
                    entry_parameter_degree='1+degree(source_symbol)-degree(target_symbol)',
                    maximum_total_base_degree=1031,constant_coefficient=P-1,
                    central_action_rank_sequence=ranks,central_nilpotence_exponent=4,
                    determinant_zero_locus_resolved=False)
    (args.output/'annihilator-circuit.json').write_text(json.dumps(definition,indent=2)+'\n')
    spread=ROOT/'research/notes/S61_AFFINE_GLOBAL_ACTION_SPREAD_20261103_v62.md'
    inputs=[Path(__file__),ROOT/'code/s61_affine_t_fourth_normal_form_20261103_v62.py',
            ROOT/'code/s61_affine_t_action_circuit_20261103_v62.py',action/'native.txt',action/'report.json',spread]
    report=dict(id='S61_AFFINE_T_DETERMINANT_ANNIHILATOR_20261103_v62',all_passed=True,prime=P,
                parameter_count=17,symbol_count=1031,full_original_equations=134,
                global_annihilator_without_cyclic_localization=True,maximum_total_base_degree=1031,
                nonzero_constant_coefficient_mod_prime=P-1,constant_coefficient_QQ=-1,
                central_action_rank_sequence=ranks,central_nilpotence_exponent=4,
                full_polynomial_coefficients_serialized=False,representation='exact cofactor/determinant circuit',
                scratch_cache_sha256=sha(args.cache),scratch_cache_committed=False,
                full_annihilator_zero_locus_resolved=False,finite_pair_rank_upper29_claimed=False,
                full_affine_exhaustion_claimed=False,
                limitations=['The global polynomial is a specified finite inverse/determinant recipe; the1031-degree coefficient expansion is not stored.',
                             'The nonzero constant follows from the central grading and holds for the generic corrected map.',
                             'A nonzero support equation does not exhaust its zero locus or exclude affine surfaces.'],
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',ranks,report['seconds'])


if __name__=='__main__':main()
