#!/usr/bin/env python3
"""One degree22 product probe on only the 36 existing accepted QQ points."""
import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import numpy as np
import sympy as s

from q_new_family_degree16_invariant_probe_20261103_v63 import ROOT, PRIME, read, sha, matrix, mod, modmat, necklaces, echelon, character
from weight1_structure import sp_basis

STEM = 'Q_NEW_FAMILY_DEGREE22_INVARIANT_PROBE_20261103_v63'


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2_200_000_000, 2_200_000_000))
    out.mkdir(parents=True, exist_ok=False)
    invpath = ROOT / 'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    p8 = ROOT / 'reports/Q_NEW_FAMILY_DEGREE8_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
    p14 = ROOT / 'reports/Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
    p18 = ROOT / 'reports/Q_NEW_FAMILY_DEGREE18_INVARIANT_EXCLUSION_20261103_v63/samples.json'
    p36 = ROOT / 'reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63/samples.json'
    info, d8, d14, d18, samples = map(read, [invpath, p8, p14, p18, p36])
    w4 = necklaces(4)
    basis8 = d8['independent_five_columns']
    assert basis8 == list(range(5))
    pairs = d14['degree14_bracket_pairs']
    selected14 = d14['degree14_basis_columns']
    bracket_pairs = [pairs[i] for i in selected14]
    assert bracket_pairs == [[1, 2], [1, 3], [2, 3]]
    words18 = d18['selected_words']
    assert len(words18) == 9 and all(len(w) == 9 for w in words18)
    chars = character(list(map(tuple, info['Cartan_weights'])), 22)
    assert chars[22]['invariant_dimension'] == 21
    omega = matrix(info['symplectic_form16'])
    gs = sp_basis()
    gram = s.Matrix([[s.trace(A * B) for B in gs] for A in gs])
    Om, Oi, Gi = map(modmat, [omega, omega.inv(), gram.inv()])
    reps = [modmat(matrix(A)) for A in info['all_ten_representation_matrices']]
    Or = [Om @ A % PRIME for A in reps]
    pairs8_14 = list(itertools.product(range(5), range(3)))
    definitions = [dict(kind='F4_times_Inv18_basis',basis_index=i,trace9_word=w) for i,w in enumerate(words18)]
    definitions += [dict(kind='Inv8_times_Inv14',basis8=i,basis14=j,trace4_word=w4[i],bracket_pair=bracket_pairs[j])
                    for i,j in pairs8_14]
    evaluations, factor_values, betas = [], [], []
    direct_gradient_check = None
    for sample_index, sample in enumerate(samples['actual_QQ_source_points']):
        beta = [mod(x) for x in sample['beta16']]
        assert mod(sample['value_pivot']) and mod(sample['value_pfaffian'])
        b = np.array(beta,dtype=np.int64).reshape(16,1)
        V = np.hstack([A @ b % PRIME for A in Or])
        mu = b.T @ V % PRIME
        mu = mu.T * pow(2,-1,PRIME) % PRIME
        co = Gi @ mu % PRIME
        theta = sum((int(c)*A for c,A in zip(co[:,0],reps)),np.zeros((16,16),dtype=np.int64)) % PRIME
        U = V @ Gi % PRIME
        H = (theta + Oi @ U % PRIME @ V.T) % PRIME
        mats = [theta,H]
        table = {():np.eye(16,dtype=np.int64)}
        for n in range(1,6):
            for word in itertools.product(range(2),repeat=n):
                table[word] = table[word[:-1]] @ mats[word[-1]] % PRIME
        def trace_word(word):
            word=tuple(word);k=len(word)//2
            return int(np.sum(table[word[:k]]*table[word[k:]].T,dtype=np.int64)%PRIME)
        gradients = {}
        for index in (1,2,3):
            word=tuple(w4[index]);K=[np.zeros((16,16),dtype=np.int64) for _ in range(2)]
            for r,letter in enumerate(word):
                K[letter]=(K[letter]+table[word[r+1:]+word[:r]])%PRIME
            t=np.array([int(np.sum((K[0]+K[1])*A.T,dtype=np.int64)%PRIME) for A in reps],dtype=np.int64).reshape(10,1)
            grad=V @ Gi % PRIME @ t % PRIME
            W=((Oi.T @ K[1].T+K[1] @ Oi)%PRIME) @ U % PRIME
            grad=(grad+sum((A.T @ W[:,j:j+1]%PRIME for j,A in enumerate(Or)),np.zeros((16,1),dtype=np.int64)))%PRIME
            gradients[index]=grad
            if sample_index==0:
                dc=Gi @ V.T % PRIME;direct=[]
                for k in range(16):
                    dt=sum((int(dc[j,k])*A for j,A in enumerate(reps)),np.zeros((16,16),dtype=np.int64))%PRIME
                    dV=np.hstack([A[:,k:k+1] for A in Or])
                    dh=(dt+Oi @ ((dV @ U.T+U @ dV.T)%PRIME))%PRIME
                    direct.append(int((np.sum(K[0]*dt.T,dtype=np.int64)+np.sum(K[1]*dh.T,dtype=np.int64))%PRIME))
                assert list(map(int,grad[:,0]))==direct
                direct_gradient_check=True
        values8=[trace_word(w4[i]) for i in basis8]
        values14=[int((gradients[i].T @ Oi % PRIME @ gradients[j])[0,0]%PRIME) for i,j in bracket_pairs]
        values18=[trace_word(w) for w in words18]
        F4=int((mu.T @ co)[0,0]%PRIME)
        evaluations.append([F4*z%PRIME for z in values18]+[values8[i]*values14[j]%PRIME for i,j in pairs8_14])
        factor_values.append(dict(F4=F4,Inv8=values8,Inv14=values14,Inv18=values18))
        betas.append(beta)
    assert len(evaluations)==36
    rr,cc,_=echelon(evaluations)
    assert len(rr)<=21
    inputs=[Path(__file__),ROOT/'code/q_new_family_degree16_invariant_probe_20261103_v63.py',ROOT/'code/weight1_structure.py',invpath,p8,p14,p18,p36]
    result=dict(id=STEM,all_passed=True,prime=PRIME,full_character_dimension22=21,
                character_records=chars,definitions=definitions,candidate_count=24,
                existing_point_count=36,new_beta_points_generated=0,
                same_accepted_QQ_samples=samples['actual_QQ_source_points'],actual_beta16_modp=betas,
                factor_values=factor_values,actual_evaluation=evaluations,
                actual_restriction_rank=len(rr),selected_rows=rr,selected_columns=cc,
                direct_dtheta_dH_gradient_check_at_same_first_point=direct_gradient_check,
                complete_Inv22_span_and_injective_restriction_if_rank21=len(rr)==21,
                no_exact_QQ_minor_in_probe=True,no_degree24_probe=True,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','candidate_count','existing_point_count','actual_restriction_rank','seconds']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output.resolve())
