#!/usr/bin/env python3
"""Certify a nonzero global affine annihilator, without inverting detK.

The global polynomial is specified by a finite determinant/resultant circuit.
Its nonzero radial specialization uses the saved exact cyclic certificate;
the whole billion-degree coefficient expansion is neither needed nor stored.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate(coefficients,value):
    answer=0
    for c in reversed(coefficients):answer=(answer*value+c)%P
    return answer


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--timeout',type=int,default=180)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    cyclic=ROOT/'reports/S61_AFFINE_CYCLIC_OPEN_20261103_v62'
    previous=json.loads((cyclic/'report.json').read_text());assert previous['all_passed']
    for name,value in previous['input_sha256'].items():assert sha(ROOT/name)==value
    for name,value in previous['evidence_sha256'].items():assert sha(cyclic/name)==value
    certificate=cyclic/'native-certificate.txt';lines=certificate.read_text().splitlines()
    chi=list(map(int,next(line for line in lines if line.startswith('CHI ')).split()[1:]))
    pt=list(map(int,next(line for line in lines if line.startswith('COORDINATE 0 ')).split()[2:]))
    g=1031;assert len(chi)==g+1 and len(pt)==g and chi[-1]==1 and pt[-1]!=0
    path=args.output/'polynomials.in';path.write_text('1031\n'+' '.join(map(str,chi))+'\n'+' '.join(map(str,pt))+'\n')
    source=ROOT/'code/s61_affine_global_resultant_annihilator_20261103_v62.cpp'
    executable=args.output/'annihilator'
    flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],capture_output=True,text=True,check=True).stdout)
    subprocess.run(['c++','-std=c++17','-O3',str(source),*flags,'-o',str(executable)],check=True)
    native=args.output/'native.txt';log=args.output/'native.log'
    with log.open('w')as stream:subprocess.run([str(executable.resolve()),str(path.resolve()),str(native.resolve())],stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=args.timeout)
    executable.unlink();out=native.read_text().splitlines()
    tokens=out[0].split();summary={k:int(v)for k,v in zip(tokens[::2],tokens[1::2])}
    assert summary==dict(G=1031,Q_DEGREE=1031,DISTINCT_RESULTANT_CHECKS=1032)
    q=list(map(int,next(line for line in out if line.startswith('Q ')).split()[1:]))
    assert len(q)==1032 and q[0]==P-1
    multiplicity=int(next(line for line in out if line.startswith('MULTIPLICITY_AT_ONE ')).split()[1]);assert multiplicity==5
    values={int(line.split()[1]):int(line.split()[2])for line in out if line.startswith('EVALUATION ')}
    assert all(evaluate(q,x)==y for x,y in values.items());assert values[2]==30737
    # Verify the lambda=2 unit identity independently using ordinary integer
    # convolution. No inverse, characteristic-polynomial, or resultant routine.
    u=np.array(list(map(int,next(line for line in out if line.startswith('BEZOUT_CHI ')).split()[1:])),dtype=np.int64)
    w=np.array(list(map(int,next(line for line in out if line.startswith('BEZOUT_2PT_MINUS_ONE ')).split()[1:])),dtype=np.int64)
    second=2*np.array(pt,dtype=np.int64);second[0]-=1;second%=P
    a=np.convolve(u,np.array(chi,dtype=np.int64));b=np.convolve(w,second)
    check=np.zeros(max(len(a),len(b)),dtype=np.int64);check[:len(a)]+=a;check[:len(b)]+=b
    expected=np.zeros_like(check);expected[0]=1;assert np.array_equal(check%P,expected)
    current=q.copy();count=0
    while len(current)>1 and evaluate(current,1)==0:
        quotient=[0]*(len(current)-1);quotient[-1]=current[-1]
        for j in range(len(quotient)-2,-1,-1):quotient[j]=(current[j+1]+quotient[j+1])%P
        assert (current[0]+quotient[0])%P==0;current=quotient;count+=1
    assert count==5
    delta=sum(range(g))-(29+2*301+3*700);assert delta==528234
    det_at_v=previous['krylov_determinant_mod_prime'];assert det_at_v==17765
    det_at_two=pow(2,delta,P)*det_at_v%P
    global_value=pow(det_at_two,2*g-1,P)*values[2]%P
    assert det_at_two==12800 and global_value==27775
    definition=dict(base='B=R_beta[Y1,...,Y17]; exact corrected global h-action as in spread note',
                    K='[e,T_h e,...,T_h^1030 e]',Delta='det K',
                    Phat_t='sum_j (adj(K))_(j,e_t) s^j',
                    chihat='Delta*s^1031-sum_j (adj(K)*T_h^1031 e)_j*s^j',
                    global_annihilator='Res_s(chihat,Phat_t-Delta), fixed degrees1031 and1030',
                    exact_global_relations=['chihat(h)=0','Phat_t(h)=Delta*t'],
                    whole_affine_consequence='In A/(t-1), the global annihilator is zero, including Delta=0.',
                    radial_norm_coefficients=q,radial_delta_degree=delta,
                    radial_norm_definition='q(lambda)=Res_s(chi_v(s),lambda*P_t,v(s)-1)',
                    independent_lambda2_bezout_coefficients=dict(chi=u.tolist(),two_pt_minus_one=w.tolist()))
    (args.output/'annihilator-circuit.json').write_text(json.dumps(definition,indent=2)+'\n')
    spread=ROOT/'research/notes/S61_AFFINE_GLOBAL_ACTION_SPREAD_20261103_v62.md'
    inputs=[Path(__file__),source,certificate,cyclic/'report.json',spread]
    report=dict(id='S61_AFFINE_GLOBAL_RESULTANT_ANNIHILATOR_20261103_v62',all_passed=True,prime=P,
                all17_bases_preserved=True,global_affine_annihilator_without_delta_localization=True,
                global_polynomial_expanded=False,global_polynomial_representation='explicit determinant/resultant circuit',
                symbol_count=g,radial_delta_degree=delta,sylvester_size=2*g-1,
                radial_norm_degree=1031,radial_norm_at_one_multiplicity=5,
                all1032_distinct_resultant_checks=True,independent_lambda2_bezout_convolution=True,
                radial_norm_at_two=values[2],delta_at_two_v=det_at_two,global_annihilator_at_two_v=global_value,
                maximum_total_base_degree=delta*(2*g-1)+g,
                generic_QQ_nonzero=True,unresolved_annihilator_zero_locus=True,
                finite_pair_generic_rank_upper29_claimed=False,full_affine_exhaustion_claimed=False,
                limitations=['Explicit coefficient specialization is frozen-beta F32003. GenericQQ uses the p-integral corrected map and its unchanged nonzero specialization.',
                             'The resultant vanishes on the full affine algebra, but its zero locus is not exhausted.',
                             'Its very high degree makes this a concrete support relation, not yet an effective rank29 bound.'],
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',delta,det_at_two,values[2],global_value,report['seconds'])


if __name__=='__main__':main()
