#!/usr/bin/env python3
"""Exact QQ selected20 minor after one bounded modular invariant probe."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import shlex
import subprocess
import time
import sympy as s
from q_new_family_degree16_invariant_probe_20261103_v63 import ROOT,PRIME,read,sha,matrix,mod
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import sp_basis


def main(output,probe):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    data=read(probe/'report.json');assert data['all_passed']
    assert data['ambient_modular_rank']==data['actual_new_map_modular_rank']==20
    for rel,digest in data['input_sha256'].items():assert sha(ROOT/rel)==digest
    invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    d8path=ROOT/'reports/Q_NEW_FAMILY_DEGREE8_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
    mappath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz'
    info,d8,payload=map(read,[invpath,d8path,mappath])
    B=matrix(info['primitive24_by16_basis']);br=list(B.T.rref()[1]);left=B[br,:].inv()
    O=matrix(info['symplectic_form16']);rep=list(map(matrix,info['all_ten_representation_matrices']))
    gs=sp_basis();gram=s.Matrix([[s.trace(A*C) for C in gs] for A in gs])
    model=S61();z,u,beta=primitive_beta(model);ones=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    import itertools
    coeff=[beta[i,j].coeff(t) for i,j in itertools.combinations(range(4),2) for t in z]
    betaU=s.Matrix([[s.expand(c).coeff(v) for v in u] for c in coeff])
    rows=data['actual_independent_rows'];cols=data['actual_independent_columns'];assert len(rows)==len(cols)==20
    points=[];betas=[]
    for row in rows:
        point=data['actual_points'][row];values=[]
        for node in payload['nodes']:
            op=node[0]
            if op=='q':v=Fraction(node[1])
            elif op=='v':v=Fraction(point[payload['input_names'][node[1]]])
            elif op=='+':v=sum((values[i] for i in node[1]),Fraction(0))
            elif op=='*':
                v=Fraction(1)
                for i in node[1]:v*=values[i]
            elif op=='inv':v=1/values[node[1]]
            else:raise AssertionError(op)
            values.append(v)
        b24=betaU*s.Matrix([s.Rational(values[payload['output_u1_through_u61'][i]]) for i in ones])
        b=left*b24[br,:];assert B*b==b24
        assert [mod(c) for c in b]==data['actual_beta16'][row]
        pivot=values[payload['normalization']['value_pivot_node']]
        pf=values[payload['normalization']['value_pfaffian_node']];assert pivot and pf
        points.append(dict(point=point,value_pivot=str(pivot),value_pfaffian=str(pf),beta16=list(map(str,b))))
        betas.append(b)
    native=ROOT/'code/q_new_family_degree16_invariant_exact_20261103_v63.cpp'
    nativeinput=output/'native-input.txt'
    with nativeinput.open('w') as f:
        def write(M):
            f.write(f'{M.rows} {M.cols}\n')
            for row in M.tolist():f.write(' '.join(map(str,row))+'\n')
        f.write('20\n')
        for M in [O,O.inv(),gram.inv(),*rep]:write(M)
        f.write(' '.join(map(str,d8['independent_five_columns']))+'\n')
        definitions=data['invariant_definitions'];f.write(str(len(definitions))+'\n')
        for definition in definitions:
            kind=definition['kind']
            if kind=='basis8_product':f.write('P '+' '.join(map(str,definition['indices']))+'\n')
            else:f.write(('T ' if kind=='trace8' else 'F ')+''.join(map(str,definition['word']))+'\n')
        f.write(' '.join(map(str,cols))+'\n')
        for b in betas:write(b)
    binary=ROOT/'work/new52-degree16-v63/exact-invariants';binary.parent.mkdir(parents=True,exist_ok=True)
    flags=shlex.split(subprocess.run(['pkg-config','--cflags','--libs','flint'],text=True,capture_output=True,check=True).stdout)
    subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True)
    evidence=output/'exact-evaluation.json'
    with (output/'native.log').open('w') as f:subprocess.run([str(binary),str(nativeinput),str(evidence)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=90)
    exact=read(evidence);assert exact['all_passed']
    assert len(exact['evaluation'])==20
    for row,native_row in zip(rows,exact['evaluation']):
        assert [mod(c) for c in native_row]==data['actual_evaluation'][row]
    assert exact['exact_selected20_minor']!='0'
    (output/'samples.json').write_text(json.dumps(dict(selected_probe_rows=rows,selected_columns=cols,
        invariant_definitions=definitions,actual_QQ_points_and_betas=points),indent=2)+'\n')
    oldreview=ROOT/'reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json'
    accepted=read(oldreview);assert accepted['all_passed'] and accepted['theta_and_H_same_conjugation_and_constant_Poisson_bracket_invariance_proved']
    inputs=[Path(__file__),native,probe/'report.json',invpath,d8path,mappath,oldreview,
            ROOT/'research/notes/Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63.md']
    report=dict(id='Q_NEW_FAMILY_DEGREE16_INVARIANT_EXCLUSION_20261103_v63',all_passed=True,
        full_invariant_dimension16=20,full_invariant_dimension17=0,full_invariant_dimension18=9,
        actual_new_image_QQ_restriction_rank16=20,
        exact20_minor_nonzero=True,exact20_minor_in='exact-evaluation.json',
        all60_QQ_invariants_at20_actual_QQ_points_exactly_evaluated=True,
        all_exact_values_reduce_to_modular_probe=True,
        universal_theta_H_conjugation_proof_reused_from_independent_review=True,
        complete_degree16_invariant_restriction_injective=True,
        defining_homogeneous_invariant_degree_lower_bound=18,
        actual_equation_or_degree18_exclusion_or_general_fiber_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},
        evidence_sha256={p.name:sha(p) for p in [nativeinput,evidence,output/'samples.json',output/'native.log']},
        seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','actual_new_image_QQ_restriction_rank16','defining_homogeneous_invariant_degree_lower_bound','seconds']},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--probe',type=Path,required=True);args=parser.parse_args();main(args.output.resolve(),args.probe.resolve())
