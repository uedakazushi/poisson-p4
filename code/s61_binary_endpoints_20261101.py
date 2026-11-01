#!/usr/bin/env python3
"""Two dominant sixteen-dimensional binary endpoint strata are locally removed."""
import argparse
import hashlib
import json
import time
from pathlib import Path
from itertools import combinations
import sympy as s
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import J
from s61_binary_endpoint_engine_20261101 import verify_endpoint


def verify():
    start=time.time();root=Path(__file__).resolve().parents[1];m=S61()
    family_path=root/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json'
    rank_path=root/'reports/S61_BINARY_UNIFORM_RANK_20261101.json'
    family=json.loads(family_path.read_text());rank=json.loads(rank_path.read_text())
    z,uu,beta=primitive_beta(m);actions=s.Matrix(family['primitive_action']).applyfunc(s.sympify)
    stabilizer=s.Matrix.hstack(*[s.Matrix(v).applyfunc(s.sympify) for v in rank['primitive_actions']])
    cases=[]
    for extra in ['u14','u49']:
        result=verify_endpoint(extra)
        constraints=[s.Symbol(v) for v in ['u12','u13','u35','u50',extra]];sub={v:0 for v in constraints}
        rows=[uu.index(v) for v in constraints]
        assert stabilizer[rows,:].subs(sub)==s.zeros(5,5)
        sample={s.Symbol(k):s.Rational(v) for k,v in result['sample'].items()}
        normal=actions[rows,:].subs(sample);columns=list(normal.to_DM().rref()[1]);assert len(columns)==5
        det=normal[:,columns].det();assert det
        result.update(beta_constraints=list(map(str,constraints)),beta_normal_action=[[str(v) for v in row] for row in normal.tolist()],
            beta_dominance_columns=columns,beta_dominance_minor=str(det),projective_stabilizer_dimension=5,
            stabilizer_preserves_slice=True,joint_image_dimension=16,first_jet_image_dimension=16,
            generic_fixed_beta_stratum_removed=True)
        cases.append(result)
    # The alternative codimension-two rank-drop slice is the same saturated
    # stratum, by a finite symplectic transformation of beta and Q.
    S=s.Matrix([[0,0,0,1],[0,0,-1,0],[0,-1,0,0],[1,0,0,0]])
    assert S.T*J*S==J and S*S==s.eye(4)
    V=s.Matrix(z);Q=s.zeros(4);Q[1,2]=z[0]*z[3];Q[2,1]=-Q[1,2]
    newQ=S*Q.subs(dict(zip(z,S*V)),simultaneous=True)*S.T
    assert (newQ+Q).applyfunc(s.expand)==s.zeros(4)
    def vector(M):return s.Matrix([s.expand(M[i,j]).coeff(z[k]) for i,j in combinations(range(4),2) for k in range(4)])
    B=vector(beta).jacobian(uu);rows=list(B.T.to_DM().rref()[1]);changed=S.T*beta.subs(dict(zip(z,S*V)),simultaneous=True)*S
    new=(B[rows,:].inv()*vector(changed)[rows,:]).applyfunc(s.expand)
    expected=[s.Symbol('u13'),s.Symbol('u12'),s.Symbol('u14')-s.Symbol('u49')/2,3*s.Symbol('u49')/2]
    assert [new[uu.index(s.Symbol(v))] for v in ['u12','u13','u35','u50']]==expected
    paths=[Path(__file__),root/'code/s61_binary_endpoint_engine_20261101.py',family_path,rank_path,
        root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py',root/'code/primitive_jet_elliptic_curve_20261101.py',root/'code/weight1_structure.py',
        root/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',root/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',root/'research/notes/S61_BINARY_UNIFORM_CUBIC_20261101.md',root/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_BINARY_ENDPOINTS_20261101',all_passed=True,prime=32003,cases=cases,
        symmetry_matrix=[[str(v) for v in row] for row in S.tolist()],symmetry_constraint_images=list(map(str,expected)),
        proof_note='research/notes/S61_BINARY_ENDPOINTS_20261101.md',
        limitations=['The two endpoint strata only; the u4=0 endpoint is not included.',
                    'No full rank-drop-locus or full affine-fiber exhaustion.'],
        input_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);args=p.parse_args();result=verify()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,cases=[{k:r[k] for k in ['extra_zero','beta_dominance_minor','quintic_initial_minor_mod_prime','relation_minor_mod_prime','augmented_homogenizer_minor_mod_prime']} for r in result['cases']],seconds=result['seconds'])))


if __name__=='__main__':main()
