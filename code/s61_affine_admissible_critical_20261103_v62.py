#!/usr/bin/env python3
"""A toric-centered admissible five-slice and global critical determinant."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring

from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS, infinitesimal_action, sp_basis
from s61_affine_fiber_macaulay import load_fiber
from s61_nonsquare_cubic_row_basis_review_20261103_v62 import determinant_mod

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def residue(value):
    value=sp.Rational(value);assert int(value.q)%P
    return int(value.p)*pow(int(value.q),-1,P)%P


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    original=ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json'
    data=json.loads(original.read_text());point_path=ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'
    record=next(item for item in json.loads(point_path.read_text())['points']if item['component_id']=='T01')
    first=[i for i,w in enumerate(data['parameter_weights'])if w==1]
    second=[i for i,w in enumerate(data['parameter_weights'])if w==2]
    center=sp.Matrix([sp.Rational(record['S61_coordinates_u1_through_u61'][i])for i in first])
    chosen_old=[35,12,26,4,24] # originalu36,u13,u27,u5,u25; all indices here0-based
    chosen=[first.index(i)for i in chosen_old]
    model=S61();z,u,beta=primitive_beta(model)
    old=sp.Matrix([[sp.expand(beta[i,j]).coeff(z[k]).coeff(v)for v in u]for(i,j),k in COORDS])
    rows=list(old.T.rref()[1]);left=old[rows,:].inv();primitive_point=old*center
    orbit=sp.Matrix.hstack(*[left*infinitesimal_action(primitive_point,A)[rows,:]for A in sp_basis()],center)
    dominance=orbit.row_join(sp.eye(16)[:,chosen]);det=sp.factor(dominance.det())
    assert det==sp.Rational(3564223982,1029)and residue(det)==9752
    coefficient_ring,*ss=ring('ss1,ss2,ss3,ss4,ss5',QQ)
    forms=[coefficient_ring(QQ(int(c.p),int(c.q)))for c in center]
    for j,i in enumerate(chosen):forms[i]+=ss[j]
    fi={i:j for j,i in enumerate(first)};qi={i:j for j,i in enumerate(second)}
    cache={():coefficient_ring.one};encoded=[];equations=[]
    for equation in data['equations']:
        row={}
        for mon,n,d in equation:
            betamon=tuple(i for i in mon if i in fi);qmon=tuple(qi[i]for i in mon if i in qi)
            if betamon not in cache:
                value=coefficient_ring.one
                for i in betamon:value*=forms[fi[i]]
                cache[betamon]=value
            row[qmon]=row.get(qmon,coefficient_ring.zero)+QQ(n,d)*cache[betamon]
        row={mon:c for mon,c in row.items()if c};assert all(len(mon)<=2 for mon in row)
        encoded.append([[list(mon),[[list(exponent),str(c)]for exponent,c in coefficient.items()]]for mon,coefficient in row.items()])
        equations.append('+'.join('('+str(coefficient).replace('**','^')+')'+''.join('*q'+str(i)for i in mon)for mon,coefficient in row.items())or'0')
    assert len(encoded)==134
    with (args.output/'restricted-original134.json.gz').open('wb')as stream:
        stream.write(gzip.compress(json.dumps(encoded,separators=(',',':')).encode(),mtime=0))
    (args.output/'admissible-slice.json').write_text(json.dumps(dict(
        original_weight1_indices_zero_based=first,center=list(map(str,center)),varying_old_u_one_based=[i+1 for i in chosen_old],
        dominance_matrix=[[str(c)for c in row]for row in dominance.tolist()],
        exact_determinant=str(det),determinant_mod32003=9752,
        existing_finite_model_frame_regular_at_center=True,
        fixed_pair_center_finiteness_on_slice_claimed=False),indent=2)+'\n')
    quintic=ROOT/'reports/S61_AFFINE_QUINTIC_OPEN_20261103_v62'
    report=json.loads((quintic/'report.json').read_text());assert report['all_passed']
    for name,value in report['input_sha256'].items():assert sha(ROOT/name)==value
    for name,value in report['evidence_sha256'].items():assert sha(quintic/name)==value
    jac=json.loads((quintic/'jacobian-certificate.json').read_text());C=np.array(jac['common_combination'],dtype=np.int64)
    F,qids,base,*_=load_fiber();assert qids==second
    five_path=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json';five=json.loads(five_path.read_text());vs=sp.symbols('v0:3')
    slopes=np.array([[residue(sp.sympify(expression).coeff(v))for v in vs]for expression in five['affine_plane_substitution']],dtype=np.int64)
    known_path=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    known=json.loads(known_path.read_text())['known_five_mod_32003'];determinants=[]
    for point in known:
        delta=slopes@np.array(point[1:],dtype=np.int64)%P;J=np.zeros((134,45),dtype=np.int64)
        for r,polynomial in enumerate(F):
            evaluation=0
            for mon,c in polynomial.items():
                co=c.numerator*pow(c.denominator,-1,P)%P
                value=co
                for i in mon:value=value*int(delta[i])%P
                evaluation=(evaluation+value)%P
                if len(mon)==1:J[r,mon[0]]=(J[r,mon[0]]+co)%P
                elif len(mon)==2:
                    a,b=mon;J[r,a]=(J[r,a]+co*delta[b])%P;J[r,b]=(J[r,b]+co*delta[a])%P
            assert evaluation==0
        determinants.append(determinant_mod((C@J%P).tolist(),P))
    assert determinants==[18384,6387,6706,13504,28620]
    (args.output/'critical-determinant.json').write_text(json.dumps(dict(
        definition='E(beta,q)=det(C*Jac_q(original134equations))',
        original_q_indices_zero_based=second,constant_row_combination=C.tolist(),
        maximum_q_degree=45,known_five_determinants_mod32003=determinants,
        all_positive_dimensional_reduced_components_contained_in_E_zero=True,
        determinant_polynomial_expanded=False,
        exact_5parameter_unit_target='(original134_on_new_slice,E)=(1) over QQ(ss1,...,ss5)',
        unit_is_sufficient_for_dimension_zero=True,unit_equivalent_to_dimension_zero=False,
        possible_remaining_points='Unknown isolated or nonreduced points can remain; E nonzero at the knownfive only.'),indent=2)+'\n')
    # Exact runnable model definition; no determinant/GB execution by this producer.
    literal='{{'+ '},{'.join(','.join(map(str,row))for row in C.tolist())+'}}'
    m2='A=QQ[ss1,ss2,ss3,ss4,ss5];K=frac A;R=K['+','.join(f'q{i}'for i in range(45))+'];\n'
    m2+='I=ideal('+','.join(equations)+');\nCC=matrix(R,'+literal+');\n'
    m2+='E=det(CC*transpose(jacobian gens I));\n-- Unexecuted exact target: gb(I+ideal(E));\n'
    (args.output/'critical-five-parameter-model.m2').write_text(m2)
    inputs=[Path(__file__),original,point_path,five_path,known_path,quintic/'report.json',
        ROOT/'code/s61_nonsquare_cubic_row_basis_review_20261103_v62.py',
        ROOT/'research/notes/S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62.md']
    answer=dict(id='S61_AFFINE_ADMISSIBLE_CRITICAL_20261103_v62',all_passed=True,
        transcendental_parameters=5,all_original_equations_retained=134,
        exact_dominance_determinant=str(det),dominance_determinant_mod32003=9752,
        finite_model_frame_nonempty_open_on_new_slice=True,
        critical_determinant_maximum_q_degree=45,
        known_five_critical_determinants_mod32003=determinants,
        global_positive_support_exclusion_proved=True,full_generic_affine_dimension_zero_claimed=False,
        critical_unit_or_GB_computed=False,unknown_isolated_points_excluded=False,
        fixed_pair_center_finiteness_on_slice_claimed=False,
        limitations=['The critical determinant cuts out a containing locus of all positive-dimensional components; it need not vanish only there.',
                    'The new five-parameter unit target is sufficient, not equivalent, to dimension0; singular isolated points may obstruct it.',
                    'No previously failed fieldGB is rerun. The determinant is saved as an exact circuit.'],
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(answer,indent=2)+'\n')
    print('PASS toric-centeredadmissible5slice;globaldegree45criticallocus;unitUNSOLVED;',answer['seconds'])


if __name__=='__main__':main()
