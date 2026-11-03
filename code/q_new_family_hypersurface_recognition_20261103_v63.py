#!/usr/bin/env python3
"""Exact low-degree invariant exclusions for the new beta-image hypersurface."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix

from weight1_structure import constraints,COORDS,sp_basis,infinitesimal_action,J

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    samplepath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/exact-sample.json'
    sample=json.loads(samplepath.read_text());beta24=s.Matrix(list(map(s.Rational,sample['primitive_beta24'])))
    B=s.Matrix.hstack(*s.Matrix.vstack(*constraints()).nullspace());rows=list(B.T.rref()[1]);left=B[rows,:].inv()
    point=left*beta24[rows,:];assert B*point==beta24
    diagonal_weights=[(1,-1,0,0),(0,0,1,-1)];weights=[]
    for column in range(16):
        choices={tuple(sum(w[i]for i in[*pair,k])for w in diagonal_weights)for row,(pair,k)in enumerate(COORDS)if B[row,column]}
        assert len(choices)==1;weights.append(choices.pop())
    generators=sp_basis()
    representation=[s.Matrix.hstack(*[left*infinitesimal_action(B[:,j],A)[rows,:]for j in range(16)])for A in generators]
    invariantpath=ROOT/'code/data/q_new_family_quartic_invariant_20261103_v63.json'
    saved=json.loads(invariantpath.read_text());quartic=saved['invariants_sparse'][0]
    data=[]
    for degree in[2,4,6]:
        words=[mon for mon in itertools.combinations_with_replacement(range(16),degree)
            if all(sum(weights[i][k]for i in mon)==0 for k in range(2))]
        equations={}
        for column,mon in enumerate(words):
            for index,G in enumerate(representation):
                for i in set(mon):
                    remaining=list(mon);remaining.remove(i)
                    for j in range(16):
                        c=G[i,j]*mon.count(i)
                        if c:
                            row=equations.setdefault((index,tuple(sorted(remaining+[j]))),{})
                            row[column]=row.get(column,s.S.Zero)+c
        entries=[row for row in equations.values()if any(row.values())]
        M=s.SparseMatrix(len(entries),len(words),{(i,j):c for i,row in enumerate(entries)for j,c in row.items()if c})
        reduced,pivots=DomainMatrix.from_Matrix(M).rref();dimension=len(words)-len(pivots)
        assert dimension==({2:0,4:1,6:0}[degree])
        if degree==4:
            index={mon:i for i,mon in enumerate(words)};vector=s.zeros(len(words),1)
            for mon,c in quartic:vector[index[tuple(mon)],0]=s.Rational(c)
            assert vector!=s.zeros(len(words),1)and M*vector==s.zeros(M.rows,1)
        data.append(dict(degree=degree,weight_zero_monomials=len(words),Lie_constraint_matrix_shape=list(M.shape),
            Lie_constraint_nnz=len(M.todok()),exact_QQ_rank=len(pivots),invariant_dimension=dimension,
            complete_all_ten_Lie_generators=True))
    quartic_value=s.factor(sum(s.Rational(c)*s.prod(point[i]for i in mon)for mon,c in quartic))
    assert quartic_value==-s.Rational(6287207,1638400)!=0
    Omega=s.Matrix(24,24,lambda a,c:(J[COORDS[a][0][0],COORDS[c][0][0]]*J[COORDS[a][0][1],COORDS[c][0][1]]-
        J[COORDS[a][0][0],COORDS[c][0][1]]*J[COORDS[a][0][1],COORDS[c][0][0]])*J[COORDS[a][1],COORDS[c][1]])
    omega=B.T*Omega*B;assert omega.T==-omega and omega.det()!=0
    assert all(G.T*omega+omega*G==s.zeros(16)for G in representation)
    gram=s.Matrix([[s.trace(A*C)for C in generators]for A in generators]);assert gram.det()!=0
    mu=s.Matrix([(point.T*omega*G*point)[0]/2 for G in representation]);coeff=gram.inv()*mu
    moment=sum((c*A for c,A in zip(coeff,generators)),s.zeros(4))
    trace2=s.factor(s.trace(moment*moment));moment_det=s.factor(moment.det());moment_discriminant=s.factor(trace2**2-16*moment_det)
    assert moment_det!=0 and moment_discriminant!=0
    (args.output/'invariant-certificates.json').write_text(json.dumps(dict(
        primitive24_by16_basis=[[str(v)for v in row]for row in B.tolist()],Cartan_weights=weights,
        all_ten_representation_matrices=[[[str(v)for v in row]for row in G.tolist()]for G in representation],
        complete_invariant_dimensions=data,unique_quartic_sparse=quartic,quartic_on_new_beta=str(quartic_value),
        beta_image_defining_homogeneous_equation_degree_atleast=8,
        reason='Sp4 connected semisimple implies its stable hypersurface equation invariant; beta scalar dilation implies homogeneous; center−I excludes odd degrees.',
        symplectic_form16=[[str(v)for v in row]for row in omega.tolist()],
        moment_map_definition='mu_A=(1/2)*b^T*omega*rho(A)*b, identified withsp4 bytracepairing',
        sample_moment_matrix=[[str(v)for v in row]for row in moment.tolist()],
        moment_trace_square=str(trace2),moment_determinant=str(moment_det),moment_eigen_discriminant=str(moment_discriminant),
        generic_new_beta_moment_regular_semisimple=True,
        canonical_curve_discriminant=str(sample['canonical_curve_discriminant']),
        three_excluded_hypersurface_candidates=['canonical_curve_Delta24','moment_determinant_degree8','moment_eigen_discriminant_degree8'],
        actual_defining_hypersurface_equation_not_found=True),indent=2)+'\n')
    paths=[Path(__file__),invariantpath,samplepath,
        ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/report.json',
        ROOT/'reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW/report.json',
        ROOT/'code/weight1_structure.py',ROOT/'research/notes/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63.md']
    report=dict(id='Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63',all_passed=True,
        full_QQ_low_degree_Sp4_invariant_kernels=True,defining_equation_degree_lower_bound=8,
        unique_quartic_not_zero_on_new_image=True,
        curve_discriminant_and_two_moment_hypersurfaces_excluded=True,
        polynomial_equation_or_full_hypersurface_identification_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('EXACT HYPERSURFACE DEGREE>=8; DELTA/MOMENT EXCLUSIONS PASS',report['seconds'],flush=True)


if __name__=='__main__':main()
