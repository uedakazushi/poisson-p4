#!/usr/bin/env python3
"""Global degree45 known-five idempotent via exact quadratic polarization.

The generic determinant/adjugate/trace is a finite recipe, not an expanded
coefficient table.  The unknown residual is retained scheme-theoretically.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp

from s61_affine_fiber_macaulay import load_fiber

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def add(poly,mon,value):
    value=poly.get(mon,Fraction(0))+value
    if value:poly[mon]=value
    else:poly.pop(mon,None)


def residue(value):
    value=Fraction(value);assert value.denominator%P
    return value.numerator*pow(value.denominator,-1,P)%P


def determinant_including_zero(matrix):
    matrix=np.array(matrix,dtype=np.int64).copy()%P;answer=1
    for column in range(len(matrix)):
        pivots=np.flatnonzero(matrix[column:,column])
        if not len(pivots):return 0
        pivot=column+int(pivots[0])
        if pivot!=column:matrix[[pivot,column]]=matrix[[column,pivot]];answer=-answer
        value=int(matrix[column,column]);answer=answer*value%P
        matrix[column]=matrix[column]*pow(value,-1,P)%P
        matrix[column+1:]=(matrix[column+1:]-matrix[column+1:,column,None]*matrix[column,None,:])%P
    return answer%P


def polarization(data):
    qi={old:new for new,old in enumerate(i for i,w in enumerate(data['parameter_weights'])if w==2)}
    statistics=[]
    for row in data['equations']:
        left={};right={};degrees=Counter()
        for mon,n,d in row:
            beta=tuple(i for i in mon if i not in qi);qmon=tuple(qi[i]for i in mon if i in qi)
            assert len(qmon)<=2;degrees[len(qmon)]+=1;c=Fraction(n,d);assert c.denominator%P
            # Keys retain ALL sixteen beta variables as original indices.
            qm=tuple((0,i)for i in qmon);pm=tuple((1,i)for i in qmon)
            add(left,(beta,qm),c);add(left,(beta,pm),-c)
            for occurrence,index in enumerate(qmon):
                others=qmon[:occurrence]+qmon[occurrence+1:]
                derivative=[((),c)]if not others else[((0,others[0]),Fraction(1,2)*c),((1,others[0]),Fraction(1,2)*c)]
                for factor,coefficient in derivative:
                    prefix=()if not factor else(factor,)
                    add(right,(beta,tuple(sorted(prefix+((0,index),)))),coefficient)
                    add(right,(beta,tuple(sorted(prefix+((1,index),)))),-coefficient)
        assert left==right
        statistics.append(dict(original_term_counts_by_q_degree=dict(sorted(degrees.items())),
            exact_polarized_nonzero_terms=len(left),all_QQ_coefficients_match=True))
    assert len(statistics)==134
    return statistics


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    original=ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json'
    data=json.loads(original.read_text());statistics=polarization(data)
    specpath=ROOT/'code/data/s61_affine_linear_block_20261103_v63.json'
    spec=json.loads(specpath.read_text());rows=spec['allowed_rows']+spec['critical_selected_remaining_rows']
    assert len(rows)==len(set(rows))==45
    F,qids,base,*_=load_fiber();assert len(qids)==45
    fivepath=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json';five=json.loads(fivepath.read_text())
    vs=sp.symbols('v0:3');slopes=np.array([[residue(sp.sympify(expr).coeff(v))for v in vs]for expr in five['affine_plane_substitution']],dtype=np.int64)
    knownpath=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    known=json.loads(knownpath.read_text())['known_five_mod_32003']
    deltas=[slopes@np.array(point[1:],dtype=np.int64)%P for point in known]
    for delta in deltas:
        for polynomial in F:
            ev=0
            for mon,c in polynomial.items():
                term=residue(c)
                for i in mon:term=term*int(delta[i])%P
                ev=(ev+term)%P
            assert ev==0
    determinants=[];midpoint_kernel_checks=0
    for i in range(5):
        out=[]
        for j in range(5):
            middle=(deltas[i]+deltas[j])*pow(2,-1,P)%P;J=np.zeros((134,45),dtype=np.int64)
            for r,polynomial in enumerate(F):
                for mon,c in polynomial.items():
                    co=residue(c)
                    if len(mon)==1:J[r,mon[0]]=(J[r,mon[0]]+co)%P
                    elif len(mon)==2:
                        a,b=mon;J[r,a]=(J[r,a]+co*middle[b])%P;J[r,b]=(J[r,b]+co*middle[a])%P
            assert not np.any(J@((deltas[j]-deltas[i])%P)%P)
            midpoint_kernel_checks+=1;out.append(determinant_including_zero(J[rows,:]))
        determinants.append(out)
    diagonal=[determinants[i][i]for i in range(5)]
    assert all(diagonal)and all(determinants[i][j]==0 for i in range(5)for j in range(5)if i!=j)
    # The row order is R10+T35 and columns here are original0..44.  The
    # sparse-block packet used a different column order; compare up to sign.
    ordered=spec['selected_columns']+spec['remaining_q_columns']
    sign=(-1)**sum(ordered[i]>ordered[j]for i in range(45)for j in range(i+1,45))
    assert diagonal==[sign*d%P for d in spec['known_five_original_block_determinants']]
    hcoeff=list(range(1,46));absolute=[(delta+np.array(list(map(residue,base)),dtype=np.int64))%P for delta in deltas]
    hvalues=[int(np.array(hcoeff,dtype=np.int64)@point)%P for point in absolute]
    assert len(set(hvalues))==5
    vandermonde=1
    for i in range(5):
        for j in range(i+1,5):vandermonde=vandermonde*(hvalues[j]-hvalues[i])%P
    norm=1
    for d in diagonal:norm=norm*d%P
    assert norm and vandermonde
    # A small exact fixture retains a nonreduced UNKNOWN residual.  It checks
    # that the construction makes no unwarranted reducedness/exhaustion claim.
    x,y=sp.symbols('x y');toyI=[x*x-x,x*y,y*y]
    toyD=sp.expand(x*((x+1)/2+y));toyG=sp.groebner(toyI,x,y,domain=sp.QQ)
    assert toyG.reduce(toyD-x)[1]==0 and toyG.reduce(toyD*toyD-toyD)[1]==0
    assert list(sp.groebner(toyI+[toyD],x,y,domain=sp.QQ).exprs)==[x,y*y]
    (args.output/'polarization134.json').write_text(json.dumps(dict(
        all_sixteen_beta_directions_retained=True,all134_exact_QQ_identities=True,
        identity='F(q)-F(p)=Jac_q(F)((q+p)/2)*(q-p)',rows=statistics),indent=2)+'\n')
    (args.output/'known-five-separator.json').write_text(json.dumps(dict(
        selected_original_rows_zero_based=rows,prime=P,midpoint25_determinants=determinants,
        all25_full134_kernel_checks=True,known_five_jacobian_determinants=diagonal,
        Jacobian_norm_mod_prime=norm,separating_h_coefficients=hcoeff,separating_h_values=hvalues,
        h_Vandermonde_mod_prime=vandermonde,h_discriminant_mod_prime=vandermonde*vandermonde%P,
        known_five_algebra_source='Actual quotient from the universal P29 finite fiber, not point counts.',
        generic_beta_units='Norm_E(d(p)) and determinant of the separating h-power basis; both beta-only.',
        generic_recipe='e(q)=Tr_E/R(det(J_selected((q+p)/2))/det(J_selected(p)))',
        maximum_q_degree=45,exact_generic_polynomial_expanded=False,
        whole_affine_algebra_split='A = E_known5 x A/(e)',
        residual_exact_ideal='I+(e)',known_five_exact_ideal='I+(1-e)',
        no_q_or_cyclic_or_H_localization=True,
        nonreduced_residual_fixture={'ideal':['x^2-x','x*y','y^2'],'known_point':[1,0],
            'idempotent':'x*((x+1)/2+y)','idempotent_mod_I':'x','residual':'QQ[y]/(y^2)'},
        residual_emptiness_claimed=False),indent=2)+'\n')
    dependencies=[ROOT/'reports/S61_AFFINE_SEVENTEEN_BASE_FIBER_20261103_v62/report.json',
        ROOT/'reports/S61_AFFINE_SEVENTEEN_BASE_LIFT_20261103_v62/report.json',
        ROOT/'reports/S61_AFFINE_LINEAR_BLOCK_20261103_v63/report.json']
    for dependency in dependencies:
        report=json.loads(dependency.read_text());assert report['all_passed']
        for path,value in report['input_sha256'].items():assert sha(ROOT/path)==value
        for path,value in report.get('evidence_sha256',{}).items():assert sha(dependency.parent/path)==value
    inputs=[Path(__file__),original,specpath,fivepath,knownpath,*dependencies,
        ROOT/'research/notes/S61_AFFINE_SEVENTEEN_BASE_LIFT_20261103_v62.md',
        ROOT/'research/notes/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63.md']
    report=dict(id='S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63',all_passed=True,
        full134_QQ_polarization_with_all16beta_verified=True,
        all25_midpoint_determinants_and_full_kernel_checks=True,
        knownfive_Jacobian_and_separating_discriminant_are_beta_units=True,
        global_degree45_idempotent_is_finite_determinant_trace_recipe=True,
        whole_affine_scheme_splits_knownfree5_and_unresolved_residual=True,
        no_q_localization_or_cyclic_or_quintic_H_assumption=True,
        actual_quotient_E_dependency_audited_by_hashes_not_recomputed=True,
        finite_over_beta_or_full_length5_or_residual_empty_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS globaldegree45 knownfive idempotent;exact134QQ+25midpoint;residualretained;',report['seconds'])


if __name__=='__main__':main()
