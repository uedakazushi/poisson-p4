#!/usr/bin/env python3
"""Independent exact check: reconstruct bivectors first, then take Jacobi.

The construction scripts expand the Jacobi equations first and then eliminate.
This auditor reverses those operations, treating parameters as extra formal
indeterminates and differentiating only the five homogeneous coordinates.
Every selected generator is compared coefficient by coefficient over Q.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys
import time

HERE=Path(__file__).resolve().parent
import poisson_p4_components_check as p

def audit(folder,prefix,constants,row_key):
    t0=time.time()
    folder=Path(folder)
    cert=json.loads((folder/(prefix+'_certificate.json')).read_text())
    path=folder/(prefix+'_equations.json')
    data=json.loads(path.read_text())
    n=data['parameter_count']
    raw={}
    for k,r in enumerate(cert['raw_free_indices']):raw[r]=[((k,),F(1))]
    for entry in data['eliminated_raw_variables']:
        raw[entry['raw_index']]=[(tuple(m),F(a,b)) for m,a,b in entry['polynomial']]
    bv=defaultdict(dict)
    for pair,c in constants.items():bv[pair][(2,0,0,0,0)+(0,)*n]=F(c)
    for r,terms in raw.items():
        pair,alpha=p.BIVECTOR_BASIS[r]
        for mon,c in terms:
            powers=[0]*n
            for k in mon:powers[k]+=1
            exps=tuple(alpha)+tuple(powers)
            bv[pair][exps]=bv[pair].get(exps,F(0))+c
    jac=p.affine_jacobi_coefficients(dict(bv))
    rows=[defaultdict(F) for _ in p.FOUR_VECTOR_BASIS]
    for triple,poly in jac.items():
        for e in range(5):
            sign=p.permutation_sign((e,)+triple)
            if not sign:continue
            quad=tuple(sorted((e,)+triple))
            for exps,c in poly.items():
                alpha=list(exps[:5]);alpha[e]+=1
                row=p.FOUR_VECTOR_INDEX[(quad,tuple(alpha))]
                mon=tuple(k for k,count in enumerate(exps[5:]) for _ in range(count))
                rows[row][mon]+=sign*c
    rows=[{m:c for m,c in q.items() if c} for q in rows]
    indices=cert[row_key]
    assert len(indices)==len(data['equations'])
    for i,eq in zip(indices,data['equations']):
        if isinstance(eq,dict):
            assert eq['original_row']==i
            eq=eq['polynomial']
        expected={tuple(m):F(a,b) for m,a,b in eq}
        assert rows[i]==expected,(prefix,i)
    erows=cert.get('equation_rows_used_for_elimination',cert.get('elimination_rows'))
    assert all(not rows[i] for i in erows)
    report={
        'system':prefix,'parameters':n,'checked_generator_rows':len(indices),
        'checked_eliminated_rows':len(erows),
        'all_coefficient_comparisons_exact_over_Q':True,
        'method':'Reconstruct all coefficient polynomials; independently expand Schouten; compare exact coefficient dictionaries.',
        'equation_json_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'seconds':round(time.time()-t0,3)
    }
    print(json.dumps(report),flush=True)
    return report

if __name__=='__main__':
    target=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE.parent/'data'
    reports=[
      audit(target,'universal_slice',{(3,4):1},'independent_original_equation_rows'),
      audit(target,'symplectic_slice',{(1,2):1,(3,4):1},'remaining_equation_rows')
    ]
    (target/'independent_reconstruction_audit.json').write_text(json.dumps(reports,indent=2)+'\n')
