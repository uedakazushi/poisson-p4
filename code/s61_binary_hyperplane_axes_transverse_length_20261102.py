#!/usr/bin/env python3
"""Certify original-cut transverse length65 at a=0 and F=0 over QQ."""
from fractions import Fraction
from pathlib import Path
import json
import time

import numpy as np
import sympy as sp

from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank,determinant_mod
from s61_binary_hyperplane_C_transverse_length_20261102 import monomials,substitute
from s61_binary_hyperplane_root_model_20261102 import (
    ROOT,CODE,P,S61,actual_binary_tangent,digest,hashes,old_cut_rows,original_point,
    output_argument,residue,source_paths,write_new,
)

QUARTICS=ROOT/'reports/S61_BINARY_HYPERPLANE_ORIGINAL_ROOT_QUARTICS_20261102.json'
NOTE=ROOT/'research/notes/S61_BINARY_HYPERPLANE_AXES_TRANSVERSE_LENGTH_20261102.md'


def rank_certificate(matrix):
    """Compute rank over QQ and exhibit a matching good-reduction minor."""
    exact=int(matrix.to_DM().rank())
    if not exact:
        assert matrix==sp.zeros(*matrix.shape)
        return {'rank_QQ':0,'minor_rows':[],'minor_columns':[],'minor_modp':1}
    reduced=np.array([[residue(c) for c in row] for row in matrix.tolist()],dtype=np.int64)
    rank,columns=elimination_rank(reduced,P)
    assert rank==exact
    _,rows=elimination_rank(reduced[:,columns].T,P)
    value=int(determinant_mod(reduced[np.ix_(rows,columns)],P))
    assert value
    return {'rank_QQ':exact,'minor_rows':rows,'minor_columns':columns,'minor_modp':value}


def check_point(model,name,a,quartic):
    beta,two,old,weights,rw,full=original_point(model,a)
    curve,plane,tangent=actual_binary_tangent(model,beta,two,old)
    cuts=old_cut_rows(old)
    determinant=tangent[cuts,:].det()
    expected={'a':sp.Rational(-52,25),'F':sp.Rational(82,125)}[name]
    assert determinant==expected and residue(determinant)
    linear=sp.Matrix([[q.get((j,),0) for j in range(45)] for q in full])
    assert linear*tangent==sp.zeros(135,3)
    raw_variables=[j for j in range(45) if j not in cuts]
    linear_ranks={}
    for w in (1,2,3,4):
        cols=[j for j in raw_variables if weights[j]==w]
        rows=[i for i,weight in enumerate(rw) if weight==w]
        audit=rank_certificate(linear.extract(rows,cols))
        audit.update(raw_rows=rows,raw_columns=cols,remaining=len(cols)-audit['rank_QQ'])
        linear_ranks[str(w)]=audit
    assert linear_ranks['3']['remaining']==linear_ranks['4']['remaining']==0
    equations=[{m:c for m,c in q.items() if not set(cuts).intersection(m)}
               if rw[i]<=3 else {} for i,q in enumerate(full)]
    variables=[j for j in raw_variables if weights[j]<=3]
    pivots=[]
    for w in (1,2,3):
        while True:
            choices=[(len(q),i,m[0]) for i,q in enumerate(equations) if rw[i]==w
                     for m in q if len(m)==1]
            if not choices:
                break
            _,row,var=min(choices)
            q=equations[row]
            pivot=q[(var,)]
            assert residue(pivot) and all(var not in m for m in q if m!=(var,))
            replacement={m:-c/pivot for m,c in q.items() if m!=(var,)}
            assert all(c.denominator%P for c in replacement.values())
            equations=[{} if k==row else substitute(f,var,replacement)
                       for k,f in enumerate(equations)]
            variables.remove(var)
            pivots.append([row,var,str(pivot),residue(pivot)])
    for i,q in enumerate(equations):
        assert all(set(m).issubset(variables) and len(m)>=2
                   and sum(weights[j] for j in m)==rw[i] for m in q)
    assert all(sum(weights[j]==w for j in variables)==linear_ranks[str(w)]['remaining']
               for w in (1,2,3))
    assert sorted(set(weights[j] for j in variables))==[1,2]
    mons={w:list(monomials(variables,weights,w)) for w in range(4)}
    pieces={}
    for w in range(4):
        index={m:i for i,m in enumerate(mons[w])}
        rows=[]
        labels=[]
        for i,q in enumerate(equations):
            if not q or rw[i]>w:
                continue
            for shift in mons[w-rw[i]]:
                row=[0]*len(index)
                for m,c in q.items():
                    row[index[tuple(sorted(m+shift))]]=sp.Rational(c.numerator,c.denominator)
                rows.append(row)
                labels.append([i,list(shift)])
        matrix=sp.Matrix(rows) if rows else sp.zeros(0,len(index))
        audit=rank_certificate(matrix)
        audit.update(columns=len(index),rows=len(rows),dimension=len(index)-audit['rank_QQ'],
                     minor_source_labels=[labels[i] for i in audit['minor_rows']])
        pieces[str(w)]=audit
    assert quartic['rank_134']==quartic['rank_135']==2814 and quartic['weight4_quotient_dimension']==0
    assert quartic['minor_det_modp']%P
    hf=[pieces[str(w)]['dimension'] for w in range(4)]+[0]
    assert hf==[1,7,23,34,0]
    return {'a':str(a),'e':2,'d':3,'cuts':cuts,'actual_transverse_determinant_QQ':str(determinant),
            'transverse_determinant_modp':residue(determinant),
            'actual_binary_tangent_45x3_QQ':[[str(v) for v in row] for row in tangent.tolist()],
            'curve_jacobian_QQ':[[str(v) for v in row] for row in curve.tolist()],
            'plane_tangent_QQ':list(map(str,plane)),'all135_Jacobian_rows_annihilate_tangent':True,
            'linear_ranks_by_weight':linear_ranks,'low_pivots':pivots,'remaining_variables':variables,
            'remaining_weights':[weights[j] for j in variables],'low_pieces':pieces,
            'weight4_minor_modp':quartic['minor_det_modp'],'weighted_Hilbert_function_QQ':hf,
            'weighted_Hilbert_function_modp':hf,'all_weights_at_least_four_zero':True,'length':sum(hf)}


def verify():
    start=time.time()
    quartic=json.loads(QUARTICS.read_text())
    assert quartic['all_passed'] and quartic['cuts']==[5,19,23]
    for name,value in quartic['input_sha256'].items():
        assert digest(ROOT/name)==value
    model=S61()
    points={name:check_point(model,name,a,quartic['points'][name])
            for name,a in [('a',Fraction(0)),('F',Fraction(-171,250))]}
    static=source_paths()+[Path(__file__),QUARTICS,NOTE,
        CODE/'s61_binary_hyperplane_C_transverse_length_20261102.py',
        CODE/'s61_binary_deeper_uniform_quadratic_20261101.py']
    return {'all_passed':True,'id':'S61_BINARY_HYPERPLANE_AXES_TRANSVERSE_LENGTH_20261102',
            'scope':'Two original-cut transverse algebras at specified rational points, with matching good reduction.',
            'points':points,'length_at_each_point':65,'prime':P,'input_sha256':hashes(static),
            'quartic_dependency_inputs':quartic['input_sha256'],'seconds':time.time()-start}


if __name__=='__main__':
    path=output_argument(__doc__)
    result=verify()
    write_new(path,result)
    print('PASS',[(k,v['weighted_Hilbert_function_QQ']) for k,v in result['points'].items()],result['seconds'])
