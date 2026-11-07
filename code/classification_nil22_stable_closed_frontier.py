#!/usr/bin/env python3
"""Finite divisor/cofactor exclusions for both whole stable magnetic sources.

The closed common-tangent incidence excludes CPLUS_EXT. A dependent-three-H
collision lemma excludes its known 24-dimensional component. Neither result
is a componenthood or MPS nonmembership assertion.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s
from classification_nil22_image_source_review import build_P,coords,mag,x,y

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic()
    h,u,v=s.symbols('h u v');base=(h,x,y)
    specs=[('stable_line_conic',s.Matrix([[3,-15],[-15,3]]),s.Matrix([s.Rational(3,2),-s.Rational(3,2)]),[h,h,x-y-12*h],[36*h*h+3*h*(x-y)+2*x*y]),
        ('stable_three_lines',s.Matrix([[-6,3],[3,-6]]),s.zeros(2,1),[h,h,x+3*h,y-3*h,x-y-3*h],[])]
    records=[]
    for name,S,V0,lines,quadrics in specs:
        P,M,D=build_P(x*y*(x-y),S,V0,-s.eye(2)/3)
        cofactor={}
        for L in lines[2:]:
            affineL=L.subs(h,1);Z=s.Matrix([-s.cancel(sum(M[i,j]*s.diff(affineL,(x,y)[j])for j in range(2))/affineL)for i in range(2)]+[0,0])
            assert all(not s.denom(z).free_symbols for z in Z)
            assert not any(z.has(*mag)for z in Z)
            linear=Z.jacobian(s.Matrix(coords));assert linear.rank()==2
            assert linear*linear==s.zeros(4)
            cofactor[str(L)]=dict(full_logarithmic_Hamiltonian=list(map(str,Z)),linear_matrix=linear.tolist(),linear_rank=2)
        plus=[]
        for outer_index,L in enumerate(lines):
            if outer_index==1:continue
            remaining=lines[:outer_index]+lines[outer_index+1:]
            if quadrics:
                assert len(remaining)==2;partitions=[(quadrics[0],s.prod(remaining))]
            else:
                assert len(remaining)==4
                partitions=[(remaining[0]*remaining[k],s.prod([remaining[j]for j in range(1,4)if j!=k]))for k in range(1,4)]
            row=s.Matrix([[s.diff(L,z)for z in base]]);K=s.Matrix.hstack(*row.nullspace());assert K.shape==(3,2)
            values=dict(zip(base,K*s.Matrix([u,v])))
            for Q1,Q2 in partitions:
                restricted=[s.expand(Q.subs(values,simultaneous=True))for Q in (Q1,Q2)]
                ranks=[s.hessian(Q,(u,v)).rank()for Q in restricted]
                coefficient_rows=s.Matrix([[s.Poly(Q,u,v).coeff_monomial(mon)for mon in [u*u,u*v,v*v]]for Q in restricted])
                common=ranks[0]<=1 and ranks[1]<=1 and coefficient_rows.rank()<=1
                assert not common
                plus.append(dict(outer_H=str(L),quadrics=list(map(str,(Q1,Q2))),restricted_quadrics=list(map(str,restricted)),restriction_ranks=ranks,proportionality_row_rank=coefficient_rows.rank(),closed_common_tangent_incidence=False))
        dependent=[]
        for indices in itertools.combinations(range(len(lines)),3):
            chosen=[lines[i]for i in indices]
            # Any unused degree-two factor is permitted; this only strengthens exclusion.
            if sum(1 for i in indices if i>=2)>len(lines)-2:raise AssertionError
            A=s.Matrix([[s.diff(L,z)for z in base]for L in chosen]);rank=A.rank()
            record=dict(selected_H=list(map(str,chosen)),linear_form_span=rank)
            if rank==3:record['excluded_by_closed_dependence']=True
            else:
                assert indices[0:2]==(0,1)
                L=chosen[2];record.update(excluded_by_collision_translation_lemma=True,limiting_double_H=str(h),simple_H=str(L),nonzero_linear_cofactor=cofactor[str(L)])
            dependent.append(record)
        records.append(dict(name=name,all_six_magnetic_parameters=list(map(str,mag)),simple_factor_full_cofactors=cofactor,
            all_CPLUS_factor_assignments=plus,all_dependent_three_H_assignments=dependent,
            entire_source_outside_CPLUS_EXT=True,entire_source_outside_dependent_3H_Q_component=True,
            scope='Uses the adopted proper common-tangent incidence and the proved regular double-plane collision lemma. No conclusion about another known component or a new component.'))
    paths=[Path(__file__),ROOT/'code/classification_nil22_image_source_review.py',ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_SOURCE_REVIEW_V2_20261107.json',ROOT/'research/classification/TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md',ROOT/'research/classification/THREE_HYPERPLANES_QUADRIC_DEPENDENT.md']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,cases=records,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6))
    with args.output.open('x')as handle:json.dump(result,handle,indent=2,default=str);handle.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))


if __name__=='__main__':main()
