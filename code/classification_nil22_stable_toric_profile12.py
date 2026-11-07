#!/usr/bin/env python3
"""One bounded structural (1,2) toric confluence for whole stable_three.

Retains the full linear r and all six regular period coefficients 0,...,4.
No other weighted profile or larger Taylor range is tried.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_nil22_image_source_review import build_P,coords,mag,x,y,q1,q2

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(4),2))


def jacobi(P,variables):
    return [s.expand(sum(P[i,l]*s.diff(P[j,k],variables[l])+P[j,l]*s.diff(P[k,i],variables[l])+P[k,l]*s.diff(P[i,j],variables[l])for l in range(4)))for i,j,k in itertools.combinations(range(4),3)]


def encode(P):return {f'{i},{j}':str(s.factor(P[i,j]))for i,j in PAIRS if P[i,j]!=0}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    h,H,Z,W,U,V,tau=s.symbols('h H Z W U V tau');r0,r1,r2=s.symbols('r0 r1 r2');a0,a1,a2=s.symbols('a0 a1 a2');variables=(Z,W,U,V);ell=1+U+V
    # Necessary full-six-magnetic bridge, not only a desired normal-form point.
    original,_,_=build_P(x*y*(x-y),s.Matrix([[-6,3],[3,-6]]),s.zeros(2,1),-s.eye(2)/3)
    old=s.Matrix([h,q1,q2,x,y]);new=s.Matrix([H,Z,W,U,V]);image=s.Matrix([3*h+x,-q2/3,(q2-q1)/3,3*h-y,3*h-x+y]);G=image.jacobian(old);assert G.det()!=0
    homogeneous=s.zeros(5)
    for i,j in PAIRS:
        entry=sum(c*h**(2-sum(e))*s.prod(z**k for z,k in zip(coords,e))for e,c in s.Poly(original[i,j],*coords).terms());homogeneous[i+1,j+1]=entry;homogeneous[j+1,i+1]=-entry
    raw=(G*homogeneous*G.T).subs(dict(zip(old,G.inv()*new)),simultaneous=True).applyfunc(s.expand)
    correction=s.Matrix([0,*[s.cancel(raw[0,j]/H)for j in range(1,5)]]);assert all(not s.denom(v).has(H)for v in correction)
    bridged=(9*(raw-new*correction.T+correction*new.T)).subs(H,1)[1:,1:].applyfunc(s.expand)
    target=s.zeros(4);target[0,2]=ell*U;target[1,3]=ell*V;target[0,1]=-V*Z+U*W
    for i,j in PAIRS:target[j,i]=-target[i,j]
    R=s.expand(bridged[0,1]-target[0,1]);assert not R.has(Z,W)
    assert all(s.expand(bridged[i,j]-target[i,j])==0 for i,j in PAIRS if(i,j)!=(0,1))
    poly=s.Poly(R,U,V);cU=poly.coeff_monomial(U);cV=poly.coeff_monomial(V);cUV=poly.coeff_monomial(U*V)
    LZ=-cV-cUV*U;LW=cU;images=s.Matrix([Z+LZ,W+LW,U,V]);J=images.jacobian(s.Matrix(variables))
    normalized=(J*bridged*J.T).subs({Z:Z-LZ,W:W-LW},simultaneous=True).applyfunc(s.expand)
    canonical_R=s.expand(normalized[0,1]-target[0,1]);assert s.Poly(canonical_R,U,V).monoms()==[(2,0),(0,2),(0,0)]
    A=s.Matrix([s.Poly(canonical_R,U,V).coeff_monomial(m)for m in (1,U*U,V*V)]);assert A.jacobian(mag).rank()==3
    target[0,1]+=a0+a1*U*U+a2*V*V;target[1,0]=-target[0,1];assert not any(jacobi(target,variables))
    # Toric affine coordinates are (U,V,L3,L4); H is its Casimir coordinate.
    r=r0+r1*U+r2*V;forms=s.Matrix([U,V,ell+tau*Z,ell+tau*Z+tau*r+tau*tau*W]);frame=forms.jacobian(s.Matrix(variables));assert s.factor(frame.det())==tau**3
    periods=s.symbols('k0:30');toric=s.zeros(4)
    for index,(i,j)in enumerate(PAIRS):
        value=sum(periods[5*index+d]*tau**d for d in range(5))*forms[i]*forms[j];toric[i,j]=value;toric[j,i]=-value
    family=(frame.inv()*toric*frame.inv().T).applyfunc(s.expand)
    assert all(s.denom(s.cancel(tau**4*v)).free_symbols<={r0,r1,r2}for v in family)
    equations=[];row_labels=[]
    for i,j in PAIRS:
        entry=s.expand(family[i,j])
        for degree in range(-4,1):
            expression=s.expand(entry.coeff(tau,degree)-(target[i,j]if degree==0 else 0))
            for exponent,coefficient in s.Poly(expression,*variables).terms():
                if coefficient!=0:
                    equations.append(coefficient);row_labels.append(dict(pair=[i,j],tau_degree=degree,spatial_exponent=list(exponent)))
    C,rhs=s.linear_eq_to_matrix(equations,(*periods,a0,a1,a2));aug=C.row_join(rhs)
    field=s.QQ.frac_field(r0,r1,r2);reduced,pivots=DomainMatrix.from_Matrix(aug).convert_to(field).rref();RR=reduced.to_Matrix()
    period_rank=sum(p<30 for p in pivots);compatibility=[];solution=dict.fromkeys(periods,s.Integer(0))
    for row,pivot in enumerate(pivots):
        if pivot<30:
            solution[periods[pivot]]=s.cancel(RR[row,33]-sum(RR[row,30+i]*a for i,a in enumerate((a0,a1,a2))))
        else:compatibility.append(s.factor(sum(RR[row,30+i]*a for i,a in enumerate((a0,a1,a2)))-RR[row,33]))
    obstruction=None
    if compatibility==[-1]:
        chosen_rows=DomainMatrix.from_Matrix(aug.T).convert_to(field).rref()[1]
        transpose=aug[list(chosen_rows),:].T;unit=s.eye(34)[:,33]
        witness_reduced,witness_pivots=DomainMatrix.from_Matrix(transpose.row_join(unit)).convert_to(field).rref();witness_matrix=witness_reduced.to_Matrix()
        assert list(witness_pivots)==list(range(len(chosen_rows)))
        weights=s.Matrix([witness_matrix[i,len(chosen_rows)]for i in range(len(chosen_rows))])
        assert (transpose*weights-unit).applyfunc(s.cancel)==s.zeros(34,1)
        denominator=s.factor(s.lcm([s.denom(value)for value in weights]))
        obstruction=dict(original_row_indices=list(chosen_rows),exact_row_weights=list(map(str,weights)),all34_coefficient_identity_verified=True,
            identity='sum(weights*original_rows)=(0,...,0,1), so 0=1',common_denominator=str(denominator),uniform_for_all_r=not denominator.free_symbols)
    progress=dict(status='generic_r_full_linear_reduction_complete',matrix_shape=list(aug.shape),period_rank=period_rank,total_augmented_rank=len(pivots),all_r_variables=list(map(str,(r0,r1,r2))),all_period_coefficients=list(map(str,periods)),
        augmented_column_order=[*map(str,periods),*map(str,(a0,a1,a2)),'RHS'],all_original_row_labels=row_labels,all_original_augmented_rows=[[str(v)for v in aug.row(i)]for i in range(aug.rows)],
        RREF_pivot_columns=list(pivots),all_nonzero_RREF_rows=[[str(v)for v in RR.row(i)]for i in range(len(pivots))],RREF_zero_row_count=aug.rows-len(pivots),explicit_linear_obstruction=obstruction,
        compatibility_conditions=list(map(str,compatibility)),generic_period_solution={str(k):str(v)for k,v in solution.items()},full_six_bridge=dict(PGL_matrix=G.tolist(),tensor_scale=9,Euler_correction=list(map(str,correction)),canonical_magnetic_parameters=list(map(str,A)),global_polynomial_fiber_translation=list(map(str,(LZ,LW)))))
    Path('/tmp/nil22_stable_toric_profile12_linear.json').write_text(json.dumps(progress,indent=2,default=str)+'\n')
    whole=not compatibility;actual=None
    if whole:
        denominator=s.prod(s.denom(value)for value in solution.values());assert not denominator.has(a0,a1,a2)
        bound=max(s.degree(denominator,z)for z in (r0,r1,r2));chosen=None
        for values in itertools.product(range(int(bound)+1),repeat=3):
            trial=dict(zip((r0,r1,r2),values))
            if denominator.subs(trial)!=0:chosen=trial;break
        assert chosen is not None
        chosen_solution={k:s.cancel(v.subs(chosen))for k,v in solution.items()}
        regular=family.subs(chosen_solution,simultaneous=True).subs(chosen,simultaneous=True).applyfunc(s.cancel)
        assert all(not s.denom(v).has(tau)for v in regular)
        assert (regular.subs(tau,0)-target).applyfunc(s.cancel)==s.zeros(4)
        assert not any(jacobi(regular,variables))
        actual=dict(chosen_r={str(k):int(v)for k,v in chosen.items()},all_six_period_polynomials=[str(s.expand(sum(chosen_solution[periods[5*i+d]]*tau**d for d in range(5))))for i in range(6)],whole_regular_family=encode(regular),endpoint=encode(target),all_tau_strict_Jacobi_zero=True,no_negative_powers=True)
    paths=[Path(__file__),ROOT/'code/classification_nil22_image_source_review.py',ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_SOURCE_REVIEW_V2_20261107.json',ROOT/'reports/CLASSIFICATION_NIL22_STABLE_PROJECTIVE_ACTION_V2_20261107.json',ROOT/'research/classification/FIVE_HYPERPLANE_SUPPORT_REVIEW.md']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,profile=[1,2],inverse_bivector_pole_bound=4,period_degree_bound=4,no_other_profile_or_larger_degree_tried=True,linear_reduction=progress,whole_three_magnetic_source_T01_covered=whole,actual_family=actual,
        scope='Positive: literal regular toric Poisson family and complete six-magnetic bridge. Negative: only this full-linear-r structural profile on the generic-r open; special r or other profiles are not inferred and will not be scanned.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as handle:json.dump(result,handle,indent=2,default=str);handle.write('\n')
    print(json.dumps(dict(all_passed=True,whole_three_magnetic_source_T01_covered=whole,compatibility_conditions=list(map(str,compatibility)),seconds=result['seconds'])))


if __name__=='__main__':main()
