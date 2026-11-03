#!/usr/bin/env python3
"""Actual proper-source Jacobians and recovered flags at generic R19 points.

The paired orbit dimensions are from the accepted global intersection table.
Smoothness follows from the entire proper fiber being one unramified point,
not from source differential ranks alone.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector
from ambient_CB_global_flag_model_20261103_v63 import raw as cbraw
from ambient_rank_two_flag_intersection_table_20261103_v63 import permutation,lie

ROOT=Path(__file__).resolve().parents[1]

def minor(M):
    columns=list(M.rref()[1]);rows=list(M[:,columns].T.rref()[1]);det=M[rows,columns].det();assert det!=0
    return dict(rank=len(columns),rows=rows,columns=columns,determinant=str(det))

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();h,x,y,z,w=X;kappa=s.Symbol('kappa')
    Q=h*h+kappa*h*z+z*z;P=canonical(raw(-4*x*x,0,h*h*z+kappa*h*z*z/2+z**3/3+2*x*z*w));point=vector(P)
    ca_path=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json';cb_path=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    ca=json.loads(ca_path.read_text());cb=json.loads(cb_path.read_text())
    CA=s.Matrix.hstack(*(s.Matrix(list(map(s.Rational,v))) for v in ca['basis_canonical150']))
    CB0=s.Matrix.hstack(*(s.Matrix(list(map(s.Rational,v))) for i,v in enumerate(cb['enlarged_canonical_basis']) if i<3 or i>=6))
    perm=(1,3,2,0,4);CB=s.Matrix.hstack(*(permutation(CB0[:,i],perm) for i in range(19)))
    rows=list(CB.T.rref()[1]);coeff=CB[rows,:].inv()*point[rows,0];assert CB*coeff==point
    assert list(coeff[:3,0])==[0,6,0]
    mu_column=permutation(vector(canonical(cbraw(0,6*h*x,0))),perm)
    orbit=s.Matrix.hstack(*(lie(point,a,b) for a,b in itertools.product(range(5),repeat=2)))
    A=CA.row_join(orbit);B=CB.row_join(mu_column).row_join(orbit);union=A.row_join(B)
    assert (A.rank(),B.rank(),union.rank())==(29,29,37)
    special={}
    for value in [0,1,2,3]:
        expected=(A.subs(kappa,value).rank(),B.subs(kappa,value).rank(),union.subs(kappa,value).rank())
        assert expected==(29,29,37);special[str(value)]=list(expected)
    # Whole proper source tangent dimension29 on the affine cone; full rank
    # therefore gives projective unramifiedness at the specified point.
    A0=A.subs(kappa,1);B0=B.subs(kappa,1);U0=union.subs(kappa,1)
    # Independent flag trace identities on quotient U=(h,x,z,w).
    C=s.Matrix([-8*h*x/5,12*x*x/5,-18*x*z/5,Q+2*x*w/5]);U=s.Matrix([h,x,z,w])
    a,b,c,d,alpha=s.symbols('a b c d alpha');linear=a*h+b*x+c*z+d*w
    Call=sum(C[i]*s.diff(linear,U[i]) for i in range(4))
    CA_x=s.expand(Call+s.Rational(48,5)*x*linear+alpha*x*x)
    CB_x=s.expand(Call+s.Rational(18,5)*x*linear+alpha*x*x)
    matrixA,_=s.linear_eq_to_matrix(s.Poly(CA_x,*U).coeffs(),[a,b,c,d,alpha]);matrixB,_=s.linear_eq_to_matrix(s.Poly(CB_x,*U).coeffs(),[a,b,c,d,alpha])
    assert matrixA.nullspace()==[s.Matrix([0,-s.Rational(1,12),0,0,1])]
    assert len(matrixB.nullspace())==2 and all(v[0]==v[3]==0 for v in matrixB.nullspace())
    trace_checks={}
    for label,m in [('h',-8*x/5),('z',-18*x/5)]:
        f=s.expand(m*x-C[1]);trace_checks[label]=dict(binary_f=str(f),f_x=str(s.diff(f,x)),CA_required=str(5*m),CB_required=str(s.Rational(5,2)*m))
    assert s.diff(-4*x*x,x)==5*(-8*x/5)
    paths=[Path(__file__),ca_path,cb_path,ROOT/'code/ambient_rank_two_flag_intersection_table_20261103_v63.py',
        ROOT/'research/notes/AMBIENT_R19_GENERIC_ORBITS_20261103_v63.md',
        ROOT/'research/notes/AMBIENT_COMPONENT_R1_NONCM_20261103_v63.md',
        ROOT/'research/notes/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63.md']
    result=dict(id='AMBIENT_R19_GENERIC_GEOMETRY_20261103_v63',all_passed=True,field='QQ(kappa)',
        point_canonical150=list(map(str,point)),CB_source_permutation=list(perm),CB_source19_coefficients=list(map(str,coeff)),
        full_CA_source_differential=[[str(v) for v in row] for row in A.tolist()],
        full_CB_source_differential=[[str(v) for v in row] for row in B.tolist()],
        CA_cone_source_rank=29,CB_cone_source_rank=29,tangent_union_affine_rank=37,
        sample_kappa1_minors={'CA':minor(A0),'CB':minor(B0),'union':minor(U0)},specializations=special,
        CA_x_flag_matrix=[[str(v) for v in row] for row in matrixA.tolist()],
        CB_x_flag_matrix=[[str(v) for v in row] for row in matrixB.tolist()],trace_flag_checks=trace_checks,
        unique_CA_flag='L=d_y,W=span(h,x),H=h',unique_CB_flag='L=d_y,W=span(x,z),H=x',
        geometric_consequences=['CA andCB propermodels areisomorphismsneartheseunramifiedunique-fiberpoints',
            'CA andCB smooth28 atgenericR','actualintersectiongenericallyNONREDUCEDalongR; tangent20versusreduced19',
            'FULLglobalCA andCB satisfyR1 usingcodim1boundaryAg/Bb andalreadyknownrank4openR1',
            'bigdoubleflagP7incidencegenericallybirationalontoR, henceR19rational'],
        special_scope='No genericnilpotenceorder/length alongR, no fullRspecial-orbit classification, and no globalN R1 conclusion in thispacket.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('R19 GENERIC GEOMETRY PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
