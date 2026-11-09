#!/usr/bin/env python3
"""Complete quadratic anchors and magnetic sources for twelve plane configurations."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import minor_certificate
ROOT=Path(__file__).resolve().parents[1]

def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(90,91));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    begun=time.monotonic();X,Y,Z=s.symbols('X Y Z');original=(X,Y,Z)
    node=Y**2*Z-X**2*(X+Z);nonsplit=Y**2*Z-X**2*(X-3*Z);cusp=Y**2*Z-X**3
    tangent=-3*X+2*Y+Z
    cases=[('node_two_flexes',nonsplit,Z,-3*X+Y+8*Z),
           ('node_tangent_flex',node,Y-X,Z),
           ('node_flex_secant',node,X,Z),
           ('node_two_node_tangents',node,Y-X,Y+X),
           ('node_flex_tangent_through_flex',node,Z,X+Z),
           ('cusp_flex_tangent',cusp,Z,tangent),
           ('cusp_flex_secant_through_flex',cusp,Z,X-Z),
           ('cusp_two_tangents_incident',cusp,tangent,12*X+16*Y-Z),
           ('cusp_flex_cusp_secant',cusp,Z,Y-X),
           ('cusp_tangent_cusp_secant',cusp,tangent,Y-X),
           ('cusp_cusp_tangent_ordinary_tangent',cusp,Y,tangent),
           ('cusp_cusp_tangent_cusp_secant',cusp,Y,Y-X)]
    h,x,y,z,w=xx=s.symbols('h x y z w');base=(h,x,y);mons=(h*h,h*x,h*y,x*x,x*y,y*y)
    records=[]
    for name,C,L1,L2 in cases:
        started=time.monotonic()
        frame=next(s.Matrix([[s.diff(L,v) for v in original]for L in (L1,L2,L3)])for L3 in original
                   if s.Matrix([[s.diff(L,v)for v in original]for L in (L1,L2,L3)]).det()!=0)
        inverse=frame.inv()*s.Matrix(base);q=s.expand(C.subs(dict(zip(original,inverse)),simultaneous=True))
        F=s.expand(h*x*q);a=s.symbols('a0:18')
        candidate=[sum(a[6*i+j]*m for j,m in enumerate(mons))for i in range(3)]
        eq=s.expand(sum(candidate[i]*s.diff(F,base[i])for i in range(3)))
        matrix,rhs=s.linear_eq_to_matrix(s.Poly(eq,*base).coeffs(),a);assert rhs==s.zeros(matrix.rows,1)
        kernel=matrix.nullspace();record=dict(case=name,original_cubic=str(C),original_lines=[str(L1),str(L2)],
            plane_frame=[[str(v)for v in row]for row in frame.tolist()],cubic=str(q),plane_quintic=str(F),
            quadratic_gradient_syzygy_dimension=len(kernel),full_linear_matrix_shape=list(matrix.shape))
        if len(kernel)!=2:
            record.update(free_exponents_2_2=False,reason='The required full quadratic syzygy space does not have dimension two.')
            records.append(record);print(json.dumps(record),flush=True);continue
        syzygies=[];anchors=[]
        for column in kernel:
            values=[s.expand(sum(column[6*i+j]*m for j,m in enumerate(mons)))for i in range(3)]
            assert s.expand(sum(values[i]*s.diff(F,base[i])for i in range(3)))==0
            multiplier=s.cancel(values[0]/h);assert multiplier.is_polynomial(*base)
            anchor=[s.expand(values[1]-x*multiplier),s.expand(values[2]-y*multiplier)]
            assert s.cancel(anchor[0]/x).is_polynomial(*base)
            syzygies.append(values);anchors.append(anchor)
        U,V=map(s.Matrix,anchors);det=s.expand(U[0]*V[1]-U[1]*V[0]);ratio=s.cancel(det/(x*q))
        if ratio==0:
            record.update(free_exponents_2_2=False,reason='All quadratic logarithmic anchors have generic rank at most one.')
            records.append(record);print(json.dumps(record),flush=True);continue
        assert not ratio.free_symbols;V=V/ratio
        Ua,Va=U.subs(h,1),V.subs(h,1)
        bracket=s.Matrix([s.expand(sum(Ua[j]*s.diff(Va[i],(x,y)[j])-Va[j]*s.diff(Ua[i],(x,y)[j])for j in range(2)))for i in range(2)])
        da=s.expand(x*q.subs(h,1));c=s.cancel((bracket[0]*Va[1]-bracket[1]*Va[0])/da);d=s.cancel((Ua[0]*bracket[1]-Ua[1]*bracket[0])/da)
        assert all(s.Poly(f,x,y).total_degree()<=1 for f in [c,d])
        ch=s.expand(h*c.subs({x:x/h,y:y/h},simultaneous=True));dh=s.expand(h*d.subs({x:x/h,y:y/h},simultaneous=True))
        params=s.symbols('r0:6');mag=sum(a*m for a,m in zip(params,mons));P=s.zeros(5)
        for ij,f in {(1,3):U[0],(2,3):U[1],(1,4):V[0],(2,4):V[1],(3,4):-ch*z-dh*w+mag}.items():P[ij]=f;P[ij[::-1]]=-f
        common=s.ilcm(1,*[s.denom(a)for f in P for a in s.Poly(f,*xx,*params).coeffs()]);P=P*common
        assert all(s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l])for l in range(5)))==0 for i,j,k in itertools.combinations(range(5),3))
        pf=s.expand(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3]);assert s.expand(pf+common**2*x*q)==0
        values=dict(zip(params,(1,2,3,5,7,11)));point=P.subs(values);directions=[P.diff(p)for p in params]
        D,Vert,Orbit,Family=ambient_matrices(point,xx,directions)
        matrices={'D':D,'vertical':Vert,'orbit_scale':Orbit,'full_source':Family}
        ranks={k:len(DomainMatrix.from_Matrix(A).convert_to(s.QQ).rref()[1])for k,A in matrices.items()}
        certificates={k:minor_certificate(A,ranks[k])for k,A in matrices.items()}
        encode=lambda A:{str(ij):str(A[ij])for ij in itertools.combinations(range(5),2)if A[ij]}
        record.update(free_exponents_2_2=True,all_quadratic_syzygies=[[str(v)for v in row]for row in syzygies],
            anchors=[list(map(str,U)),list(map(str,V))],anchor_bracket_coefficients=list(map(str,(ch,dh))),
            full_six_magnetic_source=encode(P),point_brackets=encode(point),parameter_values={str(k):str(v)for k,v in values.items()},
            all6_parameter_directions=list(map(encode,directions)),all_parameter_strict_Jacobi_zero=True,
            integer_scale=str(common),actual_degeneracy_quintic=str(-common**2*F),
            ranks_QQ=ranks,certificates=certificates,actual_dimension_lower=ranks['full_source']-25,
            ambient_tangent_upper=125-ranks['D'],smooth_component_squeeze=ranks['full_source']+ranks['D']==150,
            seconds=time.monotonic()-started)
        records.append(record)
        print(json.dumps({k:record[k]for k in ['case','ranks_QQ','actual_dimension_lower','ambient_tangent_upper','smooth_component_squeeze','seconds']}),flush=True)
    paths=['code/classification_two_lines_plane_cubic.py','code/classification_three_hyperplanes_quadric_independent.py','code/classification_rank3_mixed_component_review.py']
    report=dict(task_id='E_reduced_nonnormal',all_passed=True,cases=records,
        geometric_configuration_exhaustion_requires_separate_written_proof=True,
        all_coefficients_of_each_actual_source_retained=True,component_distinctness_requires_geometric_proof=True,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as handle:json.dump(report,handle,indent=2);handle.write('\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)
