#!/usr/bin/env python3
"""Closed line-incidence certificates for full reduced CA/CB first-jet fibers.

No frozen v62 source is changed. Exact QQ tests prove line-parabolic stability
of the raw affine16 families, their closed normalized affine10-plane fibers,
and their scheme-theoretic affine9 intersection. Proper homogeneous-bundle
descent and the unique rank3 kernel line are the separate geometric proof.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s

from additional_rank4_families_20261103_v62 import four,ROOT
from ambient_weight1_image import common_plane,primitive
from s61_model import S61,X,J,coefficients,quotient_coefficients,pushforward
from weight1_structure import sp_basis,COORDS,coefficient_map

def selector_left_inverse(matrix):
    rows=list(matrix.T.rref()[1]);selected=matrix[rows,:];assert selected.rows==selected.cols
    inverse=selected.inv();left=s.zeros(matrix.cols,matrix.rows)
    for i in range(matrix.cols):
        for j,row in enumerate(rows):left[i,row]=inverse[i,j]
    assert left*matrix==s.eye(matrix.cols)
    rest=[i for i in range(matrix.rows) if i not in rows]
    equations=s.eye(matrix.rows)[rest,:]-matrix[rest,:]*left
    assert equations*matrix==s.zeros(len(rest),matrix.cols)
    return left,equations,rows

def linear_beta_adapter(model,indices):
    columns=[];y=s.Matrix(s.symbols('y0:4'));at_zero=dict.fromkeys(y,0)
    for index in indices:
        u=[s.S.Zero]*61;u[index]=1;P=model.reconstruct(u)
        affine=s.Matrix(4,4,lambda i,j:s.expand((P[i+1,j+1]-X[i+1]*P[0,j+1]+X[j+1]*P[0,i+1])
            .subs({X[0]:1}|dict(zip(X[1:],y)))))
        linear=sum((affine.diff(v).subs(at_zero)*v for v in y),s.zeros(4))
        form=-J[1:5,1:5]*primitive(linear,y)*J[1:5,1:5]
        columns.append(s.Matrix([s.diff(form[i,j],y[k]) for (i,j),k in COORDS]))
    adapter=s.Matrix.hstack(*columns);assert adapter.rank()==16
    return adapter

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();h,x,y,z,w=X
    t,A,B,C,D,F=s.symbols('t A B C D F');base=[t,A,B,C,D,F];q=s.Symbol('q');k=s.symbols('k0:9');fiber=[q,*k]
    h2=A*x*z+B*x*w+C*z*z+D*z*w+F*w*w
    cubic=[s.prod((x,z,w)[i] for i in ids) for ids in itertools.combinations_with_replacement(range(3),3) if ids!=(0,0,0)]
    potential=h*h2+sum(a*b for a,b in zip(k,cubic));model=S61();weights=[r['weight'] for r in model.certificate['free_variable_basis']]
    beta_indices=[i for i,v in enumerate(weights) if v==1];q_indices=[i for i,v in enumerate(weights) if v==2]
    lie=sp_basis();line_action=s.Matrix.hstack(*(s.Matrix([g[i,1] for i in (0,2,3)]) for g in lie))
    null=line_action.nullspace();assert len(null)==7
    parabolic=[sum((c*g for c,g in zip(v,lie)),s.zeros(4)) for v in null]
    records={};vectors={}
    for name,p,gkind in [('CA',5*t/2,'A'),('CB',-5*t,'B')]:
        f=h*h+p*h*x+q*x*x;radial=0 if gkind=='A' else -s.diff(f,x)/2
        P=four([f,0,0,radial*z-s.diff(potential,w),radial*w+s.diff(potential,z),h*h if gkind=='A' else f])
        directions=[P.diff(v) for v in [*base,*fiber]];matrix=s.Matrix.hstack(*(s.Matrix(coefficients(R)) for R in directions))
        assert matrix.rank()==16
        assert s.Matrix(coefficients(P-h*h*J))==matrix*s.Matrix([*base,*fiber])
        left,_,_=selector_left_inverse(matrix);representations=[]
        for g in parabolic:
            full=s.zeros(5);full[1:5,1:5]=g;field=full*s.Matrix(X)
            assert full*J+J*full.T==s.zeros(5)
            action=s.Matrix.hstack(*(s.Matrix(coefficients(full*R+R*full.T-sum((field[i]*R.diff(X[i]) for i in range(5)),s.zeros(5)))) for R in directions))
            on_parameters=left*action;assert action==matrix*on_parameters
            representations.append(on_parameters)
        normalized=model.normalize(P);u=s.Matrix(normalized['parameters']);vectors[name]=u
        beta=u[beta_indices,0];Bmat=beta.jacobian(base);assert Bmat.rank()==6 and beta==Bmat*s.Matrix(base)
        bleft,bconstraint,brows=selector_left_inverse(Bmat)
        q0=u[q_indices,0].subs(dict.fromkeys(fiber,0));M=u[q_indices,0].jacobian(fiber)
        assert not set().union(*(v.free_symbols for v in M));assert M.rank()==10
        assert u[q_indices,0]==q0+M*s.Matrix(fiber)
        fleft,fconstraint,frows=selector_left_inverse(M)
        assert all(s.Poly(v,*base).total_degree()<=2 for v in q0)
        records[name]=dict(raw_direction_matrix=matrix.tolist(),parabolic_parameter_actions=[m.tolist() for m in representations],
            normalized_61_coordinates=list(map(str,u)),beta_embedding=Bmat.tolist(),beta_left_inverse=bleft.tolist(),
            beta_closed_constraints=bconstraint.tolist(),beta_selected_rows=brows,q_shift=list(map(str,q0)),
            q_direction_matrix=M.tolist(),fiber_left_inverse=fleft.tolist(),fiber_closed_constraints=fconstraint.tolist(),fiber_selected_rows=frows,
            raw_affine_family_dimension=16,normalized_base_dimension=6,normalized_fiber_dimension=10,
            connected_line_parabolic_preserves_whole_raw_family=True)
        print(name,'PARABOLIC7 + CLOSED AFFINE10 PASS',time.monotonic()-started,flush=True)
    assert records['CA']['beta_embedding']==records['CB']['beta_embedding']
    adapter=linear_beta_adapter(model,beta_indices);Bmat=s.Matrix(records['CA']['beta_embedding']);params,expected=common_plane()
    assert list(map(str,params))==list(map(str,base));assert adapter*Bmat==expected.jacobian(base)
    witness=coefficient_map(expected.subs(dict(zip(base,[2,3,5,7,11,13]))));assert witness.rank()==3
    assert witness.nullspace()==[s.Matrix([0,1,0,0])]
    MA=s.Matrix(records['CA']['q_direction_matrix']);MB=s.Matrix(records['CB']['q_direction_matrix'])
    assert MA[:,1:]==MB[:,1:];assert MA.row_join(MB).rank()==11
    target=s.Matrix(records['CB']['q_shift'])-s.Matrix(records['CA']['q_shift'])
    solve=MA.row_join(-MB[:,:1]);solution=solve.gauss_jordan_solve(target)[0]
    assert (solve*solution-target).applyfunc(s.expand)==s.zeros(45,1)
    assert solution[0]==0 and s.expand(solution[10]-25*t*t/4)==0
    assert list(map(s.expand,solution[1:10,0]))==[5*A*t/2,5*B*t/2,5*C*t/2,5*D*t/2,5*F*t/2,0,0,0,0]
    # A global section of Sp4 -> P3 over the line chart v_y!=0.
    r,s1,s2=s.symbols('r s1 s2');chart=s.Matrix([[1,r,-s2,s1],[0,1,0,0],[0,s1,1,0],[0,s2,0,1]])
    assert chart.T*J[1:5,1:5]*chart==J[1:5,1:5] and chart.det()==1
    assert chart[:,1]==s.Matrix([r,1,s1,s2]);assert chart.inv()==chart.subs({r:-r,s1:-s1,s2:-s2})
    paths=[Path(__file__),ROOT/'code/ambient_weight1_image.py',ROOT/'code/s61_model.py',ROOT/'code/weight1_structure.py',
        ROOT/'code/additional_rank4_families_20261103_v62.py',ROOT/'reports/AMBIENT_FIXED_FIRSTJET_FAMILIES_20261103_v62/report.json']
    result=dict(id='AMBIENT_FIBER_BUNDLE_20261103_v63',all_passed=True,base_parameters=list(map(str,base)),fiber_parameters=list(map(str,fiber)),
        beta_indices_zero_based=beta_indices,q_indices_zero_based=q_indices,parabolic_dimension=7,
        parabolic_lie_matrices=[g.tolist() for g in parabolic],line_chart_section=[[str(v) for v in row] for row in chart.tolist()],
        beta16_to_intrinsic24_adapter=adapter.tolist(),rank3_kernel_witness=witness.tolist(),families=records,
        fixed_line_intersection_dimension=9,intersection_parameter_difference=list(map(str,solution)),
        q_action_proof_dependency='Weights imply normalized Sp4 action is linear in q with beta-quadratic shift; inverse group action makes the q-linear part invertible.',
        geometric_proof_dependency='Connected parabolic stability, closed associated-bundle embedding into P3xS61, proper projection, unique rank3 kernel line and matching dimension19 identify full reduced components.',
        conclusions_scope='Reduced CA/CB components: all rank3-image fibers are Zariski-locally affine10; their scheme intersection has affine9 fibers. No full ambient S61 fiber or rankdrop fiber exhaustion is claimed.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.monotonic()-started)
    def serial(value):
        if isinstance(value,s.Rational):return str(value)
        raise TypeError(type(value))
    (output/'report.json').write_text(json.dumps(result,indent=2,default=serial)+'\n');print('FULL REDUCED COMPONENT INCIDENCE CERTIFICATES PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();main(args.output)
