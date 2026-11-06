#!/usr/bin/env python3
"""All four zero-residue conic/triangle sources, with all six magnetic terms.

Construct the quadratic logarithmic anchor basis from the full syzygy space,
then test one fixed point of each actual source in all150 projective variables.
No component identification is inferred from a divisor name.
"""
import argparse,hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_rank3_mixed_component_review import COLUMNS,jacobi,sparse,vector,minor_certificate

ROOT=Path(__file__).resolve().parents[1]

def ambient_matrices(P,xx,directions):
    for A in [P,*directions]:
        assert all(c.q==1 for ij in itertools.combinations(range(5),2) for c in s.Poly(A[ij],*xx).coeffs())
    ps=sparse(P,xx);assert jacobi(ps)==s.zeros(350,1)
    columns=[]
    for ij,e in COLUMNS:
        variation={key:dict(value) for key,value in ps.items()}
        variation.setdefault(ij,{})[e]=variation.get(ij,{}).get(e,0)+1
        columns.append(jacobi(variation))
    D=s.Matrix.hstack(*columns);X=s.Matrix(xx);vertical=[];orbit=[]
    for i,j in itertools.product(range(5),repeat=2):
        field=s.zeros(5,1);field[i]=xx[j]
        vertical.append(vector(sparse(X*field.T-field*X.T,xx)))
        A=s.zeros(5);A[i,j]=1
        orbit.append(vector(sparse(A*P+P*A.T-xx[j]*P.diff(xx[i]),xx)))
    V=s.Matrix.hstack(*vertical)
    M=V.row_join(s.Matrix.hstack(*orbit)).row_join(vector(ps))
    F=M.row_join(s.Matrix.hstack(*(vector(sparse(R,xx)) for R in directions)))
    assert D*F==s.zeros(350,F.cols)
    return D,V,M,F

def source(case):
    h,x,y,z,w=xx=s.symbols('h x y z w');base=(h,x,y)
    qs={'m0_t3':h*h+x*x+y*y-2*h*x-2*h*y-2*x*y,
        'm1_t2':(x+y)**2+h*y,
        'm2_t1':h*x+x*y+y*y,
        'm3_t0':h*x+h*y+x*y}
    q=qs[case];f=s.expand(h*x*y*q)
    mons=(h*h,h*x,h*y,x*x,x*y,y*y)
    coefficients=s.symbols('c0:18')
    candidate=[sum(coefficients[6*i+j]*m for j,m in enumerate(mons)) for i in range(3)]
    equation=s.expand(sum(candidate[i]*s.diff(f,base[i]) for i in range(3)))
    A,rhs=s.linear_eq_to_matrix(s.Poly(equation,*base).coeffs(),coefficients)
    assert rhs==s.zeros(A.rows,1)
    kernel=A.nullspace();assert len(kernel)==2
    syzygies=[];anchors=[]
    for column in kernel:
        values=[s.expand(sum(column[6*i+j]*m for j,m in enumerate(mons))) for i in range(3)]
        assert s.expand(sum(values[i]*s.diff(f,base[i]) for i in range(3)))==0
        logarithmic_h=s.cancel(values[0]/h)
        assert s.Poly(logarithmic_h,*base).total_degree()<=1
        anchor=[s.expand(values[1]-x*logarithmic_h),s.expand(values[2]-y*logarithmic_h)]
        assert s.cancel(anchor[0]/x).is_polynomial(*base)
        assert s.cancel(anchor[1]/y).is_polynomial(*base)
        syzygies.append(values);anchors.append(anchor)
    U,V=map(s.Matrix,anchors)
    det=s.expand(U[0]*V[1]-U[1]*V[0]);ratio=s.cancel(det/(x*y*q))
    assert ratio and not ratio.free_symbols
    V=V/ratio
    assert s.expand(U[0]*V[1]-U[1]*V[0]-x*y*q)==0
    Ua=U.subs(h,1);Va=V.subs(h,1)
    bracket=s.Matrix([s.expand(sum(Ua[j]*s.diff(Va[i],(x,y)[j])-Va[j]*s.diff(Ua[i],(x,y)[j]) for j in range(2))) for i in range(2)])
    deta=s.expand(x*y*q.subs(h,1))
    c=s.cancel((bracket[0]*Va[1]-bracket[1]*Va[0])/deta)
    d=s.cancel((Ua[0]*bracket[1]-Ua[1]*bracket[0])/deta)
    assert all(s.Poly(t,x,y).total_degree()<=1 for t in [c,d])
    assert (bracket-c*Ua-d*Va).applyfunc(s.expand)==s.zeros(2,1)
    ch=s.expand(h*c.subs({x:x/h,y:y/h},simultaneous=True))
    dh=s.expand(h*d.subs({x:x/h,y:y/h},simultaneous=True))
    params=s.symbols('r0:6');magnetic=sum(p*m for p,m in zip(params,mons))
    P=s.zeros(5)
    entries={(1,3):U[0],(2,3):U[1],(1,4):V[0],(2,4):V[1],(3,4):-ch*z-dh*w+magnetic}
    for ij,value in entries.items():P[ij]=s.expand(value);P[ij[::-1]]=-P[ij]
    common=s.ilcm(1,*[s.denom(a) for ij in itertools.combinations(range(5),2) for a in s.Poly(P[ij],*xx,*params).coeffs()])
    P=P*common
    universal=[s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l]) for l in range(5))) for i,j,k in itertools.combinations(range(5),3)]
    assert universal==[0]*10
    pf=s.expand(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3])
    assert s.expand(pf+common**2*x*y*q)==0
    meta={'case':case,'quadric':str(q),'plane_quintic':str(f),
          'full_quadratic_syzygy_matrix_shape':list(A.shape),'quadratic_syzygy_rank':16,
          'full_quadratic_syzygy_dimension':2,'two_full_syzygy_vectors':[[str(v) for v in row] for row in syzygies],
          'anchor_U':list(map(str,U)),'anchor_V':list(map(str,V)),
          'anchor_determinant':'x*y*('+str(q)+')','affine_anchor_bracket_coefficients':[str(c),str(d)],
          'all_six_magnetic_monomials':list(map(str,mons)),'source_integer_scale':str(common),
          'universal_all_ten_raw_Jacobi_zero':True,'universal_affine_Pfaffian':str(-common**2*x*y*q),
          'universal_actual_quintic':str(-common**2*h*x*y*q),
          'all_six_parameter_source':{str(ij):str(P[ij]) for ij in itertools.combinations(range(5),2) if P[ij]}}
    return P,xx,params,meta

def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(90,91))
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();records=[]
    for case in ['m0_t3','m1_t2','m2_t1','m3_t0']:
        case_begun=time.monotonic();Pgen,xx,params,record=source(case)
        values=dict(zip(params,(1,2,3,5,7,11)));P=Pgen.subs(values)
        directions=[Pgen.diff(p).subs(values) for p in params]
        D,V,M,F=ambient_matrices(P,xx,directions)
        matrices={'full350_D':D,'vertical':V,'orbit_vertical_scale':M,'actual_full_source':F}
        ranks={name:len(DomainMatrix.from_Matrix(A).convert_to(s.QQ).rref()[1]) for name,A in matrices.items()}
        certificates={name:minor_certificate(A,ranks[name]) for name,A in matrices.items()}
        lower=ranks['actual_full_source']-25;upper=125-ranks['full350_D']
        record.update(point_substitution={str(k):str(v) for k,v in values.items()},
            point_brackets={str(ij):str(P[ij]) for ij in itertools.combinations(range(5),2) if P[ij]},
            all6_actual_parameter_directions=[{str(ij):str(R[ij]) for ij in itertools.combinations(range(5),2) if R[ij]} for R in directions],
            ranks_QQ=ranks,certificates=certificates,
            all_full350_times_actual_directions_zero=True,
            projective_source_dimension_lower=lower,projective_tangent_dimension_upper=upper,
            true_transverse_tangent_dimension=upper-lower,
            smooth_full_component_proved=(lower==upper),known_component_identification_not_performed=True,
            matrix_sha256={name:hashlib.sha256(json.dumps(A.tolist(),default=str).encode()).hexdigest() for name,A in matrices.items()},
            seconds=time.monotonic()-case_begun)
        records.append(record)
        print(json.dumps({'case':case,'ranks':ranks,'lower':lower,'upper':upper,'seconds':record['seconds']}),flush=True)
    paths=['code/classification_three_hyperplanes_quadric_independent.py','code/classification_rank3_mixed_component_review.py',
           'research/classification/THREE_HYPERPLANES_QUADRIC_INDEPENDENT.md']
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,
       'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
       'source_snapshots':[{'path':p,'sha256':hashlib.sha256((ROOT/p).read_bytes()).hexdigest(),'source':(ROOT/p).read_text()} for p in paths[:2]],
       'scope':'Four complete zero-residue anchor types with all six magnetic coefficients. One fixed full-ambient test per type; no new-component identification or whole classification inferred without the recorded dimension and containment proofs.',
       'global_exhaustion_proved':False,'seconds':time.monotonic()-begun}
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
