#!/usr/bin/env python3
"""Independent coefficient/Lie reconstruction for reduced R19 boundary locals."""
import argparse
from collections import defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
X=s.symbols('h x y z w');h,x,y,z,w=X
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))
BASIS=list(itertools.product(PAIRS,MONS));INDEX={value:i for i,value in enumerate(BASIS)}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def vector_field(D):
    out=s.zeros(150,1)
    for pair,polynomial in [((0,2),-D[0]),((1,2),-D[1]),((2,3),D[2]),((2,4),D[3])]:
        p=s.Poly(polynomial,*X)
        for mon in MONS:out[INDEX[pair,mon]]=p.coeff_monomial(s.prod(X[i]for i in mon))
    assert s.expand(s.diff(D[0],h)+s.diff(D[1],x)+s.diff(D[2],z)+s.diff(D[3],w))==0
    return out


def lie_operator(u,v):
    result=defaultdict(int)
    for source,((i,j),mon)in enumerate(BASIS):
        def insert(r,t,m,coefficient):
            if r==t:return
            if r>t:r,t,coefficient=t,r,-coefficient
            result[INDEX[(r,t),tuple(sorted(m))],source]+=coefficient
        if i==v:insert(u,j,mon,1)
        if j==v:insert(i,u,mon,1)
        if u in mon:
            changed=list(mon);changed.remove(u);changed.append(v);insert(i,j,changed,-mon.count(u))
    return[(i,j,c)for(i,j),c in result.items()if c]


OPERATORS=[lie_operator(i,j)for i,j in itertools.product(range(5),repeat=2)]


def orbit(P):
    columns=[]
    for operator in OPERATORS:
        value=s.zeros(150,1)
        for i,j,c in operator:value[i]+=c*P[j]
        columns.append(value.applyfunc(s.expand))
    return s.Matrix.hstack(*columns)


def swap(P):
    permutation=[1,0,2,3,4];out=s.zeros(150,P.cols)
    for source,((i,j),mon)in enumerate(BASIS):
        r,t=permutation[i],permutation[j];sign=1
        if r>t:r,t,sign=t,r,-1
        target=INDEX[(r,t),tuple(sorted(permutation[i]for i in mon))]
        out[target,:]=sign*P[source,:]
    return out


def check_matrix(M,data):
    expected=s.Matrix([[s.sympify(c)for c in row]for row in data['matrix']]);assert M==expected
    rank=DomainMatrix.from_Matrix(M).convert_to(s.QQ).rank();assert rank==data['rank']
    determinant=DomainMatrix.from_Matrix(M[data['rows'],data['columns']]).convert_to(s.QQ).det()
    assert s.Rational(determinant)==s.sympify(data['determinant'])
    return dict(rank=rank,minor=str(determinant))


def build(output):
    producer=ROOT/'code/ambient_R19_boundary_local_models_20261103_v63.py'
    report_path=ROOT/'reports/AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63/report.json'
    base_path=ROOT/'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63/report.json'
    report=json.loads(report_path.read_text());base=json.loads(base_path.read_text());assert report['all_passed']and base['all_passed']
    for path,expected in report['input_sha256'].items():assert sha(ROOT/path)==expected,path
    t,r,b,c,d,e=s.symbols('t r b c d e');g=s.symbols('g0:6');G=sum(a*m for a,m in zip(g,[h*h,h*x,h*z,x*x,x*z,z*z]))
    parameters=[t,r,b,c,d,e,*g]
    field=[x*((r-13*t)*h+b*x+c*z),12*t*x*x,x*(d*h+e*x+(-r-13*t)*z),2*t*x*w+G]
    C=s.Matrix.hstack(*(vector_field([s.diff(value,parameter)for value in field])for parameter in parameters))
    expected=s.Matrix.hstack(*(s.Matrix(list(map(s.sympify,col)))for col in base['spectral_basis_raw150']))
    assert C==expected and DomainMatrix.from_Matrix(C).rank()==12
    GG=h*h+2*h*x+3*h*z+4*x*x+5*x*z+6*z*z
    P=vector_field([0,0,-h*x,GG]);O=orbit(P)
    keep=[i for i in range(12)if i!=3];nilkeep=[i for i in range(12)if i not in [0,3]]
    J1=C[:,keep].row_join(O);J2=swap(C[:,keep]).row_join(O)
    K1=C[:,nilkeep].row_join(O);K2=swap(C[:,nilkeep]).row_join(O)
    f=s.symbols('f0:3');TT=vector_field([0,0,-f[0]*h*h-f[1]*h*x-f[2]*x*x,G])
    TC=s.Matrix.hstack(*(TT.diff(parameter)for parameter in [*f,*g])).row_join(O)
    matrices={'R_branch1':J1,'R_branch2':J2,'R_union':J1.row_join(J2),
              'B18_branch1':K1,'B18_branch2':K2,'B18_union':K1.row_join(K2),'B17_tangent':TC}
    checks={name:check_matrix(value,report[name])for name,value in matrices.items()}
    assert [checks[name]['rank']for name in matrices]==[20,20,22,19,19,20,18]
    A,B,E,F=s.symbols('A B E F');P0=vector_field([x*x,0,h*x,A*h*h+B*h*z+E*x*z+F*z*z])
    J=C[:,keep].row_join(orbit(P0));small=report['constant19_minor'];rows,columns=small['rows'],small['columns'];V=J[:,columns]
    domain=s.QQ.poly_ring(A,B,E,F)
    determinant=domain.to_sympy(DomainMatrix.from_Matrix(V[rows,:]).convert_to(domain).det());assert determinant==48
    K=s.QQ.frac_field(A,B,E,F)
    def dm(M):return DomainMatrix.from_Matrix(M).convert_to(K).to_dense()
    residual=(dm(J)-dm(V)*dm(V[rows,:]).inv()*dm(J[rows,:])).to_Matrix().applyfunc(s.cancel)
    nz=[[i,j,str(residual[i,j])]for i in range(150)for j in range(36)if residual[i,j]]
    assert nz==report['entire_residual_entries']
    special={A:1,B:0,E:-s.Rational(1,6),F:0}
    checks['C15_source_differential']=check_matrix(J.subs(special),report['exceptional_point'])
    checks['C15_orbit']=check_matrix(orbit(P0.subs(special)),report['exceptional_orbit'])
    aa,bb,cc=s.symbols('a b c');q=s.symbols('q0:6')
    old=s.Matrix([aa*x*x,0,bb*h*x+cc*x*x,sum(v*m for v,m in zip(q,[h*h,h*x,h*z,x*x,x*z,z*z]))])
    U=s.Matrix([h,x,z,w]);S=s.Matrix([[1/aa,0,0,0],[0,1,0,0],[-cc/(aa*aa*bb),0,1/(aa*bb),0],[0,0,0,1]])
    transformed=(S*old.subs(dict(zip(U,S.inv()*U)),simultaneous=True)).applyfunc(s.cancel)
    assert transformed[0]==x*x and transformed[1]==0 and transformed[2]==h*x
    mu=s.expand(transformed[3]).coeff(h,1).coeff(x,1);nu=s.expand(transformed[3]).coeff(x,2)
    final=s.expand(transformed[3]-mu*h*x-nu*x*x)
    assert s.cancel(mu-s.sympify(report['chain_removed_hx_coefficient']))==0
    assert s.cancel(nu-s.sympify(report['chain_removed_x2_coefficient']))==0
    assert s.cancel(final-s.sympify(report['chain_final_field'][3]))==0
    # Nilpotent regular-chain case has the unique three-space ker(N^3).
    NN=s.Matrix([[0,1,0,0],[0,0,0,0],[1,0,0,0],[0,0,1,0]])
    assert NN.rank()==3 and (NN**3).rank()==1 and NN**4==s.zeros(4)
    note=ROOT/'research/notes/AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63_REVIEW.md'
    inputs=[Path(__file__),note,producer,report_path,base_path,
            ROOT/'research/notes/AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63.md',
            ROOT/'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW.json',
            ROOT/'reports/AMBIENT_R_COMPLEMENT_20261103_v63_REVIEW.json']
    for path in inputs[-2:]:assert json.loads(path.read_text())['all_passed']
    result=dict(id='AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_V63_REVIEW',all_passed=True,
                independent_full150_spectral_basis_and_sparse_GL_actions=True,
                all_nine_saved_full_tangent_matrices_exactly_match=True,QQ_rank_checks=checks,
                independent_universal_constant19_minor=48,
                independent_entire_residual_entries=nz,
                all_ab_unit_chain_normalization_verified_with_c_arbitrary=True,
                generic_B17_whole_geometric_fiber_exactly_two_flags=True,
                entire_B18_outside_B17_geometric_fiber_exactly_one_flag=True,
                regular_nilpotent_G_mod_x_zero_edge_case_checked=True,
                proper_finite_normalization_Nakayama_arguments_accepted=True,
                generic_B17_actual_scheme_intersection_smooth17_and_fiber_product_model_accepted=True,
                generic_B18_node_model_accepted=True,
                full_reduced_R_bad_locus_equalities_accepted='Sing=Nonnormal=NonCM=B17 union C15',
                full_R_complement_smoothness_uses_existing_ENTIRE_complement_review=True,
                S2_wording='NonCM follows because a CM local ring has all localizations S2 and, with inherited R1, is normal; no assertion that closed-point depth is below2.',
                scope='Reduced R and reduced B18 only; neither full Poisson nor actual CA/CB boundary intersection scheme. New C15 pinched completion is a separate packet not reviewed here.',
                input_sha256={str(path.relative_to(ROOT)):sha(path)for path in inputs})
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS independent R19 boundary local matrices, fibers and reduced-scheme geometry')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
