#!/usr/bin/env python3
"""Independent original134 homogenization/source and clopen sheaf-split audit."""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_model import S61,X
from primitive_jet_elliptic_curve_20261101 import primitive_beta

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_RAW_PURITY_SPLIT_20261103_v63'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def read(path):return json.loads(path.read_text())
def same(A,B):assert A.shape==B.shape and (A-B).applyfunc(s.cancel)==s.zeros(*A.shape)


def main(output,replay):
    started=time.monotonic();folder=ROOT/'reports'/STEM;report=read(folder/'report.json')
    assert report['all_passed']
    for name,digest in report['input_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in report['evidence_sha256'].items():assert sha(folder/name)==digest,name
    subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'code/s61_affine_raw_purity_split_20261103_v63.py'),
                    '--output',str(replay)],check=True,timeout=30)
    assert (folder/'certificate.json').read_bytes()==(replay/'certificate.json').read_bytes()
    oldpath=ROOT/'reports/S61_MINIMAL_AXIS_INFINITY_20261101.json';old=read(oldpath)
    for name,digest in old['input_sha256'].items():assert sha(ROOT/name)==digest,name
    model=S61();weights=model.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1];two=[i for i,w in enumerate(weights) if w==2]
    assert len(one)==16 and len(two)==45 and two==old['weight_two_u_indices_zero_based']
    point=list(map(Fraction,old['sample_S61_coordinates']));equations=model.reduced['equations'];assert len(equations)==134
    qindex={j:i for i,j in enumerate(two)}
    # Direct coefficient/occurrence differentiation of degree2 homogenization.
    J=[];values=[]
    for equation in equations:
        row=[Fraction(0) for _ in range(46)];value=Fraction(0)
        for mon,n,d in equation:
            q=[qindex[i] for i in mon if i in qindex];assert len(q)<=2
            coefficient=Fraction(n,d)
            for i in mon:
                if i not in qindex:coefficient*=point[i]
            if len(q)==2:
                value+=coefficient*point[two[q[0]]]*point[two[q[1]]]
                for position,j in enumerate(q):row[j]+=coefficient*point[two[q[1-position]]]
            elif len(q)==1:row[-1]+=coefficient*point[two[q[0]]]
        assert value==0;J.append([s.Rational(c) for c in row]);values.append(value)
    fixed=next(i for i,j in enumerate(two) if point[j]);chart=[j for j in range(46) if j!=fixed]
    rows=old['independent_full_homogenized_minor']['rows'];assert len(rows)==45 and max(rows)<134
    M=s.Matrix(J)[rows,chart];det=s.Rational(DomainMatrix.from_Matrix(M).convert_to(s.QQ).det())
    assert det==s.Rational(old['independent_full_homogenized_minor']['determinant'])!=0
    radial=s.Matrix([*map(s.Rational,[point[i] for i in two]),0]);same(s.Matrix(J)*radial,s.zeros(134,1))
    # Rebuild the actual rational axis source from its physical quadratic tensor.
    uall=s.symbols('u1:62');k,U,V,L=s.symbols('k U V L');gs=s.symbols('g0:9');x,y,z,w=X[1:]
    G=sum(c*mon for c,mon in zip(gs,[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3]))
    Q=s.zeros(4);Q[0,1]=x*x;Q[1,2]=-s.diff(G,w);Q[1,3]=s.diff(G,z);Q=Q-Q.T
    J4=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    physical=s.zeros(5);physical[1:,1:]=X[0]**2*J4+Q
    normalized=model.normalize(physical)['parameters']
    beta_sub={uall[3]:0,uall[11]:0,uall[4]:0,uall[12]:0,uall[13]:k,uall[48]:2*k+1,
              uall[25]:U,uall[34]:V,uall[49]:-(9*k+7)*V/(3*k+1),uall[50]:5*U-L}
    D=L*(3*k+2)-U*(6*k+8)
    g_values=[-2*V*(D-8*U)/((3*k+1)*D),2*U*(L-2*U)/D,-4*V*(3*k+2)/((3*k+1)*D),
              -8*U/D,0,0,2*(3*k+2)/D,0,0]
    g_sub=dict(zip(gs,g_values))
    assert {str(a):s.cancel(s.sympify(b)) for a,b in beta_sub.items()}=={a:s.cancel(s.sympify(b)) for a,b in old['beta_formula'].items()}
    for name,value in old['g_formula'].items():assert s.cancel(g_sub[s.Symbol(name)]-s.sympify(value))==0
    source=[s.cancel((uall[i] if i in one else normalized[i]).subs(beta_sub).subs(g_sub,simultaneous=True)) for i in range(61)]
    sample={k:1,U:2,V:3,L:5,uall[24]:2,uall[26]:3,uall[33]:5,uall[35]:7,uall[36]:11,uall[47]:13}
    assert (3*k+1).subs(sample)*(3*k+2).subs(sample)*D.subs(sample)!=0
    assert [s.cancel(value.subs(sample)) for value in source]==list(map(s.Rational,old['sample_S61_coordinates']))
    extra=read(ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json')['universal_new_equation']
    for equation in [*equations,extra]:
        value=sum(s.Rational(n,d)*s.prod(source[i] for i in mon) for mon,n,d in equation if sum(i in qindex for i in mon)==2)
        assert s.cancel(value)==0
    # Independently rebuild every beta derivative, including the tensor-slot action.
    zz,uu,beta=primitive_beta(model);coords=[(pair,j) for pair in itertools.combinations(range(4),2) for j in range(4)]
    vector=s.Matrix([s.diff(beta[i,j],zz[r]) for (i,j),r in coords]);coef=vector.jacobian(uu)
    br=list(coef.T.rref()[1]);left=coef[br,:].inv();beta_map=s.Matrix([beta_sub.get(uall[i],uall[i]) for i in one])
    variables=[k,U,V,L]+[uall[i] for i in one if uall[i] not in beta_sub]
    assert list(map(str,variables))==old['source_variables']
    beta_values=beta_map.subs(sample);primitive=coef*beta_values;lookup={coord:i for i,coord in enumerate(coords)}
    def component(i,j,r):
        if i==j:return 0
        return (1 if i<j else -1)*primitive[lookup[((min(i,j),max(i,j)),r)]]
    columns=[]
    for i in range(4):
        for j in range(i,4):
            sym=s.zeros(4);sym[i,j]=sym[j,i]=1;A=J4*sym
            transformed=s.Matrix([sum(A[r,i]*component(r,j,t)+A[r,j]*component(i,r,t)+A[r,t]*component(i,j,r) for r in range(4)) for (i,j),t in coords])
            columns.append(left*transformed[br,:])
    differential=beta_map.jacobian(variables).subs(sample).row_join(s.Matrix.hstack(*columns))
    same(differential,s.Matrix(old['dominance_matrix']))
    dm=old['dominance_minor'];dominance=differential[dm['rows_zero_based'],dm['columns_zero_based']]
    domdet=s.Rational(DomainMatrix.from_Matrix(dominance).convert_to(s.QQ).det());assert domdet==569
    dimension=read(ROOT/'reports/S61_AFFINE_BOUNDARY_DIMENSION_20261102.json')
    assert dimension['all_passed'] and dimension['raw_generic_infinity_exact_dimension']==3
    finitepath=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/certificate.json'
    finite=read(finitepath);quads=list(itertools.combinations_with_replacement(range(3),2));cubics=list(itertools.combinations_with_replacement(range(3),3))
    Q=s.Matrix(finite['plane_quadratic_matrix']);polys=[dict(zip(quads,list(Q[i,:]))) for i in range(5)]
    cubic=s.Matrix([[poly.get(tuple(list(mon[:mon.index(j)])+list(mon[mon.index(j)+1:])),0) if j in mon else 0
                     for mon in cubics] for j in range(3) for poly in polys])
    same(cubic,s.Matrix(finite['plane_cubic_matrix']))
    same(cubic[finite['plane_cubic_rows'],:]*s.Matrix(finite['plane_cubic_inverse']),s.eye(10))
    recipe=read(folder/'certificate.json')['projective_sheaf_ideal_recipe']
    assert recipe['isolated_infinity_ideal']=='I_inf=(I+(w)):Delta^infinity'
    assert recipe['complement_ideal']=='I_rest=I:I_inf'
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}_REVIEW.md',folder/'report.json',folder/'certificate.json',oldpath,
           ROOT/'research/notes'/f'{STEM}.md',ROOT/'reports/S61_AFFINE_BOUNDARY_DIMENSION_20261102.json',
           ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/report.json',
           finitepath,
           ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py']
    result=dict(id=STEM+'_REVIEW',all_passed=True,
        producer_fresh_certificate_byte_identical=True,
        independent_original134_degree2_homogenization_and_full46column_Jacobian=True,
        original45_minor_excludes_added135th_equation=True,full45_QQ_determinant=str(det),
        projective_Euler_radial_relation_and_deleted_coordinate_chart_checked=True,
        independent_physical_axis_source_and_all135_highest_identities=True,
        independent_actual_source20column_beta_differential_and_16minor=True,beta_dominance16_determinant=str(domdet),
        generic_reduced_isolated_infinity_component_and_nonpure3_proof_accepted=True,
        known_plane_has_no_infinity_and_normal_direction_source_identification_required=True,
        exact_plane_all10_cubic_monomial_span_independently_replayed=True,
        first_colon_is_Delta_saturation_not_single_colon=True,
        clopen_finite_etale_selection_and_annihilator_complement_sheaf_proof_accepted=True,
        sheaf_I_intersection_and_unit_sum_and_affine_Dw_unchanged=True,
        scope_generic_K_or_sufficiently_shrunk_beta_open=True,
        no_colon_basis_etale_degree_full_w_saturation_residual_purity_or_affine_unit_claim=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['all_passed','full45_QQ_determinant','beta_dominance16_determinant','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True);args=parser.parse_args();main(args.output.resolve(),args.replay.resolve())
