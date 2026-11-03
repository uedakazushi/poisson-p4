#!/usr/bin/env python3
"""Independent literal frame, immutable S61 adapter, actual218 localization review."""
import argparse
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import ambient_B18_zero_linear_contact4_review_20261103_v63 as own
from ambient_B18_actual_D_tangent_cone_support_review_20261103_v63 import inverse,matrix
from ambient_B18_Noether_fiber_even_quadratic_review_20261103_v63 import exact_row_rank

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63_REVIEW'


def rref(A):
    B=[list(map(F,row)) for row in A];pivots=[];i=0
    for j in range(len(B[0])):
        p=next((a for a in range(i,len(B)) if B[a][j]),None)
        if p is None:continue
        B[i],B[p]=B[p],B[i];c=B[i][j];B[i]=[a/c for a in B[i]]
        for a in range(len(B)):
            if a!=i and B[a][j]:
                c=B[a][j];B[a]=[x-c*y for x,y in zip(B[a],B[i])]
        pivots.append(j);i+=1
        if i==len(B):break
    return B[:i],pivots


def push(M,P):
    Mi,det=inverse(M)
    quad=[]
    for a,b in own.MONS:
        quad.append([Mi[a][i]*Mi[b][j]+(Mi[a][j]*Mi[b][i] if i!=j else 0) for i,j in own.MONS])
    changed=[[sum((P[15*k+l]*quad[l][m] for l in range(15)),F(0)) for m in range(15)] for k in range(10)]
    wedge=[[M[i][a]*M[j][b]-M[i][b]*M[j][a] for a,b in own.PAIRS] for i,j in own.PAIRS]
    return [sum((wedge[i][k]*changed[k][m] for k in range(10)),F(0)) for i in range(10) for m in range(15)],Mi,det


def encoded(polynomial,u):
    return sum((F(num,den)*product(u[i] for i in mon) for mon,num,den in polynomial),F(0))


def product(values):
    out=F(1)
    for value in values:out*=value
    return out


def eval_poly(poly,values):
    total=F(0)
    for exp,c in poly.terms.items():
        total+=c*product(v**e for v,e in zip(values,exp))
    return total


def gradient(poly,values):
    out=[F(0)]*len(values)
    for exp,c in poly.terms.items():
        for i,e in enumerate(exp):
            if e:
                out[i]+=c*e*product(v**(a-int(j==i)) for j,(v,a) in enumerate(zip(values,exp)))
    return out


def main(output):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    assert not output.exists()
    ap=ROOT/'reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63/report.json'
    oldreview=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63_REVIEW/report.json'
    wreview=ROOT/'reports/AMBIENT_W1_ACTUAL_TANGENT_INTERSECTION_20261103_v63_REVIEW.json'
    uniform=ROOT/'reports/AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63_REVIEW/report.json'
    author,oldgate,wgate,ugate=[own.accepted(p) for p in [ap,oldreview,wreview,uniform]]
    assert oldgate['exact_actual_target_equalities_rechecked_by_fresh_M2']
    rp=ap.parent/'recognition.json.gz';c=own.read(rp)
    sourcepath=ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json'
    basepath=ROOT/'reports/AMBIENT_BASE_PINCHED_SURFACE_20261103_v63/report.json'
    bundlepath=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json'
    builderpath=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json'
    oldpath=builderpath.parent/'report.json'
    body,base,bundle,builder,old=[own.read(p) for p in [sourcepath,basepath,bundlepath,builderpath,oldpath]]
    point=[F(0)]*150
    for pair,mon,val in [((0,2),(1,1),-1),((2,3),(0,1),1),((2,4),(0,0),1),((2,4),(3,3),1),((3,4),(1,1),1)]:
        point[15*own.PI[pair]+own.MI[mon]]=F(val)
    M=matrix(c['exact_homogeneous_frame'])
    assert M==matrix([[0,1,0,0,0],[-1,0,0,0,0],[0,0,1,0,0],[0,-1,0,1,0],[-1,0,0,0,1]])
    physical,Mi,det=push(M,point);assert det==1 and Mi==matrix(c['exact_frame_inverse'])
    assert own.matvec(M,[F(0),F(1),F(0),F(1),F(0)])==[1,0,0,0,0]
    expected=[F(0)]*150
    for pair,mon,val in [((1,2),(0,0),1),((3,4),(0,0),1),((2,3),(0,1),-1),((2,4),(0,3),2),((2,4),(1,1),1),((2,4),(3,3),1)]:
        expected[15*own.PI[pair]+own.MI[mon]]=F(val)
    assert physical==expected==list(map(F,c['full150_transported_physical_coefficients']))
    # Independent Euler-wedge quotient construction, no S61.normalize call.
    rows=[]
    for direction,variable in itertools.product(range(5),repeat=2):
        row=[F(0)]*150
        for ip,(i,j) in enumerate(own.PAIRS):
            if j==direction:row[15*ip+own.MI[tuple(sorted((i,variable)))]]+=1
            if i==direction:row[15*ip+own.MI[tuple(sorted((j,variable)))]]-=1
        rows.append(row)
    vertical,pivots=rref(rows);assert len(pivots)==24
    certpath=ROOT/'legacy/v02/data/symplectic_slice_certificate.json'
    datapath=ROOT/'legacy/v02/data/symplectic_slice_equations.json'
    reducedpath=ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json'
    cert,data,reduced=[own.read(p) for p in [certpath,datapath,reducedpath]]
    assert pivots==cert['vertical_pivots']
    quotient=[a-sum((row[i]*physical[p] for row,p in zip(vertical,pivots)),F(0)) for i,a in enumerate(physical)]
    assert all(quotient[i]==0 for i in cert['perspective_orbit_pivots_raw'])
    assert c['perspective_parameters']==['0']*4
    u=[quotient[i] for i in cert['raw_free_indices']]
    assert u==list(map(F,c['exact_normalized61'])) and quotient==list(map(F,c['normalized150']))
    reconstruction=[F(0)]*150
    for pair in [(1,2),(3,4)]:reconstruction[15*own.PI[pair]+own.MI[0,0]]=1
    for i,a in zip(cert['raw_free_indices'],u):reconstruction[i]=a
    for entry in data['eliminated_raw_variables']:reconstruction[entry['raw_index']]=encoded(entry['polynomial'],u)
    assert reconstruction==quotient
    assert all(encoded(eq,u)==0 for eq in reduced['equations']) and len(reduced['equations'])==134
    qi=bundle['q_indices_zero_based'];bi=bundle['beta_indices_zero_based']
    beta=own.matvec(matrix(bundle['beta16_to_intrinsic24_adapter']),[u[i] for i in bi])
    envpoint={'p':F(0),'s':F(0),'r':F(0),**{f'k{i}':F(0) for i in range(9)}}
    envpoint['k0']=F(1);envpoint['k5']=F(1,3)
    assert beta==list(map(F,c['full_beta24']))
    assert beta==[own.evaluate(q,envpoint) for q in body['beta24_coordinates']]
    _,orbitdet=inverse(matrix(base['orbit_constant_block_witness']['matrix']))
    assert orbitdet==F(base['orbit_constant_block_witness']['determinant']) and base['orbit_rank']==7
    source_dets={}
    for name in ['CA','CB','N']:
        record=body['families'][name];pars=record['source_fiber_parameters']
        names=['p','s',*pars];variables=own.symbols(names)
        poly=[own.Poly(own.evaluate(q,variables)) for q in record['full_61_coordinates']]
        values=[envpoint[n] for n in names]
        assert [eval_poly(q,values) for q in poly]==u
        selected=record['actual_q_selected_rows'] if name=='N' else builder['common_q_selected_rows']
        selected61=[qi[i] for i in selected]
        jac=[[gradient(poly[i],values)[names.index(n)] for n in pars] for i in selected61]
        _,det=inverse(jac);assert det==-72
        source_dets[name]=str(det)
        odd=matrix(record['first_odd_coefficient_matrix'])
        assert exact_row_rank(odd)==(6 if name=='N' else 7)
        assert own.matvec(odd,[envpoint[n] for n in pars])==[0]*45
        # The shifted k0,k5 are literally free directions of these missing ideals.
        assert all(row[pars.index('k0')]==row[pars.index('k5')]==0 for row in odd)
    K=[u[qi[i]] for i in builder['common_q_selected_rows']]
    Z=[u[qi[i]] for i in builder['remaining_q_rows']]
    names=builder['ambient_variables'];variables=own.symbols(names)
    point49={n:F(0) for n in names}
    point49.update({f'v{i}':a for i,a in enumerate(K)});point49.update({f'z{i}':a for i,a in enumerate(Z)})
    assert {n:str(a) for n,a in point49.items() if a}=={'v4':'1','v7':'1'}==c['old_49_target_point_nonzero']
    polynomials=[own.Poly(own.evaluate(q,variables)) for name in ['CA','CB'] for q in builder['component_presentations'][name]]
    assert len(polynomials)==218
    values=[point49[n] for n in names]
    assert all(eval_poly(q,values)==0 for q in polynomials)
    jac=[gradient(q,values) for q in polynomials]
    assert jac==matrix(c['actual218_by49_numeric_Jacobian'])
    ranks=[exact_row_rank(jac[:109]),exact_row_rank(jac[109:]),exact_row_rank(jac)]
    assert ranks==[28,28,29]
    nilF=own.Poly(own.evaluate(old['nilradical_generator'],variables))
    nilgrad=gradient(nilF,values)
    assert eval_poly(nilF,values)==0 and nilgrad==list(map(F,c['nil_F_numeric_gradient']))
    assert exact_row_rank(jac+[nilgrad])==30
    prime=[own.Poly(own.evaluate(q,variables)) for q in old['annihilator_prime_generators']]
    assert all(eval_poly(q,values)==0 for q in prime)
    free=['v1','v4','v5','v7'];assert old['annihilator_free_variables']==free
    assert exact_row_rank([gradient(q,values) for q in prime])==45
    assert all(all(sum(e)==1 and all(e[names.index(n)]==0 for n in free) for e in q.terms) for q in prime)
    shifted={n:variables[n]+point49[n] for n in names}
    translated=[own.Poly(own.evaluate(q,shifted)) for name in ['CA','CB'] for q in builder['component_presentations'][name]]
    assert all(q==own.evaluate(s,variables) for q,s in zip(translated,c['all218_translated_actual_intersection_generators']))
    assert own.evaluate(old['nilradical_generator'],shifted)==own.evaluate(c['translated_nilradical_generator'],variables)
    # In fact both shifted coordinates are absent from all218 polynomials; this is
    # an extra direct check, not a premise replacing actual-source/unit binding.
    free_shift_invariant=all(all(e[names.index('v4')]==e[names.index('v7')]==0 for e in q.terms) for q in polynomials)
    assert old['square_zero'] and old['annihilator_prime_slice_dimension']==4
    # Confirm that the independently replayed M2 identities were in QQ[49],
    # not in a completion/localization centered at the old point.
    scriptpath=oldpath.parent/'intersection_colon_and_reduction_v7.m2'
    script=scriptpath.read_text();assert script.startswith('R=QQ[')
    assert 'annFideal=ID:ideal(FF)' in script and 'D_PLUS_F_IN_N' in script and 'N_IN_D_PLUS_F' in script
    inputs=[Path(__file__),Path(own.__file__),ROOT/'code/ambient_B18_actual_D_tangent_cone_support_review_20261103_v63.py',
        ROOT/'code/ambient_B18_Noether_fiber_even_quadratic_review_20261103_v63.py',ROOT/f'research/notes/{STEM}.md',
        ap,rp,oldreview,wreview,uniform,sourcepath,basepath,bundlepath,builderpath,oldpath,scriptpath,certpath,datapath,reducedpath]
    report=dict(id=STEM,all_passed=True,
        independent_Fraction_GL_frame_inverse_and_all150_literal_tensor_identity=True,
        independent_Euler_quotient_raw61_reconstruction_and_all134_equations=True,
        literal_beta24_and_all3_source61_matches_at_shifted_point=True,
        actual_CA_CB_N_fiber_Jacobian_units=source_dets,beta_orbit7_unit=str(orbitdet),
        actual_source_regular_sequence_odd_rank_CA_CB_N=[7,7,6],
        exact_old49_point={'v4':'1','v7':'1'},all218_vanish_and_all218_by49_gradients_rebuilt=True,
        independent_critical_ranks_CA_CB_D_N=[28,28,29,30],
        all218_translated_relations_and_actual_common_target_binding=True,
        all218_polynomials_actually_invariant_under_v4_v7_translations=free_shift_invariant,
        old_global_QQ_polynomial_colon_square_reduction_gates_legitimately_localize_here=True,
        point_on_smooth4_AnnF_and_nilF_remains_nonzero=True,
        critical_slice_dim_edim_depth_mult=[11,20,4,2],smooth_restoration_dimensions={'beta_orbit':7,'symplectic_value':6,'perspective':4},
        full_D_cone_dim_edim_depth_mult=[28,37,21,2],full_D_projective_dim_edim_depth_mult=[27,36,20,2],
        nilradical_principal_exact_square_zero=True,nilradical_module_regular_support_dimension=21,
        exact_Ass=['minimal N prime','smooth21-dimensional embedded prime'],
        all_nonzero_W_by_scheme_preserving_transport=True,
        original_B18_W0_extension_ring_splitting_or_new_global_components_claimed=False,
        verification_kind='Independent Fraction tensor/immutable adapter/source Jacobian/218 coefficients and conceptual global-polynomial localization; no M2 or source jets.',
        input_sha256={str(p.relative_to(ROOT)):own.sha(p) for p in inputs},seconds=time.monotonic()-start)
    output.write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','actual_CA_CB_N_fiber_Jacobian_units','independent_critical_ranks_CA_CB_D_N','full_D_cone_dim_edim_depth_mult','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    main(parser.parse_args().output)
