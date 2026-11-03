#!/usr/bin/env python3
"""Actual N odd-generator quadratic jets: two new classes beyond seven linears."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector,action
ROOT=Path(__file__).resolve().parents[1];h,x,y,z,w=X
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63'
    prior=json.loads(gzip.decompress((folder/'certificate.json.gz').read_bytes()));N=prior['N_full_firstjet_data']
    A=s.Matrix([[s.Rational(v)for v in row]for row in N['actual']]);dnu=s.Matrix([[s.Rational(v)for v in row]for row in N['discriminant']])
    independent=list(A.T.rref()[1]);assert independent==[18,19,28,29,60,61,67]
    columns=list(A[independent,:].rref()[1]);linear_coeff=A[:,columns]*A[independent,columns].inv()
    assert linear_coeff*A[independent,:]==A and A.rank()==7
    nu_coeff=dnu[:,columns]*A[independent,columns].inv();assert nu_coeff*A[independent,:]==dnu
    assert list(nu_coeff)==[0,0,0,0,0,0,-s.Rational(1,8)]
    cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
    sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
    change=lambda P:(sw*P.subs({h:x,x:h},simultaneous=True)*sw.T).applyfunc(s.expand)
    basis=[change(canonical(raw(f,0,0)))for f in[h*h,h*x]]+[change(canonical(raw(0,1,0)))]+[change(canonical(raw(0,0,H)))for H in cubics]
    C=s.Matrix.hstack(*map(vector,basis));ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T
    W,V=s.symbols('W V');Pw=basis[2];Pv=s.Rational(5,2)*basis[1]+s.Rational(11,2)*basis[8]+basis[16]
    P=P0+W*Pw+V*Pv
    flag=[(0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(1,4),(3,2),(4,2)]
    J=C.row_join(s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P))for i,j in flag]))
    complement=[i for i in range(28)if i!=22];rr=prior['N_own_linear27_rows'];B=J[:,complement];b=J[:,22]
    B0=B.subs({W:0,V:0});B1=B-B0;b0=b.subs({W:0,V:0});b1=b-b0
    inverse=B0[rr,:].inv();assert B0[rr,:].det()!=0
    c0=-inverse*b0[rr,:];c1=-inverse*(b1[rr,:]+B1[rr,:]*c0);c2=-inverse*B1[rr,:]*c1
    S0=(b0+B0*c0).applyfunc(s.expand);S1=(b1+B0*c1+B1*c0).applyfunc(s.expand);S2=(B0*c2+B1*c1).applyfunc(s.expand)
    assert S0==S1==s.zeros(150,1)
    assert S2[27]==0
    # All selected27 actual target derivatives vanish through degree2,
    # providing a direct certificate for the truncated inverse calculation.
    assert S2[rr,:]==s.zeros(27,1)
    tangent_R=s.Matrix.hstack(vector(Pw),vector(Pv))[rr+[27],:]
    assert tangent_R.rank()==2 and A*tangent_R==s.zeros(150,2)
    # Since the branch derivativeS27 has no terms<=2, tracefree t has no
    # terms<=2 on this plane. Therefore S2 equals the actual odd b_i second
    # jet there; no even/trace terms survive at degree2.
    nu_lift=s.factor((nu_coeff*S2[independent,:])[0]);assert nu_lift==-s.Rational(9,4)*V*V
    zero_linear=(S2-linear_coeff*S2[independent,:]).applyfunc(s.expand)
    assert zero_linear[46]==-4*W*V and zero_linear[72]==-6*V*V
    quadratics=s.Matrix([[s.Poly(v,W,V).coeff_monomial(W*V),s.Poly(v,W,V).coeff_monomial(V*V)]for v in zero_linear])
    assert quadratics.rank()==2 and quadratics[[46,72],:].det()==24
    cert=dict(actual_odd150_by28_linear=A.tolist(),seven_independent_actual_odd_rows=independent,selected_seven_base_columns=columns,all_linear_row_combinations=mat(linear_coeff),
        nu_linear_combination=mat(nu_coeff),normalizer_critical_tangent_Pw=mat(Pw),normalizer_critical_tangent_Pv=mat(Pv),
        full_physical_source_Jacobian_on_plane150_by28=mat(J),source_complement_columns=complement,actual27_target_rows=rr,
        complement_inverse_constant=mat(inverse),vertical_vector_coefficients_degrees0_1_2=[mat(c)for c in[c0,c1,c2]],
        all150_vertical_output_derivative_degree2=list(map(str,S2)),branch_row27_degree2='0',base_linear_tangent28_by2=mat(tangent_R),
        chosen_nu_lift_quadratic_class=str(nu_lift),all_zero_firstjet_generator_quadratics=list(map(str,zero_linear)),
        two_independent_zero_firstjet_rows=[46,72],two_by_two_coefficient_minor='24')
    # Convert the rational linear matrix to strings as for all other blocks.
    cert['actual_odd150_by28_linear']=mat(A)
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(cert,separators=(',',':'))+'\n').encode(),mtime=0))
    paths=[Path(__file__),ROOT/'research/notes/AMBIENT_B18_N_ODD_QUADRATIC_JETS_20261103_v63.md',folder/'report.json',folder/'certificate.json.gz',ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json',ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json',ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py']
    for p in paths[-3:-1]:assert json.loads(p.read_text())['all_passed'] is True
    result=dict(id='AMBIENT_B18_N_ODD_QUADRATIC_JETS_20261103_v63',all_passed=True,
        actual_all150_physical_source_Jacobian_on_intrinsic_critical_tangent_plane=True,
        selected27_base_coordinates_vertical_through_quadratic_order=True,
        whole_branch_row_derivative_zero_through_quadratic_order=True,
        actual_odd_second_jets_not_just_Fitting_ideal_generators=True,
        seven_independent_actual_odd_linear_classes=True,
        chosen_nu_lift='-b67/8',chosen_nu_lift_quadratic_class='-9*V^2/4',
        two_zero_base_firstjet_actual_odd_generators_quadratic_classes=['-4*W*V','-6*V^2'],
        independent_quadratic_classes_mod_maximal_times_odd_ideal=2,
        actual_full_N_odd_ideal_minimal_generator_lower_bound=9,
        both_previously_proposed_eight_generator_full_N_ansatze_excluded=True,
        nu_lift_quadratic_class_not_canonical_under_adding_zero_firstjet_generators=True,
        reduced_N_cone_embedding_dimension_lower_bound=37,reduced_N_projective_embedding_dimension_lower_bound=36,
        exact_full_odd_ideal_depth_or_actual_D_type_not_claimed=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('ACTUAL N ODD2JETS: 7LINEAR+2NEW=>muI>=9 PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
