#!/usr/bin/env python3
"""All actual CA/CB odd2jets and the missingUV class: nu is NOT inI."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector,action
from ambient_CB_global_flag_model_20261103_v63 import raw as rawB
ROOT=Path(__file__).resolve().parents[1];h,x,y,z,w=X
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63';prior=json.loads(gzip.decompress((folder/'certificate.json.gz').read_bytes()))
    W,U,V=s.symbols('W U V');coordinates=[W,U,V];zero=dict.fromkeys(coordinates,0)
    cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
    sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
    change=lambda P:(sw*P.subs({h:x,x:h},simultaneous=True)*sw.T).applyfunc(s.expand)
    AA=[change(canonical(raw(f,0,0)))for f in[h*h,h*x,x*x]]+[change(canonical(raw(0,1,0)))]+[change(canonical(raw(0,0,H)))for H in cubics]
    FF=[change(canonical(rawB(f,0,0)))for f in[h*h,h*x,x*x]];GG=[change(canonical(rawB(0,f,0)))for f in[h*h,h*x,x*x]];HH=[change(canonical(rawB(0,0,H)))for H in cubics]
    flag=[(0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(1,4),(3,2),(4,2)]
    ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T
    free=[1,2,3,9,10,14,15,17,18,19];Y=prior['common_Y28'];branch=61;records={}
    def homogeneous(M,degree):
        return M.applyfunc(lambda f:s.expand(sum(c*s.prod(q**i for q,i in zip(coordinates,e))for e,c in s.Poly(f,*coordinates).terms()if sum(e)==degree)))
    for name in ['CA','CB']:
        earlier=prior['CA_CB_full_firstjet_data'][name];A=s.Matrix(earlier['actual']);dnu=s.Matrix(earlier['discriminant'])
        rr=list(A.T.rref()[1]);cc=list(A[rr,:].rref()[1]);linear=A[:,cc]*A[rr,cc].inv();assert linear*A[rr,:]==A and len(rr)==7
        nu_linear=dnu[:,cc]*A[rr,cc].inv();assert nu_linear*A[rr,:]==dnu
        if name=='CA':
            s3=(3*U+V)/10;s7=(7*U-V)/20;coeff=[(V-U)/10,0,W,s3,0,0,0,s7,0,0]
            P=P0+sum((q*AA[i]for q,i in zip(coeff,free)),s.zeros(5));C=s.Matrix.hstack(*map(vector,AA))
        else:
            s3=(U+V)/2;s7=(U-V)/2;coeff=[-7*U/9+V,0,W,s3,0,0,0,s7,0,0]
            alpha=list(map(s.sympify,[-1,0,0,0,0,0,0,-1,0,0,0,1,0,0,0,0,s.Rational(1,3),0,0,0]))
            for i,q in zip(free,coeff):alpha[i]+=q
            q=alpha[3]
            P=sum((alpha[i]*FF[i]for i in range(3)),s.zeros(5))-q*sum((alpha[i]*GG[i]for i in range(3)),s.zeros(5))+sum((alpha[i+4]*HH[i]for i in range(16)),s.zeros(5))
            C=s.Matrix.hstack(*[vector(FF[i]-q*GG[i])for i in range(3)],-sum((alpha[i]*vector(GG[i])for i in range(3)),s.zeros(150,1)),*map(vector,HH))
        J=C.row_join(s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P))for i,j in flag]))
        comp=[i for i in range(29)if i!=23];B=J[:,comp];b=J[:,23]
        B0=B.subs(zero);B1=homogeneous(B,1);B2=homogeneous(B,2);b0=b.subs(zero);b1=homogeneous(b,1);b2=homogeneous(b,2)
        inverse=B0[Y,:].inv();c0=-inverse*b0[Y,:];c1=-inverse*(b1[Y,:]+B1[Y,:]*c0);c2=-inverse*(b2[Y,:]+B1[Y,:]*c1+B2[Y,:]*c0)
        S0=(b0+B0*c0).applyfunc(s.expand);S1=(b1+B0*c1+B1*c0).applyfunc(s.expand);S2=(b2+B0*c2+B1*c1+B2*c0).applyfunc(s.expand)
        assert S0==S1==s.zeros(150,1) and S2[Y,:]==s.zeros(28,1) and S2[61]==0
        tangentR=vector(P).jacobian(coordinates).subs(zero)[Y+[61],:]
        assert tangentR.rank()==3 and A*tangentR==s.zeros(150,3)
        residual=(S2-linear*S2[rr,:]).applyfunc(s.expand)
        monomials=[W*U,W*V,U*V,U*U,V*V,W*W]
        quad=s.Matrix([[s.Poly(f,*coordinates).coeff_monomial(m)for m in monomials]for f in residual])
        assert quad.rank()==2 and quad[:,2:]==s.zeros(150,4)
        basisrows=list(quad.T.rref()[1]);assert basisrows==[46,95]
        minor=s.factor(quad[basisrows,:2].det());assert minor!=0
        lift=s.factor((nu_linear*S2[rr,:])[0]);expected=-U*V/75 if name=='CA'else-s.Rational(16,9)*U*V
        assert lift==expected
        records[name]=dict(full_source_P=mat(P),variable_fiber150_by20=mat(C),full_physical_Jacobian150_by29=mat(J),critical_source_free_indices=free,critical_linear_source_coefficients=list(map(str,coeff)),
            source_complement_columns=comp,selected_actual28_Y=Y,constant_inverse=mat(inverse),vertical_coefficients012=[mat(c)for c in[c0,c1,c2]],
            all150_vertical_derivatives2=list(map(str,S2)),branch_row61_all012_zero=True,base_critical_tangent29_by3=mat(tangentR),
            actual_seven_odd_basis_rows=rr,selected_actual_base_columns=cc,all_actual_linear_combinations=mat(linear),
            all_zero_R_firstjet_odd_quadratics=list(map(str,residual)),quadratic_coefficient_order=list(map(str,monomials)),quadratic_coefficient_matrix=mat(quad),
            independent_quadratic_rows=basisrows,independent_WU_WV_minor=str(minor),chosen_nu_lift_combination=mat(nu_linear),chosen_nu_lift_quadratic=str(lift),
            exact_zero_firstjet_quadratic_span=['W*U','W*V'],nu_not_in_actual_I_by_missing_UV=True)
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(records,separators=(',',':'))+'\n').encode(),mtime=0))
    paths=[Path(__file__),ROOT/'research/notes/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63.md',folder/'report.json',folder/'certificate.json.gz',ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json',ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json',ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.json',ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py']
    for path in paths[4:7]:assert json.loads(path.read_text())['all_passed'] is True
    result=dict(id='AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63',all_passed=True,
        full_physical_CA_and_CB_sources_and_entire_fiber_flag_Jacobians=True,
        CB_all_ratio_times_f_mixed_and_quadratic_Jacobian_terms_retained=True,
        source_vertical_neumann_inverse_through2_all150_checked=True,
        actual_branch_row61_derivative_zero_through_degree2_both=True,
        actual_critical_tangent_maps_to_kernel_of_all7_odd_linears=True,
        actual_zero_R_firstjet_odd_quadratic_ranks={'CA':2,'CB':2},
        exact_quadratic_span_both=['W*U','W*V'],missing_UV_class=True,
        chosen_nu_lift_quadratic={'CA':'-U*V/75','CB':'-16*U*V/9'},
        exact_discriminant_not_in_actual_odd_ideal_CA_and_CB=True,
        tracefree_t_cubed_not_in_actual_CA_or_CB_image=True,
        full_odd_ideal_eq_J_Nakayama_strategy_excluded_for_CA_CB=True,
        full_odd_ideal_cubic_class_generator_number_depth_or_actual_D_not_claimed=True,
        R19_not_used_as_individual_CA_CB_conductor_branch=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('ACTUAL CA CB ODD2 DEFECT: spanWU/WV, nu NOTinI PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
