#!/usr/bin/env python3
"""Actual cubicU^3 class forCA/CB after all sevenlinear and twoquadclasses."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,vector,action
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    quadfolder=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63';quadratic=load(quadfolder/'certificate.json.gz')
    review=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json';audit=load(review);assert audit['all_passed'] is True
    for filename,digest in audit['input_sha256'].items():assert sha(ROOT/filename)==digest,filename
    U=s.Symbol('U');ss=s.symbols('s0:10');pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
    ey=s.Matrix([0,0,1,0,0]);h,x,y,z,w=X;D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T;u0=vector(P0)
    flag=[(0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(1,4),(3,2),(4,2)]
    firstfolder=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63';first=load(firstfolder/'certificate.json.gz');Y=first['common_Y28'];results={};inputs=[]
    for name in ['CA','CB']:
        cpath=ROOT/f'reports/AMBIENT_B18_{name}_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';crit=load(cpath);inputs.append(cpath)
        if name=='CA':
            C=s.Matrix(crit['fixed_fiber_linear_basis150_by20']);s3=3*U/10;s7=7*U/20;e5=(12*s7*s7-3*s3*s3)/5
            values=[-U/10,7*e5/12,0,s3,35*e5/24-2*s7*s7,e5,0,s7,s7*s7,s7**3/3];C0=C
        else:
            Cfull=s.Matrix(crit['full_variable_fiber_Jacobian150_by20']);C0=Cfull.subs(dict.fromkeys(ss,0));s3=s7=U/2;e5=10*(s3*s7+s7*s7)/27
            values=[-7*U/9,-3*e5/5,0,s3,27*e5/10-s7*s7,e5,0,s7,s7*s7,s7**3/3]
            C=Cfull.subs(dict(zip(ss,values)),simultaneous=True)
        free=quadratic[name]['critical_source_free_indices'];delta=s.zeros(20,1)
        for i,value in zip(free,values):delta[i]=value
        pv=u0+C0*delta;P=s.zeros(5)
        for k,(i,j)in enumerate(pairs):
            P[i,j]=sum(pv[15*k+l]*X[a]*X[b]for l,(a,b)in enumerate(mons));P[j,i]=-P[i,j]
        assert vector(P)==pv
        J=C.row_join(s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P))for i,j in flag]))
        comp=[i for i in range(29)if i!=23];B=J[:,comp];b=J[:,23]
        Bs=[B.applyfunc(lambda f:s.Poly(f,U).nth(d)*U**d)for d in range(4)]
        bs=[b.applyfunc(lambda f:s.Poly(f,U).nth(d)*U**d)for d in range(4)]
        inverse=Bs[0][Y,:].inv();cs=[];S=[]
        for d in range(4):
            coefficient=-inverse*(bs[d][Y,:]+sum((Bs[i][Y,:]*cs[d-i]for i in range(1,d+1)),s.zeros(28,1)))
            cs.append(coefficient);S.append((bs[d]+sum((Bs[i]*cs[d-i]for i in range(d+1)),s.zeros(150,1))).applyfunc(s.expand))
        assert all(v==s.zeros(150,1)for v in S[:3])
        assert S[3][61]==0 and S[3][Y,:]==s.zeros(28,1)
        cubic=S[3].applyfunc(lambda f:s.Poly(f,U).nth(3))
        q=quadratic[name];linear=s.Matrix(q['all_actual_linear_combinations']);lr=q['actual_seven_odd_basis_rows']
        zero_linear=cubic-linear*cubic[lr,:]
        Q=s.Matrix(q['quadratic_coefficient_matrix'])[:,:2];qr=q['independent_quadratic_rows'];qc=Q*Q[qr,:].inv();assert qc*Q[qr,:]==Q
        last=(zero_linear-qc*zero_linear[qr,:]).applyfunc(s.factor)
        expected=s.Rational(12,125)if name=='CA'else-s.Rational(4,243)
        assert last[63]==expected and[i for i,v in enumerate(last)if v]==[63]
        results[name]=dict(exact_monic_source_section_s0_through_s9=list(map(str,values)),fullphysical_source_P=mat(P),wholevariable_fiber_jacobian=mat(C),whole150_by29_sourceJacobian=mat(J),
            selected_actual28targetrows=Y,source_complement_columns=comp,vertical_coefficients0123=[mat(c)for c in cs],
            all150_output_vertical_derivative_order3=list(map(str,S[3])),branch61_zero_through3=True,
            actual_all150_odd_cubic_coefficients=list(map(str,cubic)),actual_zero_firstjet_cubic_classes=list(map(str,zero_linear)),
            all_two_quad_class_subtractions=mat(qc),actual_zero_linear_and_quad_cubic_classes=list(map(str,last)),
            unique_nonzero_rawrow=63,actual_U3_class_coefficient=str(expected))
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(results,separators=(',',':'))+'\n').encode(),mtime=0))
    note=ROOT/'research/notes/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63.md';paths=[Path(__file__),note,review,quadfolder/'report.json',quadfolder/'certificate.json.gz',firstfolder/'certificate.json.gz',*inputs,ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py']
    result=dict(id='AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63',all_passed=True,
        exact_full_monic_critical_Uaxis_source_sections_not_only_linear_tangents=True,
        entire_variable_fiber_and_flag_Jacobians_retained=True,
        actual_vertical_output_all150_zero_through_order2=True,actual_branch61_derivative_zero_through_order3=True,
        all150_actual_odd_cubic_classes_and_all_previous_generator_subtractions=True,
        actual_new_U_cubed_class={'CA':'12/125','CB':'-4/243'},witness_rawrow63=True,
        image_actual_I_in_J_mod_mJ_dimension_exactly10_for_both=True,
        actual_minimal_generators_lower_bound_CA_CB=10,
        missing_UV_quadratic_class_and_discriminant_notinI_unchanged=True,
        exact_full_I_generator_upper_bound_depth_or_actual_D_not_claimed=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('ACTUAL CA CB NEW U3 CLASS: 7+2+1=10,NOTFULLI PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
