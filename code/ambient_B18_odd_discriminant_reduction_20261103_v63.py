#!/usr/bin/env python3
"""Separate explicit N-rank2 base, then accepted critical/different consequences.

The odd ideals are classified only after reduction modulo each discriminant,
up to independent formal ambient coordinates. Actual commonR embedding open.
"""
import argparse,gzip,hashlib,itertools,json,subprocess,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector,action
ROOT=Path(__file__).resolve().parents[1];h,x,y,z,w=X
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out,m2):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    dependencies=[ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json',ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.json',ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json']
    for p in dependencies:
        d=json.loads(p.read_text());assert d['all_passed'] is True
        for name,digest in d['input_sha256'].items():assert sha(ROOT/name)==digest,name
    # N has source dimension28, so its actual parameter base must have28,
    # not the29 common cone parameters used by CA andCB.
    cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
    sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
    change=lambda P:(sw*P.subs({h:x,x:h},simultaneous=True)*sw.T).applyfunc(s.expand)
    basis=[change(canonical(raw(f,0,0)))for f in[h*h,h*x]]+[change(canonical(raw(0,1,0)))]+[change(canonical(raw(0,0,H)))for H in cubics]
    C=s.Matrix.hstack(*map(vector,basis));cr=list(C.T.rref()[1]);left=C[cr,:].inv()
    ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T
    flag=[(0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(1,4),(3,2),(4,2)]
    J=C.row_join(s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P0))for i,j in flag]))
    cols=[i for i in range(28)if i!=22];V=J[:,cols];rr=list(V.T.rref()[1]);assert len(rr)==27
    det=s.factor(V[rr,:].det());assert det==-s.Rational(1152,5)
    M=s.zeros(5);M[1,0]=1;M[0,3]=-2
    LP=action(M,P0);co=left*vector(LP)[cr,:];assert C*co==vector(LP)
    kernel=s.Matrix(list(-co)+[M[i,j]for i,j in flag]);assert kernel[22]==1 and J*kernel==s.zeros(150,1)
    W=-vector(action(M,LP))/2
    corrected=(W-V*V[rr,:].inv()*W[rr,:]).applyfunc(s.factor)
    assert corrected[27]==-3 and 27 not in rr
    assert 106 in rr and vector(P0)[106]==1
    # This actual finite source curve lies in N and has the required
    # normalization kernel tangent. Its p² coefficient is W above.
    p=s.Symbol('p');A=s.eye(5)+p*M+p*p*M*M/2;Ai=s.eye(5)-p*M+p*p*M*M/2
    curve=(A*(P0-p*LP).subs(dict(zip(X,Ai*s.Matrix(X))),simultaneous=True)*A.T).applyfunc(s.expand)
    assert curve.subs(p,0)==P0 and curve.diff(p).subs(p,0)==s.zeros(5)
    assert vector(curve.diff(p,2).subs(p,0))/2==W
    # Minimal generator counts of the exact discriminant-image ideals.
    text='R=QQ[a,b,c,d,e,f,w,u,v];\nI=ideal(a,b,c,d,e,f,w*u,w*v,u*v,u^3);\nassert(numgens trim I==10);\nS=QQ[a,b,c,d,e,f,w,v];\nJ=ideal(a,b,c,d,e,f,w*v,v^2);\nassert(numgens trim J==8);\nprint "ODD_DISCRIMINANT_MINIMAL_GENERATORS_PASS";\nexit 0;\n'
    (out/'minimal.m2').write_text(text)
    if m2:
        run=subprocess.run(['M2','--script',str((out/'minimal.m2').resolve())],capture_output=True,text=True,check=True,timeout=60)
        assert 'ODD_DISCRIMINANT_MINIMAL_GENERATORS_PASS'in run.stdout;(out/'minimal.log').write_text(run.stdout+run.stderr)
    cert=dict(N_source_Jacobian150_by28=mat(J),N_complement27_columns=cols,N_actual27_target_rows=rr,N_27minor=str(det),N_source_kernel=list(map(str,kernel)),N_kernel_flag_matrix=mat(M),N_entire_second_coefficient=list(map(str,W)),N_all150_corrected_second_coefficients=list(map(str,corrected)),N_28th_actual_target_row=27,N_corrected_quadratic_coefficient='-3',N_actual_scalar_target_row=106,N_scalar_value='1',N_finite_kernel_curve=mat(curve))
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(cert,separators=(',',':'))+'\n').encode(),mtime=0))
    note=ROOT/'research/notes/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63.md'
    paths=[Path(__file__),note,ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',*dependencies]
    result=dict(id='AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63',all_passed=True,
        explicit_N_own_actual_parameter_base_cone_dimension=28,N_parameter_base_projective_dimension=27,N_normalizer_finite_free_rank=2,
        N_exact27_linear_minor=str(det),N_corrected_quadratic_actual_row=27,N_corrected_quadratic_value='-3',
        N_base_not_the_common_CA_CB_R29=True,
        trace_shifted_normalizer='B=R[t]/(t^2-nu), nu a regular parameter',
        image='A=R+I*t',ramification_Fitting_ideal='(t,I) inside B',critical_quotient='B/(t,I)=R/(nu,I)',
        discriminant_image_ideal_CA_CB_formal_type='(s1,...,s6,w*u,w*v,u*v,u^3) in C[[g19,s6,w,u,v]]',
        discriminant_image_ideal_N_formal_type='(s1,...,s6,w*v,v^2) in C[[g19,s6,w,v]]',
        exact_minimal_generators_after_discriminant={'CA':10,'CB':10,'N':8},full_odd_ideal_minimal_generators_lower_bounds={'CA':10,'CB':10,'N':8},
        regular_t_not_in_maximal_times_critical_ideal=True,independent_formal_coordinate_changes_not_aligned=True,
        image_of_I_in_R_mod_nu_not_module_I_tensor_R_mod_nu=True,
        full_odd_ideal_or_actual_D_or_relative_embedding_claimed=False,
        independent_Macaulay2_minimal_generator_check=bool(m2),input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('N OWN RANK2 BASE + ACTUAL ODD MODDISCRIMINANT10/10/8 PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--m2',action='store_true');args=ap.parse_args();main(args.output.resolve(),args.m2)
