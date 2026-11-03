#!/usr/bin/env python3
"""Independent physical N19/fullflag9 and actual odd-quadratic Neumann audit."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_N_ODD_QUADRATIC_JETS_20261103_v63'
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()
def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return json.load(f)
def mat(rows):return s.Matrix([[s.sympify(c) for c in row] for row in rows])
def same(A,B):assert A.shape==B.shape and (A-B).applyfunc(s.expand)==s.zeros(*A.shape)


def main(output,replay):
    started=time.monotonic();folder=ROOT/'reports'/STEM;report=read(folder/'report.json');assert report['all_passed']
    for name,digest in report['input_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in report['evidence_sha256'].items():assert sha(folder/name)==digest,name
    subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'code/ambient_B18_N_odd_quadratic_jets_20261103_v63.py'),
                    '--output',str(replay)],check=True,timeout=30)
    assert (folder/'certificate.json.gz').read_bytes()==(replay/'certificate.json.gz').read_bytes()
    cert=read(folder/'certificate.json.gz')
    priorpath=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz'
    prior=read(priorpath);A=mat(prior['N_full_firstjet_data']['actual']);dnu=mat(prior['N_full_firstjet_data']['discriminant'])
    same(A,mat(cert['actual_odd150_by28_linear']))
    basisrows=list(A.T.rref()[1]);assert basisrows==cert['seven_independent_actual_odd_rows']==[18,19,28,29,60,61,67]
    columns=list(A[basisrows,:].rref()[1]);inverse7=A[basisrows,columns].inv()
    linear=A[:,columns]*inverse7;same(linear,mat(cert['all_linear_row_combinations']));same(linear*A[basisrows,:],A)
    nu_coeff=dnu[:,columns]*inverse7;same(nu_coeff,mat(cert['nu_linear_combination']));same(nu_coeff*A[basisrows,:],dnu)
    assert list(nu_coeff)==[0,0,0,0,0,0,-s.Rational(1,8)]
    h,x,y,z,w=s.symbols('h x y z w');X=s.Matrix([h,x,y,z,w]);W,V=s.symbols('W V')
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
    def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[a]*X[b]) for i,j in pairs for a,b in mons])
    def canonical(P):
        d=s.Matrix([sum(s.diff(P[i,j],X[i]) for i in range(5)) for j in range(5)])
        return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
    def raw(f,g,H):
        P=s.zeros(5);P[1,2]=f;P[3,4]=g*h*h;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z)
        return P-P.T
    sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
    def change(P):return (sw*P.subs({h:x,x:h},simultaneous=True)*sw.T).applyfunc(s.expand)
    def lie(M,P):return (M*P+P*M.T-sum(((M*X)[i]*P.diff(X[i]) for i in range(5)),s.zeros(5))).applyfunc(s.expand)
    cubic=[s.prod((h,x,z,w)[i] for i in mon) for mon in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in mon)]
    basis=[change(canonical(raw(f,0,0))) for f in [h*h,h*x]]+[change(canonical(raw(0,1,0)))]+[change(canonical(raw(0,0,H))) for H in cubic]
    C=s.Matrix.hstack(*map(vec,basis));cr=list(C.T.rref()[1]);left=C[cr,:].inv()
    ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T
    Pw=basis[2];Pv=s.Rational(5,2)*basis[1]+s.Rational(11,2)*basis[8]+basis[16]
    same(Pw,mat(cert['normalizer_critical_tangent_Pw']));same(Pv,mat(cert['normalizer_critical_tangent_Pv']))
    P=P0+W*Pw+V*Pv
    flag=[(0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(1,4),(3,2),(4,2)]
    generators=[s.eye(5)[:,i]*s.eye(5)[j,:] for i,j in flag]
    J=C.row_join(s.Matrix.hstack(*[vec(lie(M,P)) for M in generators]))
    same(J,mat(cert['full_physical_source_Jacobian_on_plane150_by28']))
    criticalpath=ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';critical=read(criticalpath)
    same(C,mat(critical['fiber19_basis']));assert flag==list(map(tuple,critical['flag_positions']))
    parabolic=[(i,j) for i,j in itertools.product(range(5),repeat=2) if (i,j) not in flag]
    orbit=[]
    for i,j in parabolic:
        M=s.eye(5)[:,i]*s.eye(5)[j,:]
        actions=s.Matrix.hstack(*[vec(lie(M,Q)) for Q in basis])
        same(C*left*actions[cr,:],actions)
        orbit.append(left*vec(lie(M,P0))[cr,:])
    orbit=s.Matrix.hstack(*orbit);same(orbit,mat(critical['parabolic_orbit']))
    det=orbit[critical['fixed_rows'],critical['orbit10_columns']].det()
    assert det==s.sympify(critical['orbit10_minor'])!=0
    assert critical['free_indices']==[1,2,8,9,13,14,16,17,18]
    def flags(M):
        L=M[:,2];L[2]=0;ell=-M[1,:].T;ell[1]=0
        wedge=[0 if (i,j)==(0,1) else -(M[0,i]*int(j==1)-M[0,j]*int(i==1))-(int(i==0)*M[1,j]-int(j==0)*M[1,i]) for i,j in pairs]
        return s.Matrix(list(L)+list(ell)+wedge)
    H=s.Matrix.hstack(*map(flags,generators));assert H.rank()==9
    augmented=J.subs({W:0,V:0}).col_join(s.zeros(20,19).row_join(H));assert augmented.rank()==28
    rr=cert['actual27_target_rows'];complement=[i for i in range(28) if i!=22];assert complement==cert['source_complement_columns']
    B=J[:,complement];b=J[:,22];B0=B.subs({W:0,V:0});B1=B-B0;b0=b.subs({W:0,V:0});b1=b-b0
    T0=B0[rr,:];T1=B1[rr,:];R0=T0.inv();R1=-R0*T1*R0;R2=R0*T1*R0*T1*R0
    same(R0,mat(cert['complement_inverse_constant']));same(R0*T0,s.eye(27))
    same(R0*T1+R1*T0,s.zeros(27));same(R1*T1+R2*T0,s.zeros(27))
    c0=-R0*b0[rr,:];c1=-(R0*b1[rr,:]+R1*b0[rr,:]);c2=-(R1*b1[rr,:]+R2*b0[rr,:])
    for c,expected in zip([c0,c1,c2],cert['vertical_vector_coefficients_degrees0_1_2']):same(c,mat(expected))
    S0=(b0+B0*c0).applyfunc(s.expand);S1=(b1+B0*c1+B1*c0).applyfunc(s.expand);S2=(B0*c2+B1*c1).applyfunc(s.expand)
    same(S0,s.zeros(150,1));same(S1,s.zeros(150,1));same(S2[rr,:],s.zeros(27,1));assert S2[27]==0
    same(S2,s.Matrix(list(map(s.sympify,cert['all150_vertical_output_derivative_degree2']))))
    tangent=s.Matrix.hstack(vec(Pw),vec(Pv))[rr+[27],:];assert tangent.rank()==2;same(A*tangent,s.zeros(150,2))
    residual=(S2-linear*S2[basisrows,:]).applyfunc(s.expand)
    same(residual,s.Matrix(list(map(s.sympify,cert['all_zero_firstjet_generator_quadratics']))))
    assert residual[46]==-4*W*V and residual[72]==-6*V**2
    coefficient=s.Matrix([[s.Poly(v,W,V).coeff_monomial(mon) for mon in [W**2,W*V,V**2]] for v in residual])
    assert coefficient.rank()==2 and coefficient[[46,72],[1,2]].det()==24
    assert s.expand((nu_coeff*S2[basisrows,:])[0])==-s.Rational(9,4)*V**2
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}_REVIEW.md',folder/'report.json',folder/'certificate.json.gz',priorpath,criticalpath,
           ROOT/'research/notes'/f'{STEM}.md',ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json']
    result=dict(id=STEM+'_REVIEW',all_passed=True,
        author_fresh_certificate_byte_identical=True,
        independent_N19_physical_basis_and_full150_by28_source_Jacobian=True,
        all16_by19_parabolic_columns_stable_and_10plus9_etale_gauge_checked=True,
        full9_flag_chart_and_augmented_source28_checked=True,
        independent_truncated_inverse2_all150_Schur_coefficients=True,
        whole_branch_derivative_zero_through2_and_trace_correction_proof_accepted=True,
        seven_actual_linear_odd_generators_and_critical_tangent_rank2_checked=True,
        zero_full_R_firstjet_two_odd_classes_WV_V2_minor=24,
        actual_odd_minimal_generator_lower_bound9=True,
        core_flag_tracefree_coordinates_differ_by_unit_center1_zero_firstjet2classes_invariant=True,
        firstjet_neither_full_ideal_nor_pair_arc_surrogate=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','actual_odd_minimal_generator_lower_bound9','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True);args=parser.parse_args();main(args.output.resolve(),args.replay.resolve())
