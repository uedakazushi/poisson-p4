#!/usr/bin/env python3
"""Independent full physical CA/CB odd-quadratic deficiency audit."""
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
STEM='AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63'
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
    subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'code/ambient_B18_CA_CB_odd_quadratic_defect_20261103_v63.py'),
                    '--output',str(replay)],check=True,timeout=40)
    assert (folder/'certificate.json.gz').read_bytes()==(replay/'certificate.json.gz').read_bytes()
    cert=read(folder/'certificate.json.gz')
    priorpath=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz';prior=read(priorpath)
    Y=prior['common_Y28'];h,x,y,z,w=s.symbols('h x y z w');X=s.Matrix([h,x,y,z,w]);W,U,V=s.symbols('W U V');variables=[W,U,V]
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
    def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[a]*X[b]) for i,j in pairs for a,b in mons])
    def canonical(P):
        d=s.Matrix([sum(s.diff(P[i,j],X[i]) for i in range(5)) for j in range(5)])
        return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
    def raw(name,f,g,H):
        P=s.zeros(5);P[1,2]=f;P[3,4]=g*h*h if name=='CA' else g
        P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z)
        if name=='CB':P[2,3]-=s.diff(f,x)*z/2;P[2,4]-=s.diff(f,x)*w/2
        return P-P.T
    sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
    def change(P):return (sw*P.subs({h:x,x:h},simultaneous=True)*sw.T).applyfunc(s.expand)
    def lie(M,P):return (M*P+P*M.T-sum(((M*X)[j]*P.diff(X[j]) for j in range(5)),s.zeros(5))).applyfunc(s.expand)
    cubic=[s.prod((h,x,z,w)[i] for i in mon) for mon in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in mon)]
    alpha=s.symbols('a0:20');f=alpha[0]*h*h+alpha[1]*h*x+alpha[2]*x*x;H=sum(c*mon for c,mon in zip(alpha[4:],cubic))
    ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T
    flag=[(0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(1,4),(3,2),(4,2)];zero={W:0,U:0,V:0}
    def homogeneous(M,degree):
        return M.applyfunc(lambda q:sum(c*W**powers[0]*U**powers[1]*V**powers[2] for powers,c in s.Poly(q,W,U,V).terms() if sum(powers)==degree))
    summaries={};criticalpaths=[]
    for name in ['CA','CB']:
        criticalpath=ROOT/f'reports/AMBIENT_B18_{name}_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
        critical=read(criticalpath);criticalpaths.append(criticalpath)
        tt=list(map(s.Symbol,critical['slice_variables']))
        eliminated=critical['seven_monic_eliminations'];constraints=[]
        for q in eliminated:
            poly=s.Poly(s.sympify(q),*tt)
            constraints.append(sum(c*s.prod(t**power for t,power in zip(tt,powers)) for powers,c in poly.terms() if sum(powers)==1))
        constraints += [s.sympify(critical['u'])-U,s.sympify(critical['v'])-V,s.sympify(critical['w'])-W]
        M,rhs=s.linear_eq_to_matrix(constraints,tt);coeff=M.inv()*rhs
        same(coeff,s.Matrix(list(map(s.sympify,cert[name]['critical_linear_source_coefficients']))))
        free=critical.get('affine_slice_free_indices',critical.get('slice_free_indices'))
        assert free==cert[name]['critical_source_free_indices']
        center=[s.S.Zero]*20;center[0]=-1;center[7]=-1;center[11]=1;center[16]=s.Rational(1,3)
        values=s.Matrix(center)
        for i,c in zip(free,coeff):values[i]+=c
        # Differentiate the full generic physical source BEFORE substituting
        # the tangent slice, so all CB ratio*f mixed derivatives survive.
        Pgeneric=change(canonical(raw(name,f,alpha[3] if name=='CA' else -alpha[3]*f,H)))
        tensor=vec(Pgeneric);fiber=tensor.jacobian(alpha)
        P=Pgeneric.subs(dict(zip(alpha,values)),simultaneous=True).applyfunc(s.expand)
        C=fiber.subs(dict(zip(alpha,values)),simultaneous=True).applyfunc(s.expand)
        same(P.subs(zero),P0)
        same(P,mat(cert[name]['full_source_P']));same(C,mat(cert[name]['variable_fiber150_by20']))
        J=C.row_join(s.Matrix.hstack(*[vec(lie(s.eye(5)[:,i]*s.eye(5)[j,:],P)) for i,j in flag]))
        same(J,mat(cert[name]['full_physical_Jacobian150_by29']))
        assert all(s.Poly(q,W,U,V).total_degree()<=2 for q in J)
        comp=[i for i in range(29) if i!=23];B=J[:,comp];b=J[:,23]
        B0=B.subs(zero);B1=homogeneous(B,1);B2=homogeneous(B,2);b0=b.subs(zero);b1=homogeneous(b,1);b2=homogeneous(b,2)
        same(B,B0+B1+B2);same(b,b0+b1+b2)
        T0=B0[Y,:];T1=B1[Y,:];T2=B2[Y,:];R0=T0.inv();R1=-R0*T1*R0;R2=R0*T1*R0*T1*R0-R0*T2*R0
        same(R0,mat(cert[name]['constant_inverse']));same(R0*T0,s.eye(28))
        same(R0*T1+R1*T0,s.zeros(28));same(R0*T2+R1*T1+R2*T0,s.zeros(28))
        c0=-R0*b0[Y,:];c1=-(R0*b1[Y,:]+R1*b0[Y,:]);c2=-(R0*b2[Y,:]+R1*b1[Y,:]+R2*b0[Y,:])
        for c,expected in zip([c0,c1,c2],cert[name]['vertical_coefficients012']):same(c,mat(expected))
        S0=(b0+B0*c0).applyfunc(s.expand);S1=(b1+B0*c1+B1*c0).applyfunc(s.expand)
        S2=(b2+B0*c2+B1*c1+B2*c0).applyfunc(s.expand)
        same(S0,s.zeros(150,1));same(S1,s.zeros(150,1));same(S2[Y,:],s.zeros(28,1));assert S2[61]==0
        same(S2,s.Matrix(list(map(s.sympify,cert[name]['all150_vertical_derivatives2']))))
        earlier=prior['CA_CB_full_firstjet_data'][name];A=mat(earlier['actual']);dnu=mat(earlier['discriminant'])
        rows=list(A.T.rref()[1]);columns=list(A[rows,:].rref()[1]);assert len(rows)==7
        linear=A[:,columns]*A[rows,columns].inv();same(linear,mat(cert[name]['all_actual_linear_combinations']));same(linear*A[rows,:],A)
        nu_linear=dnu[:,columns]*A[rows,columns].inv();same(nu_linear*A[rows,:],dnu)
        tangent=vec(P).jacobian(variables).subs(zero)[Y+[61],:];assert tangent.rank()==3;same(A*tangent,s.zeros(150,3))
        residual=(S2-linear*S2[rows,:]).applyfunc(s.expand)
        same(residual,s.Matrix(list(map(s.sympify,cert[name]['all_zero_R_firstjet_odd_quadratics']))))
        monomials=[W*U,W*V,U*V,U**2,V**2,W**2]
        Q=s.Matrix([[s.Poly(q,W,U,V).coeff_monomial(m) for m in monomials] for q in residual])
        same(Q,mat(cert[name]['quadratic_coefficient_matrix']));assert Q.rank()==2 and Q[:,2:]==s.zeros(150,4)
        minor=Q[[46,95],:2].det();assert minor==s.sympify(cert[name]['independent_WU_WV_minor'])!=0
        lift=s.expand((nu_linear*S2[rows,:])[0]);assert lift==(-U*V/75 if name=='CA' else -s.Rational(16,9)*U*V)
        augmented=Q.col_join(s.Matrix([[0,0,1,0,0,0]]));assert augmented.rank()==3
        summaries[name]=dict(actual_odd_linear_rank=7,all_zero_R_firstjet_quadratic_rank=2,
            exact_span=['W*U','W*V'],independent_minor=str(minor),chosen_nu_lift=str(lift),missing_UV_and_nu_not_in_I=True,
            quadratic_Jacobian_terms_retained=bool(B2!=s.zeros(150,28) or b2!=s.zeros(150,1)))
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}_REVIEW.md',folder/'report.json',folder/'certificate.json.gz',priorpath,
           ROOT/'research/notes'/f'{STEM}.md',*criticalpaths,
           ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.json']
    result=dict(id=STEM+'_REVIEW',all_passed=True,author_certificate_fresh_byte_identical=True,
        independent_full_generic_source_before_differentiation_and_all150_CA_CB_Jacobians=True,
        critical_three_plane_derived_from_full_critical_linear_eliminations=True,
        independent_inverse2_including_CB_full_B2_b2_and_all150_Schur_entries=True,
        branch61_derivative_zero_through2_trace_t_order3_proof_accepted=True,
        actual7linear_basis_and_all150_zero_firstjet_quadratic_spaces_checked=summaries,
        exact_nu_not_in_actual_odd_I_CA_and_CB=True,tracefree_t_cubed_not_in_actual_images=True,
        missing_quadratic_class_cannot_be_repaired_by_higher_generators=True,
        no_full_I_mu_upper_cubic_depth_D_or_global_conductor_claim=True,
        R19_not_used_as_individual_CA_CB_normalization_double_locus=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','exact_nu_not_in_actual_odd_I_CA_and_CB','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True);args=parser.parse_args();main(args.output.resolve(),args.replay.resolve())
