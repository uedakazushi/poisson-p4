#!/usr/bin/env python3
"""Independent physical source, inverse3 and actual cubic odd-class audit."""
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
STEM='AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63'
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
    subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'code/ambient_B18_CA_CB_odd_cubic_class_20261103_v63.py'),
                    '--output',str(replay)],check=True,timeout=40)
    assert (folder/'certificate.json.gz').read_bytes()==(replay/'certificate.json.gz').read_bytes()
    cert=read(folder/'certificate.json.gz')
    quadpath=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/certificate.json.gz';quad=read(quadpath)
    priorpath=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz';prior=read(priorpath);Y=prior['common_Y28']
    h,x,y,z,w=s.symbols('h x y z w');X=s.Matrix([h,x,y,z,w]);U=s.Symbol('U')
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
    flag=[(0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(1,4),(3,2),(4,2)];summaries={};criticalpaths=[]
    for name in ['CA','CB']:
        criticalpath=ROOT/f'reports/AMBIENT_B18_{name}_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';critical=read(criticalpath);criticalpaths.append(criticalpath)
        tt=list(map(s.Symbol,critical['slice_variables']))
        equations=[*map(s.sympify,critical['seven_monic_eliminations']),s.sympify(critical['u'])-U,s.sympify(critical['v']),s.sympify(critical['w'])]
        solutions=s.solve(equations,tt,dict=True);assert len(solutions)==1
        section=s.Matrix([solutions[0][t] for t in tt]);same(section,s.Matrix(list(map(s.sympify,cert[name]['exact_monic_source_section_s0_through_s9']))))
        for equation in equations:assert s.expand(equation.subs(solutions[0],simultaneous=True))==0
        polynomials=critical.get('all_nonzero_critical_polynomials',critical.get('all_nonzero_critical_generators'))
        restriction=[s.expand(s.sympify(q).subs(solutions[0],simultaneous=True)) for q in polynomials]
        assert s.Poly(s.gcd_list(restriction),U).monic().as_expr()==U**3
        center=[s.S.Zero]*20;center[0]=-1;center[7]=-1;center[11]=1;center[16]=s.Rational(1,3)
        values=s.Matrix(center)
        for i,c in zip(quad[name]['critical_source_free_indices'],section):values[i]+=c
        assert values[3]==0
        Pgeneric=change(canonical(raw(name,f,alpha[3] if name=='CA' else -alpha[3]*f,H)))
        Cgeneric=vec(Pgeneric).jacobian(alpha)
        P=Pgeneric.subs(dict(zip(alpha,values)),simultaneous=True).applyfunc(s.expand)
        C=Cgeneric.subs(dict(zip(alpha,values)),simultaneous=True).applyfunc(s.expand)
        same(P.subs(U,0),P0);same(P,mat(cert[name]['fullphysical_source_P']));same(C,mat(cert[name]['wholevariable_fiber_jacobian']))
        J=C.row_join(s.Matrix.hstack(*[vec(lie(s.eye(5)[:,i]*s.eye(5)[j,:],P)) for i,j in flag]))
        same(J,mat(cert[name]['whole150_by29_sourceJacobian']))
        comp=[i for i in range(29) if i!=23];B=J[:,comp];b=J[:,23]
        Bs=[B.applyfunc(lambda q:s.Poly(q,U).nth(d)*U**d) for d in range(4)]
        bs=[b.applyfunc(lambda q:s.Poly(q,U).nth(d)*U**d) for d in range(4)]
        same(B,sum(Bs,s.zeros(150,28)));same(b,sum(bs,s.zeros(150,1)))
        T=[M[Y,:] for M in Bs];R0=T[0].inv()
        R1=-R0*T[1]*R0
        R2=R0*T[1]*R0*T[1]*R0-R0*T[2]*R0
        R3=-R0*T[3]*R0+R0*T[1]*R0*T[2]*R0+R0*T[2]*R0*T[1]*R0-R0*T[1]*R0*T[1]*R0*T[1]*R0
        Rs=[R0,R1,R2,R3]
        for d in range(4):same(sum((Rs[i]*T[d-i] for i in range(d+1)),s.zeros(28)),s.eye(28) if d==0 else s.zeros(28))
        cs=[-sum((Rs[i]*bs[d-i][Y,:] for i in range(d+1)),s.zeros(28,1)) for d in range(4)]
        for c,expected in zip(cs,cert[name]['vertical_coefficients0123']):same(c,mat(expected))
        S=[(bs[d]+sum((Bs[i]*cs[d-i] for i in range(d+1)),s.zeros(150,1))).applyfunc(s.expand) for d in range(4)]
        assert all(M==s.zeros(150,1) for M in S[:3]);same(S[3][Y,:],s.zeros(28,1));assert S[3][61]==0
        same(S[3],s.Matrix(list(map(s.sympify,cert[name]['all150_output_vertical_derivative_order3']))))
        cubic_coeff=S[3].applyfunc(lambda q:s.Poly(q,U).nth(3));same(cubic_coeff,s.Matrix(list(map(s.sympify,cert[name]['actual_all150_odd_cubic_coefficients']))))
        linear=mat(quad[name]['all_actual_linear_combinations']);rows=quad[name]['actual_seven_odd_basis_rows']
        first=cubic_coeff-linear*cubic_coeff[rows,:]
        Q=mat(quad[name]['quadratic_coefficient_matrix'])[:,:2];qr=quad[name]['independent_quadratic_rows'];QC=Q*Q[qr,:].inv()
        same(QC,mat(cert[name]['all_two_quad_class_subtractions']));same(QC*Q[qr,:],Q)
        last=(first-QC*first[qr,:]).applyfunc(s.factor)
        same(first,s.Matrix(list(map(s.sympify,cert[name]['actual_zero_firstjet_cubic_classes']))))
        same(last,s.Matrix(list(map(s.sympify,cert[name]['actual_zero_linear_and_quad_cubic_classes']))))
        assert [i for i,value in enumerate(last) if value]==[63]
        expected=s.Rational(12,125) if name=='CA' else -s.Rational(4,243)
        assert last[63]==expected
        summaries[name]=dict(full_monic_section_and_all_source_Jacobian_terms_checked=True,
            independent_inverse3_and_all150_S0_S1_S2_zero=True,branch61_zero_through3=True,
            all_linear_and_quadratic_class_subtractions_checked=True,witness_row63_coefficient=str(expected))
    previous=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json';assert read(previous)['all_passed']
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}_REVIEW.md',folder/'report.json',folder/'certificate.json.gz',quadpath,priorpath,previous,
           ROOT/'research/notes'/f'{STEM}.md',*criticalpaths]
    result=dict(id=STEM+'_REVIEW',all_passed=True,author_fresh_certificate_byte_identical=True,
        independent_full_original_physical_CA_CB_sources_and_variable_fiber_derivatives=True,
        exact_all7_nonlinear_monic_eliminations_and_Uaxis_not_linear_tangent=True,
        independent_noncommuting_inverse3_and_all150_source_output_coefficients=summaries,
        tracefree_t_order4_and_actual_odd_cubic_interpretation_accepted=True,
        all_J_generators_order_atleast3_on_section_so_mJ_orderatleast4=True,
        new_cubic_class_independent_of_previous7linear_and2quadratic=True,
        image_I_in_J_mod_mJ_dimension_exactly10=True,
        actual_I_minimal_generators_lower_bound10=True,
        no_I_intersect_mJ_equals_mI_assumption=True,
        no_mu_upper_full_I_depth_analytic_Uaxis_type_or_actual_D_claim=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','image_I_in_J_mod_mJ_dimension_exactly10','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True);args=parser.parse_args();main(args.output.resolve(),args.replay.resolve())
