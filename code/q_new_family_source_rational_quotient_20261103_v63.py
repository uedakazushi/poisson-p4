#!/usr/bin/env python3
"""Exact five first integrals, rational fiber parameter, and symbolic five-source circuit."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

from q_new_family_six_parameter_covariants_20261103_v63 import DAG

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63'


def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def rows(M):return [[str(x) for x in row] for row in M.tolist()]


def cancel(expr):return s.cancel(expr).doit()


def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    kernelpath=ROOT/'reports/Q_NEW_FAMILY_SOURCE_KERNEL_20261103_v63/kernel-certificate.json.gz'
    sixpath=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json'
    circuitpath=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/source-six-circuit.json.gz'
    kernel,six,old=map(read,[kernelpath,sixpath,circuitpath])
    u,v,c,X,W,B=s.symbols('u v c X W B');A=s.Integer(256)
    oldvars=[u,v,c,X,W,B]
    source=kernel['source_six_direction']
    velocity=[s.sympify(source['u']),s.sympify(source['v']),c*s.sympify(source['c_d_dc']),
              s.sympify(source['X']),s.sympify(source['W']),s.sympify(source['B'])]
    xi=384*X;zeta=3072*W-1
    rr=v*(xi+zeta)+A;ss=v*(xi-zeta)+A
    def old_derivative(expr):return s.expand(v/s.Integer(3072)*sum(w*s.diff(expr,z) for z,w in zip(oldvars,velocity)))
    r,z=s.symbols('r s')
    newvars=[u,v,r,z,c,B]
    new_velocity=[0,v*(r*r+z*z)+16*u*(r-z)**2,r*(r*r-A*A),z*(z*z-A*A),
                  c*(-r*r-z*z+2*A*A-32*u*(r-z)**2/v),
                  B*(-r*r-z*z+2*A*A-48*u*(r-z)**2/v)]
    images=[u,v,rr,ss,c,B]
    assert all(cancel(old_derivative(image)-w.subs({r:rr,z:ss},simultaneous=True))==0
               for image,w in zip(images,map(s.sympify,new_velocity)))
    inverse_X=(r+z-2*A)/(768*v)
    inverse_W=(1+(r-z)/(2*v))/3072
    assert cancel(rr.subs({X:inverse_X,W:inverse_W},simultaneous=True)-r)==0
    assert cancel(ss.subs({X:inverse_X,W:inverse_W},simultaneous=True)-z)==0
    zr,zs=1-A*A/r**2,1-A*A/z**2
    Vtilde=v-8*u*(r-z)**2/(r*z)
    invariants=[u,zr/zs,Vtilde/(r*z*zr),c*v*v/(r*z*zr*zr),B*v**3/((r*z)**2*zr**3)]
    def derivative(expr):return cancel(sum(w*s.diff(expr,z0) for z0,w in zip(newvars,new_velocity)))
    assert all(derivative(I)==0 for I in invariants)
    fibre=(r/z-1)/(r/A-1)
    fibre_derivative=s.factor(derivative(fibre))
    expected=A*r*(r-z)*(z-A)/(z*(r-A))
    assert cancel(fibre_derivative-expected)==0 and fibre_derivative!=0
    U,K,Iv,Ic,Ib,t=s.symbols('U K Iv Ic Ib t')
    den=1-K+K*t*t
    r_param=A*(1-2*(1-K+K*t)/den)
    ratio=1-2*t*(1-K+K*t)/den
    s_param=r_param/ratio
    assert cancel((1-K)*(r_param/A)**2+K*(r_param/s_param)**2-1)==0
    assert cancel((1-A*A/r_param**2)/(1-A*A/s_param**2)-K)==0
    assert cancel((r_param/s_param-1)/(r_param/A-1)-t)==0
    assert cancel(r_param.subs({K:invariants[1],t:fibre},simultaneous=True)-r)==0
    assert cancel(s_param.subs({K:invariants[1],t:fibre},simultaneous=True)-z)==0
    Zr=1-A*A/r_param**2
    v_param=Iv*r_param*s_param*Zr+8*U*(r_param-s_param)**2/(r_param*s_param)
    c_param=Ic*r_param*s_param*Zr**2/v_param**2
    B_param=Ib*(r_param*s_param)**2*Zr**3/v_param**3
    # Remaining inverse identities are checked before expanding the conic formulas.
    vp=Iv*r*z*zr+8*U*(r-z)**2/(r*z)
    cp=Ic*r*z*zr**2/vp**2;Bp=Ib*(r*z)**2*zr**3/vp**3
    targets=[U,K,Iv,Ic,Ib]
    assert cancel((vp-8*U*(r-z)**2/(r*z))/(r*z*zr)-Iv)==0
    assert cancel(cp*vp**2/(r*z*zr**2)-Ic)==0
    assert cancel(Bp*vp**3/((r*z)**2*zr**3)-Ib)==0
    back={U:u,Iv:invariants[2],Ic:invariants[3],Ib:invariants[4]}
    assert cancel(vp.subs(back,simultaneous=True)-v)==0
    assert cancel(cp.subs(back,simultaneous=True)-c)==0
    assert cancel(Bp.subs(back,simultaneous=True)-B)==0
    # t=2 is a symbolic rational section, with five freely varying parameters.
    section_condition=v*(xi*xi-zeta*zeta)+A*xi
    assert cancel((rr*ss-A*(rr+ss)/2)/v-section_condition)==0
    v_section=-A*xi/(xi*xi-zeta*zeta)
    assert cancel(section_condition.subs(v,v_section))==0
    r_section=rr.subs(v,v_section);s_section=ss.subs(v,v_section)
    assert cancel((r_section/s_section-1)/(r_section/A-1)-2)==0
    # Compose all exact old source outputs without evaluating at any numerical point.
    dag=DAG(['u','c','X','W','B'])
    replacement=[dag.expr(z0) for z0 in [u,v_section,c,X,W,B]]
    mapping=[]
    for op,arg in old['nodes']:
        if op=='q':node=dag.q(s.Rational(arg))
        elif op=='v':node=replacement[arg]
        elif op=='+':node=dag.add(*(mapping[j] for j in arg))
        elif op=='*':node=dag.mul(*(mapping[j] for j in arg))
        elif op=='inv':node=dag.inv(mapping[arg])
        else:raise ValueError(op)
        mapping.append(node)
    def remap(obj):
        if isinstance(obj,int):return mapping[obj]
        if isinstance(obj,list):return [remap(z0) for z0 in obj]
        return {key:remap(z0) for key,z0 in obj.items()}
    composed=dict(input_names=dag.names,nodes=dag.nodes,outputs=remap(old['outputs']),
                  source_v=str(v_section),generic_symbolic_section='t=2',
                  coefficient_field='QQ',no_point_evaluation_performed=True,
                  domain='u*c*X*B*(xi^2-zeta^2)*(32u+v_section) !=0, plus old m-domain when m output is used')
    with gzip.open(out/'source-five-circuit.json.gz','wt') as stream:json.dump(composed,stream);stream.write('\n')
    # A known finite full-jet equivalence, not an assertion that all ambiguities are known.
    tauold={u:-u,v:-v,c:-c,X:-X,W:W,B:B}
    rawC,rawL=s.Matrix(six['raw_value4']),s.Matrix(six['raw_linear4'])
    qvec=s.Matrix(s.symbols('q0:4'));G=s.diag(-1,1,-1,1)
    assert (G*(c*rawC)*G.T-(c*rawC).subs(tauold,simultaneous=True)).applyfunc(cancel)==s.zeros(4)
    transL=(c*rawL).subs(dict(zip(qvec,G*qvec)),simultaneous=True)
    assert (G*transL*G.T-(c*rawL).subs(tauold,simultaneous=True)).applyfunc(s.expand)==s.zeros(4)
    taunew={u:-u,v:-v,r:z,z:r,c:-c,B:B}
    tau_targets=[-invariants[0],1/invariants[1],-invariants[1]*invariants[2],
                 -invariants[1]**2*invariants[3],-invariants[1]**3*invariants[4]]
    assert all(cancel(I.subs(taunew,simultaneous=True)-T0)==0 for I,T0 in zip(invariants,tau_targets))
    assert cancel(fibre.subs(taunew,simultaneous=True)-fibre/(fibre-1))==0
    assert cancel(v_section.subs(X,-X)+v_section)==0
    cert=dict(original_source_variables=list(map(str,oldvars)),original_kernel_source_velocity=list(map(str,velocity)),
              derivation_rescaling='delta=(v/3072)*verified polynomial kernel source derivation',
              normalized_coordinates=dict(xi=str(xi),zeta=str(zeta),r=str(rr),s=str(ss)),
              inverse_X=str(inverse_X),inverse_W=str(inverse_W),
              normalized_variables=list(map(str,newvars)),normalized_derivation=list(map(str,new_velocity)),
              five_integral_names=['U','K','Iv','Ic','Ib'],five_first_integrals=list(map(str,invariants)),
              fibre_parameter_t=str(fibre),delta_t=str(fibre_derivative),
              inverse_denominator=str(den),inverse_r=str(r_param),inverse_s=str(s_param),
              inverse_v=str(v_param),inverse_c=str(c_param),inverse_B=str(B_param),
              exact_source_field='QQ(u,v,c,X,W,B)=QQ(U,K,Iv,Ic,Ib,t)',
              derivation_constant_field='QQ(U,K,Iv,Ic,Ib)',
              section_t2_condition=str(section_condition),section_v=str(v_section),
              known_finite_source_involution=['-u','-v','-c','-X','W','B'],
              known_finite_integral_involution=['-U','1/K','-K*Iv','-K^2*Ic','-K^3*Ib'],
              known_fixed_field_on_section=['X^2','u/X','c/X','W','B'],
              actual_NEW25_rational_orbit_invariant_field_finite_index_not_computed=True)
    (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    files=[Path(__file__),kernelpath,sixpath,circuitpath,ROOT/'code/q_new_family_six_parameter_covariants_20261103_v63.py']
    report=dict(id=STEM,all_passed=True,all_five_rational_first_integrals_verified=True,
                complete_rational_source_constant_field_proved=True,explicit_rational_inverse_with_t=True,
                symbolic_t2_section_gives_free_five_parameter_circuit=True,five_source_circuit_nodes=len(dag.nodes),
                known_order2_full_rawjet_equivalence_and_fixed_field_verified=True,
                actual_NEW25_invariant_field_equals_source_constant_field_claimed=False,
                actual_NEW25_field_finite_index_unknown=True,
                no_new_beta_point_no_high_degree_fit_no_GB=True,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in files},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','five_source_circuit_nodes','seconds']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output.resolve())
