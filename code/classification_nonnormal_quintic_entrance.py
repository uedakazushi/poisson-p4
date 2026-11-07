#!/usr/bin/env python3
"""Finite characteristic-zero checks for the nonnormal irreducible quintic entrance.

Human geometric proofs accompany the output.  The optional Macaulay2 job is
over QQ(ell), not a parameter sample, and computes the T40 zero-surface ideal.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import subprocess
import time

from classification_quadratic_cofactor_review import Ring, PAIRS


def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c


def binomial_shift(a):
    v=[F(1)]
    for i in range(4):v=mul(v,[F(a-i),F(1)])
    return [c/24 for c in v]


def finite_checks():
    p=[a-b for a,b in zip(binomial_shift(4),binomial_shift(-1))]
    assert p==[0,F(25,6),0,F(5,6),0]
    numerical=[dict(sectional_genus=g,weighted_delta_surface_degree=6-g,
                    normalization_conductor_degree=12-2*g) for g in range(6)]
    R=Ring(0)
    def jacobi(entries):
        P=[[{} for j in range(4)] for i in range(4)]
        for (i,j),v in zip(PAIRS,entries):P[i][j],P[j][i]=v,R.scale(v,-1)
        return [R.add(*(R.mul(P[a][l],R.derivative(P[b][c],l)) for a,b,c in [(i,j,k),(j,k,i),(k,i,j)] for l in range(4))) for i,j,k in it.combinations(range(4),3)]
    def pf(p):return R.add(R.mul(p[0],p[5]),R.scale(R.mul(p[1],p[4]),-1),R.mul(p[2],p[3]))
    # Primary MP Lemma3.4.1: (x,y,z,w)=(h1,h2,e1,e2).
    linear=list(map(R.parse,['0','2*z','-w','-z','2*w','0']))
    assert all(not j for j in jacobi(linear))
    assert pf(linear)==R.parse('-3*z*w')
    # Local counterexample to automatic logarithmic residue on normalization.
    # This is a local Poisson germ, not asserted to be a global P4 example.
    cusp=list(map(R.parse,['1','0','0','0','0','w**2-z**3']))
    assert all(not j for j in jacobi(cusp))
    assert pf(cusp)==R.parse('w**2-z**3')
    return dict(quintic_hypersurface_Hilbert_polynomial_coefficients=list(map(str,p)),
                all_nonnormal_numeric_possibilities=numerical,
                T40_generic_conormal_linearization_Pfaffian='-3*e1*e2',
                T40_generic_transverse_delta=1,T40_zero_scroll_degree=5,
                T40_sectional_genus=1,T40_normalization_conductor_degree=10,
                local_cusp_Poisson_counterexample=True,
                cusp_normalization='u=a^2,v=a^3; residue du/(2v)=da/a^2=-d(1/a)',
                cusp_residue_pole_order=2,cusp_global_P4_existence_claimed=False)


def ordinary_double_checks():
    """Cleared log-symplectic inverse and the complete RR degree bounds."""
    R=Ring(6)
    lam,a,b,c,d,e=[R.variable(f'n{i}') for i in range(6)]
    z={}
    neg=lambda v:R.scale(v,-1)
    W=[[z,lam,a,b],[neg(lam),z,c,d],
       [neg(a),neg(c),z,e],[neg(b),neg(d),neg(e),z]]
    C=[[z,neg(e),d,neg(c)],[e,z,neg(b),a],
       [neg(d),b,z,neg(lam)],[c,neg(a),lam,z]]
    pf=R.add(R.mul(lam,e),neg(R.mul(a,d)),R.mul(b,c))
    for U,V in [(W,C),(C,W)]:
        for i,j in it.product(range(4),repeat=2):
            assert R.add(*(R.mul(U[i][k],V[k][j]) for k in range(4)))==(pf if i==j else {})
    # Actual bivector is diag(u,v,1,1) C diag(u,v,1,1)/pf.
    normal=[R.variable('x'),R.variable('y'),R.constant(1),R.constant(1)]
    P=[[R.mul(normal[i],C[i][j],normal[j]) for j in range(4)] for i in range(4)]
    on_double=lambda poly:{f:v for f,v in poly.items() if f[R.names.index('x')]==f[R.names.index('y')]==0}
    assert all(on_double(P[i][j])==(neg(lam) if (i,j)==(2,3) else lam if (i,j)==(3,2) else {}) for i,j in it.product(range(4),repeat=2))
    zero_lam=lambda poly:{f:v for f,v in poly.items() if not f[0]}
    coframe=R.add(R.mul(a,d),neg(R.mul(b,c)))
    assert zero_lam(pf)==neg(coframe)
    # For d<=5, d(k+1)>=10k forces k<=d/(10-d)<=1.
    ruled=[]
    for degree in range(1,6):
        for fiber_degree in range(1,degree//(10-degree)+1):
            if degree*(fiber_degree+1)>=10*fiber_degree:
                ruled.append([degree,fiber_degree])
    assert ruled==[[5,1]]
    assert all(F(degree,2)<5 for degree in range(1,6))
    # (1-t)^5 * [1+5t/(1-t)^3].
    numerator=mul(mul([1,-1],[1,-1]),mul([1,-1],mul([1,-1],[1,-1])))
    correction=[0]+[5*v for v in mul([1,-1],[1,-1])]
    for i,v in enumerate(correction):numerator[i]+=v
    assert numerator==[1,0,0,-5,5,-1]
    return dict(logarithmic_matrix_upper_entries=['lambda','a','b','c','d','e'],
                inverse_numerator=[[R.encode(v) for v in row] for row in C],
                both_cleared_inverse_products_all_sixteen_entries_passed=True,
                pfaffian='lambda*e-a*d+b*c',
                actual_bivector='diag(u,v,1,1)*C*diag(u,v,1,1)/pfaffian',
                double_surface_tangent_bivector_numerator='-lambda',
                zero_lambda_residue_coframe_determinant='a*d-b*c=-pfaffian',
                complete_positive_integer_RR_degree_solutions=ruled,
                free_quotient_h0_bound_strictly_less_than_five=True,
                scroll_coordinate_ring_Hilbert_numerator=list(map(int,numerator)),
                all_geometric_hypotheses_proved_by_computation=False)


def m2_input():
    x=[f'x{i}' for i in range(5)]
    B=[['0']*5 for _ in range(5)]
    a=['(2*ell^2+ell^(-3))/5','ell','-ell^(-1)']
    b=['(-ell^2+2*ell^(-3))/5','-ell^(-2)','1']
    for i in range(5):
        for step,co in [(1,a),(2,b)]:
            j=(i+step)%5
            extra=(f'{x[(i+2)%5]}*{x[(i+4)%5]}',f'{x[(i+3)%5]}^2') if step==1 else (f'{x[(i+3)%5]}*{x[(i+4)%5]}',f'{x[(i+1)%5]}^2')
            value=f'({co[0]})*{x[i]}*{x[j]}+({co[1]})*{extra[0]}+({co[2]})*{extra[1]}'
            B[i][j]=f'({value})';B[j][i]=f'-({value})'
    lines=['T=QQ[ell0];','K=frac T;','ell=sub(ell0,K);','R=K[x0,x1,x2,x3,x4];']
    for i,j,k in it.combinations(range(5),3):
        jac='+'.join(f'({B[a][l]})*diff({x[l]},({B[b][c]}))' for a,b,c in [(i,j,k),(j,k,i),(k,i,j)] for l in range(5))
        lines.append(f'assert(({jac})==0);')
    generators=[f'{x[i]}*({B[j][k]})-{x[j]}*({B[i][k]})+{x[k]}*({B[i][j]})' for i,j,k in it.combinations(range(5),3)]
    lines.extend(['J=trim ideal('+','.join(generators)+');','J=saturate J;','C=res(coker gens J);',
                  'assert(codim J==2);','assert(degree(R/J)==5);','assert(numgens J==5);',
                  'assert(rank C_0==1 and rank C_1==5 and rank C_2==5 and rank C_3==1);',
                  'assert(degrees C_1=={{3},{3},{3},{3},{3}});',
                  'assert(degrees C_2=={{4},{4},{4},{4},{4}});','assert(degrees C_3=={{5}});',
                  'assert(rank(C.dd_1*C.dd_2)==0 and rank(C.dd_2*C.dd_3)==0);',
                  'assert(codim minors(4,C.dd_2)>=2);','assert(codim minors(1,C.dd_3)>=3);',
                  'bad=saturate(J+minors(2,jacobian gens J));','assert(bad==ideal(1_R));',
                  '<<"GENERIC_FIELD_QQ_ell_PASS"<<endl;',
                  '<<"ZERO_SURFACE_BETTI "<<toString betti C<<endl;',
                  '<<"IDEAL_RESOLUTION 0->R(-5)->R(-4)^5->R(-3)^5->I->0"<<endl;',
                  '<<"BUCHSBAUM_EISENBUD_GRADES_PASS"<<endl;','<<"PROJECTIVE_SMOOTHNESS_PASS"<<endl;',
                  '<<"ZERO_CUBICS "<<toString gens J<<endl;','<<"D2 "<<toString C.dd_2<<endl;',
                  '<<"D3 "<<toString C.dd_3<<endl;','exit 0;'])
    return '\n'.join(lines)+'\n'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--m2',type=Path)
    ap.add_argument('--ordinary-double',action='store_true')
    ap.add_argument('--workdir',type=Path,default=Path('work/nonnormal-quintic-entrance'))
    args=ap.parse_args();start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    result=dict(task_id='E_reduced_nonnormal',finite_checks=ordinary_double_checks() if args.ordinary_double else finite_checks(),
                mathematical_proof='IRREDUCIBLE_NONNORMAL_ORDINARY_DOUBLE_RECOGNITION.md' if args.ordinary_double else 'IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE.md',
                global_exhaustion_proved=False)
    if args.m2:
        args.workdir.mkdir(parents=True,exist_ok=True)
        source=m2_input();input_path=args.workdir/'generic-q52-zero.m2';input_path.write_text(source)
        try:
            run=subprocess.run([str(args.m2),'--script',str(input_path)],capture_output=True,text=True,timeout=110)
        except subprocess.TimeoutExpired as exc:
            decode=lambda value:value.decode(errors='replace') if isinstance(value,bytes) else (value or '')
            run=subprocess.CompletedProcess(args=[],returncode=124,stdout=decode(exc.stdout),stderr=decode(exc.stderr)+'\nFinite110second timeout; geometric entrance proofs remain independent of this optional calculation.\n')
        log=args.workdir/'generic-q52-zero.log';log.write_text(run.stdout+run.stderr)
        result['generic_q52_zero_surface_computation']=dict(returncode=run.returncode,
            input_m2_bytes=source,input_m2_sha256=hashlib.sha256(source.encode()).hexdigest(),
            output_log_bytes=run.stdout+run.stderr,
            exact_parameter_field='QQ(ell), not a sampled elliptic parameter',
            all_assertions_passed=run.returncode==0 and 'PROJECTIVE_SMOOTHNESS_PASS' in run.stdout,
            independent_CAS_replay_done=False)
    result['all_passed']=result.get('generic_q52_zero_surface_computation',{}).get('all_assertions_passed',True)
    result.update(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_bytes=Path(__file__).read_text(),
                  seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('all_passed','seconds','maximum_RSS_KiB')}))
    if not result['all_passed']:raise SystemExit(1)


if __name__=='__main__':main()
