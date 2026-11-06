#!/usr/bin/env python3
"""Universal three-parameter MPS degeneration covering every m1 magnetic term."""
import argparse,hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def encode(A):return {str(ij):str(A[ij])for ij in itertools.combinations(range(5),2)if A[ij]}
def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();h,x,y,z,w=xx=s.symbols('h x y z w');t,alpha,beta,gamma=s.symbols('t alpha beta gamma')
    u=x+y;Z=t*z+2*u;W=t*(w-z/2)+h+u;p=t*alpha;q=s.Rational(1,2)+t*beta
    D1=s.Matrix([0,x,0,(p*Z-2*x)/t,((1-p)*W+p*Z/2-2*x)/t])
    D2=s.Matrix([0,0,y,(q*Z-2*y)/t,(-q*W+q*Z/2-2*y)/t])
    Q=s.expand(Z*W/2-h*x)
    P=t*(D1*D2.T-D2*D1.T);P[3,4]-=(1+t*gamma)*Q/t;P[4,3]=-P[3,4]
    P=P.applyfunc(lambda f:s.expand(s.cancel(f)))
    assert all(s.Poly(P[ij],*xx,t,alpha,beta,gamma).total_degree()>=0 for ij in itertools.combinations(range(5),2))
    for i,j,k in itertools.combinations(range(5),3):
        value=s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l])for l in range(5)))
        assert value==0
    pf=s.expand(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3])
    assert s.expand(pf+(1+t*gamma)*x*y*Q)==0
    limit=P.subs(t,0);q0=(x+y)**2+h*y
    U=s.Matrix([x*(x-y),2*x*y]);V=s.Matrix([-x*(h+4*y)/2,-y*(h-x+y)])
    assert all(s.expand(a-b)==0 for a,b in [(limit[1,3],U[0]),(limit[2,3],U[1]),(limit[1,4],V[0]),(limit[2,4],V[1])])
    Ralpha=-2*h*y-4*x*y-4*y*y
    Rbeta=2*x*x-2*y*y-2*h*y
    Rgamma=-q0
    R=alpha*Ralpha+beta*Rbeta+gamma*Rgamma
    assert s.expand(limit[3,4]+(h+y)*z+x*w-R)==0
    c=h+y;d=x;base=(h,x,y);cob=[]
    for phi in base:
        cob.append(s.expand(V[0]*s.diff(phi,x)+V[1]*s.diff(phi,y)+c*phi))
    for psi in base:
        cob.append(s.expand(d*psi-U[0]*s.diff(psi,x)-U[1]*s.diff(psi,y)))
    mons=[h*h,h*x,h*y,x*x,x*y,y*y]
    cobmatrix=s.Matrix([[s.Poly(f,*base).coeff_monomial(m)for f in cob]for m in mons])
    assert cobmatrix.rank()==3
    six=[Ralpha,Rbeta,Rgamma,h*h+h*y,h*x,x*y]
    matrix=s.Matrix([[s.Poly(f,*base).coeff_monomial(m)for f in six]for m in mons])
    determinant=matrix.det();assert determinant!=0
    # The coordinate change is invertible for t nonzero, and the canonical
    # quadratic ZW/2-hx has rank four, with no specialization of alpha,beta,gamma.
    T=s.Matrix([h,x,y,Z,W]).jacobian(xx);assert s.factor(T.det())==t*t
    H,X,Y,ZZ,WW=s.symbols('H X Y ZZ WW');canonicalQ=ZZ*WW/2-H*X
    assert (s.hessian(canonicalQ,(H,X,Y,ZZ,WW))/2).rank()==4
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,case='m1_t2',
        all_parameters=['t','alpha','beta','gamma'],coordinate_change=[str(f)for f in [h,x,y,Z,W]],
        coordinate_change_determinant='t**2',canonical_weights_D1=['0','1','0','p','1-p'],
        canonical_weights_D2=['0','0','1','q','-q'],weight_substitution={'p':str(p),'q':str(q)},
        full_polynomial_family=encode(P),all_ten_universal_raw_Jacobi_zero=True,
        no_negative_t_powers=True,actual_quintic=str(-h*(1+t*gamma)*x*y*Q),
        t_nonzero_quadric_rank=4,rank4_open_condition='t*(1+t*gamma) != 0',
        limit_source=encode(limit),limit_quadric=str(q0),limit_magnetic_terms=list(map(str,[Ralpha,Rbeta,Rgamma])),
        all_six_translation_coboundaries=list(map(str,cob)),translation_coboundary_rank=3,
        full_magnetic_span_six_generators=list(map(str,six)),magnetic_coefficient_order=list(map(str,mons)),
        full_magnetic_span_matrix=[[str(v)for v in matrix.row(i)]for i in range(6)],
        full_magnetic_span_determinant=str(determinant),all_six_magnetic_parameters_reached=True,
        classification_implication='For t*(1+t*gamma) nonzero the family lies in the proved nonzero-residue rank4-quadric MPS support. Every magnetic R is a limit up to an explicit projective fiber translation. Closedness of the finite MPS union covers the entire m1 source.',
        global_exhaustion_proved=False,input_sha256={'code/classification_three_hyperplanes_quadric_m1_mps.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({'all_passed':True,'magnetic_determinant':str(determinant),'seconds':result['seconds']}))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)
