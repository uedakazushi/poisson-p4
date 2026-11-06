#!/usr/bin/env python3
"""Universal rank-three-quadric MPS degenerations for both remaining sources."""
import argparse,hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def encode(A):return {str(ij):str(A[ij])for ij in itertools.combinations(range(5),2)if A[ij]}
def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();h,x,y,z,w=xx=s.symbols('h x y z w');t,alpha,beta,gamma=s.symbols('t alpha beta gamma')
    records=[]
    cases=[dict(case='m2_t1',A=h+y,B=x+y-h,sign=1,p0=-1,r0=2,b0=-2,
                quadric=h*x+x*y+y*y,U=s.Matrix([x*(2*h+y),y*(h+y)]),
                V=s.Matrix([x*(2*h-2*x-3*y),y*(h-y)]),c=y-2*h,d=2*h+y,
                Ra=-x*y-2*y*y,Rb=h*x+x*y,Rg=-2*(h*x+x*y+y*y),cobbasis=[h*h,h*y,x*x+x*y]),
           dict(case='m3_t0',A=h+x,B=h+y,sign=-1,p0=1,r0=-1,b0=1,
                quadric=h*x+h*y+x*y,U=s.Matrix([-x*(h+x),-h*y]),
                V=s.Matrix([h*x,y*(h+y)]),c=-h,d=-h,
                Ra=-h*y-x*y,Rb=h*x+x*y,Rg=h*x+h*y+x*y,cobbasis=[h*h,x*x,y*y])]
    for case in cases:
        started=time.monotonic();A,B=case['A'],case['B'];Z=A+t*z;W=B+t*w
        p=case['p0']+t*alpha;r=case['r0']+t*beta
        D1=s.Matrix([0,x,0,(p*Z-x*s.diff(A,x))/t,(-p*W-x*s.diff(B,x))/t])
        D2=s.Matrix([0,0,y,(r*Z-y*s.diff(A,y))/t,(-r*W-y*s.diff(B,y))/t])
        Q=s.expand(Z*W+case['sign']*h*h)
        P=t*(D1*D2.T-D2*D1.T);P[3,4]+=case['b0']*(1+t*gamma)*Q/t;P[4,3]=-P[3,4]
        P=P.applyfunc(lambda f:s.expand(s.cancel(f)))
        assert all(s.Poly(P[ij],*xx,t,alpha,beta,gamma).total_degree()>=0 for ij in itertools.combinations(range(5),2))
        for i,j,k in itertools.combinations(range(5),3):
            assert s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l])for l in range(5)))==0
        pf=s.expand(P[1,2]*P[3,4]-P[1,3]*P[2,4]+P[1,4]*P[2,3])
        assert s.expand(pf-case['b0']*(1+t*gamma)*x*y*Q)==0
        limit=P.subs(t,0);U,V=case['U'],case['V'];c,d=case['c'],case['d']
        assert all(s.expand(a-b)==0 for a,b in [(limit[1,3],U[0]),(limit[2,3],U[1]),(limit[1,4],V[0]),(limit[2,4],V[1])])
        R=alpha*case['Ra']+beta*case['Rb']+gamma*case['Rg']
        assert s.expand(limit[3,4]+c*z+d*w-R)==0
        bracket=V.jacobian([x,y])*U-U.jacobian([x,y])*V
        assert (bracket-c*U-d*V).applyfunc(s.expand)==s.zeros(2,1)
        assert s.expand(Q.subs(t,0)-case['quadric'])==0
        base=(h,x,y);mons=[h*h,h*x,h*y,x*x,x*y,y*y]
        cob=[s.expand(V[0]*s.diff(phi,x)+V[1]*s.diff(phi,y)+c*phi)for phi in base]
        cob += [s.expand(d*psi-U[0]*s.diff(psi,x)-U[1]*s.diff(psi,y))for psi in base]
        cobmatrix=s.Matrix([[s.Poly(f,*base).coeff_monomial(m)for f in cob]for m in mons]);assert cobmatrix.rank()==3
        cobbasis=s.Matrix([[s.Poly(f,*base).coeff_monomial(m)for f in case['cobbasis']]for m in mons])
        assert cobmatrix.row_join(cobbasis).rank()==3
        six=[case['Ra'],case['Rb'],case['Rg'],*case['cobbasis']]
        matrix=s.Matrix([[s.Poly(f,*base).coeff_monomial(m)for f in six]for m in mons]);determinant=matrix.det();assert determinant
        T=s.Matrix([h,x,y,Z,W]).jacobian(xx);assert s.factor(T.det())==t*t
        H,X,Y,ZZ,WW=s.symbols('H X Y ZZ WW');canonicalQ=ZZ*WW+case['sign']*H*H
        assert (s.hessian(canonicalQ,(H,X,Y,ZZ,WW))/2).rank()==3
        records.append(dict(case=case['case'],all_parameters=['t','alpha','beta','gamma'],
            coordinate_change=[str(f)for f in [h,x,y,Z,W]],coordinate_change_determinant='t**2',
            canonical_weights_D1=['0','1','0','p','-p'],canonical_weights_D2=['0','0','1','r','-r'],
            weight_substitution={'p':str(p),'r':str(r)},canonical_quadric=str(canonicalQ),vertical_coefficient=str(case['b0']*t*(1+t*gamma)),
            full_polynomial_family=encode(P),all_ten_universal_raw_Jacobi_zero=True,no_negative_t_powers=True,
            actual_quintic=str(case['b0']*h*(1+t*gamma)*x*y*Q),t_nonzero_quadric_rank=3,
            rank4_open_condition='t*(1+t*gamma) != 0',limit_source=encode(limit),limit_quadric=str(case['quadric']),
            limit_magnetic_terms=list(map(str,[case['Ra'],case['Rb'],case['Rg']])),
            anchor_bracket_coefficients=list(map(str,[c,d])),all_six_translation_coboundaries=list(map(str,cob)),
            translation_coboundary_rank=3,full_magnetic_span_six_generators=list(map(str,six)),
            magnetic_coefficient_order=list(map(str,mons)),full_magnetic_span_matrix=[[str(v)for v in matrix.row(i)]for i in range(6)],
            full_magnetic_span_determinant=str(determinant),all_six_magnetic_parameters_reached=True,seconds=time.monotonic()-started))
        print(json.dumps(dict(case=case['case'],all_passed=True,magnetic_determinant=str(determinant),seconds=time.monotonic()-started)),flush=True)
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,cases=records,
        classification_implication='For t*(1+t*gamma) nonzero each family lies in the proved nonzero-residue rank3-quadric MPS support. Every magnetic R is a limit up to a projective fiber translation. Closedness of the finite MPS union covers both entire sources.',
        global_exhaustion_proved=False,input_sha256={'code/classification_three_hyperplanes_quadric_m23_mps.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)
