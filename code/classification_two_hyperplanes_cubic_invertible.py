#!/usr/bin/env python3
"""All four nonscalar invertible spectral branches: complete W and SNC sources."""
import argparse,hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def brackets(lam,A,B,C,xx):
    h,x,u,v,w=xx;P=s.zeros(5)
    for ij,f in {(1,2):x*lam[0]*u,(1,3):x*lam[1]*v,(1,4):x*lam[2]*w,
                 (2,3):C,(2,4):-B,(3,4):A}.items():P[ij]=f;P[ij[::-1]]=-f
    return P
def jacobi(P,xx):
    return {ijk:s.expand(sum(P[i,l]*s.diff(P[j,k],xx[l])+P[j,l]*s.diff(P[k,i],xx[l])+P[k,l]*s.diff(P[i,j],xx[l])for l in range(5)))
            for ijk in itertools.combinations(range(5),3)for i,j,k in [ijk]}
def encode(P):return {str(ij):str(P[ij])for ij in itertools.combinations(range(5),2)if P[ij]}
def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();h,x,u,v,w=xx=s.symbols('h x u v w');a,b,c,d,e,f,g=s.symbols('a b c d e f g')
    ah,ax,axx,bh,bx,bxx,lh,lx=s.symbols('ah ax axx bh bx bxx lh lx')
    A2=ah*h*h+ax*h*x+axx*x*x;B2=bh*h*h+bx*h*x+bxx*x*x;L=lh*h+lx*x
    generic=[dict(case='2_1_minus1',lam=(2,1,-1),A=A2+a*v*w,B=v*L+b*u*w,C=c*u*v,
                  expected=(c-a)*L*v*v+(b-c)*A2*u-x*(2*u*s.diff(A2,x)+v*v*s.diff(L,x))),
             dict(case='3_1_minus1',lam=(3,1,-1),A=A2+a*v*w,B=b*v*v+d*u*w,C=c*u*v,
                  expected=b*(c-a)*v**3+(d-c)*A2*u-3*x*u*s.diff(A2,x)),
             dict(case='3_2_1',lam=(3,2,1),A=u*L+a*v*w,B=b*v*v+c*u*w,C=d*u*v,
                  expected=b*(d-a)*v**3+(c-d)*L*u*u-3*x*u*u*s.diff(L,x)),
             dict(case='1_1_minus1',lam=(1,1,-1),A=A2+w*(a*u+b*v),B=B2+w*(c*u+d*v),C=e*u*u+f*u*v+g*v*v,
                  expected=w*((e*(b+c)-a*f)*u*u+2*(d*e-a*g)*u*v+(d*f-(b+c)*g)*v*v)
                    +u*((c-f)*A2+(2*e-a)*B2-x*s.diff(A2,x))
                    +v*((d-2*g)*A2+(f-b)*B2-x*s.diff(B2,x)))]
    full_records=[]
    mons=[xx[i]*xx[j]for i,j in itertools.combinations_with_replacement(range(5),2)]
    for entry in generic:
        lam=entry['lam'];P=brackets(lam,entry['A'],entry['B'],entry['C'],xx);J=jacobi(P,xx)
        assert all(value==0 for key,value in J.items()if key!=(2,3,4))
        assert s.expand(J[(2,3,4)]-entry['expected'])==0
        weight=(0,0,*lam);weight_basis=[]
        for ij in [(2,3),(2,4),(3,4)]:
            actual=set(e for e,value in s.Poly(P[ij],*xx).terms()if value)
            expected=[]
            for mon in mons:
                exponent=s.Poly(mon,*xx).monoms()[0]
                if sum(a*b for a,b in zip(exponent,weight))==weight[ij[0]]+weight[ij[1]]:expected.append(exponent)
            assert actual==set(expected)
            # Every coefficient displayed in each W position is independent.
            parameters=sorted(set().union(*(value.free_symbols for _,value in s.Poly(P[ij],*xx).terms())),key=str)
            matrix,_=s.linear_eq_to_matrix([value for _,value in s.Poly(P[ij],*xx).terms()],parameters)
            assert matrix.rank()==len(expected)
            weight_basis.append(dict(pair=list(ij),monomials=[str(s.prod(v**power for v,power in zip(xx,exponent)))for exponent in expected]))
        F=s.expand(lam[0]*u*entry['A']+lam[1]*v*entry['B']+lam[2]*w*entry['C'])
        full_records.append(dict(case=entry['case'],spectrum=list(lam),complete_weight_basis=weight_basis,
            invariant_vertical_quadratic_space_dimension=sum(len(z['monomials'])for z in weight_basis),
            full_W=encode(P),only_possible_nonzero_Jacobi_coefficient=str(J[(2,3,4)]),
            expected_complete_Jacobi_equation=str(s.expand(entry['expected'])),actual_cubic=str(F)))
    amplitudeA,amplitudeB,tau=s.symbols('amplitudeA amplitudeB tau');families=[]
    def add(case,indices,lam,A,B,C,snc_substitution):
        P=brackets(lam,A,B,C,xx);assert all(value==0 for value in jacobi(P,xx).values())
        F=s.expand(lam[0]*u*A+lam[1]*v*B+lam[2]*w*C)
        snc=s.factor(F.subs(snc_substitution));ratio=s.cancel(snc/(u*v*w))
        assert ratio and not ratio.free_symbols
        families.append(dict(case=case,indices=indices,full_source=encode(P),all_ten_raw_Jacobi_zero=True,
            actual_cubic=str(F),SNC_substitution={str(k):str(val)for k,val in snc_substitution.items()},
            SNC_actual_quintic=str(h*x*snc),nonzero_SNC_coefficient=str(ratio)))
    for j in range(3):
        for ell in range(2):
            AA=amplitudeA*h**(2-j)*x**j;LL=amplitudeB*h**(1-ell)*x**ell
            add('2_1_minus1',[j,ell],(2,1,-1),AA+a*v*w,v*LL+(a+ell+2*j)*u*w,(a+ell)*u*v,
                {amplitudeA:0,amplitudeB:0,a:1})
    for j in range(3):
        AA=amplitudeA*h**(2-j)*x**j
        add('3_1_minus1',[j],(3,1,-1),AA+a*v*w,amplitudeB*v*v+(a+3*j)*u*w,a*u*v,
            {amplitudeA:0,amplitudeB:0,a:1})
    for j in range(2):
        LL=amplitudeA*h**(1-j)*x**j
        add('3_2_1',[j],(3,2,1),u*LL+a*v*w,amplitudeB*v*v+(a+3*j)*u*w,a*u*v,
            {amplitudeA:0,amplitudeB:0,a:1})
    for j,k in itertools.product(range(3),repeat=2):
        AA=amplitudeA*h**(2-j)*x**j;BB=amplitudeB*h**(2-k)*x**k
        add('1_1_minus1_binary_rank2',[j,k],(1,1,-1),AA+(tau-k)*v*w,BB+(tau+j)*u*w,tau*u*v,
            {amplitudeA:0,amplitudeB:0,tau:4})
    q11,q12,q22=s.symbols('q11 q12 q22')
    for j in range(3):
        AA=amplitudeA*h**(2-j)*x**j;BB=amplitudeB*h**(2-j)*x**j
        add('1_1_minus1_scalar_M_all_binary_quadrics',[j],(1,1,-1),
            AA+w*(2*tau*q11*u+(tau*q12-j)*v),
            BB+w*((tau*q12+j)*u+2*tau*q22*v),tau*(q11*u*u+q12*u*v+q22*v*v),
            {amplitudeA:0,amplitudeB:0,tau:1,q11:0,q12:1,q22:0})
    inputs=['code/classification_two_hyperplanes_cubic_invertible.py']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,
        complete_invariant_W_models=full_records,all_parametric_SNC_families=families,
        family_count=len(families),all_original_quadratic_vertical_coefficients_accounted_for=True,
        classification_implication='Together with the explicit spectral and binary-matrix classification, the affine sources above cover every invertible-B nonscalar branch. Every irreducible source contains an SNC point and therefore lies in the proved finite MPS union. Scalar B is covered separately.',
        spectral_necessity_and_branch_exhaustion_require_the_written_proof=True,global_exhaustion_proved=False,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in inputs},
        producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps(dict(all_passed=True,source_count=len(families),seconds=result['seconds'])))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)
