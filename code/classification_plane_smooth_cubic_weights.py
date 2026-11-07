#!/usr/bin/env python3
"""Exact full normal-weight arithmetic for a plane with smooth cubic restriction.

All ten cubic coefficients, all30 weight1 and all9 weight2 coefficients are
symbolic. Geometry (connection rigidity and smooth-cubic injectivity) is
separate; this checks every Schouten/Euler coefficient used in the proof.
"""
from pathlib import Path
from itertools import combinations,combinations_with_replacement
import hashlib,json,resource,signal,time
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(combinations(range(5),2));TRIPLES=list(combinations(range(5),3));QUADS=list(combinations(range(5),4))


def monomials(n,d):
    result=[]
    for terms in combinations_with_replacement(range(n),d):
        e=[0]*n
        for i in terms:e[i]+=1
        result.append(tuple(e))
    return sorted(result,reverse=True)


def wedge(first,second):return first*second.T-second*first.T


def bilinear(P,Q,x):
    return {tri:s.expand(sum(P[a,l]*s.diff(Q[b,c],x[l])+Q[a,l]*s.diff(P[b,c],x[l])for a,b,c in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(5)))for tri in TRIPLES}


def projective(strict,x):
    return {quad:s.expand(sum((-1)**i*x[quad[i]]*strict[quad[:i]+quad[i+1:]]for i in range(4)))for quad in QUADS}


def all350(jac,x):
    return [str(s.Poly(jac[quad],*x).coeff_monomial(e))for quad in QUADS for e in monomials(5,4)]


def main():
    destination=ROOT/'reports/CLASSIFICATION_PLANE_SMOOTH_CUBIC_NORMAL_WEIGHTS_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    start=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    x=s.Matrix(s.symbols('x0 x1 x2 z w'));base=s.Matrix([*x[:3],0,0]);z,w=x[3:];a=s.Symbol('a');coeff=s.symbols('f0:10')
    f=sum(c*s.prod(x[i]**e[i]for i in range(3))for c,e in zip(coeff,monomials(3,3)))
    def gradient(g):
        P=s.zeros(5)
        for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:P[i,j]=s.diff(g,x[k]);P[j,i]=-P[i,j]
        return P
    Pf=gradient(f)
    def linear(name):
        symbols=s.symbols(name+'0:9');P=s.zeros(5)
        for k,(i,j)in enumerate(combinations(range(3),2)):
            P[i,j]=sum(symbols[3*k+l]*x[l]for l in range(3));P[j,i]=-P[i,j]
        return P,symbols
    L,lp=linear('l');M,mp=linear('m');vectors={name:s.Matrix([*s.symbols(name.lower()+'0:3'),0,0])for name in ['B','C','D','E']}
    B,C,D,E=[vectors[name]for name in ['B','C','D','E']];ez=s.Matrix([0,0,0,1,0]);ew=s.Matrix([0,0,0,0,1])
    P1=z*L+w*M+z*z*wedge(D,ew)+z*w*(wedge(B,ez)+wedge(E,ew))+w*w*wedge(C,ez)
    def U(V):return Pf.applyfunc(lambda value:-sum(V[i]*s.diff(value,x[i])for i in range(3)))
    Uv={name:U(V)for name,V in vectors.items()}
    # Smoothness supplies injectivity of V->V(f). The identity is polynomial for every cubic.
    for name,V in vectors.items():
        assert s.expand(x[0]*Uv[name][1,2]-x[1]*Uv[name][0,2]+x[2]*Uv[name][0,1]+2*sum(V[i]*s.diff(f,x[i])for i in range(3)))==0
    base_strict_L=bilinear(Pf,L,x)[0,1,2];base_strict_M=bilinear(Pf,M,x)[0,1,2]
    # Curl criterion has no numerical cubic specialization.
    ell=s.Matrix([L[1,2],-L[0,2],L[0,1]]);curl=s.Matrix([s.diff(ell[2],x[1])-s.diff(ell[1],x[2]),s.diff(ell[0],x[2])-s.diff(ell[2],x[0]),s.diff(ell[1],x[0])-s.diff(ell[0],x[1])])
    curl_dot=s.expand(sum(curl[i]*s.diff(f,x[i])for i in range(3)))
    assert s.expand(base_strict_L+curl_dot)==0
    q2=s.symbols('q0:6');g=sum(c*s.prod(x[i]**e[i]for i in range(3))for c,e in zip(q2,monomials(3,2)))
    h2=s.symbols('h0:6');h=sum(c*s.prod(x[i]**e[i]for i in range(3))for c,e in zip(h2,monomials(3,2)))
    records=[]
    for name,K in [('two_roots',a*z*w),('double_root',a*z*z),('zero',s.Integer(0))]:
        P0=Pf+K*wedge(ez,ew);raw=bilinear(P0,P1,x);jac=projective(raw,x)
        if name=='two_roots':normal=a*z**3*D+a*w**3*C
        elif name=='double_root':normal=-a*z**3*E+a*z*z*w*B+2*a*z*w*w*C
        else:normal=s.zeros(5,1)
        Tz=K*M+z*w*Uv['B']+w*w*Uv['C'];Tw=-K*L+z*z*Uv['D']+z*w*Uv['E']
        expected={tri:s.Integer(0)for tri in TRIPLES};expected[0,1,2]=z*base_strict_L+w*base_strict_M
        for i,j in combinations(range(3),2):expected[i,j,3]=Tz[i,j];expected[i,j,4]=Tw[i,j]
        for i in range(3):expected[i,3,4]=normal[i]
        assert all(s.expand(raw[tri]-expected[tri])==0 for tri in TRIPLES)
        forcing=wedge(base,normal)+z*Tw-w*Tz
        for i,j in combinations(range(3),2):assert s.expand(jac[i,j,3,4]-forcing[i,j])==0
        if name=='two_roots':closed_L=Uv['E']/a;closed_M=-Uv['B']/a;closed_B=B;closed_C=s.zeros(5,1);closed_D=s.zeros(5,1);closed_E=E
        elif name=='double_root':closed_L=Uv['D']/a-wedge(base,E);closed_M=Uv['E']/a;closed_B=s.zeros(5,1);closed_C=s.zeros(5,1);closed_D=D;closed_E=E
        else:closed_L=gradient(g);closed_M=gradient(h);closed_B=closed_C=closed_D=closed_E=s.zeros(5,1)
        kernel=z*closed_L+w*closed_M+z*z*wedge(closed_D,ew)+z*w*(wedge(closed_B,ez)+wedge(closed_E,ew))+w*w*wedge(closed_C,ez)
        assert all(s.expand(value)==0 for value in projective(bilinear(P0,kernel,x),x).values())
        if K!=0:
            Y=s.Matrix([*s.symbols('y0:3'),0,0]);Z=s.Matrix([*s.symbols('v0:3'),0,0]);field=z*Y+w*Z
            orbit=s.zeros(5)
            for i,j in PAIRS:
                orbit[i,j]=s.expand(-sum(field[k]*s.diff(P0[i,j],x[k])for k in range(5))+sum(P0[k,j]*s.diff(field[i],x[k])+P0[i,k]*s.diff(field[j],x[k])for k in range(5)));orbit[j,i]=-orbit[i,j]
            first=s.Matrix([*[s.Poly(orbit[i,3],z,w).coeff_monomial(z*z)for i in range(3)],0,0]);second=s.Matrix([*[s.Poly(orbit[i,4],z,w).coeff_monomial(w*w)for i in range(3)],0,0])
            canonical_orbit=orbit+wedge(x,z*first+w*second)
            replacement={}
            if name=='two_roots':replacement.update({B[i]:-a*Z[i]for i in range(3)});replacement.update({E[i]:a*Y[i]for i in range(3)})
            else:replacement.update({D[i]:a*Y[i]for i in range(3)});replacement.update({E[i]:a*Z[i]for i in range(3)})
            assert all(s.expand(value)==0 for value in canonical_orbit-kernel.subs(replacement,simultaneous=True))
        t2=s.symbols('r0:9');P2=s.zeros(5)
        for k,(i,j)in enumerate(combinations(range(3),2)):
            P2[i,j]=t2[3*k]*z*z+t2[3*k+1]*z*w+t2[3*k+2]*w*w;P2[j,i]=-P2[i,j]
        weight2=projective(bilinear(P0,P2,x),x)
        for i,j in combinations(range(3),2):assert s.expand(weight2[i,j,3,4]+2*K*P2[i,j])==0
        actual=s.expand(sum((-1)**i*x[i]*(P0[j,k]*P0[l,m]-P0[j,l]*P0[k,m]+P0[j,m]*P0[k,l])for i in range(5)for j,k,l,m in [[t for t in range(5)if t!=i]]))
        assert s.expand(actual-3*f*K)==0
        records.append({'normal_quadratic':name,'all350_full_weight1_Jacobi':all350(jac,x),'complete30_source':{str(i)+','+str(j):str(s.expand(P1[i,j]))for i,j in PAIRS if P1[i,j]!=0},
            'strict_full_firstweight_formula_verified':True,'all_base2_normal2_forcing_coefficients':{str(i)+','+str(j):str(s.expand(forcing[i,j]))for i,j in combinations(range(3),2)},
            'all350_kernel_substitution_zero':True,'complete_kernel_representative':{str(i)+','+str(j):str(s.expand(kernel[i,j]))for i,j in PAIRS if kernel[i,j]!=0},
            'actual_shear_orbit_equals_full_kernel':K!=0,'all9_weight2_injectivity_formula_verified':True,'all350_weight2_Jacobi':all350(weight2,x),'actual_D_equals_three_f_K':True})
        print(name,'all350/full30/kernel/shear/weight2 exact PASS',flush=True)
    result={'task_id':'G','all_passed':True,'all10_cubic_coefficients_symbolic':True,'smoothness_used_only_in_geometric_constant_direction_injectivity':True,
        'nonzero_a_localization_only_for_two_nonzero_binary_quadratic_types':True,'cubic':str(f),'normal_filtration_complete30_and9':True,'cases':records,
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'seconds':time.monotonic()-start,
        'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'scope':'Exact full Schouten/Euler arithmetic with arbitrary ternary cubic coefficients. All-smooth injectivity, conormal scalar rigidity, source completeness and finite shear geometry are separate proofs in the note; no numerical Hesse specialization or generic rank argument.'}
    destination.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
