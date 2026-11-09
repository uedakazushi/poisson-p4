#!/usr/bin/env python3
"""Exact identities for the cofactor reconstruction and magnetic residue theorem.

Uses complete frozen original sources, not their producers or PASS receipts.
Universal first jets and full polynomial families are checked without sampling.
Full150 linear rank completeness and componenthood remain old proof inputs.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def wedge(a,b):
    out={}
    for ka,va in a.items():
        for kb,vb in b.items():
            word=ka+kb
            if len(set(word))!=len(word):
                continue
            sign=(-1)**sum(word[i]>word[j] for i in range(len(word)) for j in range(i+1,len(word)))
            key=tuple(sorted(word));out[key]=s.expand(out.get(key,0)+sign*va*vb)
    return out


def contract(form,vector):
    out={}
    for key,value in form.items():
        for at,k in enumerate(key):
            target=key[:at]+key[at+1:]
            out[target]=s.expand(out.get(target,0)+(-1)**at*vector[k]*value)
    return out


def dual_multivector(coefficients,n):
    volume={tuple(range(n)):s.Integer(1)}
    result={}
    for key,value in coefficients.items():
        term=volume
        for k in key:
            vector=[s.Integer(i==k) for i in range(n)]
            term=contract(term,vector)
        for target,coefficient in term.items():
            result[target]=s.expand(result.get(target,0)+value*coefficient)
    return result


def exterior(form,derivative,n):
    result={}
    for key,value in form.items():
        term=wedge({(k,):derivative(value,k) for k in range(n)},{key:1})
        for target,coefficient in term.items():
            result[target]=s.expand(result.get(target,0)+coefficient)
    return result


def skew(values,n):
    p=s.zeros(n)
    for (i,j),value in zip(itertools.combinations(range(n),2),values):
        p[i,j]=value;p[j,i]=-value
    return p


def div(p,x):
    return s.Matrix([sum(s.diff(p[i,j],x[i]) for i in range(len(x))) for j in range(len(x))])


def jac(p,x):
    return {key:s.expand(sum(p[i,l]*s.diff(p[j,k],x[l])+p[j,l]*s.diff(p[k,i],x[l])
        +p[k,l]*s.diff(p[i,j],x[l]) for l in range(len(x))))
        for key in itertools.combinations(range(len(x)),3) for i,j,k in [key]}


def pfive(p):
    out=[]
    for at in range(5):
        i,j,k,l=[i for i in range(5) if i!=at]
        out.append(s.expand((-1)**at*(p[i,j]*p[k,l]-p[i,k]*p[j,l]+p[i,l]*p[j,k])))
    return s.Matrix(out)


def run():
    names=[]
    def check(name,values):
        for value in values:
            if s.cancel(s.expand(value))!=0:
                raise AssertionError(f'{name}: {s.expand(value)}')
        names.append(name)
    n=5;pairs=list(itertools.combinations(range(n),2))
    pp=s.symbols('p0:10');p=skew(pp,n)
    dp={v:s.symbols(f'd{k}_0:5') for k,v in enumerate(pp)}
    def derivative(f,k):
        return sum(s.diff(f,v)*dp[v][k] for v in pp)
    z=s.Matrix([sum(derivative(p[i,j],i) for i in range(n)) for j in range(n)])
    alpha=dual_multivector(dict(zip(pairs,pp)),n)
    beta={(i,):v for i,v in enumerate(pfive(p))}
    jj={key:s.expand(sum(p[i,l]*derivative(p[j,k],l)+p[j,l]*derivative(p[k,i],l)
        +p[k,l]*derivative(p[i,j],l) for l in range(n)))
        for key in itertools.combinations(range(n),3) for i,j,k in [key]}
    d_beta=exterior(beta,derivative,n);iz_alpha=contract(alpha,z);dual_j=dual_multivector(jj,n)
    check('universal_five_dimensional_BV_identity',
        [d_beta.get(key,0)+iz_alpha.get(key,0)-dual_j.get(key,0) for key in pairs])
    x=s.symbols('x0:5');coeff=s.symbols('c0:150');mons=[x[i]*x[j] for i in range(n) for j in range(i,n)]
    generic=skew([sum(coeff[15*k+j]*m for j,m in enumerate(mons)) for k in range(10)],n)
    E=s.Matrix(x)
    def canonical(a):
        zz=div(a,x)
        return a-(E*zz.T-zz*E.T)/5
    canonical_generic=canonical(generic)
    check('full150_canonical_divergence',div(canonical_generic,x))
    check('full150_canonical_idempotence',canonical(canonical_generic)-canonical_generic)
    lin=s.Matrix(5,5,s.symbols('a0:25'));Y=lin*E
    gauge=E*Y.T-Y*E.T
    check('full25_Euler_gauge_kernel',canonical(gauge))

    # All31 fixed-cofactor original equations; do not infer their linear rank anew.
    path=ROOT/'code/data/quadric_cofactor_sources.json'
    frozen=json.loads(path.read_text());source_records=[]
    assert len(frozen['cases'])==31
    local=dict(zip(map(str,x),x))
    for entry in frozen['cases']:
        spec=entry['case'];ns=s.symbols('n0:'+str(entry['source_dimension']))
        loc=local|dict(zip(map(str,ns),ns))
        q=s.sympify(spec['quadric'],locals=loc)
        xx=s.Matrix([s.sympify(v,locals=loc) for v in spec['cofactor_components']]);N=xx.jacobian(x)
        po=s.zeros(5)
        for key,value in entry['brackets'].items():
            i,j=map(int,key);po[i,j]=s.sympify(value,locals=loc);po[j,i]=-po[i,j]
        check(spec['name']+'_full_source_linear_equations',list(div(po,x))+
            list(po.T*s.Matrix([s.diff(q,v) for v in x])-q*xx)+
            list(po.applyfunc(lambda f:sum(xx[i]*s.diff(f,x[i]) for i in range(5)))-N*po-po*N.T))
        raw=jac(po,x)
        actual={s.srepr(s.expand(v)) for f in raw.values() for v in s.Poly(f,*x).coeffs() if v!=0}
        reported={s.srepr(s.expand(s.sympify(v,locals=loc))) for v in entry['residual_jacobi_coefficients']}
        assert actual==reported,spec['name']
        potential=s.expand((E.T*pfive(po))[0]/5)
        check(spec['name']+'_potential_before_Jacobi',[potential-s.sympify(entry['potential_before_Jacobi'],locals=loc)])
        source_records.append({'case':spec['name'],'quadric_rank':spec['quadric_rank'],
             'source_dimension':len(ns),'all_source_linear_equations':'PASS','all_residual_Jacobi_coefficients':'PASS',
             'source_complete_linear_rank_independently_recomputed':False})

    a,b,c,r,ss=s.symbols('a b c r ss');base=(a,b,c);coords=base+(r,ss)
    qs=[225*a*a+24*a*b+1875*a*c+b*b+1000*b*c-3125*c*c,
        -3*a*a+3*a*c+b*b-c*c,6*a*a-4*a*b+a*c+b*b-4*b*c,
        3*a*a-3*a*b+b*b-b*c,2*a*a+3*a*c+b*b+c*c,7*a*a+4*a*b+8*a*c+b*b+8*b*c]
    cubics=[b*b*c-a**3]*4+[b*b*c-a*a*(a+c)]*2
    examples=[('conic_cubic_'+str(i+1),[q,cc]) for i,(q,cc) in enumerate(zip(qs,cubics))]
    examples += [('two_lines_cubic',[c,-a+b,b*b*c-a*a*(a+c)]),
                 ('integral_quintic',[a**5+a*a*b**3+b**4*c])]
    examples_out=[];m2=[a*a,a*b,a*c,b*b,b*c,c*c]
    for label,factors in examples:
        F=s.prod(factors);dF=s.Matrix([s.diff(F,v) for v in base])
        un=s.symbols('u0:18')
        uv=s.Matrix([sum(un[6*i+j]*mon for j,mon in enumerate(m2)) for i in range(3)])
        mat,_=s.linear_eq_to_matrix(s.Poly(s.expand((uv.T*dF)[0]),*base).coeffs(),un)
        kernel=mat.nullspace();assert len(kernel)==2,label
        U=uv.subs(dict(zip(un,kernel[0])));V=uv.subs(dict(zip(un,kernel[1])))
        kap=s.cancel(s.Matrix.hstack(s.Matrix(base),U,V).det()/F)
        assert kap!=0 and not kap.free_symbols,(label,kap)
        du=sum(s.diff(U[i],base[i]) for i in range(3));dv=sum(s.diff(V[i],base[i]) for i in range(3))
        ell=s.symbols('l0:3');mm=s.symbols('m0:3')
        L=sum(ell[i]*base[i] for i in range(3));M=sum(mm[i]*base[i] for i in range(3))
        delta=s.expand(sum(s.diff(L*V[i]-M*U[i],base[i]) for i in range(3)))
        delta_matrix,_=s.linear_eq_to_matrix([s.Poly(delta,*base).coeff_monomial(mon) for mon in m2],ell+mm)
        rank=delta_matrix.rank();assert rank==7-len(factors),(label,rank)
        residue_vectors=[]
        for fi in factors[:-1]:
            flast=factors[-1];di=s.Poly(fi,*base).total_degree();dj=s.Poly(flast,*base).total_degree()
            eta=s.Matrix([s.diff(fi,v)/fi-s.Rational(di,dj)*s.diff(flast,v)/flast for v in base])
            lu=s.cancel((eta.T*U)[0]);mv=s.cancel((eta.T*V)[0])
            assert s.Poly(lu,*base).total_degree()<=1 and s.Poly(mv,*base).total_degree()<=1
            vector=s.Matrix([s.expand(lu).coeff(v) for v in base]+[s.expand(mv).coeff(v) for v in base])
            check(label+'_residue_shear_kernel',delta_matrix*vector)
            residue_vectors.append(vector)
            W=s.simplify(lu*V-mv*U)
            theta=dual_multivector({(i,):W[i] for i in range(3)},3)
            contracted=contract(theta,base)
            check(label+'_residue_map_inverse',[contracted.get((i,),0)-kap*F*eta[i] for i in range(3)])
        assert (s.Matrix.hstack(*residue_vectors).rank() if residue_vectors else 0)==len(factors)-1
        rr=s.symbols('R0:6');R=sum(v*mon for v,mon in zip(rr,m2))
        pi=s.zeros(5)
        for i in range(3):
            pi[i,3]=U[i];pi[3,i]=-U[i];pi[i,4]=V[i];pi[4,i]=-V[i]
        pi[3,4]=-dv*r+du*ss+R;pi[4,3]=-pi[3,4]
        check(label+'_whole_magnetic_source',list(div(pi,coords))+list(jac(pi,coords).values())+
              list(pfive(pi)+kap*s.Matrix([s.diff(F,v) for v in coords])/5))
        change=s.Matrix([a,b,c,r+L,ss+M]);dchange=change.jacobian(coords)
        inverse={r:r-L,ss:ss-M}
        pushed=(dchange*pi*dchange.T).subs(inverse,simultaneous=True)
        target=pi.copy();target[3,4]+=delta;target[4,3]-=delta
        check(label+'_all_shear_forward_and_inverse',pushed-target)
        examples_out.append({'name':label,'factor_count':len(factors),'shear_rank':rank,
           'residue_kernel_rank':len(residue_vectors),'magnetic_quotient_dimension':6-rank,
           'all_six_magnetic_parameters_retained':True,'Saito_determinant_ratio':str(kap)})
    return {'date':'2026-11-09','version':'v69-quadric-residue-reconstruction','status':'PASS',
        'scope':'Universal BV identity, full150 canonical projection, all31 original coefficient sources and residual Jacobi, and residue kernels/full magnetic sources in eight configurations.',
        'checks':names,'quadric_sources':source_records,'plane_examples':examples_out,
        'frozen_source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'old_producers_imported':False,'sampling_used':False,'classification147_verified':False,
        'limits':['Full150 linear rank completeness is an unchanged original proof input.',
                  'No ambient componenthood or nonlinear all-source recognition is inferred from these checks.']}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'quadric_sources':len(result['quadric_sources']),
        'plane_examples':result['plane_examples'],'output':str(args.output)},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
