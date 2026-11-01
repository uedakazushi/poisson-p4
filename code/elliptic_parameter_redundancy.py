#!/usr/bin/env python3
"""Exact normalizer action on both elliptic parameter curves.

Arithmetic is in QQ[zeta]/Phi_5 and its rational function field.  The upper
bound on all projective equivalences uses the geometric proof in
research/notes/ELLIPTIC_PARAMETER_REDUNDANCY_20261101.md, not this computation.
The resulting generic representation degree3000 is not either W16 fiber degree.
"""
import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import sympy as s
from itertools import combinations, combinations_with_replacement
from sympy.polys.domains import QQ
import json

ROOT = Path(__file__).resolve().parents[1]

z = s.Symbol('zeta')
K = QQ.algebraic_field((s.Poly(sum(z**i for i in range(5)), z), z))
z5 = K.unit
zero, one = K.zero, K.one
PAIRS=list(combinations(range(5),2))
MONS=list(combinations_with_replacement(range(5),2))
EXTRACT=[(0,1,0,1),(0,1,2,4),(0,1,3,3),(0,2,0,2),(0,2,3,4),(0,2,1,1)]

def template(b):
    out={}
    for i in range(5):
        for delta,terms in [(1,[(0,i,(i+1)%5),(1,(i+2)%5,(i+4)%5),(2,(i+3)%5,(i+3)%5)]),(2,[(3,i,(i+2)%5),(4,(i+3)%5,(i+4)%5),(5,(i+1)%5,(i+1)%5)])]:
            j=(i+delta)%5
            for b0,r,t in terms:
                if b0!=b:continue
                a,c=sorted((i,j));u,v=sorted((r,t))
                out[a,c,u,v]=one if i<j else -one
    return out

def push_basis(A,B,b):
    ans={}
    for (r,t,u,v),coef in template(b).items():
        for i,j in PAIRS:
            wedge=A[i][r]*A[j][t]-A[i][t]*A[j][r]
            if not wedge:continue
            for a,c in MONS:
                subst=B[u][a]*B[v][c]
                if a!=c:subst+=B[u][c]*B[v][a]
                if subst:
                    key=i,j,a,c
                    ans[key]=ans.get(key,zero)+coef*wedge*subst
    return {k:v for k,v in ans.items() if v}

def action(A,B):
    matrix=[[zero for _ in range(6)] for _ in range(6)]
    for b in range(6):
        out=push_basis(A,B,b)
        for i,key in enumerate(EXTRACT):matrix[i][b]=out.get(key,zero)
        recon={}
        for i in range(6):
            for key,val in template(i).items():recon[key]=recon.get(key,zero)+val*matrix[i][b]
        assert out=={k:v for k,v in recon.items() if v}
    return matrix

def mul(A,B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(B))),zero) for j in range(len(B[0]))] for i in range(len(A))]

F=K.frac_field('lam');lam=F.gens[0]
def coeff(k,l):
    if k==1:return [(-3*l**2+l**-3)/5,-2/l,l**-2,(-l**2-3*l**-3)/5,F(2),l]
    return [(2*l**2+l**-3)/5,l,-1/l,(-l**2+2*l**-3)/5,-l**-2,F(1)]

fourier=[[z5**(i*j%5) for j in range(5)] for i in range(5)]
invfourier=[[z5**(-i*j%5)/5 for j in range(5)] for i in range(5)]
phase=[[z5**(i*i%5) if i==j else zero for j in range(5)] for i in range(5)]
invphase=[[z5**(-i*i%5) if i==j else zero for j in range(5)] for i in range(5)]
cycle=[[one if i==(j+1)%5 else zero for j in range(5)] for i in range(5)]
invcycle=[[one if j==(i+1)%5 else zero for j in range(5)] for i in range(5)]
diagonal=[[z5**i if i==j else zero for j in range(5)] for i in range(5)]
invdiagonal=[[z5**(-i%5) if i==j else zero for j in range(5)] for i in range(5)]
reflection=[[one if i==(-j)%5 else zero for j in range(5)] for i in range(5)]

def serial_matrix(M):
    return [[str(K.to_sympy(t)) for t in row] for row in M]

def canonical(M):
    flat=tuple(v for row in M for v in row)
    pivot=next(v for v in flat if v)
    return tuple(v/pivot for v in flat)

def enumerate_mobius(generators):
    identity=((one,zero),(zero,one))
    found={canonical(identity): ''}
    queue=[identity]
    for old in queue:
        word=found[canonical(old)]
        for name,new in generators.items():
            product=mul(new,old)
            key=canonical(product)
            if key not in found:
                assert len(found)<60
                found[key]=name+word
                queue.append((key[:2],key[2:]))
    assert len(found)==60
    return [{'word':word,'matrix':serial_matrix((key[:2],key[2:]))}
            for key,word in found.items()]

def enumerate_sl2():
    gens={'F':(0,1,4,0),'Q':(1,2,0,1)}
    def product(a,b):
        return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2))%5
                     for i in range(2) for j in range(2))
    found={(1,0,0,1):''}
    queue=[(1,0,0,1)]
    for old in queue:
        for name,new in gens.items():
            key=product(new,old)
            if key not in found:
                found[key]=name+found[old]
                queue.append(key)
    assert len(found)==120
    assert all((a*d-b*c)%5==1 for a,b,c,d in found)
    return [{'matrix':list(key),'word':word} for key,word in found.items()]

def projection_rank():
    """The six homogeneous templates remain independent on projective space."""
    x=s.symbols('x0:5')
    columns=[]
    for b in range(6):
        B=s.zeros(5)
        for (i,j,r,t),value in template(b).items():
            B[i,j]+=s.Integer(1 if value==one else -1)*x[r]*x[t]
            B[j,i]=-B[i,j]
        column={}
        for i,j in combinations(range(1,5),2):
            affine=(B[i,j]-x[i]*B[0,j]+x[j]*B[0,i]).subs(x[0],1)
            for mon,c in s.Poly(affine,*x[1:]).terms():column[(i,j,*mon)]=c
        columns.append(column)
    keys=sorted(set().union(*(set(c) for c in columns)))
    matrix=s.Matrix([[c.get(key,0) for c in columns] for key in keys])
    pivot=list(matrix.T.rref()[1])
    assert len(pivot)==6
    minor=matrix[pivot,:]
    return {'rank':6,'minor_rows':[list(keys[i]) for i in pivot],
            'minor':[[str(v) for v in row] for row in minor.tolist()],
            'determinant':str(minor.det())}

def certify():
    phi=-z5**2-z5**3
    assert phi**2-phi-one==zero
    mobius={1:{'F':[[-phi,one],[one,phi]],'Q':[[z5**3,zero],[zero,one]]},
            2:{'F':[[phi-one,one],[one,one-phi]],'Q':[[z5**4,zero],[zero,one]]}}
    gammas={1:(-7*z5**3+z5**2-6*z5-3)/25,
            2:(z5**3+7*z5**2+8*z5+4)/25}
    result={}
    for name,A,B in [('F',fourier,invfourier),('Q',phase,invphase),
                     ('C',cycle,invcycle),('D',diagonal,invdiagonal),
                     ('H',reflection,reflection)]:
        assert mul(A,B)==[[one if i==j else zero for j in range(5)] for i in range(5)]
        M=action(A,B)
        data={'projective_matrix':serial_matrix(A),'inverse_matrix':serial_matrix(B),
              'action6':serial_matrix(M),'families':{}}
        for k in (1,2):
            cs=coeff(k,lam)
            cp=[sum((F.convert(M[i][j],K)*cs[j] for j in range(6)),F.zero) for i in range(6)]
            lp=-cp[1]/(2*cp[2]) if k==1 else cp[1]/cp[5]
            scalar=cp[2]*lp**2 if k==1 else cp[5]
            assert all(not (u-scalar*v) for u,v in zip(cp,coeff(k,lp)))
            if name in ('F','Q'):
                (a,b),(c,d)=mobius[k][name]
                a,b,c,d=[F.convert(t,K) for t in (a,b,c,d)]
                expected=(a*lam+b)/(c*lam+d)
                assert not (lp-expected)
                if name=='F':
                    expected_scalar=F.convert(gammas[k],K)*(a*lam+b)**3*(c*lam+d)**2/lam**3
                else:expected_scalar=F.convert(z5**(4 if k==1 else 2),K)
                assert not (scalar-expected_scalar)
            else:
                assert not (lp-lam)
                assert not (scalar-(-1 if name=='H' else 1))
            data['families'][str(k)]={'lambda':str(F.to_sympy(lp)),
                'scalar':str(F.to_sympy(scalar)),'all_six_rational_identities':True}
        result[name]=data
    # Conjugation conventions: columns describe images of D,C respectively.
    assert canonical(mul(mul(fourier,diagonal),invfourier))==canonical(invcycle)
    assert canonical(mul(mul(fourier,cycle),invfourier))==canonical(diagonal)
    assert canonical(mul(mul(phase,diagonal),invphase))==canonical(diagonal)
    assert canonical(mul(mul(phase,cycle),invphase))==canonical(mul(mul(diagonal,diagonal),cycle))
    assert canonical(mul(fourier,fourier))==canonical(reflection)
    families={}
    for k in (1,2):
        cs=coeff(k,lam)
        # Multiplication by 5*lambda^3 gives a basis of all degree<=5 polynomials.
        matrix=s.Matrix([[K.to_sympy((5*lam**3*c).numer.get((i,),zero)
                                     /(5*lam**3*c).denom[(0,)]) for i in range(6)] for c in cs])
        assert matrix.det()!=0
        families[str(k)]={'normal_quintic_coefficient_matrix':[[str(v) for v in row] for row in matrix.tolist()],
            'coefficient_determinant':str(matrix.det()),'mobius_group':enumerate_mobius(mobius[k]),
            'mobius_group_order':60,'generic_lambda_stabilizer_order':50,
            'generic_full_parameter_redundancy':3000,
            'Fourier_scalar_gamma':str(K.to_sympy(gammas[k])),
            'injective_parameter_recovery':'-a2/(2*a3)' if k==1 else 'a2/b3'}
    return {'id':'ELLIPTIC_PARAMETER_REDUNDANCY_20261101','all_passed':True,
        'timestamp_utc':datetime.now(timezone.utc).isoformat(),
        'field':'QQ(zeta5)(lambda), Phi5(zeta5)=0',
        'pushforward_convention':'A_* B(x) = A B(A^-1 x) A^T',
        'coefficient_order':['a1','a2','a3','b1','b2','b3'],
        'generator_actions':result,'SL2_F5_elements':enumerate_sl2(),
        'normalizer_order':3000,'projective_template_independence':projection_rank(),
        'families':families,
        'sources':[
            {'url':'https://arxiv.org/pdf/1001.4422','location':'v2 equations16,18,19; equation21 sign caveat'},
            {'url':'https://arxiv.org/html/math/0610403v2','location':'Lemma7.9'},
            {'url':'https://arxiv.org/html/1110.3520v1','location':'Theorem3.6'},
            {'url':'https://arxiv.org/html/2207.07770v2','location':'Proposition2.3.1 and its n=4 proof'},
            {'url':'https://arxiv.org/html/2310.18759v2','location':'Lemma3.2.1(i), paragraph following its proof, Lemma3.2.3'}],
        'scope':['Generator identities and finite group orders are exact arithmetic certificates.',
                 'Exclusion of equivalences outside the normalizer is the separate geometric proof in the accompanying note.',
                 '3000 is the generic degree of the parameterization of either elliptic S61 component.',
                 'If d_k denotes its W16 degree and N_k the five-variable model length, N_k=3000*d_k. Neither d_k is computed here.'],
        'source_sha256':{'code/elliptic_parameter_redundancy.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=certify()
    with args.output.open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps({'all_passed':True,'normalizer_order':3000,'rho_1':3000,'rho_2':3000}))

if __name__=='__main__':main()
