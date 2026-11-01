#!/usr/bin/env python3
"""Finite remaining normal forms for generic-beta reduced R0 infinity.

The proof note supplies the orbit and incidence-dimension arguments. This
certificate verifies three exact radicals, the normalizing transformations,
and the rank-six evaluation map. It does not exhaust all Q45 infinity strata.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as s
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import COORDS,J,minor_witness

ROOT=Path(__file__).resolve().parents[1]


def run():
    source=ROOT/'reports/S61_R0_INFINITY_20261101.json'
    data=json.loads(source.read_text());g=list(map(s.Symbol,data['g_variables']))
    raw=list(map(s.sympify,data['leading_equations_70_through_133']))+[s.sympify(data['leading_G6'])]
    allu=s.symbols('u1:62');cases=[]
    for label,kval,Kval,Lval,rank_one in [
            ('resonant_plus',s.Rational(1,3),1,0,False),
            ('resonant_minus',s.Rational(-4,3),0,1,False),
            ('rank_one_nonzero_scalar',s.S.One,0,0,True)]:
        substitution={allu[3]:0,allu[11]:0,allu[4]:-int(rank_one),allu[12]:0,
            allu[13]:kval,allu[48]:2*kval+int(not rank_one),
            allu[25]:0,allu[34]:0,allu[49]:Kval,allu[50]:-Lval}
        equations=[s.expand(f.subs(substitution)) for f in raw]
        faces=([{1,4,8},{0,1}] if rank_one else [{0,2,5},{1,4,8},{0,1}])
        used=set().union(*faces)
        generators=[g[i] for i in range(9) if i not in used]+[
            g[i]*g[j] for i,j in combinations(sorted(used),2) if not any({i,j}<=f for f in faces)]
        for f in equations:
            for mon,co in s.Poly(f,*g).terms():
                if co:
                    assert any(all(a<=b for a,b in zip(s.Poly(gen,*g).monoms()[0],mon,strict=True))
                               for gen in generators)
        ideal=s.groebner(equations,*g,order='grevlex',domain=s.QQ)
        powers=[];quotients=[]
        for generator in generators:
            for power in range(1,5):
                q,remainder=ideal.reduce(generator**power)
                if remainder==0:
                    assert s.expand(sum(a*b.as_expr() for a,b in zip(q,ideal.polys,strict=True))
                                    -generator**power)==0
                    powers.append(power);quotients.append(list(map(str,q)));break
            else:raise AssertionError(f'Unproved generator {generator}')
        cases.append({'label':label,'beta_substitution':{str(a):str(b) for a,b in substitution.items()},
            'equations':list(map(str,equations)),'radical_generators':list(map(str,generators)),
            'coordinate_faces':list(map(sorted,faces)),'groebner_basis':list(map(lambda p:str(p.as_expr()),ideal.polys)),
            'generator_powers_in_original_ideal':powers,'power_membership_quotients':quotients})

    z,u,beta=primitive_beta(S61());v=s.Matrix([0,1,0,0]);at_v=beta.subs(dict(zip(z,v,strict=True)))
    alpha=s.expand(at_v[0,1]);assert alpha==-3*allu[48]/5
    a,b,c=s.symbols('a b c');vec=s.Matrix([0,a,b,c])
    q=s.expand((vec.T*beta.subs(dict(zip(z,vec,strict=True)))*v)[0])
    assert s.expand(q-a*b*allu[11]+a*c*allu[3]-b*b*allu[12]-b*c*(2*allu[13]-allu[48])+c*c*allu[4])==0
    evaluation=s.Matrix([allu[3],allu[11],allu[12],2*allu[13]-allu[48],-allu[4],alpha]).jacobian(u)
    evaluation_minor=minor_witness(evaluation);assert evaluation_minor['rank']==6
    coefficients=s.Matrix([s.diff(beta[i,j],z[t]) for (i,j),t in COORDS]);basis=coefficients.jacobian(u)
    rows=list(basis.T.rref()[1]);left=basis[rows,:].inv()
    k,K,L,A,B,t=s.symbols('k K L A B t',nonzero=True)
    T=s.Matrix([[1,0,0,0],[0,1,A,-B],[B,0,1,0],[A,0,0,1]])
    assert T.T*J*T==J
    transforms=[]
    for rank_one in [False,True]:
        sub={allu[3]:0,allu[11]:0,allu[4]:-int(rank_one),allu[12]:0,
             allu[13]:k,allu[48]:2*k+int(not rank_one),
             allu[25]:0,allu[34]:0,allu[49]:K,allu[50]:-L}
        def pullback(matrix):
            pulled=(matrix.T*beta.subs(sub).subs(dict(zip(z,matrix*s.Matrix(z))),simultaneous=True)*matrix).applyfunc(s.expand)
            coeff=s.Matrix([s.diff(pulled[i,j],z[h]) for (i,j),h in COORDS])
            return (left*coeff[rows,:]).applyfunc(s.cancel)
        newu=pullback(T.inv())
        Kp=s.factor(5*newu[u.index(allu[34])]+newu[u.index(allu[49])])
        Lp=s.factor(5*newu[u.index(allu[25])]-newu[u.index(allu[50])])
        if rank_one:
            assert s.expand(Kp-(3*A*k+K))==0
            assert s.expand(Lp-(-5*A+3*B*k+L))==0
            removal={A:-K/(3*k),B:(-5*K-3*L*k)/(9*k*k)}
            assert s.cancel(Kp.subs(removal,simultaneous=True))==0
            assert s.cancel(Lp.subs(removal,simultaneous=True))==0
            torus=s.diag(t*t,1/(t*t),1/t,t)
        else:
            assert s.expand(Kp-((3*k-1)*A+K))==0
            assert s.expand(Lp-((3*k+4)*B+L))==0
            removal={A:-K/(3*k-1),B:-L/(3*k+4)}
            torus=s.diag(1,1,1/t,t)
        assert torus.T*J*torus==J
        torus_u=pullback(torus)
        assert torus_u[u.index(allu[13])]==(k/t**2 if rank_one else k)
        assert torus_u[u.index(allu[4])]==-int(rank_one)
        assert s.cancel(5*torus_u[u.index(allu[34])]+torus_u[u.index(allu[49])]-K/t)==0
        assert s.cancel(5*torus_u[u.index(allu[25])]-torus_u[u.index(allu[50])]-L*t)==0
        transforms.append({'rank_one_conic':rank_one,'K_prime':str(Kp),'L_prime':str(Lp),
            'translation_parameters':{str(aa):str(bb) for aa,bb in removal.items()},
            'translated_u':{str(aa):str(bb) for aa,bb in zip(u,newu,strict=True)},
            'torus':[[str(value) for value in row] for row in torus.tolist()],
            'torus_u':{str(aa):str(bb) for aa,bb in zip(u,torus_u,strict=True)}})
    inputs=[Path(__file__),source,ROOT/'reports/S61_R0_UNIFORM_VERTEX_20261101.json',
            ROOT/'code/primitive_jet_elliptic_curve_20261101.py']
    return {'id':'S61_R0_GENERIC_EXHAUSTION_20261101','all_passed':True,'field':'QQ',
        'fixed_normal_form_radicals':cases,'tau_evaluation':['u4','u12'],
        'q_on_v_perp':str(q),'scalar_alpha':str(alpha),
        'six_jet_evaluation_matrix':[[str(value) for value in row] for row in evaluation.tolist()],
        'six_jet_evaluation_minor':evaluation_minor,'normalizing_transformations':transforms,
        'excluded_incidence_dimension_bounds':{'tau_zero_q_zero':14,'tau_zero_q_rank_one_alpha_zero':15},
        'scope':'Together with the open-recognition and uniform-vertex certificates, exhausts the '
                'REDUCED R0 part of UNSATURATED infinity for sufficiently general beta. '
                'Does not exhaust other Q45 strata, determine nonreduced R0, or bound the affine S61 fiber.',
        'proof_note':'research/notes/S61_R0_GENERIC_EXHAUSTION_20261101.md',
        'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();report=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'fixed_radicals':3,'evaluation_rank':6,'output':str(args.output)}))


if __name__=='__main__':main()
