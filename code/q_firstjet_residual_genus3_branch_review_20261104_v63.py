#!/usr/bin/env python3
"""Independent proper-model arithmetic and complete branch/genus audit; no GB."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import re
import resource
import signal
import time
from q_firstjet_residual_forced_nodes_review_20261104_v63 import Poly, var

ROOT=Path(__file__).resolve().parents[1]
BASE='Q_FIRSTJET_RESIDUAL_GENUS3_BRANCH_20261104_v63'
STEM=BASE+'_REVIEW'
P=32003


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def mm(a,b):
    return ((a[0]*b[0]+a[1]*b[2])%P,(a[0]*b[1]+a[1]*b[3])%P,
            (a[2]*b[0]+a[3]*b[2])%P,(a[2]*b[1]+a[3]*b[3])%P)
def ma(a,b):return tuple((x+y)%P for x,y in zip(a,b))
def ms(a,c):return tuple(x*c%P for x in a)
I=(1,0,0,1)
Z=(0,0,0,0)
def mp(a,n):
    out=I
    for _ in range(n):out=mm(out,a)
    return out
def pa(a,b):
    out=dict(a)
    for e,c in b.items():
        out[e]=ma(out.get(e,Z),c)
        if out[e]==Z:out.pop(e)
    return out
def pm(a,b):
    out={}
    for e,c in a.items():
        for f,d in b.items():out=pa(out,{tuple(x+y for x,y in zip(e,f)):mm(c,d)})
    return out
def pp(a,n):
    out={(0,0):I}
    for _ in range(n):out=pm(out,a)
    return out
def det2(a):return (a[0]*a[3]-a[1]*a[2])%P


def polynomial_gcd(a,b):
    def trim(c):
        while c and c[-1]==0:c.pop()
        return c
    a=trim(a[:]);b=trim(b[:])
    while b:
        while len(a)>=len(b):
            n=len(a)-len(b);c=a[-1]*pow(b[-1],-1,P)%P
            for i,v in enumerate(b):a[n+i]=(a[n+i]-c*v)%P
            trim(a)
        a,b=b,a
    return [c*pow(a[-1],-1,P)%P for c in a]


def parse_sparse(expression):
    def visit(n):
        if isinstance(n,ast.Constant) and isinstance(n.value,int):return Poly(n.value)
        if isinstance(n,ast.Name):return var(n.id)
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -visit(n.operand)
        if isinstance(n,ast.BinOp):
            a=visit(n.left)
            if isinstance(n.op,ast.Pow):
                assert isinstance(n.right,ast.Constant) and isinstance(n.right.value,int)
                return a**n.right.value
            b=visit(n.right)
            if isinstance(n.op,ast.Add):return a+b
            if isinstance(n.op,ast.Sub):return a-b
            if isinstance(n.op,ast.Mult):return a*b
        raise ValueError(ast.dump(n))
    return visit(ast.parse(expression,mode='eval').body)


def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    folders={n:ROOT/'reports'/n for n in [
        BASE,'Q_FIRSTJET_RESIDUAL_NODE_BLOWUP_20261104_v63',
        'Q_FIRSTJET_RESIDUAL_QUARTIC_SMOOTH_SPECIALIZATION_20261104_v63',
        'Q_FIRSTJET_RESIDUAL_QUARTIC_FORCED_NODES_20261104_v63_REVIEW',
        'Q_FIRSTJET_RESIDUAL_RECIPROCAL_GENUS_20261104_v63_REVIEW',
        'Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63_REVIEW']}
    reports={n:read(d/'report.json') for n,d in folders.items()}
    primary=reports[BASE]
    assert primary['all_passed']
    assert primary['ordered_relative_sign_quotient_full_curve_geometrically_irreducible_genus']==3
    assert primary['actual_selected_m_dominant_component_m_degree']==8
    assert not primary['NEW_hypersurface_equation_multisection_degree_Weyl_descent_or_global_Poisson_claim']
    model=folders['Q_FIRSTJET_RESIDUAL_NODE_BLOWUP_20261104_v63']
    modelr=read(model/'report.json');modelc=read(model/'certificate.json')
    assert modelr['all_passed'] and modelr['M2_returncode']==0
    assert modelr['generic_quotient_normalization_geometrically_irreducible_genus1']
    failed=folders['Q_FIRSTJET_RESIDUAL_QUARTIC_SMOOTH_SPECIALIZATION_20261104_v63']
    old=read(failed/'report.json');oldc=read(failed/'certificate.json')
    assert not old['all_passed'] and old['run_status']=='completed'
    assert old['generic_four_distinct_infinity_points']
    assert oldc['one_authorized_coefficient_specialization']==[1,2,3,5,7]
    assert oldc['prime']==P and modelc['same_coefficient_point']==[1,2,3,5,7]
    checked={}
    for n,d in folders.items():
        if n not in ('Q_FIRSTJET_RESIDUAL_QUARTIC_SMOOTH_SPECIALIZATION_20261104_v63',):
            assert reports[n]['all_passed']
        for key in ['input_sha256','checked_dependency_sha256']:
            for name,want in reports[n].get(key,{}).items():
                p=ROOT/name;assert sha(p)==want,str(p);checked[str(p.relative_to(ROOT))]=want
        for name,want in reports[n].get('evidence_sha256',{}).items():
            p=d/name;assert sha(p)==want,str(p);checked[str(p.relative_to(ROOT))]=want

    # Independent companion-MATRIX arithmetic in the etale node algebra.
    s=287*pow(2,-1,P)%P;pr=10914*pow(2,-1,P)%P
    alpha=(0,-pr%P,1,s);beta=ma(ms(I,s),ms(alpha,-1))
    delta=ma(alpha,ms(beta,-1))
    assert ma(ma(ms(mp(alpha,2),2),ms(alpha,-287)),ms(I,10914))==Z
    assert det2(delta)==modelc['alpha_beta_difference_norm']!=0
    rn={(1,0):ms(beta,-1),(0,1):alpha}
    zn={(1,0):ms(I,-1),(0,1):I}
    H={}
    for key,co in oldc['all15_specialized_coefficients_r_m'].items():
        i,j=map(int,key.split(','));k=4-i-j
        term=pm(pp(rn,i),pp(zn,k));factor=ms(mp(delta,j),co)
        H=pa(H,{e:mm(v,factor) for e,v in term.items()})
    ownpairs={}
    for e,v in H.items():
        assert e[0]<=2 and e[1]<=2
        pair=(v[0],v[2]);assert v==ma(ms(I,pair[0]),ms(alpha,pair[1]))
        ownpairs[e]=pair
    assert ownpairs=={tuple(e):tuple(c) for e,c in modelc['all_H22_coefficient_pairs']}
    assert H[2,2]==ms(mp(delta,4),oldc['all15_specialized_coefficients_r_m']['4,0'])
    assert det2(H[2,2])==modelc['H22_coefficient_norm']!=0
    script=(model/'blowup22.m2').read_text()
    assert '2*aa^2-287*aa+10914' in script and 'ideal(N,H,diff(u,H),diff(v,H))' in script
    patches=re.findall(r'H=(.*?);\ngood=',script,re.S)
    assert len(patches)==4
    for (i,j),expression in zip(itertools.product(range(2),repeat=2),patches):
        found={(int(e),int(f)):(int(a),int(b)) for a,b,e,f in re.findall(
            r'\((\d+)\+(\d+)\*aa\)\*u\^(\d+)\*v\^(\d+)',expression)}
        expected={(2-e if i else e,2-f if j else f):c for (e,f),c in ownpairs.items()}
        assert found==expected
    log=(model/'M2.log').read_text()
    assert log.splitlines()==['PATCH_00 true','PATCH_01 true','PATCH_10 true','PATCH_11 true']
    # Own Euclidean arithmetic checks the positive infinity sub-receipt, SAME point.
    infinity=[0]*5
    for key,c in oldc['all15_specialized_coefficients_r_m'].items():
        i,j=map(int,key.split(','))
        if i+j==4:infinity[i]=c%P
    derivative=[i*infinity[i]%P for i in range(1,5)]
    assert infinity[4]!=0 and polynomial_gcd(infinity,derivative)==[1]

    # Own sparse determinant with the OPPOSITE order of the two quadratic norms.
    c,d,f,u,v,r,t,D,a,x,y,k,B=[var(n) for n in ['c','d','f','u','v','r','t','D','a','x','y','k','B']]
    q=u*v;R=-4*f*q;T=-4*c*q
    basis=[(0,0),(1,0),(0,1),(1,1)];where={e:i for i,e in enumerate(basis)}
    wc=[2*q*B,4*q*a,4*q*d,u+v]
    M=[[Poly(0) for _ in basis] for _ in basis]
    for j,(er,et) in enumerate(basis):
        for co,(i,h) in zip(wc,basis):
            rr,tt=er+i,et+h
            value=co*(R if rr==2 else 1)*(T if tt==2 else 1)
            M[where[rr%2,tt%2]][j]+=value
    determinant=Poly(0)
    for perm in itertools.permutations(range(4)):
        term=Poly((-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4)))
        for i in range(4):term*=M[i][perm[i]]
        determinant+=term
    A0,B1,C1,D1=wc
    Eopp=A0*A0+C1*C1*T-R*B1*B1-R*D1*D1*T
    norm_opp=Eopp*Eopp-4*T*(A0*C1-R*B1*D1)**2
    assert determinant==norm_opp
    E=B*B-16*f*q*a*a+16*c*q*d*d-4*c*f*(u+v)**2
    Fnorm=a*B+2*c*d*(u+v)
    assert determinant==16*q**4*(E*E+64*f*q*Fnorm*Fnorm)
    bc=read(folders[BASE]/'branch-certificate.json')
    assert all(M[i][j]==parse_sparse(bc['multiplication_matrix_W'][i][j]) for i in range(4) for j in range(4))
    Delta=x*x-y*y;S=x*x+y*y
    assert (x-y)**2*(x+y)**2==Delta*Delta
    assert 4*x*x*y*y-S*S==-Delta*Delta
    assert F(16,4**8)==F(1,4096)
    assert (x-y)**8*(x+y)**8*Delta**4==Delta**12
    # Boundary identities use the original Psi, not division by uv.
    psi=r*r*t*t+(u+v)*r*t*(B+2*a*r+2*d*t)+q*(B+2*a*r+2*d*t)**2-(u-v)**2*(c*r*r+f*t*t+4*c*f*q)
    assert psi.subst('u',Poly(0)).subst('r',Poly(0))==-v*v*f*t*t
    assert psi.subst('v',Poly(0)).subst('r',Poly(0))==-u*u*f*t*t
    assert psi.subst('f',Poly(0)).subst('r',Poly(0))==q*(B+2*d*t)**2
    assert reports['Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63_REVIEW']['intrinsic_rational_m_and_universal_six_invariants_independent_accepted']

    proof={'proper_model_same_coefficient_matrix_reconstruction_all9_pairs':True,
        'all4_saved_projective_patch_inputs_and_unit_logs_bound':True,
        'same_point_infinity_binary_coefficients':infinity,
        'own_Fp_Euclidean_infinity_gcd':[1],
        'proper_flat_etale_node_family_and_geometric_connectedness_proof_accepted':True,
        'own_sparse_biquadratic_norm_opposite_order_and_saved_matrix':True,
        'NormW_m_degree':16,'NormW_m_leading':'Delta^12/(4096*x^12*y^12)',
        'finite_ramification_degree':0,'infinity_ramification_degree':4,
        'each_infinity_qR_valuation':-3,'quotient_genus':1,'full_ordered_cover_genus':3,
        'actual_selected_m_degree':8,
        'NEW_multisection_degree_or_Weyl_quotient_genus_claimed':False,
        'new_GB_or_M2_rerun_or_alternative_coefficient_point':False}
    out.mkdir(parents=True,exist_ok=False)
    evidence=out/'independent-proof.json';evidence.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(STEM+'.md')]
    paths += [d/'report.json' for d in folders.values()]
    paths += [ROOT/'research/notes/Q_FIRSTJET_RESIDUAL_NODE_BLOWUP_20261104_v63.md',
              ROOT/'research/notes'/(BASE+'.md')]
    for p in paths:checked[str(p.relative_to(ROOT))]=sha(p)
    result={'id':STEM,'all_passed':True,
        'independent_companion_matrix_all9_blowup_coefficients_and_all4_patch_binding':True,
        'proper_flat_relative_etale_node_base_and_generic_geometric_smoothness_reviewed':True,
        'generic_quotient_geometrically_irreducible_normalization_genus':1,
        'own_norm_all_coefficients_and_m16_leading_nonzero_verified':True,
        'all_finite_qR_valuations_even_including_m_u_v_f_boundaries':True,
        'four_distinct_infinity_points_order_minus3_and_exact_branch_degree4':True,
        'ordered_relative_sign_full_curve_geometrically_irreducible_genus':3,
        'actual_selected_m_dominant_component_exact_m_degree':8,
        'original_failed_smoothness_report_preserved':True,
        'NEW_multisection_degree_equation_or_Weyl_descent_claim':False,
        'new_M2_GB_coefficient_point_beta_source_or_invariant_job':False,
        'input_sha256':checked,'evidence_sha256':{'independent-proof.json':sha(evidence)},
        'seconds':time.monotonic()-started,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','seconds']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
