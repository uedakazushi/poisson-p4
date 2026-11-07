#!/usr/bin/env python3
"""Independent finite review of all Ga9 quartic sources.

No producer imports. Symbolic coefficient operations use SymPy; full ambient
matrices use the prior independent Fraction implementation. No jets or
component membership are computed here.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

import sympy as s
import classification_translated_minimal_quartic_review as exact
from classification_quadratic_cofactor_review import linear_system, complete_solution

ROOT=Path(__file__).resolve().parents[1]
PREFIX='reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL211_'
SUFFIX='_20261107.json'
n=s.symbols('n0:28')
h,x,y,z,w,u,t=s.symbols('h x y z w u t')
oldcoords=(x,y,z,w)
basecoords=(u,y,w)
allnames={str(a):a for a in (*n,h,x,y,z,w,u,t)}
CHECKS=[]
CERTIFICATES=[]


def read(tag,key):
    return json.loads((ROOT/(PREFIX+tag+SUFFIX)).read_text())[key]


def parse(value):
    # In particular E is a coefficient symbol, never Euler's number.
    names={name:s.Symbol(name)for name in ['E','H','A','D','G','C','L','N','W','r','b','alpha','p','q','k','D1','D2','D3','D4','d0',*(f'm{i}'for i in range(6))]}
    return s.sympify(value,locals={**allnames,**names})


def zero(value,label):
    value=s.cancel(value)
    assert value==0,(label,str(value))
    CHECKS.append(label)


def eqs(polynomial,coords=basecoords,minimum=1):
    return [v for e,v in s.Poly(s.expand(polynomial),*coords).terms()if sum(e)>=minimum]


def subst(value,bindings):
    return s.expand(value.subs(bindings,simultaneous=True))


def bindings(saved):
    return {parse(k):parse(v)for k,v in saved.items()}


def matrix(entries,dimension=4):
    P=s.zeros(dimension)
    for (i,j),v in zip(it.combinations(range(dimension),2),entries):
        P[i,j],P[j,i]=v,-v
    return P


def bracket(P,f,g,coords):
    return s.expand(sum(s.diff(f,coords[i])*P[i,j]*s.diff(g,coords[j])for i,j in it.product(range(len(coords)),repeat=2)))


def jacobi(P,coords):
    return [s.expand(sum(P[a,l]*s.diff(P[b,c],coords[l])for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]for l in range(len(coords))))for i,j,k in it.combinations(range(len(coords)),3)]


def affine_change(P,G):
    return (G*P*G.T).subs(dict(zip(oldcoords,G.inv()*s.Matrix(oldcoords))),simultaneous=True).applyfunc(s.expand)


def force(poly,targets,label,coords=basecoords,minimum=1,allowed=()):
    """Prove necessity without localizing unnamed coefficients.

    Select original coefficients linear in the target variables which vanish
    after the proposed substitution. A full square minor whose determinant
    is a unit in the explicitly allowed localization proves all relations.
    Extra equations are retained for the subsequent reduction.
    """
    variables=list(targets)
    selected=[]
    for c in eqs(poly,coords,minimum):
        if s.Poly(c,*variables).total_degree()<=1 and s.cancel(c.subs(targets,simultaneous=True))==0:
            if any(c.has(v)for v in variables):selected.append(c)
    A,b=s.linear_eq_to_matrix(selected,variables)
    assert A.cols==len(variables)
    columns=A.T.rref()[1]
    assert len(columns)==len(variables),(label,'rank',len(columns),len(variables))
    square=A[list(columns),:]
    determinant=s.factor(square.det())
    numerator,denominator=s.fraction(determinant)
    for element in [numerator,denominator]:
        factors=s.factor_list(element)[1]
        assert all(any(s.cancel(f-v)==0 or s.cancel(f+v)==0 for v in allowed)for f,power in factors),(label,'unannounced localization',determinant)
    zero((A*s.Matrix([targets[v]for v in variables])-b).norm()**2,label+' all selected coefficients')
    CERTIFICATES.append(dict(label=label,variables=list(map(str,variables)),coefficient_rows=list(map(str,selected)),
        square_row_indices=list(columns),determinant=str(determinant),inverted=list(map(str,allowed))))
    return targets


def geometry(row):
    P=matrix([parse(v)for v in row['brackets']])
    # Independently recover every one of the 90 coefficients, then the full
    # 110-row divergence/Lie linear problem and its 28-dimensional solution.
    affine_monomials=[s.Integer(1),x,y,z,w]
    mons=[affine_monomials[i]*affine_monomials[j]for i,j in it.combinations_with_replacement(range(5),2)]
    coefficients=[s.Poly(P[i,j],*oldcoords).coeff_monomial(m)for i,j in it.combinations(range(4),2)for m in mons]
    for a,b in zip(coefficients,row['all_ninety_coefficients']):zero(a-parse(b),'original90')
    A=[[0]*4 for _ in range(4)];A[0][1]=1
    linear,rhs=linear_system(A,[0,0,1,0])
    free,pivots,solution=complete_solution(linear,rhs)
    assert len(pivots)==62 and len(free)==28
    for r,b in zip(linear,rhs):zero(sum(v*coefficients[j]for j,v in r.items())-b,'full110 source equations')
    assert s.Matrix(coefficients).jacobian(n).rank()==28
    c=n[18]+n[19]*y+n[20]*w-n[25]*y*y/3
    phi=n[0]*w+n[1]*w*w/2+n[2]*w**3/3-n[11]*y*w+(n[23]-n[13])*y*y*w/2-n[14]*y*w*w+n[26]*y**3*w/3+n[27]*y*y*w*w/2-n[9]*y-(n[10]-n[21])*y*y/2-(n[12]-n[22])*y**3/3+n[24]*y**4/4
    Phi=u*c+phi
    magnetic=n[3]+n[4]*y+n[5]*w+n[6]*y*y+n[7]*y*w+n[8]*w*w
    V=s.Matrix([(n[16]-n[23]-1-n[26]*y/3-n[27]*w)*u-magnetic,-n[15]-n[16]*y-n[17]*w-n[26]*y*y/3-n[27]*y*w,n[21]+n[22]*y+n[23]*w+n[24]*y*y+n[26]*y*w+n[27]*w*w-n[25]*u])
    Q=s.expand(Phi+sum(V[i]*s.diff(Phi,a)for i,a in enumerate(basecoords)))
    T=(z,x-y*z,y,w)
    changed=matrix([subst(bracket(P,a,b,oldcoords),{x:u+y*t,z:t})for a,b in it.combinations(T,2)])
    expected=matrix([*V,s.diff(Phi,w),-s.diff(Phi,y),s.diff(Phi,u)])
    for a,b in zip(changed,expected):zero(a-b,'full28 straightening')
    assert s.Matrix(T).jacobian(oldcoords).det()==1
    assert s.Matrix((x-y*z,y,w,z)).jacobian(oldcoords).det()==-1
    zero(sum(s.diff(V[i],a)for i,a in enumerate(basecoords))+1,'div V')
    for a,b in zip(jacobi(changed,(t,u,y,w)),[s.diff(Q,w),-s.diff(Q,y),s.diff(Q,u),0]):zero(a-b,'all Jacobi derivatives Q')
    pf=P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2]
    zero(subst(pf,{x:u+y*t,z:t})-(Q-Phi),'original actual Pfaffian')
    zero(s.Poly(Q,u,y,w).coeff_monomial(u*u)+n[20]*n[25],'first exhaustive c equation')
    return P,Phi,V,Q


def normalization_review():
    # General affine centralizer. In the displayed formula beta is the old
    # y value sent to0; shear and translation are in the old y coordinate.
    lam,mu,beta,nu,a,Y,W=s.symbols('lambda mu beta nu a Y W')
    c=n[18]+n[19]*y+n[20]*w-n[25]*y*y/3
    transported=lam*mu*c.subs({y:Y/lam+beta,w:(W-a-nu*(Y/lam+beta))/mu},simultaneous=True)
    zero(transported.subs({n[25]:0,lam:1/n[20],mu:1,beta:0,nu:n[19]/n[20],a:n[18]/n[20]},simultaneous=True)-W,'c=w explicit full affine normalization')
    zero(transported.subs({n[25]:0,n[20]:0,lam:1,mu:1/n[19],beta:-n[18]/n[19],nu:0,a:0},simultaneous=True)-Y,'c=y explicit affine normalization')
    zero(transported.subs({n[25]:0,n[20]:0,n[19]:0,lam:1,mu:1/n[18],beta:0,nu:0,a:0},simultaneous=True)-1,'constant c explicit affine normalization')
    lead,middle,constant=s.symbols('lead middle constant')
    quadratic=lam*mu*(lead*(Y/lam+beta)**2+middle*(Y/lam+beta)+constant)
    centered=s.cancel(quadratic.subs({beta:-middle/(2*lead),mu:lam/lead},simultaneous=True))
    zero(centered-(Y*Y+lam*lam*(constant-middle*middle/(4*lead))/lead),'quadratic c centered normalized formula')
    # Entire c=0 top vector field, transformed by actual affine coordinates.
    Vy=-n[26]*y*y/3-n[27]*w*y
    Vw=n[24]*y*y+n[26]*w*y+n[27]*w*w
    changed=[(lam*Vy).subs({y:Y/lam,w:(W-nu*Y/lam)/mu},simultaneous=True),
             (mu*Vw+nu*Vy).subs({y:Y/lam,w:(W-nu*Y/lam)/mu},simultaneous=True)]
    for value,expected in zip(changed,[-Y*W,W*W+(n[27]*n[24]-2*n[26]**2/9)*Y*Y]):
        zero(value.subs({lam:1,mu:n[27],nu:n[26]/3},simultaneous=True)-expected,'n27 all-top normalization')
    for value,expected in zip(changed,[-Y*Y,3*Y*W]):
        zero(value.subs({n[27]:0,lam:n[26]/3,mu:1,nu:3*n[24]/(4*n[26])},simultaneous=True)-expected,'n26 all-top normalization')
    return dict(affine_centralizer='(x_new,y_new,z_new,w_new)=(lambda*(x-beta*z),lambda*(y-beta),z,mu*w+nu*y+a)',
        determinant='lambda^2*mu',potential_c_transport='lambda*mu*c(y_new/lambda+beta,(w_new-a-nu*(y_new/lambda+beta))/mu)',
        quadratic_center='beta=-middle/(2*lead), mu=lambda/lead, gamma=constant-middle^2/(4*lead)',
        distinct_root_choice='lambda^2=-lead/gamma, gamma!=0',double_root_choice='lambda=1, gamma=0',
        distinct_undo_dimension=5,distinct_undo_coordinates=['three arbitrary quadratic c coefficients','n23 shear','n26 shear'],
        constant_undo_dimension=1,zero_c_undo_dimension=0,all_c_zero_and_discriminant_complements_retained=True)


def reductions(P,Phi,V,Q):
    summaries={}
    cbase={n[18]:0,n[19]:0,n[20]:0,n[25]:0}
    # Three entirely empty normalized c cases; no scalar of V is inverted.
    a={**cbase,n[20]:1}
    target={n[i]:s.Integer(0)for i in [26,21,22,24,16]}
    force(s.diff(subst(Q,a),u),target,'c=w entire u equation',minimum=0)
    a.update(target)
    D=subst(Phi,a).subs({u:0,w:0})
    zero(subst(Q,a).subs(w,0)-(D-n[15]*s.diff(D,y)),'c=w restriction ODE')
    a.update({n[9]:0,n[10]:0,n[12]:0})
    zero(subst(Phi,a).subs(w,0),'c=w Pf factor after degree induction')
    a={**cbase,n[19]:1};target={n[i]:s.Integer(0)for i in [15,23,17,26,27]}
    force(s.diff(subst(Q,a),u),target,'c=y entire u equation',minimum=0)
    a.update(target);B=subst(Phi,a).subs({u:0,y:0})
    zero(subst(Q,a).subs(y,0)-(B+n[21]*s.diff(B,w)),'c=y restriction ODE')
    a.update({n[0]:0,n[1]:0,n[2]:0});zero(subst(Phi,a).subs(y,0),'c=y Pf factor')
    a={**cbase,n[25]:-3}
    target={n[2]:0,n[0]:0,n[1]:0,n[13]:(-2*n[16]+n[23])/3,n[14]:-n[17]/3,n[11]:-2*n[15]/3}
    force(s.diff(subst(Q,a),u),target,'c=y2 entire u equation',minimum=0)
    a.update(target);restriction=subst(Q,a).subs(y,0)
    zero(restriction+(n[15]+n[17]*w)*(n[17]*w*w/3+2*n[15]*w/3-n[9]),'double root squares')
    a.update({n[17]:0,n[15]:0,n[14]:0,n[11]:0});zero(subst(Phi,a).subs(y,0),'double root Pf factor')
    summaries['empty_c']=['w','y','y2']

    # All four distinct-root flags, and every old coefficient under undo.
    a={**cbase,n[25]:-3,n[18]:-1}
    target={n[2]:0,n[13]:(-2*n[16]+n[23])/3,n[14]:-n[17]/3,n[1]:-n[27]/3,n[0]:(n[16]-n[23])/3,n[11]:n[26]/9-2*n[15]/3}
    force(s.diff(subst(Q,a),u),target,'distinct c entire u equation',minimum=0)
    a.update(target);R=subst(Q,a)
    for rho in [-1,1]:
        zero(s.Poly(R,w).coeff_monomial(w**3).subs(y,rho)+(n[17]+rho*n[27])**2/3,'root cubic square')
    R=R.subs({n[17]:0,n[27]:0})
    for rho in [-1,1]:
        beta=n[15]+rho*n[16]+n[26]/3
        zero(s.Poly(R,w).coeff_monomial(w).subs(y,rho)-2*beta*(rho-beta)/3,'root linear flag')
    saved=read('SHORT1_DISTINCT_C','nil_short_translation_distinct_c')
    for item in saved['all_four_root_flags']:
        positive,negative=item['flags'];a0=s.Rational(positive+negative,2);b0=s.Rational(positive-negative,2)
        flag={n[16]:a0,n[15]:b0-n[26]/3,n[17]:0,n[27]:0}
        original=P.applyfunc(lambda v:subst(subst(v,a),flag))
        G=s.eye(4);G[0,3]=(a0+n[23])/3;G[2,3]=-n[26]/3
        assert G.det()==1
        changed=affine_change(original,G)
        new=bindings(item['transformed_parameters'])
        for old,newentry in zip(changed,P.applyfunc(lambda v:subst(v,new))):zero(old-newentry,'flag full90 affine undo')
        canonical=bindings(item['canonical_parameters']);ph=subst(Phi,canonical);q=subst(Q,canonical)
        D=s.expand(ph-(y*y-1)*u-((positive-negative)*y+positive+negative)*w/3)
        assert not D.has(u,w)
        for rho,active in [(1,positive),(-1,negative)]:
            if not active:zero(q.subs(y,rho)-D.subs(y,rho),'inactive root all-parameter factor')
    canonical=bindings(saved['all_four_root_flags'][-1]['canonical_parameters'])
    unknown=[n[j]for j in [3,4,5,6,7,8,9,10,12,21,22,24]]
    M,rhs=s.linear_eq_to_matrix(eqs(subst(Q,canonical)),unknown)
    assert M.rank()==7 and M.cols==12
    full=bindings(saved['canonical_five_parameter_source'])
    zero((M*s.Matrix([full[v]for v in unknown])-rhs).norm()**2,'whole canonical five source complete linear system')
    assert s.Matrix([full[v]for v in unknown]).jacobian(s.symbols('D1 D2 D3 D4 d0')).rank()==5
    summaries['distinct_flags']=['00 empty','01 empty','10 empty','11 complete source5']

    # c=1: u equations, all quartic and cubic necessities, then A2 exactly.
    a={**cbase,n[18]:1}
    target={n[16]:n[23],n[26]:0,n[27]:0}
    force(s.diff(subst(Q,a),u),target,'constant c full u',minimum=0);a.update(target)
    target={n[2]:0,n[13]:s.Rational(1,2)-n[23],n[17]:-2*n[14]}
    force(subst(Q,a),target,'constant c quartic',minimum=4,allowed=(n[24],));a.update(target)
    zero(s.Poly(subst(Q,a),u,y,w).coeff_monomial(w**3)+2*n[14]**2,'constant c cubic square')
    a.update({n[14]:0,n[17]:0})
    target={n[1]:(n[23]-1)*(n[23]-s.Rational(1,4))/n[24],
            n[11]:-n[15]+((1-3*n[23])*(n[22]-n[12])/3+n[22]*(n[23]-s.Rational(1,4)))/n[24]}
    force(subst(Q,a),target,'constant c full cubic',minimum=3,allowed=(n[24],))
    constant=read('SHORT1_CONSTANT_C','nil_short_translation_constant_c')
    summaries['constant_source_parameters']=constant['parameters']

    # c=0 nonquadruple leading cases: all coefficients, including n24=0.
    a={**cbase,n[27]:1,n[26]:0}
    target={n[2]:0,n[17]:-n[14],n[13]:2*n[16]-n[23]-1,n[12]:n[24]*n[14]}
    force(subst(Q,a),target,'zero c n27 quartic',minimum=4);a.update(target)
    target={n[1]:n[14]**2,n[15]:-n[14]*n[16]}
    force(subst(Q,a),target,'zero c n27 cubic',minimum=3);a.update(target)
    ph=subst(Phi,a).subs(y,n[14]);q=subst(Q,a).subs(y,n[14])
    assert s.Poly(ph,w).degree()<=1
    zero(s.Poly(q,w).coeff_monomial(w*w)-s.diff(ph,w),'n27 root factor coefficient')
    zero(subst(V[1],a).subs(y,n[14]),'n27 root Vy')
    a={**cbase,n[27]:0,n[26]:3,n[24]:0}
    target={n[2]:0,n[14]:-3*n[17]/5,n[13]:3*n[23]-6*n[16]+2,n[12]:0}
    force(subst(Q,a),target,'zero c n26 quartic',minimum=4);a.update(target)
    zero(s.Poly(subst(Q,a),u,y,w).coeff_monomial(w**3)+3*n[17]**2/5,'n26 cubic square')
    a.update({n[17]:0,n[14]:0})
    zero(s.Poly(subst(Q,a),u,y,w).coeff_monomial(y*w*w)-3*n[1],'n26 residual n1')
    a[n[1]]=0;ph=subst(Phi,a);Ac=s.diff(ph,w);D=ph.subs(w,0);L=y*y+n[16]*y+n[15]
    assert s.Poly(Ac,y).degree()==3 and s.Poly(Ac,y).LC()==1
    zero(s.diff(subst(Q,a),w)+L*s.diff(Ac,y)-(3*y+n[23]+1)*Ac,'n26 full A root equation')
    zero(subst(Q,a).subs(w,0)+L*s.diff(D,y)-(n[22]*y+n[21])*Ac-D,'n26 full D root equation')

    a={**cbase,n[26]:0,n[27]:0}
    target={n[2]:0,n[13]:n[23]-2*n[16]+s.Rational(1,2),n[17]:-2*n[14]}
    force(subst(Q,a),target,'zero c quadruple quartic',minimum=4,allowed=(n[24],));a.update(target)
    zero(s.Poly(subst(Q,a),u,y,w).coeff_monomial(w**3)+2*n[14]**2,'quadruple cubic square')
    a.update({n[14]:0,n[17]:0})
    q=subst(Q,a)
    zero(s.Poly(q,u,y,w).coeff_monomial(w*w)-n[1]*(n[23]+s.Rational(1,2)),'Nnonzero forces d')
    zero(s.Poly(q,u,y,w).coeff_monomial(y*y*w)-n[24]*n[1]-(n[16]-s.Rational(1,4))*(n[23]-2*n[16]+1),'N exact relation')
    # General completed-square equation before the finite odd-degree split.
    b,alpha,r,N,E,H,p,q0,k,W=s.symbols('b alpha r N E H p q k W')
    B=b-s.Rational(1,4);A=B*y*y+E*y+H;D1,D2,D3=s.symbols('D1 D2 D3');D=r*y**4/4+D3*y**3+D2*y*y+D1*y
    R=k+A*A/(2*N)-D
    Vw=r*y*y+q0*y-w/2+p
    T=s.expand(Vw+(-alpha-b*y)*s.diff(A,y)/N+W/2).subs(w,W-A/N)
    T=s.expand(T)
    zero(T-(r*y*y+q0*y+p+A/(2*N)-(alpha+b*y)*s.diff(A,y)/N),'completed W translation')
    zero(s.Poly(T,y).coeff_monomial(y*y).subs(N,2*B*B/r),'completed W highest coefficient')
    assert s.degree(s.cancel(R.subs(N,2*B*B/r)),y)<=3
    formalPhi=N*w*w/2+A*w+D
    formalQ=formalPhi+(-alpha-b*y)*s.diff(formalPhi,y)+Vw*s.diff(formalPhi,w)
    zero(formalQ.subs(w,W-A/N)-(k+(alpha+b*y)*s.diff(R,y)-R+N*W*T),'full completed square Q identity')
    # The two coefficients T determine exactly the recorded p,q, with no
    # division by alpha,b,E,H. Only r,N have already been inverted.
    target={q0:-((s.Rational(1,2)-b)*E-2*B*alpha)/N,p:-(H/2-E*alpha)/N}
    force(T,target,'completed square p q',coords=(y,),minimum=0,allowed=(N,))
    # N=0, constant A=L: complete residual D coefficients linearly.
    L=s.Symbol('L');D=r*y**4/4+D3*y**3+D2*y*y+D1*y
    linearQ=L*w+D+(-alpha-y/4)*s.diff(D,y)+(r*y*y+q0*y-w+p)*L
    target={D3:4*alpha*r,D2:24*alpha*alpha*r-2*r*L,D1:8*alpha*(24*alpha*alpha*r-2*r*L)/3-4*q0*L/3}
    force(linearQ,target,'linear w complete D',coords=(y,w),minimum=1)
    summaries['zero_source_count']=3
    summaries['ODE_m_values']=[1,3]
    return summaries


def source_families(P,Phi,V,Q):
    saved5=read('SHORT1_DISTINCT_C','nil_short_translation_distinct_c')
    saved10=read('SHORT1_CONSTANT_C','nil_short_translation_constant_c')
    saved11=read('SHORT1_ZERO_C','nil_short_translation_zero_c')['all_three_complete_families']
    specs=[dict(name='distinct_quadratic_c',parameters=saved5['free_source_parameters'],complete_source_substitutions=saved5['canonical_five_parameter_source'],constant_k=saved5['constant_k'],actual_pfaffian=saved5['actual_pfaffian']),dict(name='constant_c',constant_k='k',**saved10),*saved11]
    results=[]
    universal=[]
    for item in specs:
        params=[parse(a)for a in item['parameters']]
        mapping=bindings(item['complete_source_substitutions'])
        assert len(mapping)==28
        ph=subst(Phi,mapping);q=subst(Q,mapping);kap=parse(item['constant_k'])
        zero(q-kap,item['name']+' universal eigenpotential')
        family=P.applyfunc(lambda value:s.cancel(value.subs(mapping,simultaneous=True)))
        for value in jacobi(family,oldcoords):zero(value,item['name']+' original full Jacobi')
        actual=family[0,1]*family[2,3]-family[0,2]*family[1,3]+family[0,3]*family[1,2]
        zero(actual-(kap-ph).subs(u,x-y*z),item['name']+' actual Pf')
        zero(actual-parse(item['actual_pfaffian']).subs(u,x-y*z),item['name']+' frozen Pf')
        if 'universal_brackets'in item:
            for (i,j),text in zip(it.combinations(range(4),2),item['universal_brackets']):zero(family[i,j]-parse(text),item['name']+' every frozen source coefficient')
        if len(params)==10:
            zero(s.diff(ph,u)-1,'constant c truly1')
            assert s.Poly(s.cancel(subst(V[0],mapping)+u),y,w).total_degree()<=2
        results.append(dict(name=item['name'],parameter_count=len(params),irreducible_domain=True,
            full28_substitution_checked=True,all_original_Jacobi_zero=True,actual_Pf_checked=True))
        universal.append((item,family,params,mapping))
    return results,universal


def to_fraction(P,scale=1):
    upper={}
    for i,j in it.combinations(range(5),2):
        v={e+(0,0):Fraction(c)*Fraction(scale)for e,c in s.Poly(P[i,j],h,x,y,z,w).terms()if c}
        if v:upper[i,j]=v
    return exact.arithmetic.antisymmetric(upper)


def homogeneous(P):
    result=s.zeros(5)
    for i,j in it.combinations(range(4),2):
        entry=sum(c*h**(2-sum(e))*s.prod(a**b for a,b in zip(oldcoords,e))for e,c in s.Poly(P[i,j],*oldcoords).terms())
        assert all(sum(e)<=2 for e,c in s.Poly(P[i,j],*oldcoords).terms())
        result[i+1,j+1],result[j+1,i+1]=entry,-entry
    return result


def ambient_and_normalforms(universal):
    results=[]
    for item,family,params,mapping in universal[2:]:
        tag='ZERO_C_'+item['name'].upper()
        frozen=read(tag+'_COMPONENT','source_point')
        point={parse(k):parse(v)for k,v in frozen['point']['family_parameter_point'].items()}
        H=homogeneous(family)
        scale=parse(frozen['integer_scale'])
        P=H.subs(point).applyfunc(s.cancel)*scale
        assert s.Matrix([mapping[v]for v in n]).jacobian(params).subs(point).rank()==11
        D,G,rows=exact.matrices(to_fraction(P))
        derivatives=[exact.vector(to_fraction(H.diff(a).subs(point).applyfunc(s.cancel),scale))for a in params]
        F=[row+[col[i]for col in derivatives]for i,row in enumerate(G)]
        assert all(sum(v*column[j]for j,v in row)==0 for row in rows for column in derivatives)
        ranks={};certificates={}
        for name,M,key in [('D',D,'derivative'),('G',G,'orbit_vertical_scale'),('F',F,'actual_full_family')]:
            certificate=frozen[key]
            rank,cols=exact.arithmetic.rank(M)
            assert rank==certificate['rank']
            ranks[name]=rank
            row_indices=exact.arithmetic.rank([[row[j]for row in M]for j in cols])[1]
            determinant=exact.minor(M,row_indices,cols)
            assert determinant
            certificates[name]=dict(rank=rank,rows=row_indices,columns=cols,determinant=str(determinant),
                full_matrix_sha256=hashlib.sha256(json.dumps(M,default=str,separators=(',',':')).encode()).hexdigest())
        assert ranks['G']==ranks['F']==48
        assert exact.arithmetic.rank([row[:25]for row in G])[0]==24
        normal=read(tag+'_NORMALFORM','source_magnetic_normalform')
        T=s.Matrix([[parse(v)for v in row]for row in normal['complete_projective_transform']])
        assert T.det()!=0
        oldpoint=family.subs(point).applyfunc(s.cancel)
        G4=T[1:,1:];shift=T[1:,0]
        transformed=(G4*oldpoint*G4.T).subs(dict(zip(oldcoords,G4.inv()*(s.Matrix(oldcoords)-shift))),simultaneous=True)
        canonical_map=bindings(normal['canonical_source_substitutions'])
        canonical=matrix([subst(parse(v),canonical_map)for v in SOURCE['brackets']])
        for a,b in zip(transformed,canonical):zero(a-parse(normal['tensor_scale'])*b,item['name']+' actual affine canonical conjugacy')
        Hcanonical=homogeneous(canonical)*parse(normal['canonical_integer_scale'])
        for i,j in it.combinations(range(5),2):zero(Hcanonical[i,j]-parse(normal['canonical_brackets'].get(f'{i},{j}','0')),item['name']+' canonical all150')
        _,canonical_G,_=exact.matrices(to_fraction(Hcanonical))
        assert exact.arithmetic.rank(canonical_G)[0]==48
        results.append(dict(name=item['name'],full_ambient_D_rank=ranks['D'],actual_vertical_rank=24,
            actual_GL5_vertical_scale_rank=48,actual_whole_source_rank=48,projective_orbit_dimension=23,
            whole_source_saturation_upper_bound=23,whole_source_in_canonical_orbit_closure=True,
            independent_full_matrix_certificates=certificates,
            all_points_in_single_actual_orbit_claimed=False,ambient_componenthood_claimed=False))
    return results


def main():
    global SOURCE
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    started=time.monotonic()
    source_path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'
    SOURCE=json.loads(source_path.read_text())['full_nilpotent_sources'][8]
    P,Phi,V,Q=geometry(SOURCE)
    normalizations=normalization_review()
    reduction=reductions(P,Phi,V,Q)
    sources,universal=source_families(P,Phi,V,Q)
    orbits=ambient_and_normalforms(universal)
    # Full affine cofactor orbit, directly in all20 affine Lie directions.
    coords=s.Matrix(oldcoords);Z=s.Matrix([y,0,1,0]);columns=[]
    for a,b in it.product(range(4),repeat=2):
        Y=s.zeros(4,1);Y[a]=coords[b]
        value=Z.jacobian(coords)*Y-Y.jacobian(coords)*Z
        columns.append([s.Poly(entry,*oldcoords).coeff_monomial(mon)for entry in value for mon in [1,*oldcoords]])
    for a in range(4):
        Y=s.zeros(4,1);Y[a]=1;value=Z.jacobian(coords)*Y
        columns.append([s.Poly(entry,*oldcoords).coeff_monomial(mon)for entry in value for mon in [1,*oldcoords]])
    assert s.Matrix.hstack(*map(s.Matrix,columns)).rank()==9
    source_bytes=Path(__file__).read_bytes()
    inputs=[source_path,*[ROOT/(PREFIX+tag+SUFFIX)for tag in ['SHORT1_GEOMETRY','SHORT1_EMPTY_C','SHORT1_DISTINCT_C','SHORT1_CONSTANT_C','SHORT1_ZERO_C',*[f'ZERO_C_{kind}_{mode}'for kind in ['QUADRATIC_W_POWER1','QUADRATIC_W_POWER3','LINEAR_W']for mode in ['COMPONENT','NORMALFORM']]]]]
    report=dict(task='E_reduced_nonnormal',all_passed=True,reviewer_source_sha256=hashlib.sha256(source_bytes).hexdigest(),reviewer_source_bytes=source_bytes.decode(),
        input_hashes=[dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest())for p in inputs],
        checked_identities=len(CHECKS),full28_to_all90_linear_rank=62,full110_linear_nullity=28,
        coordinates_ordered=dict(t_u_y_w_determinant=1,u_y_w_t_determinant=-1),finite_reduction=reduction,
        full_affine_normalizations=normalizations,
        exact_linear_necessity_certificates=CERTIFICATES,complete_five_sources=sources,zero_c_orbits=orbits,
        affine_cofactor_orbit_dimension=9,fixed_cofactor_dimension_bounds=dict(distinct_quadratic_c=10,constant_c=11,zero_c=11),
        projective_saturation_dimension_bounds=dict(distinct_quadratic_c=22,constant_c=23,zero_c=23),
        whole_Ga9_projective_dimension_upper_bound=23,G_dimension25_to27_generic_dominance_excluded=True,
        component_membership_claimed=False,elapsed_seconds=round(time.monotonic()-started,6),peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(report,indent=2,default=str)+'\n')
    print(json.dumps({k:report[k]for k in ['all_passed','checked_identities','whole_Ga9_projective_dimension_upper_bound','elapsed_seconds','peak_rss_kib']}))


if __name__=='__main__':main()
