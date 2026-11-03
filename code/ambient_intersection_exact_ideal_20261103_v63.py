#!/usr/bin/env python3
"""Exact finite image ideals for CA/CB in their common critical-beta slice.

The common target fiber selection has determinant -72 identically. After
solving it, each map is a polynomial even+odd expression in p. The image
T+J*p, p²=a has its full presentation by module syzygies and pair products.
The generated M2 file compares the ACTUAL target ideals, not normalizations.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]

def even_odd(f,p,a):
    poly=s.Poly(s.expand(f),p);even=0;odd=0
    for (n,),coefficient in poly.terms():
        if n%2:odd+=coefficient*a**(n//2)
        else:even+=coefficient*a**(n//2)
    return s.expand(even),s.expand(odd)

def a_s_decompose(f,a,sv):
    U=s.cancel(f.subs(sv,0)/a);V=s.cancel((f-a*U)/sv)
    assert s.expand(a*U+sv*V-f)==0
    assert not s.denom(U).free_symbols and not s.denom(V).free_symbols
    return s.expand(U),s.expand(V)

def m2(f):return str(s.expand(f)).replace('**','^')

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    path=ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json';body=json.loads(path.read_text())
    base=json.loads((ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json').read_text());qi=base['q_indices_zero_based'];bi=base['beta_indices_zero_based']
    p,old_s=s.symbols('p s');sv,a,b,c=s.symbols('ss aa bb cc');K=s.Matrix(s.symbols('v0:10'));z=s.Matrix(s.symbols('z0:35'))
    all_q_selected=body['families']['CA']['actual_q_selected_rows'];remaining=[i for i in range(45) if i not in all_q_selected]
    assert all_q_selected==body['families']['CB']['actual_q_selected_rows']
    families={}
    for name in ['CA','CB']:
        record=body['families'][name];pars=list(map(s.Symbol,record['source_fiber_parameters']));u=s.Matrix(list(map(s.sympify,record['full_61_coordinates'])))
        q=u[qi,0];M=q.jacobian(pars);offset=q.subs(dict.fromkeys(pars,0));L=M[all_q_selected,:]
        assert s.factor(L.det())==-72
        fiber=(L.inv()*(K-offset[all_q_selected,0])).applyfunc(s.expand)
        solved=q.subs(dict(zip(pars,list(fiber))),simultaneous=True).subs(old_s,sv).applyfunc(s.expand)
        assert solved[all_q_selected,0]==K
        eo=[even_odd(f,p,a) for f in solved[remaining,0]];E=s.Matrix([v[0] for v in eo]);O=s.Matrix([v[1] for v in eo]);O0=O.subs({a:0,sv:0})
        R=O0.jacobian(K);assert O0==R*K and R.rank()==7
        rows=list(R.T.rref()[1]);j7=O0[rows,0];Jmat=R[rows,:];cols=list(Jmat.rref()[1]);C=R[:,cols]*Jmat[:,cols].inv();assert C*j7==O0
        UV=[a_s_decompose(f,a,sv) for f in O-C*j7];U=s.Matrix([v[0] for v in UV]);V=s.Matrix([v[1] for v in UV])
        W=z[rows,0]-E[rows,0]-b*U[rows,0]-c*V[rows,0];odd=s.Matrix([b,c,*W]);j=s.Matrix([a,sv,*j7])
        linear=(z-E-b*U-c*V-C*W).applyfunc(s.expand)
        equations=[v for v in linear if v!=0]
        equations += [s.expand(j[i]*odd[k]-j[k]*odd[i]) for i in range(9) for k in range(i)]
        equations += [s.expand(odd[i]*odd[k]-a*j[i]*j[k]) for i in range(9) for k in range(i+1)]
        subst={a:p*p,b:p**3,c:sv*p}|dict(zip(z,list(solved[remaining,0])))
        assert all(s.expand(f.subs(subst,simultaneous=True))==0 for f in equations)
        families[name]=dict(equations=equations,source_fiber_solution=fiber,solved=solved,E=E,O=O,j7=j7,C=C,W=W,rows=rows,cols=cols,Jmat=Jmat,
            solved_61=u.subs(dict(zip(pars,list(fiber))),simultaneous=True).subs(old_s,sv))
        print(name,'EXACT IMAGE PRESENTATION',len(equations),'seconds',time.monotonic()-started,flush=True)
    # Common49 ambient target coordinates: s,a,b,c,K10,q_remaining35.
    target_u=s.zeros(61,1)
    for index,kk in zip(all_q_selected,K):target_u[qi[index]]=kk
    for index,zz in zip(remaining,z):target_u[qi[index]]=zz
    beta=families['CA']['solved_61'][bi,0]
    for index,f in zip(bi,beta):
        ev,od=even_odd(f,p,a);U,V=a_s_decompose(od,a,sv);target_u[index]=s.expand(ev+b*U+c*V)
    F=s.expand(target_u[0]-4*target_u[34]+target_u[51])
    fa=families['CA'];Fe,Fo=even_odd(fa['solved_61'][0]-4*fa['solved_61'][34]+fa['solved_61'][51],p,a)
    E0=Fe.subs({a:0,sv:0});d=E0.jacobian(K) if isinstance(E0,s.MatrixBase) else s.Matrix([E0]).jacobian(K)
    coeff=d[:,fa['cols']]*fa['Jmat'][:,fa['cols']].inv();res=s.expand(Fe-(coeff*fa['j7'])[0]);U,V=a_s_decompose(res,a,sv)
    G=s.expand(b*U+c*V+(coeff*fa['W'])[0]+a*Fo)
    substA={a:p*p,b:p**3,c:sv*p}|dict(zip(z,list(fa['solved'][remaining,0])))
    assert s.expand((G-p*F).subs(substA,simultaneous=True))==0
    # The common nonnormal stratum has only k0,k2,k5 free at p=s=0.
    nr=body['families']['N'];npars=list(map(s.Symbol,nr['source_fiber_parameters']));nu=s.Matrix(list(map(s.sympify,nr['full_61_coordinates'])))
    nq=nu[qi,0];NM=nq.subs({p:0,old_s:0}).jacobian(npars);Nsel=nr['actual_q_selected_rows'];Nleft=NM[Nsel,:].inv()
    target_nk=s.Matrix([K[all_q_selected.index(i)] for i in Nsel]);residue_fiber=Nleft*target_nk
    nfree=s.zeros(9,1)
    for i in [0,2,5]:nfree[i]=residue_fiber[i]
    nonnormal=nu.subs({p:0,old_s:0}|dict(zip(npars,list(nfree))),simultaneous=True)
    annihilator=[s.expand(v) for v in target_u-nonnormal if s.expand(v)!=0]
    # F has coefficient1 at q_remaining corresponding to raw61 index51.
    star=z[remaining.index(qi.index(51))];assert s.diff(F,star)==1
    zero_F={star:s.expand(star-F)}
    IN=[*fa['equations'],F,G]
    section_generators=[s.expand(f.subs(zero_F,simultaneous=True)) for f in IN]
    split=[f for f in section_generators if f!=0]+[F*F]+[s.expand(F*h) for h in annihilator]
    names=list(map(str,[*z,b,c,*K,sv,a]));text=['R=QQ['+','.join(names)+',MonomialOrder=>GRevLex];']
    def ideal(name,equations):return name+'=ideal('+','.join(m2(f) for f in equations if f!=0)+');'
    text += [ideal('IA',fa['equations']),ideal('IB',families['CB']['equations']),ideal('ISPLIT',split),
        'ID=IA+IB;','t0=cpuTime();','GD=gb ID;','print("GB_D_SECONDS "|toString(cpuTime()-t0));',
        'FF='+m2(F)+';','GG='+m2(G)+';',
        'print("F_REMAINDER "|toString(FF % GD));','print("G_REMAINDER "|toString(GG % GD));',
        'print("F_SQUARE_REMAINDER "|toString(FF^2 % GD));',
        'print("SPLIT_IN_D "|toString(all(flatten entries gens ISPLIT,f->f % GD==0)));',
        't1=cpuTime();','GS=gb ISPLIT;','print("GB_SPLIT_SECONDS "|toString(cpuTime()-t1));',
        'print("D_IN_SPLIT "|toString(all(flatten entries gens ID,f->f % GS==0)));',
        'print("DIM_D "|toString dim ID);','print("EXACT_INTERSECTION_CHECK_FINISHED");']
    (output/'intersection.m2').write_text('\n'.join(text)+'\n')
    paths=[Path(__file__),path,ROOT/'reports/AMBIENT_WHOLE_REDUCED_INTERSECTION_20261103_v63/report.json',
        ROOT/'reports/AMBIENT_INTERSECTION_TANGENT_EXCESS_20261103_v63/report.json']
    report=dict(id='AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63',builder_passed=True,field='QQ',ambient_variables=names,
        common_q_selected_rows=all_q_selected,remaining_q_rows=remaining,component_generator_counts={k:len(v['equations']) for k,v in families.items()},
        actual_nilpotent_candidate=str(F),candidate_p_multiple=str(G),F_zero_coordinate=str(star),
        candidate_annihilator=list(map(str,annihilator)),
        component_presentations={k:list(map(str,v['equations'])) for k,v in families.items()},split_presentation=list(map(str,split)),
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'builder.json').write_text(json.dumps(report,indent=2)+'\n');print('EXACT IDEAL BUILDER SAVED',report['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
