#!/usr/bin/env python3
"""Independent complete-chart, pivot, and exception audit for minimal infinity."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as s
from s61_model import S61, X, J as J5

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json'


def run():
    report=json.loads(INPUT.read_text());model=S61()
    g=s.symbols('g0:9');k,K,L,U,V,h=s.symbols('k K L U V h')
    u=s.symbols('u1:62');x,y,z,w=X[1:]
    weights=model.reduced['parameter_weights']
    one=[i for i,v in enumerate(weights) if v==1]
    two=[i for i,v in enumerate(weights) if v==2]
    base={u[3]:0,u[11]:0,u[4]:0,u[12]:0,u[13]:k,u[48]:2*k+1,
          u[25]:U,u[34]:V,u[49]:K-5*V,u[50]:5*U-L}
    rankone=base|{u[4]:-1,u[13]:1,u[48]:2,u[49]:-5*V,u[50]:5*U}
    field=s.QQ.frac_field(K,L,U,V)
    high=[g[5]*g[6],g[5]*g[7],g[6]*g[8],g[7]*g[8],9*g[5]*g[8]+g[6]*g[7]]
    highgb=s.groebner(high,*reversed(g[5:]),domain=s.QQ)
    assert all(highgb.reduce((a*b)**2)[1]==0 for a,b in combinations(g[5:],2))
    H=(3*k+1)*K-(6*k-2)*V;D=(3*k+2)*L-(6*k+8)*U
    assert s.gcd(3*k+1,6*k-2)==s.gcd(3*k+2,6*k+8)==1
    results=[]
    for family in report['families']:
        label=family['label']
        potential=sum(c*m for c,m in zip(g,[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3]))
        Q=s.zeros(4)
        Q[0,1]=x*x if label=='M0' else x*x/2
        Q[1,2]=-s.diff(potential,w)-(x*z if label=='M1' else 0)
        Q[1,3]=s.diff(potential,z)
        Q[2,3]=x*x/2 if label=='M1' else 0
        Q=Q-Q.T
        B=X[0]**2*J5;B[1:5,1:5]+=Q
        q=model.normalize(B)['parameters']
        values=[u[i] if i in one else q[i] for i in range(61)]
        full=[]
        for equation in model.reduced['equations']:
            full.append(s.expand(s.Add(*[s.Rational(n,d)*s.prod(values[i] for i in mon)
                                         for mon,n,d in equation if sum(i in two for i in mon)==2])))
        F=[s.expand(v.subs(base)) for v in full]
        Fr=[s.expand(v.subs(rankone)) for v in full]
        assert all(s.expand(f-s.sympify(v))==0 for f,v in zip(F,family['vertex_equations']))
        assert all(s.expand(f-s.sympify(v))==0 for f,v in zip(Fr,family['rankone_equations']))
        assert [s.expand(-F[100]/6),s.expand(-F[101]/18),s.expand(F[74]/18),
                s.expand(F[75]/6),s.expand(-F[102]/3)]==high

        # A uniform triangular proof, stronger than reading denominators of
        # a function-field Groebner solution: all three pivots depend only on k.
        zero=dict.fromkeys([g[i] for i in [2,4,5,6,7,8]],0)
        baseline=[s.factor(f.subs(zero)) for f in F]
        assert {i for i,f in enumerate(baseline) if f!=0}=={70,84,87}
        A,bv=s.linear_eq_to_matrix([baseline[87]],[g[3]])
        value3=s.factor(bv[0]/A[0,0])
        sub={g[3]:value3}
        C,cv=s.linear_eq_to_matrix([baseline[84].subs(sub),baseline[70].subs(sub)],[g[0],g[1]])
        assert not s.factor(C.det()).has(K,L,U,V)
        solution={**zero,**sub,**dict(zip([g[0],g[1]],C.inv()*cv))}
        assert all(s.cancel(f.subs(solution,simultaneous=True))==0 for f in F)
        pivots=[s.factor(A[0,0]),s.factor(C[0,0]),s.factor(C[1,1])]
        singular=set()
        for p in pivots:
            numerator,denominator=s.fraction(p)
            singular.update(s.solve(numerator,k));singular.update(s.solve(denominator,k))
        old_exceptions={-s.Rational(1,2),-s.Rational(7,9),-s.Rational(2,9)}
        free_values={-s.Rational(1,3):0,-s.Rational(2,3):1} if label=='M0' else {
            -s.Rational(4,9):0,-s.Rational(5,9):1}
        assert singular==old_exceptions|set(free_values)
        for kk in old_exceptions:
            gb=s.groebner([f.subs(k,kk) for f in baseline],g[3],g[1],g[0],domain=field)
            assert gb.contains(s.S.One)
        released=[]
        for kk,axis in free_values.items():
            gb=s.groebner([f.subs(k,kk) for f in baseline],g[3],g[1],g[0],domain=field)
            assert len(gb.polys)==2 and all(not f.has(g[axis]) for f in gb)
            released.append({'k':str(kk),'free_coordinate':str(g[axis]),
                             'groebner_basis':list(map(str,gb))})
        saved_free={(s.Rational(v['k']),v['parameter']) for v in family['baseline_free_lines']}
        assert saved_free=={(kk,str(g[axis])) for kk,axis in free_values.items()}

        # Recompute each saturated graph from the independently reconstructed
        # restrictions, rather than trusting its reported shape alone.
        tails=[]
        for tail in family['rational_tails']:
            axis=tail['axis'];kk=s.Rational(tail['k'])
            zeros=[i for i in [5,6,7,8] if i!=axis]
            if axis in [2,4]:zeros=[5,6,7,8,6-axis]
            subzero=dict.fromkeys([g[i] for i in zeros],0)
            fs=[f.subs(k,kk).subs(subzero) for f in F]+[g[i] for i in zeros]+[h*g[axis]-1]
            variables=[h,*[v for v in g if v!=g[axis]],g[axis]]
            gb=s.groebner(fs,*variables,domain=field,order='lex')
            graph=[s.sympify(v) for v in tail['groebner_basis']]
            graphgb=s.groebner(graph,*variables,domain=field,order='lex')
            assert all(gb.reduce(v)[1]==0 for v in graph)
            assert all(graphgb.reduce(v)[1]==0 for v in fs)
            sol={s.sympify(a):s.sympify(b) for a,b in tail['coordinate_formulas'].items()}
            assert all(s.cancel(f.subs(k,kk).subs(sol,simultaneous=True))==0 for f in F)
            assert all(not s.denom(s.factor(v)).has(g[axis]) for v in sol.values())
            tails.append({'axis':axis,'k':str(kk),'both_ideal_containments':True})

        axes=[]
        for axis in [6,7]:
            other=4 if axis==6 else 2
            zeroaxis=dict.fromkeys([g[i] for i in [5,6,7,8] if i!=axis],0)
            Fa=[s.factor(f.subs(zeroaxis)) for f in F]
            Fz=[s.factor(f.subs(g[other],0)) for f in Fa]
            rows=[91,92] if axis==6 else [72,73]
            unknown=[g[2],g[3]] if axis==6 else [g[3],g[4]]
            A,av=s.linear_eq_to_matrix([Fz[i] for i in rows],unknown)
            firstdet=s.factor(A.det())
            assert firstdet==-3*g[axis]**2*(3*k-1)*(3*k+4)/25
            sol=dict(zip(unknown,A.inv()*av))
            rows2=[86,87] if axis==6 else [71,87]
            C,cv=s.linear_eq_to_matrix([Fz[i].subs(sol) for i in rows2],[g[0],g[1]])
            seconddet=s.factor(C.det())
            expected=-6*g[axis]**2*(2*k+1)*(3*k+4 if axis==6 else 3*k-1)/25
            assert s.expand(seconddet-expected)==0
            sol2=dict(zip([g[0],g[1]],C.inv()*cv))
            residual=[s.factor(v.subs(sol).subs(sol2)) for v in Fz]
            residual=[s.fraction(v)[0] for v in residual if v!=0]
            item=next(v for v in family['isolated_axis_charts'] if v['axis']==axis)
            target=[(k-s.Rational(item['tail_k']))*s.sympify(item['base_divisor']),
                    (k-s.Rational(item['tail_k']))*s.sympify(item['axis_equation'])]
            var=[k,V,U,L,K,g[axis]]
            bg=s.groebner(residual,*var,domain=s.QQ,order='grevlex')
            tg=s.groebner(target,*var,domain=s.QQ,order='grevlex')
            assert all(bg.reduce(v)[1]==0 for v in target)
            assert all(tg.reduce(v)[1]==0 for v in residual)
            exceptions=[v for v in family['axis_exceptional_k'] if v['axis']==axis]
            assert {s.Rational(v['k']) for v in exceptions}=={
                s.Integer(-3 if axis==6 else 2),s.Rational(1,3),-s.Rational(4,3),-s.Rational(1,2)}
            for entry in exceptions:
                fs=[f.subs(k,s.Rational(entry['k'])) for f in Fa]+[h*g[axis]-1]+list(zeroaxis)
                gb=s.groebner(fs,h,*g,domain=field,order='grevlex')
                assert gb.contains(s.S.One)
            axes.append({'axis':axis,'first_pivot_determinant':str(firstdet),
                         'second_pivot_determinant':str(seconddet),'residual_ideal_equality':True,
                         'exceptional_unit_ideals':len(exceptions)})
        rankgb=s.groebner(Fr,*g,domain=s.QQ.frac_field(U,V),order='grevlex')
        saved=[s.sympify(v) for v in family['rankone_groebner']]
        savedgb=s.groebner(saved,*g,domain=s.QQ.frac_field(U,V),order='grevlex')
        assert all(rankgb.reduce(v)[1]==0 for v in saved)
        assert all(savedgb.reduce(v)[1]==0 for v in Fr)
        results.append({'label':label,'all134_restrictions_reconstructed':True,
                        'uniform_baseline_pivots':list(map(str,pivots)),
                        'baseline_free_lines':released,'rational_graphs':tails,
                        'axis_pivots':axes,'rankone_ideal_equality':True})
    paths=[Path(__file__),INPUT,ROOT/'code/s61_minimal_vertex_exhaustion_20261101.py',
           ROOT/'research/notes/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.md',
           ROOT/'research/notes/S61_R0_GENERIC_EXHAUSTION_20261101.md']
    return {'id':'S61_MINIMAL_VERTEX_REVIEW_20261101','all_passed':True,'field':'QQ and declared nuisance function fields',
            'families':results,'H_D_independent_on_every_k_fiber':True,
            'scope':'Independent restrictions, saturated graph and residual ideal equalities, uniform pivot denominators, exceptional unit ideals; smooth normalized-presentation and historical normal-form dependencies are reviewed in the note.',
            'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();result=run()
    with args.output.open('x') as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps({'all_passed':True,'output':str(args.output)}))
