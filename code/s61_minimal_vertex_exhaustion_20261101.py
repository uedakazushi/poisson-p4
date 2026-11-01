#!/usr/bin/env python3
"""Finite reduced M0/M1 vertex case split over a sufficiently general first jet.

The output distinguishes generic-base charts, divisor charts, and codimension-
two base exceptions. It does not classify saturated local rings at intersections.
"""
from __future__ import annotations
import argparse,hashlib,json
from itertools import combinations
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J5
from s61_nil22_infinity_20261101 import tensor,beta_matrix,sm
ROOT=Path(__file__).resolve().parents[1]
def graph_basis(F,g,axis,kvalue,k,field):
    h=s.Symbol('h');zero=[j for j in[5,6,7,8] if j!=axis]
    if axis in[2,4]:zero=[5,6,7,8,6-axis]
    sub={k:kvalue,**dict.fromkeys([g[j] for j in zero],0)}
    fs=[f.subs(sub) for f in F]+[h*g[axis]-1]+[g[j] for j in zero]
    variables=[h]+[gi for j,gi in enumerate(g) if j!=axis]+[g[axis]]
    gb=s.groebner(fs,*variables,domain=field,order='lex');assert len(gb.polys)==9 and h*g[axis]-1 in gb
    solution={}
    for gi in g:
        if gi==g[axis]:continue
        candidates=[f for f in gb if s.diff(f,gi)==1 and not (f-gi).has(*[v for v in variables if v not in[gi,g[axis]]])]
        assert len(candidates)==1;solution[gi]=s.factor(gi-candidates[0])
    assert all(s.cancel(f.subs(k,kvalue).subs(solution,simultaneous=True))==0 for f in F)
    return {'axis':axis,'k':str(kvalue),'groebner_basis':list(map(str,gb)),'coordinate_formulas':{str(a):str(b) for a,b in solution.items()},'parameter':str(g[axis]),'graph_extends_to_parameter_zero':True}
def roots_of(poly,k):
    poly=s.factor(poly);factors=s.factor_list(poly,k)[1];assert all(s.degree(f,k)==1 for f,e in factors)
    return sorted({s.solve(f,k)[0] for f,e in factors}),str(poly)
def run():
    m=S61();x,y,z,w=X[1:];u=s.symbols('u1:62');g=s.symbols('g0:9');k,K,L,U,V,h=s.symbols('k K L U V h');weights=m.reduced['parameter_weights'];one=[i for i,v in enumerate(weights) if v==1];two=[i for i,v in enumerate(weights) if v==2]
    G=sum(a*b for a,b in zip(g,[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3],strict=True));base={u[3]:0,u[11]:0,u[4]:0,u[12]:0,u[13]:k,u[48]:2*k+1,u[25]:U,u[34]:V,u[49]:K-5*V,u[50]:5*U-L}
    field=s.QQ.frac_field(K,L,U,V);extra_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';extra=json.loads(extra_path.read_text())['universal_new_equation'];families=[]
    high=[g[5]*g[6],g[5]*g[7],g[6]*g[8],g[7]*g[8],9*g[5]*g[8]+g[6]*g[7]];high_gb=s.groebner(high,*g[5:]);assert all(high_gb.reduce((a*b)**2)[1]==0 for a,b in combinations(g[5:],2))
    for label in['M0','M1']:
        Q=tensor([x*x,0,0,-s.diff(G,w),s.diff(G,z),0]) if label=='M0' else tensor([x*x/2,0,0,-x*z-s.diff(G,w),s.diff(G,z),x*x/2])
        M=beta_matrix(m,Q);F=[s.factor(f.subs(base)) for f in M*s.Matrix([u[i] for i in one])]
        assert [s.expand(-F[100]/6),s.expand(-F[101]/18),s.expand(F[74]/18),s.expand(F[75]/6),s.expand(-F[102]/3)]==high
        B=X[0]**2*J5;B[1:,1:]+=Q;normalized=m.normalize(B)['parameters'];values=[s.sympify(u[i] if i in one else normalized[i]).subs(base) for i in range(61)]
        G6=s.factor(sum(s.Rational(n,d)*s.prod(values[i] for i in mon) for mon,n,d in extra if sum(i in two for i in mon)==2));assert s.expand(G6-9*(2*k+1)*g[6]*g[7])==0,(label,str(G6))
        zero_high=dict.fromkeys(g[5:],0);F0=[s.factor(f.subs(zero_high)) for f in F];assert s.expand(g[2]*F0[71]+g[4]*F0[86]-s.Rational(6,5)*g[2]*g[4])==0
        baseline=[f.subs({g[2]:0,g[4]:0}) for f in F0]+[g[i] for i in[2,4,5,6,7,8]];baseline_gb=s.groebner(baseline,*g,domain=s.QQ.frac_field(k,K,L,U,V),order='lex');assert len(baseline_gb.polys)==9
        baseline_solution=s.solve(list(baseline_gb),list(g),dict=True)[0]
        # A canceled numerator in the generic solution can conceal a free
        # specialized coordinate. Enumerate zeros of actual triangular pivots,
        # not merely poles of the simplified solution.
        assert [i for i,f in enumerate(baseline[:134]) if f!=0]==[70,84,87]
        cpivot=s.factor(s.diff(baseline[87],g[3]));csol=s.solve(baseline[87],g[3])[0]
        apivot=s.factor(s.diff(baseline[84].subs(g[3],csol),g[0]));bpivot=s.factor(s.diff(baseline[70].subs(g[3],csol),g[1]))
        pivot_numerators=s.prod(s.fraction(p)[0] for p in [cpivot,apivot,bpivot]);bad_baseline=sorted({s.solve(f,k)[0] for f,e in s.factor_list(pivot_numerators,k)[1]})
        excluded_baseline=[-s.Rational(7,9),-s.Rational(1,2),-s.Rational(2,9)]
        free_baseline={-s.Rational(1,3):0,-s.Rational(2,3):1} if label=='M0' else {-s.Rational(4,9):0,-s.Rational(5,9):1}
        assert set(bad_baseline)==set(excluded_baseline)|set(free_baseline)
        baseline_ex=[];baseline_tails=[]
        for kk in bad_baseline:
            gb=s.groebner([f.subs(k,kk) for f in baseline],*g,domain=field,order='grevlex')
            if kk in excluded_baseline:
                assert list(gb)==[1];baseline_ex.append({'k':str(kk),'groebner_basis':['1']})
            else:
                axis=free_baseline[kk];assert len(gb.polys)==8
                sol=s.solve(list(gb),[gi for j,gi in enumerate(g) if j!=axis],dict=True)[0]
                assert all(s.expand(f.subs(k,kk).subs(sol,simultaneous=True))==0 for f in F)
                assert all(s.denom(s.factor(v)).is_Rational for v in sol.values())
                baseline_tails.append({'axis':axis,'k':str(kk),'parameter':str(g[axis]),'groebner_basis':list(map(str,gb)),'coordinate_formulas':{str(a):str(b) for a,b in sol.items()},'all_nuisance_parameters_valid':True})
        tails=[];case_proofs=[]
        for axis,row in[(2,86),(4,71)]:
            constraint=s.cancel(F0[row].subs(g[6-axis],0)/g[axis]);resultant=s.resultant(constraint,F0[87].subs(g[6-axis],0),g[3]);ks,factor=roots_of(resultant,k)
            expected={('M0',2):[-s.Rational(11,12),s.Rational(1,6)],('M0',4):[-s.Rational(7,6),-s.Rational(1,12)],('M1',2):[-s.Rational(11,12),-s.Rational(1,3)],('M1',4):[-s.Rational(2,3),-s.Rational(1,12)]}[label,axis];assert ks==expected
            case_proofs.append({'axis':axis,'resultant':factor,'allowed_k':list(map(str,ks))})
            for kk in ks:tails.append(graph_basis(F,g,axis,kk,k,field))
        for axis,other,cr,pr,lr in[(5,4,91,92,71),(8,2,73,72,86)]:
            sub=dict.fromkeys([g[i] for i in[5,6,7,8] if i!=axis],0);Fa=[s.factor(f.subs(sub)) for f in F]
            independent=[s.cancel(Fa[cr]/g[axis]),s.cancel(Fa[pr]/g[axis]),Fa[lr],h*g[other]-1]
            gb=s.groebner(independent,h,k,g[3],g[other],order='lex');assert list(gb)==[1]
            resultant=s.resultant(s.cancel(Fa[cr]/g[axis]),Fa[87].subs(g[other],0),g[3]);ks,factor=roots_of(resultant,k)
            assert ks==([-1] if axis==5 else[0]) if label=='M0' else ks==[-1,0]
            case_proofs.append({'axis':axis,'other_quadratic_axis':other,'other_nonzero_groebner':['1'],'resultant':factor,'allowed_k':list(map(str,ks))})
            for kk in ks:tails.append(graph_basis(F,g,axis,kk,k,field))
        isolated=[];axis_exceptions=[]
        for axis in[6,7]:
            other=4 if axis==6 else 2;sub=dict.fromkeys([g[i] for i in[5,6,7,8] if i!=axis],0);Fa=[s.factor(f.subs(sub)) for f in F]
            forcing_row=72 if axis==6 else 92;forcing=s.cancel(Fa[forcing_row]/(g[axis]*g[other]));expected_forcing=s.Rational(12,5)*(k+3 if axis==6 else k-2);assert s.expand(forcing-expected_forcing)==0
            Fz=[s.factor(f.subs(g[other],0)) for f in Fa];primary=[91,92] if axis==6 else[72,73];unknown=[g[2],g[3]] if axis==6 else[g[3],g[4]];sol=s.solve([Fz[i] for i in primary],unknown,dict=True)[0];Fz=[s.factor(f.subs(sol)) for f in Fz];sol2=s.solve([Fz[i] for i in([86,87] if axis==6 else[71,87])],[g[0],g[1]],dict=True)[0];residual={i:s.factor(f.subs(sol2)) for i,f in enumerate(Fz) if s.factor(f.subs(sol2))!=0};assert set(residual)=={70,84}
            allowed_denominator=g[axis]*(2*k+1)*(3*k-1)*(3*k+4)
            assert all(s.denom(s.cancel(allowed_denominator**8/s.denom(s.factor(value)))).is_Rational for value in list(sol.values())+list(sol2.values()))
            H=(3*k+1)*K-(6*k-2)*V;D=(3*k+2)*L-(6*k+8)*U;tailk=-s.Rational(2,3) if axis==6 else-s.Rational(1,3)
            base_eq=H if axis==6 else D;denom=D if axis==6 else H
            numerator=(2*(3*k+2) if axis==6 else 2*(3*k+1)) if label=='M0' else(3*k*(3*k+1)/2 if axis==6 else-3*(k+1)*(3*k+2)/2)
            cleared=[s.fraction(s.factor(v))[0] for v in residual.values()];target=[(k-tailk)*base_eq,(k-tailk)*(denom*g[axis]-numerator)]
            assert list(s.groebner(cleared,g[axis],K,L,U,V,k,order='grevlex',domain=s.QQ))==list(s.groebner(target,g[axis],K,L,U,V,k,order='grevlex',domain=s.QQ)),(label,axis,cleared,target)
            isolated.append({'axis':axis,'forced_other_zero_away_from':str(-3 if axis==6 else 2),'rational_solution':{str(a):str(b) for a,b in{g[other]:0,**sol,**sol2}.items()},'residuals':{str(i):str(v) for i,v in residual.items()},'base_divisor':str(base_eq),'axis_equation':str(denom*g[axis]-numerator),'tail_k':str(tailk)})
            tails.append(graph_basis(F,g,axis,tailk,k,field))
            exceptions=[-3 if axis==6 else 2,s.Rational(1,3),-s.Rational(4,3),-s.Rational(1,2)]
            for kk in exceptions:
                gb=s.groebner([f.subs(k,kk) for f in Fa]+[h*g[axis]-1]+[g[i] for i in[5,6,7,8] if i!=axis],h,*g,domain=field,order='grevlex');assert list(gb)==[1];axis_exceptions.append({'axis':axis,'k':str(kk),'groebner_basis':['1']})
        rankone={u[3]:0,u[11]:0,u[4]:-1,u[12]:0,u[13]:1,u[48]:2,u[25]:U,u[34]:V,u[49]:-5*V,u[50]:5*U};Fr=[s.factor(f.subs(rankone)) for f in M*s.Matrix([u[i] for i in one])];rfield=s.QQ.frac_field(U,V);rgb=s.groebner(Fr,*g,domain=rfield,order='grevlex');assert len(rgb.polys)==(8 if label=='M0' else 9)
        rs=s.solve(list(rgb),list(g[:-1] if label=='M0' else g),dict=True)[0];assert all(s.cancel(f.subs(rs,simultaneous=True))==0 for f in Fr)
        rankone_values=[s.sympify(u[i] if i in one else normalized[i]).subs(rankone) for i in range(61)]
        rankone_G6=s.expand(sum(s.Rational(n,d)*s.prod(rankone_values[i] for i in mon) for mon,n,d in extra if sum(i in two for i in mon)==2))
        assert s.expand(rankone_G6+9*(5*g[5]*g[6]+6*g[5]*g[7]-4*g[6]**2))==0
        assert rgb.reduce(rankone_G6)[1]==0
        families.append({'label':label,'vertex_equations':list(map(str,F)),'G6':str(G6),'baseline_groebner':list(map(str,baseline_gb)),'baseline_solution':{str(a):str(b) for a,b in baseline_solution.items()},'baseline_pivots':{'g3':str(cpivot),'g0':str(apivot),'g1':str(bpivot)},'baseline_pivot_roots':list(map(str,bad_baseline)),'baseline_exceptional_k':baseline_ex,'baseline_free_lines':baseline_tails,'case_resultants':case_proofs,'rational_tails':tails,'isolated_axis_charts':isolated,'axis_exceptional_k':axis_exceptions,'rankone_G6':str(rankone_G6),'rankone_equations':list(map(str,Fr)),'rankone_groebner':list(map(str,rgb)),'rankone_solution':{str(a):str(b) for a,b in rs.items()}})
        print(json.dumps({'component':label,'tails':len(tails),'additional_baseline_lines':len(baseline_tails),'isolated_charts':2,'rankone_dimension':1 if label=='M0' else 0}),flush=True)
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_infinity_flags.py',m.source/'symplectic_slice_reduced_equations.json',extra_path,ROOT/'research/notes/S61_R0_GENERIC_EXHAUSTION_20261101.md',ROOT/'reports/S61_MINIMAL_INFINITY_VERTEX_NECESSITY_20261101.json']
    return {'id':'S61_MINIMAL_VERTEX_EXHAUSTION_20261101','all_passed':True,'field':'QQ with explicitly declared parameter function fields','highest_cubic_radical_generators':list(map(str,[a*b for a,b in combinations(g[5:],2)])),'families':families,'scope':'Complete reduced M0/M1 infinity case coverage for sufficiently general beta after excluding codimension>=2 vertex-base exceptions. Does not determine saturation or local rings at rational-tail intersections or special points.','proof_note':'research/notes/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.md','input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'output':str(args.output)}))
if __name__=='__main__':main()
