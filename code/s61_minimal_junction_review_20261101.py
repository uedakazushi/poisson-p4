#!/usr/bin/env python3
"""Independent raw135 series and actual-branch audit at eight minimal junctions."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)

def standards(gb):
    assert gb.is_zero_dimensional
    leads=[p.LM(order=gb.order).exponents for p in gb.polys];basis=[];counts=[]
    for degree in range(30):
        cur=[]
        for inds in combinations_with_replacement(range(len(gb.gens)),degree):
            mon=tuple(inds.count(i) for i in range(len(gb.gens)))
            if not any(all(a<=b for a,b in zip(l,mon)) for l in leads):cur.append(mon)
        counts.append(len(cur));basis+=cur
        if not cur:return basis,counts
    raise AssertionError('standard monomial enumeration not exhausted')

def run():
    started=time.time();bp=ROOT/'reports/S61_MINIMAL_JUNCTION_BRANCHES_20261101.json';up=ROOT/'reports/S61_MINIMAL_JUNCTION_SATURATION_20261101.json'
    lower=json.loads(bp.read_text());upper=json.loads(up.read_text());assert lower['all_passed'] and upper['all_passed']
    m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];index={i:j for j,i in enumerate(two)}
    eqs=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    def polynomials(beta):
        out=[]
        for eq in eqs:
            f={}
            for mon,a,b in eq:
                mm=tuple(index[i] for i in mon if i in index);cc=s.Rational(a,b)*s.prod(beta[one.index(i)] for i in mon if i in one);f[mm]=f.get(mm,s.S.Zero)+cc
            out.append({mm:s.expand(cc) for mm,cc in f.items() if cc})
        return out
    fam=json.loads((ROOT/'reports/S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101.json').read_text());vertex=json.loads((ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json').read_text())
    g=s.symbols('g0:9');k,K,L,UU,VV=s.symbols('k K L U V');b=s.symbols('b0:6');pars=[k,K,L,UU,VV,*b]
    bv=s.Matrix([0,0,0,0,k,b[0],UU,b[1],b[2],VV,b[3],b[4],b[5],2*k+1,K-5*VV,5*UU-L])
    zz,bu,form=primitive_beta(m);vbeta=s.Matrix([form[i,j].diff(zz[h]) for (i,j),h in COORDS]);basis=vbeta.jacobian(bu);rr=list(basis.T.rref()[1]);binv=basis[rr,:].inv();actions=s.Matrix.hstack(*[binv*infinitesimal_action(vbeta,H)[rr,:] for H in sp_basis()])
    out=[]
    for low,high in zip(lower['cases'],upper['cases']):
        timer=time.time();assert (low['kind'],low['k'])==(high['kind'],high['k']);kind=low['kind'];kk=s.Rational(low['k']);ld=low['local'];point=list(map(s.Rational,ld['point61']));q=s.Matrix([point[i] for i in two]);beta=[point[i] for i in one];F=polynomials(beta)
        cols=ld['chart_columns'];fixed=next(i for i in range(45) if i not in cols);A=s.zeros(135,45)
        for row,f in enumerate(F):
            assert sum(cc*s.prod(q[i] for i in mm) for mm,cc in f.items() if len(mm)==2)==0
            for mm,cc in f.items():
                if len(mm)==2:
                    for j,i in enumerate(cols):A[row,j]+=cc*sum(q[mm[1-h]] for h in range(2) if mm[h]==i)
                elif len(mm)==1:A[row,44]+=cc*q[mm[0]]
        assert A==M(ld['jacobian']);T=M(ld['tangent_basis']);n=T.cols;assert A*T==s.zeros(135,n) and T.rank()==n
        rows=high['jacobian_pivot_rows'];piv=high['jacobian_pivot_columns'];rank=A.to_DM().rank();assert rank==len(rows)==45-n
        U=s.zeros(45,rank);inverse=A.extract(rows,piv).inv(method='DM')
        for j,i in enumerate(piv):U[i,:]=inverse[j,:]
        assert A[rows,:]*U==s.eye(rank)
        # Work in the full46 homogeneous coordinates and multiply graded coefficient arrays.
        E={1:s.zeros(46,n)}
        for j,i in enumerate(cols):E[1][i,:]=T[j,:]
        E[1][45,:]=T[44,:]
        mons={d:list(combinations_with_replacement(range(n),d)) for d in range(1,high['initial_degree_bound']+1)}
        residual={}
        for degree in range(2,high['initial_degree_bound']+1):
            pos={mm:i for i,mm in enumerate(mons[degree])};raw=s.zeros(135,len(pos))
            for row,f in enumerate(F):
                for mm,cc in f.items():
                    i,j=(*mm,*([45]*(2-len(mm))))
                    if i==fixed or j==fixed:continue
                    for aa in range(1,degree):
                        bb=degree-aa
                        left=[(h,c) for h,c in enumerate(E[aa][i,:]) if c];right=[(h,c) for h,c in enumerate(E[bb][j,:]) if c]
                        for h,c in left:
                            for l,d in right:raw[row,pos[tuple(sorted(mons[aa][h]+mons[bb][l]))]]+=cc*c*d
            lift=-U*raw[rows,:];residual[degree]=raw+A*lift;assert residual[degree][rows,:]==s.zeros(rank,len(pos))
            E[degree]=s.zeros(46,len(pos))
            for j,i in enumerate(cols):E[degree][i,:]=lift[j,:]
            E[degree][45,:]=lift[44,:]
        # Build the entire truncated row ideal through the chosen degree. This
        # differs from the main cubic left-kernel elimination and includes all multiples.
        degree_bound=high['initial_degree_bound'];mmall=sum([mons[d] for d in range(2,degree_bound+1)],[]);ix={mm:i for i,mm in enumerate(mmall)}
        truncated=s.Matrix.hstack(*[residual[d] for d in range(2,degree_bound+1)]);rr,pv=truncated.to_DM().rref();base=rr.to_Matrix()[:len(pv),:];multiples=[]
        for row in base.tolist():
            for d in range(1,degree_bound-1):
                for fac in mons[d]:
                    val=[s.S.Zero]*len(mmall)
                    for h,mm in enumerate(mmall):
                        if row[h] and len(mm)+d<=degree_bound:val[ix[tuple(sorted(mm+fac))]]+=row[h]
                    if any(val):multiples.append(val)
        full=base.col_join(s.Matrix(multiples)) if multiples else base;rr,pv=full.to_DM().rref();rr=rr.to_Matrix();v=s.symbols(f'v0:{n}');forms=[];initial_ranks={}
        for row,pc in enumerate(pv):
            d=len(mmall[pc]);initial_ranks[d]=initial_ranks.get(d,0)+1
            forms.append(sum(rr[row,j]*s.prod(v[i] for i in mm) for j,mm in enumerate(mmall) if len(mm)==d))
        coeff=high['parameter_coefficients'];assert coeff[0]==1
        substitute=-sum(coeff[i]*v[i] for i in range(1,n))
        # Reversed variable order checks the same ideals with different pivots.
        gb0=s.groebner([f.subs(v[0],substitute) for f in forms],*reversed(v[1:]),order='grevlex');basis0,count0=standards(gb0)
        gb1=s.groebner([f.subs(v[0],1+substitute) for f in forms],*reversed(v[1:]),order='grevlex');basis1,count1=standards(gb1)
        assert len(basis1)==high['generic_rank']==low['length_lower']
        # Main source: exact differential reconstructed from the baseline formulas.
        block=next(c for c in vertex['families'] if c['label']==kind);qexpr=M((fam if kind=='M0' else fam['M1_family'])['q_coordinates']);sol={s.sympify(a):s.sympify(val) for a,val in block['baseline_solution'].items()};qmain=qexpr.subs(sol,simultaneous=True)
        sample={k:kk,K:19,L:5,UU:2,VV:3,**dict(zip(b,[2,3,5,7,11,13]))};assert bv.subs(sample)==s.Matrix(beta)
        canonical_c=-2/(3*(2*k+1)) if kind=='M0' else -(3*k+1)/(3*(2*k+1));assert 2*kk+1!=0 and s.diff(canonical_c,k).subs(k,kk)!=0
        assert qmain.subs(sample)==q and all(s.denom(s.cancel(c)).subs(sample)!=0 for c in qmain)
        BD=bv.jacobian(pars).row_join(actions.subs(dict(zip(bu,beta))));assert BD==M(low['beta_differential']) and BD.rank()==16 and BD[:,1:].rank()==16
        DQ=qmain.jacobian(pars).subs(sample);DQ-=q*DQ[fixed,:]/q[fixed]
        x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);shape=[sol[gi].subs(sample) for gi in g];G=sum(c*f for c,f in zip(shape,[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3]));Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-s.diff(G,w)-(x*z if kind=='M1' else 0);Q[1,3]=s.diff(G,z);Q-=Q.T;orbit=[]
        for H in sp_basis():
            dQ=sum(((H*Y)[i]*Q.diff(Y[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T;P=X[0]**2*J;P[1:,1:]+=dQ;pp=m.normalize(P)['parameters'];dq=s.Matrix([pp[i] for i in two]);dq-=q*dq[fixed]/q[fixed];orbit.append(dq[cols,:])
        QD=DQ[cols,:].row_join(s.Matrix.hstack(*orbit));assert QD==M(low['main_q_differential'])
        tangents=[];t=s.Symbol('t');lower_sum=0
        for branch in low['branches']:
            direction=M(branch['tangent']);assert A*T*direction==s.zeros(135,1)
            parameter=sum(c*d for c,d in zip(coeff,direction));assert parameter!=0;tangents.append(direction);lower_sum+=branch['lower']
            if branch['kind']=='main':
                sd=M(branch['source_direction']);assert BD*sd==s.zeros(16,1) and sd[0]!=0 and sd[0]==s.Rational(branch['dk']);assert T*direction==(QD*sd).col_join(s.zeros(1,1))
            else:
                curve=M(branch['q_curve']);center=s.Rational(branch['center']);assert curve.subs(t,center)==q
                assert all(s.expand(sum(cc*s.prod(curve[i] for i in mm) for mm,cc in f.items() if len(mm)==2))==0 for f in F)
                dq=curve.diff(t).subs(t,center);dq-=q*dq[fixed]/q[fixed];assert T*direction==dq[cols,:].col_join(s.zeros(1,1))
                if branch['kind']=='free' or branch['lower']==3:
                    path=ROOT/'reports'/('S61_MINIMAL_BASELINE_THICKENINGS_20261101.json' if branch['kind']=='free' else 'S61_MINIMAL_CUBIC_TAIL_20261101.json');known=json.loads(path.read_text());known=known['cases'][kind]['point']['local_data'] if branch['kind']=='free' else known['point']['local_data'];pcheck=point.copy()
                    for j,i in enumerate(two):pcheck[i]=curve[j].subs(t,1)
                    assert list(map(str,pcheck))==known['point61']
        assert s.Matrix.hstack(*tangents).rank()==len(tangents) and lower_sum==len(basis1)
        record=dict(kind=kind,k=str(kk),tangent_dimension=n,initial_degree_bound=degree_bound,initial_ranks=initial_ranks,zero_section_standard_degrees=count0,localized_standard_degrees=count1,generic_rank=len(basis1),branch_number=len(tangents),all_actual_branch_differentials_verified=True,seconds=time.time()-timer);out.append(record);print(json.dumps(record),flush=True)
    paths=[Path(__file__),bp,up,ROOT/'code/s61_minimal_junction_branches_20261101.py',ROOT/'code/s61_minimal_junction_saturation_20261101.py',ROOT/'research/notes/S61_MINIMAL_JUNCTION_BRANCHES_20261101.md',ROOT/'research/notes/S61_MINIMAL_JUNCTION_SATURATION_20261101.md',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json',ROOT/'reports/S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101.json',ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json',ROOT/'reports/S61_MINIMAL_CUBIC_TAIL_20261101.json',ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json',m.source/'symplectic_slice_reduced_equations.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    return dict(id='S61_MINIMAL_JUNCTION_REVIEW_20261101',all_passed=True,cases=out,accepted_pointwise_saturation_removal=True,limitations=['Universal generic rank and known curve formal-ring results remain explicit prior dependencies.','No whole-stratum generic-beta removal or complete local multiplication is asserted.'],seconds=time.time()-started,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args();result=run()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'cases':len(result['cases']),'seconds':result['seconds']}))
