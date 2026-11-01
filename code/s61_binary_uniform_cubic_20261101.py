#!/usr/bin/env python3
"""A universal cubic rank bound from a supported lift and full Bianchi syzygies."""
import argparse
import hashlib
import json
import math
import sys
import time
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61, X
from s61_flag_transverse_obstruction_20261101 import specialize
from elliptic_weight1_dominance import projective_affine


def full_raw_jacobi_maps(root, model):
    """Reconstruct FULL old Jacobi coefficients, not their q-leading parts."""
    sys.path.insert(0, str(root/'legacy/v02/code'))
    import reduce_symplectic_bianchi as old
    source=root/'legacy/v02/data'
    cert=json.loads((source/'symplectic_slice_certificate.json').read_text())
    data=json.loads((source/'symplectic_slice_equations.json').read_text())
    raw=cert['raw_slice_indices'];free=cert['raw_free_indices']
    expr_raw={k:{(i,):old.F(1)} for i,k in enumerate(free)}
    expr_raw.update({v['raw_index']:old.decode(v['polynomial']) for v in data['eliminated_raw_variables']})
    expr={i:expr_raw[k] for i,k in enumerate(raw)}
    allrows=[old.substitute(v,expr) for v in old.s.equations_from_raw(raw)]
    gens=[old.decode(v) for v in data['equations']]
    weights=data['parameter_weights']
    degree=lambda q:sum(weights[i] for i in next(iter(q)))
    i4=[i for i,q in enumerate(gens) if degree(q)==4]
    i5=[i for i,q in enumerate(gens) if degree(q)==5]
    assert len(i4)==70 and len(i5)==84
    assert [old.decode(v) for v in model.reduced['equations']]==[gens[i] for i in model.reduced['old_generator_indices']]
    bases=[old.LinearCoordinates([gens[i] for i in ids]) for ids in [i4,i5]]
    maps=[[],[]];lower_count=0
    for triple in combinations(range(1,5),3):
        for alpha in old.p.homogeneous_monomials(4):
            target=allrows[old.p.FOUR_VECTOR_INDEX[((0,)+triple,alpha)]]
            yd=4-alpha[0]
            if yd<3:
                assert not target
                lower_count+=1
                continue
            k=yd-3;coords=bases[k].coordinates(target);pred={}
            ids=[i4,i5][k]
            for i,c in coords.items():old.sparse_add(pred,gens[ids[i]],c)
            assert pred==target
            maps[k].append((tuple(i-1 for i in triple),tuple(alpha[1:]),coords))
    assert lower_count==60 and [len(v) for v in maps]==[80,140]
    return maps,i4,i5,lower_count


def bianchi_constants(root,model,u,fixed,raw):
    maps,i4,i5,lower_count=raw
    oldmap={v:i for i,v in enumerate(model.reduced['old_generator_indices'])}
    identities=json.loads((root/'legacy/v02/data/symplectic_bianchi_identities.json').read_text())
    eliminated={v['eliminated_old_generator']:v for v in identities['generator_redundancy_identities']}
    e=s.symbols('e0:134');y=s.symbols('y0:4')
    def generator(old):
        if old in oldmap:return e[oldmap[old]]
        a=eliminated[old]
        return sum(s.Rational(n,d)*e[oldmap[j]] for j,n,d in a['constant_coefficients_on_retained_old_generators'])+sum(s.Rational(n,d)*u[v]*e[oldmap[j]] for j,v,n,d in a['linear_coefficients_on_weight4_old_generators'])
    T=[]
    for entries,ids in zip(maps,[i4,i5]):
        comp={tr:s.S.Zero for tr in combinations(range(4),3)}
        for triple,mon,coords in entries:
            comp[triple]+=s.prod(v**a for v,a in zip(y,mon))*sum(s.Rational(c)*generator(ids[i]) for i,c in coords.items())
        T.append({tr:s.expand(v.subs(fixed)) for tr,v in comp.items()})
    values=[s.S.Zero]*61
    for i,w in enumerate(model.reduced['parameter_weights']):
        if w==1:values[i]=u[i].subs(fixed)
    p0=model.reconstruct(values);values[29]=1;p1=model.reconstruct(values)
    P=projective_affine(p1-p0,X,y).applyfunc(s.expand);Q=s.zeros(4);C=s.zeros(4)
    for i,j in combinations(range(4),2):
        for mon,c in s.Poly(P[i,j],*y).terms():
            term=c*s.prod(v**a for v,a in zip(y,mon))
            if sum(mon)==2:Q[i,j]+=term
            elif sum(mon)==3:C[i,j]+=term
            else:assert c==0
        Q[j,i]=-Q[i,j];C[j,i]=-C[i,j]
    assert Q[1,2]==y[0]*y[3] and sum(bool(Q[i,j]) for i,j in combinations(range(4),2))==1
    def bracket(B,comp):
        def tri(indices):
            if len(set(indices))<3:return 0
            sign=(-1)**sum(indices[i]>indices[j] for i in range(3) for j in range(i+1,3))
            return sign*comp[tuple(sorted(indices))]
        result=0
        for i in range(4):
            others=tuple(k for k in range(4) if k!=i)
            for h in range(4):
                if B[i,h]:result+=(-1)**i*B[i,h]*s.diff(tri(others),y[h])
        for i,k in combinations(range(4),2):
            others=tuple(a for a in range(4) if a not in (i,k))
            for h in range(4):
                if B[i,k].has(y[h]):result+=(-1)**(i+k)*s.diff(B[i,k],y[h])*tri((h,*others))
        return s.expand(result)
    expressions=[bracket(Q,T[0]),bracket(C,T[0])+bracket(Q,T[1]),bracket(C,T[1])]
    rows=[];labels=[]
    for degree,expr in zip([4,5,6],expressions):
        for mon,c in s.Poly(expr,*y).terms():
            if c==0:continue
            assert sum(mon)==degree
            rows.append([s.expand(c).coeff(v) for v in e]+[s.S.Zero]);labels.append([degree,list(mon)])
    return s.Matrix(rows),labels,C,lower_count


def verify():
    start=time.time();root=Path(__file__).resolve().parents[1]
    data_path=root/'reports/S61_BINARY_QUADRATIC_PROJECTION_DATA_20261101.json'
    lift_path=root/'reports/S61_BINARY_CUBIC_LIFT_DATA_20261101.json'
    data=json.loads(data_path.read_text());lift=json.loads(lift_path.read_text())
    R,*gens=ring(','.join(data['variables']),s.QQ);FF=R.to_field()
    symbols=s.symbols(' '.join(data['variables']));local=dict(zip(data['variables'],symbols))
    pol=lambda v:R.from_expr(s.sympify(v,locals=local))
    frac=lambda v:FF.from_expr(s.sympify(v,locals=local))
    projection=[[frac(v) for v in row] for row in data['projection']]
    selected=data['jacobian_rows'];assert selected==lift['jacobian_rows']
    normal=[2,4,5,6,7,8,9];parameters=[0,1,3]
    T=[[frac(row[j]) for j in normal] for row in data['tangent']]
    for j in range(7):
        den=R.one
        for row in T:den=den.lcm(row[j].denom)
        for row in T:row[j]=row[j].numer*den.exquo(row[j].denom)
    assert all(v==0 for v in T[44])
    model=S61();u=s.symbols('u1:62')
    fixedbeta={s.Symbol(k):s.Rational(v) for k,v in data['fixed_beta'].items()}
    beta={i:u[i].subs(fixedbeta) for i,w in enumerate(model.reduced['parameter_weights']) if w==1}
    beta[11]=beta[12]=s.S.Zero
    _,equations=specialize(model,beta,True)
    chart_fixed=data['chart_fixed_q'];qcols=[i for i in range(45) if i!=chart_fixed]
    J=[[R.zero]*45 for _ in range(135)];terms=[]
    for r,f in enumerate(equations):
        row=[]
        for mon,c in f.items():
            c=pol(c)
            if len(mon)==2:
                i,j=mon
                if i!=chart_fixed and j==chart_fixed:J[r][qcols.index(i)]+=c
                if j!=chart_fixed and i==chart_fixed:J[r][qcols.index(j)]+=c
            elif len(mon)==1 and mon[0]==chart_fixed:J[r][-1]+=c
            if chart_fixed not in mon:row.append(tuple([qcols.index(k) for k in mon]+[44]*(2-len(mon)))+(c,))
        terms.append(row)
    Pnum=[];Pden=[]
    for r in range(135):
        den=R.one
        for v in projection[r]:den=den.lcm(v.denom)
        nums=[v.numer*den.exquo(v.denom) for v in projection[r]];Pnum.append(nums);Pden.append(den)
        for j in range(45):assert den*J[r][j]==sum((v*J[selected[k]][j] for k,v in enumerate(nums) if v),R.zero)
    for i,r in enumerate(selected):assert all(projection[r][j]==FF(int(i==j)) for j in range(35))
    for r in range(135):
        for a in range(7):assert sum((c*T[j][a] for j,c in enumerate(J[r]) if c),R.zero)==0
    mon2=list(combinations_with_replacement(range(7),2));ix={m:i for i,m in enumerate(mon2)}
    pair_cache={}
    def product(i,j):
        key=tuple(sorted((i,j)))
        if key not in pair_cache:
            result=[R.zero]*28
            for a,v in enumerate(T[i]):
                if not v:continue
                for b,w in enumerate(T[j]):
                    if w:result[ix[tuple(sorted((a,b)))]]+=v*w
            pair_cache[key]=result
        return pair_cache[key]
    Q=[]
    for row in terms:
        values=[R.zero]*28
        for i,j,c in row:
            if i==44 or j==44:continue
            for n,v in enumerate(product(i,j)):
                if v:values[n]+=c*v
        Q.append(values)
    A=[];nz=[]
    for r in range(135):
        row=[Pden[r]*Q[r][j]-sum((v*Q[selected[k]][j] for k,v in enumerate(Pnum[r]) if v),R.zero) for j in range(28)]
        A.append(row)
        if any(row):nz.append(r)
    assert nz==[19,60,117,123,133]
    support=lift['support_columns'];assert len(support)==19 and 44 not in support
    L=[]
    for rr in lift['annihilator']:
        values=[frac(v) for v in rr];den=R.one
        for v in values:den=den.lcm(v.denom)
        L.append([v.numer*den.exquo(v.denom) for v in values])
    assert len(L)==16
    for row in L:
        for j in support:assert sum((c*J[selected[k]][j] for k,c in enumerate(row) if c),R.zero)==0
        for n in range(28):assert sum((c*Q[selected[k]][n] for k,c in enumerate(row) if c),R.zero)==0
    # These annihilator identities establish existence of a correct quadratic
    # implicit lift supported in the19 columns, without expanding that lift.
    B=[]
    for row in terms:
        values=[[R.zero]*19 for _ in range(7)]
        for i,j,c in row:
            for a in range(7):
                if j in support and T[i][a]:values[a][support.index(j)]+=c*T[i][a]
                if i in support and T[j][a]:values[a][support.index(i)]+=c*T[j][a]
        B.append(values)
    cubic_rows=[]
    for r in range(135):
        if r in nz:continue
        nonzero=False
        for a in range(7):
            for j in range(19):
                value=Pden[r]*B[r][a][j]-sum((c*B[selected[k]][a][j] for k,c in enumerate(Pnum[r]) if c),R.zero)
                nonzero=nonzero or bool(value)
        if nonzero:cubic_rows.append(r)
    assert cubic_rows==[18,20,23,39,40,43,59,61,68,90,92,99,102,113,115,116,120,122,129,132,134]
    sample={u[i]:s.Rational(j+2) for j,i in enumerate([i for i,w in enumerate(model.reduced['parameter_weights']) if w==1])}
    sample.update(fixedbeta);sample[u[11]]=sample[u[12]]=0
    vals=[s.QQ.from_sympy(sample[v]) for v in symbols];powers={}
    def ev(f):
        value=s.QQ.zero
        for mon,c in f.items():
            if mon not in powers:powers[mon]=math.prod(v**a for v,a in zip(vals,mon))
            value+=c*powers[mon]
        return s.QQ.to_sympy(value)
    # Check that the three chosen parameter coordinates restrict independently
    # to the common quadratic radical, which equals the smooth binary tangent
    # space on the constant-rank open used in the proof.
    Tf=s.Matrix(data['tangent']).subs(sample)
    assert Tf[-1,:]==s.zeros(1,10)
    Tf=[[s.QQ.from_sympy(v) for v in row] for row in Tf.tolist()]
    fullmon=list(combinations_with_replacement(range(10),2));fullindex={v:i for i,v in enumerate(fullmon)};numeric_pairs={}
    def numeric_product(i,j):
        key=tuple(sorted((i,j)))
        if key not in numeric_pairs:
            values=[s.QQ.zero]*55
            for a,v in enumerate(Tf[i]):
                if not v:continue
                for b,w in enumerate(Tf[j]):
                    if w:values[fullindex[tuple(sorted((a,b)))]]+=v*w
            numeric_pairs[key]=values
        return numeric_pairs[key]
    Qfull={}
    for r in set(selected+nz):
        row=[s.QQ.zero]*55
        for i,j,c in terms[r]:
            if i==44 or j==44:continue
            c=s.QQ.from_sympy(ev(c))
            if c:
                for n,v in enumerate(numeric_product(i,j)):
                    if v:row[n]+=c*v
        Qfull[r]=row
    hessians=[]
    for r in nz:
        ps=[s.QQ.from_sympy(ev(v.numer))/s.QQ.from_sympy(ev(v.denom)) for v in projection[r]]
        row=[Qfull[r][j]-sum((c*Qfull[selected[k]][j] for k,c in enumerate(ps) if c),s.QQ.zero) for j in range(55)]
        h=s.zeros(10)
        for c,(i,j) in zip(row,fullmon):
            h[i,j]+=s.QQ.to_sympy(c);h[j,i]+=s.QQ.to_sympy(c)
        hessians.extend(h.tolist())
    radical=s.Matrix.hstack(*s.Matrix(hessians).nullspace())
    assert radical.shape==(10,3)
    parameter_minor=radical[parameters,:].det();assert parameter_minor!=0
    Ls=s.Matrix([[ev(v) for v in row] for row in L]);Js=s.Matrix([[ev(J[r][j]) for j in support] for r in selected])
    assert len(Ls.to_DM().rref()[1])==16 and len(Js.to_DM().rref()[1])==19
    Aq=s.Matrix([[ev(v) for v in A[r]] for r in nz]);assert len(Aq.to_DM().rref()[1])==5
    raw=full_raw_jacobi_maps(root,model)
    M,labels,C,lower_count=bianchi_constants(root,model,u,fixedbeta|{u[11]:0,u[12]:0},raw)
    assert M[:,nz]==s.zeros(M.rows,5)
    Ms=M[:,cubic_rows].subs(sample);piv=list(Ms.to_DM().rref()[1]);br=list(Ms[:,piv].T.to_DM().rref()[1]);assert len(piv)==6
    for r in br:
        coeff=[pol(v) for v in M[r,:]]
        for j in range(45):assert sum((c*J[k][j] for k,c in enumerate(coeff) if c),R.zero)==0
    minor=Ms.extract(br,piv).det();assert minor!=0
    paths=[Path(__file__),data_path,lift_path,root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py',root/'code/elliptic_weight1_dominance.py',root/'reports/S61_BINARY_UNIFORM_RANK_20261101.json',root/'reports/S61_BINARY_UNIFORM_QUADRATIC_20261101.json',root/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    paths += [root/'legacy/v02'/v for v in ['code/reduce_symplectic_bianchi.py','code/build_symplectic_slice.py','code/build_universal_slice.py','data/symplectic_slice_certificate.json','data/symplectic_slice_equations.json','data/symplectic_bianchi_identities.json']]
    return dict(id='S61_BINARY_UNIFORM_CUBIC_20261101',all_passed=True,field='QQ(nine dominant slice coordinates)',normal_tangent_columns=normal,parameter_columns=parameters,common_quadratic_radical_dimension_at_sample=3,parameter_projection_minor=str(parameter_minor),quadratic_rank=5,lift_support_columns=support,annihilator_rank=16,supported_jacobian_rank=19,lift_existence_identities_checked=True,potential_cubic_rows=cubic_rows,other109_zero_quadratic_cubic_rows_zero=True,full_raw_J4_coefficients_verified=80,full_raw_J5_coefficients_verified=140,lower_affine_Jacobi_coefficients_identically_zero=lower_count,selected_Bianchi_labels=[labels[r] for r in br],selected_Bianchi_multipliers=[[str(v) for v in M[r,:]] for r in br],Bianchi_minor_columns=piv,Bianchi_minor=str(minor),independent_cubic_relations_mod_linear_quadrics=6,cubic_relation_rank_upper_bound=50,normal_Hilbert_function_lower_bound_through_degree3=[1,7,23,34],normal_length_lower_bound=65,proof_note='research/notes/S61_BINARY_UNIFORM_CUBIC_20261101.md',limitations=['A generic normal cubic rank and length lower bound; use the proof note for parameter and formal-ring quantifiers.','No all-beta or full affine-fiber exhaustion assertion.'],source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();report=verify()
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({k:report[k] for k in ['all_passed','cubic_relation_rank_upper_bound','normal_length_lower_bound','seconds']}))


if __name__=='__main__':main()
