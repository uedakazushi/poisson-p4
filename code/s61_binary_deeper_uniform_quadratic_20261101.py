#!/usr/bin/env python3
"""Universal quadratic rank24 and proper-image rank drop in the deeper chart."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import numpy as np
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
PIVOT_ROWS=[28,38,19,68,60,18,61,103,131,39,40,20,133,23,84,86,115,129,117,122,114,85,99,120]
PIVOT_COLS=[37,42,23,43,24,22,28,29,35,40,15,26,18,7,25,38,20,31,27,5,11,36,13,33]
FACTOR_POINTS=[(8727,1),(0,31998),(0,17735),(0,217),(1,9757),(16629,1),(1,16899),(0,29224)]


def build_quadratics(data):
    start=time.time()
    R,*gens=ring(','.join(data['variables']),s.QQ);FF=R.to_field();pol=lambda x:R.from_expr(s.sympify(x));frac=lambda x:FF.from_expr(s.sympify(x));z=R.zero
    T=[[frac(v) for v in row] for row in data['T']];U=[[frac(v) for v in row] for row in data['U']]
    Dt=R.one;Du=R.one
    for row in T:
     for v in row:Dt=Dt.lcm(v.denom)
    for row in U:
     for v in row:Du=Du.lcm(v.denom)
    print('den T',s.factor(Dt.as_expr()),'den U',s.factor(Du.as_expr()),time.time()-start,flush=True)
    normal=data['normal']
    dens=[]
    for j in normal:
     dd=R.one
     for row in T:dd=dd.lcm(row[j].denom)
     dens.append(dd)
    Tn=[[row[j].numer*dens[k].exquo(row[j].denom) for k,j in enumerate(normal)] for row in T]
    Un=[[v.numer*Du.exquo(v.denom) for v in row] for row in U]
    assert all(v==0 for v in Tn[-1]);JJ=[[pol(v) for v in row] for row in data['J']];rows=data['pivot_rows'];fixed=data['fixed'];cols=data['cols']
    m=S61();beta={int(k[1:])-1:s.sympify(v) for k,v in data['beta'].items()};two,F=specialize(m,beta,True)
    ms=list(combinations_with_replacement(range(9),2));ix={mon:i for i,mon in enumerate(ms)};cache={}
    def prod(i,j):
     key=tuple(sorted((i,j)))
     if key not in cache:
      out=[z]*45
      for a,v in enumerate(Tn[i]):
       if not v:continue
       for b,w in enumerate(Tn[j]):
        if w:out[ix[tuple(sorted((a,b)))]]+=v*w
      cache[key]=out
     return cache[key]
    Q=[]
    for f in F:
     out=[z]*45
     for mon,c in f.items():
      if fixed in mon:continue
      ij=[cols.index(v) for v in mon]+[44]*(2-len(mon));i,j=ij
      if i==44 or j==44:continue
      c=pol(c)
      for a,v in enumerate(prod(i,j)):
       if v:out[a]+=c*v
     Q.append(out)
    print('Q built',sum(len(v) for row in Q for v in row),time.time()-start,flush=True)
    W=[[-sum((v*Q[rows[k]][j] for k,v in enumerate(row) if v),z) for j in range(45)] for row in Un]
    A=[[Du*Q[i][j]+sum((v*W[k][j] for k,v in enumerate(row) if v),z) for j in range(45)] for i,row in enumerate(JJ)]
    nonzero=[i for i,row in enumerate(A) if any(row)]
    print('A built',nonzero,sum(len(v) for row in A for v in row),time.time()-start,flush=True)
    Ag=[];content=[]
    for row in A:
     g=z
     for v in row:g=g.gcd(v)
     content.append(g)
     Ag.append([v.exquo(g) if g else z for v in row])
    print('content cleared',sum(len(v) for row in Ag for v in row),time.time()-start,flush=True)
    return R,dict(Tden=dens,Uden=Du,Tn=Tn,Un=Un,Q=Q,W=W,A=A,Ag=Ag,contents=content,nonzero=nonzero)


def elimination_rank(M,p):
    M=M.copy();rr=0;piv=[]
    for col in range(M.shape[1]):
        ids=np.flatnonzero(M[rr:,col])
        if len(ids)==0:continue
        i=rr+int(ids[0]);M[[rr,i]]=M[[i,rr]];M[rr]=M[rr]*pow(int(M[rr,col]),-1,p)%p
        ids=np.flatnonzero(M[rr+1:,col])+rr+1
        if len(ids):M[ids]=(M[ids]-M[ids,col,None]*M[rr,None,:])%p
        piv.append(col);rr+=1
        if rr==M.shape[0]:break
    return rr,piv


def determinant_mod(M,p):
    M=M.copy();det=1
    for col in range(len(M)):
        ids=np.flatnonzero(M[col:,col]);assert len(ids)
        i=col+int(ids[0])
        if i!=col:M[[col,i]]=M[[i,col]];det=-det
        v=int(M[col,col]);det=det*v%p;M[col]=M[col]*pow(v,-1,p)%p
        ids=np.flatnonzero(M[col+1:,col])+col+1
        if len(ids):M[ids]=(M[ids]-M[ids,col,None]*M[col,None,:])%p
    return det%p


def verify():
    start=time.time();path=ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json'
    old=json.loads(path.read_text());assert old['all_passed'];data=old['linear']
    R,quad=build_quadratics(data);e,d=s.symbols('e d')
    K=s.QQ.frac_field('e','d');R5,*v=ring('p,q,r,h,k',K);gm=(0,)*5
    def convert(poly):
        ans={}
        for mon,c in poly.items():
            term=K.convert(c,R.domain)*K.gens[0]**mon[0]*K.gens[1]**mon[1]
            key=mon[2:];ans[key]=ans.get(key,K.zero)+term
        return R5.from_dict(ans)
    roworder=PIVOT_ROWS+[i for i in quad['nonzero'] if i not in PIVOT_ROWS]
    M=[[convert(v) for v in quad['Ag'][i]] for i in roworder]
    det=K.one;pivots=[]
    for rr,j in enumerate(PIVOT_COLS):
        pivot=M[rr][j];assert len(pivot)==1 and gm in pivot
        c=pivot[gm];assert c;det*=c;pivots.append(str(K.to_sympy(c)))
        M[rr]=[v.quo_ground(c) for v in M[rr]]
        for i in range(rr+1,len(M)):
            h=M[i][j]
            if h:M[i]=[x-h*y for x,y in zip(M[i],M[rr])]
        print('quadratic pivot',rr+1,'seconds',time.time()-start,flush=True)
    assert all(not v for row in M[24:] for v in row)
    # Every row content is a polynomial in the two core variables. The minor
    # of the original A differs by exactly the selected24 row contents.
    fac=s.factor_list(K.to_sympy(det))[1]
    for i in PIVOT_ROWS:
        value=quad['contents'][i].as_expr();assert value.free_symbols<={e,d}
        fac+=s.factor_list(value)[1]
    factors={str(s.Poly(f,e,d).monic().as_expr()):f for f,power in fac}
    known=[2*e-1,10*e-1,15*e+8,24*e-7,90*d*e*e-93*d*e+4*d+600*e-400]
    mask=s.prod(known)*(33*e-4);p=32003
    def mod(c):
        c=s.Rational(c);return int(c.p)%p*pow(int(c.q),-1,p)%p
    matrices={};matrix_minors={}
    for ev,dv in FACTOR_POINTS:
        assert int(s.Poly(mask,e,d).eval({e:ev,d:dv}))%p
        array=[]
        for row in quad['A']:
            values=[]
            for poly in row:
                val=sum(mod(c)*pow(ev,mon[0],p)*pow(dv,mon[1],p)
                    for mon,c in poly.items() if not any(mon[2:]))%p
                values.append(val)
            array.append(values)
        mat=np.array(array,dtype=np.int64);rk,cols=elimination_rank(mat,p);assert rk==24
        _,rows=elimination_rank(mat[:,cols].T,p);assert len(rows)==24
        determinant=determinant_mod(mat[np.ix_(rows,cols)],p);assert determinant
        matrices[(ev,dv)]=mat
        matrix_minors[(ev,dv)]=dict(sample_core=[ev,dv],sample_other_five=[0]*5,
            rows=rows,columns=cols,determinant_mod_prime=determinant)
    tests=[];excluded=[]
    for f in factors.values():
        if any(s.rem(f,g,e,d)==0 for g in known):
            excluded.append(str(f));continue
        poly=s.Poly(f,e,d)
        matches=[]
        for ev,dv in FACTOR_POINTS:
            value=sum(mod(c)*pow(ev,mon[0],p)*pow(dv,mon[1],p) for mon,c in poly.terms())%p
            if value==0:matches.append((ev,dv))
        assert matches,str(f)
        ev,dv=matches[0];tests.append(dict(factor=str(f),total_degree=poly.total_degree(),
            factor_vanishes_mod_prime=True,**matrix_minors[(ev,dv)]))
    assert len(factors)==13 and len(tests)==8 and len(excluded)==5
    paths=[Path(__file__),path,ROOT/'code/s61_model.py',
        ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
        ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',
        ROOT/'research/notes/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.md']
    return dict(id='S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101',all_passed=True,
        coefficient_ring='QQ(e,d)[p,q,r,h,k]',normal_variable_count=9,
        universal_quadratic_rank_upper_bound=24,generic_quadratic_rank=24,
        gaussian_pivot_rows=PIVOT_ROWS,gaussian_pivot_columns=PIVOT_COLS,gaussian_pivots=pivots,
        all_remaining_rows_zero=True,
        T_column_denominators=[str(v.as_expr()) for v in quad['Tden']],
        U_denominator=str(quad['Uden'].as_expr()),
        row_contents=[str(v.as_expr()) for v in quad['contents']],
        primitive_quadratic_rows={str(i):[str(v.as_expr()) for v in quad['Ag'][i]] for i in quad['nonzero']},
        irreducible_minor_factors=list(map(str,factors.values())),
        previously_removed_factors=excluded,prime=p,factor_nondivisibility_witnesses=tests,
        residual_rank_drop_slice_codimension_lower_bound=2,
        residual_rank_drop_slice_dimension_upper_bound=5,
        residual_rank_drop_Sp4_saturation_dimension_upper_bound=15,
        all_remaining_general_beta_deeper_points_quadratic_rank_24=True,
        proof_note='research/notes/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.md',
        limitations=['This proves only the universal quadratic rank on the remaining deeper chart.',
            'It does not prove the higher-order colength65 at every point.',
            'Full infinity saturation and affine five-point exhaustion remain open.'],
        input_sha256={str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest() for q in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();data=verify()
    with args.output.open('x') as stream:json.dump(data,stream,indent=2);stream.write('\n')
    print(json.dumps({k:data[k] for k in ['all_passed','generic_quadratic_rank','seconds']}))
