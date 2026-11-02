#!/usr/bin/env python3
"""Exact no-divisor rank-drop certificate for six weighted quadrics."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations,combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
root=Path(__file__).resolve().parents[1]
NOTE=root/'research/notes/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.md'
def verify():
    begin=time.time()
    data=json.loads((root/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json').read_text())
    old=json.loads((root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json').read_text())['linear']
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();pol=lambda x:R.from_expr(s.sympify(x));frac=lambda x:K.from_expr(s.sympify(x))
    T=[[frac(row[j]) for j in data['normal34']] for row in data['T34']]
    U=[[frac(v) for v in row] for row in data['U34']]
    td=[]
    for j in range(8):
        den=R.one
        for row in T:den=den.lcm(row[j].denom)
        td.append(den)
    ud=R.one
    for row in U:
        for v in row:ud=ud.lcm(v.denom)
    Tn=[[v.numer*td[j].exquo(v.denom) for j,v in enumerate(row)] for row in T]
    Un=[[v.numer*ud.exquo(v.denom) for v in row] for row in U]
    m=S61();beta={int(k[1:])-1:s.sympify(v) for k,v in data['beta'].items()};two,F=specialize(m,beta,True)
    fixed=old['fixed'];cols=old['cols'];J=[[R.zero for j in range(45)] for i in range(135)];H=[]
    for i,f in enumerate(F):
        row=[]
        for mon,c in f.items():
            cc=pol(c)
            if len(mon)==2:
                x,y=mon
                if x!=fixed and y==fixed:J[i][cols.index(x)]+=cc
                if y!=fixed and x==fixed:J[i][cols.index(y)]+=cc
            elif len(mon)==1 and mon[0]==fixed:J[i][44]+=cc
            if fixed not in mon:
                ij=[cols.index(j) for j in mon]+[44]*(2-len(mon))
                if cc:row.append((ij,cc))
        H.append(row)

    den=ud
    start=time.time();zero=R.zero;light=[0,1,2,3,4,6,7];monoms=list(combinations_with_replacement(range(7),2));lookup={m:i for i,m in enumerate(monoms)};cache={}
    def prod(i,j):
     key=tuple(sorted((i,j)))
     if key not in cache:
      row=[zero]*28
      for k,kk in enumerate(light):
       for l,ll in enumerate(light):
        if Tn[i][kk] and Tn[j][ll]:row[lookup[tuple(sorted((k,l)))]]+=Tn[i][kk]*Tn[j][ll]
      cache[key]=row
     return cache[key]
    Q=[]
    for terms in H:
     row=[zero]*28
     for (i,j),c in terms:
      for n,p in enumerate(prod(i,j)):
       if p:row[n]+=c*p
     Q.append(row)
    W=[[-sum((v*Q[data['selected34_rows'][j]][m] for j,v in enumerate(row) if v),zero) for m in range(28)] for row in Un]
    A=[[den*Q[i][m]+sum((v*W[j][m] for j,v in enumerate(row) if v),zero) for m in range(28)] for i,row in enumerate(J)]
    rows=[];known0=[a,500*a+90*e**2+3*e-24,50*a*d+18*d*e**2-9*d*e-2*d+60*e-40]
    for i,row in enumerate(A):
     if not any(row):continue
     assert data['equation_weights'][i]==2
     g=zero
     for p in row:g=g.gcd(p)
     rem=g
     for f in known0:
      while rem.rem(f)==0:rem=rem.exquo(f)
     assert max(map(sum,rem))==0,(i,rem)
     rr=[p.exquo(g) for p in row];lc=next(p.LC for p in rr if p);rr=[p.quo_ground(lc) for p in rr]
     encoded={}
     for mon,p in zip(monoms,rr):
      if not p:continue
      power=[0]*8
      for j in mon:power[light[j]]+=1
      encoded[tuple(power)]=p
     rows.append((i,2,encoded))
    assert len(rows)==6
    print('DIRECT_QUADRATICS',len(rows),sum(len(p) for _,_,r in rows for p in r.values()),'seconds',time.time()-start,flush=True)
    K=R.to_field();known=[a,500*a+90*e**2+3*e-24,50*a*d+18*d*e**2-9*d*e-2*d+60*e-40,100*a+18*e**2-9*e-2]
    def unit(poly):
        if not poly:return False
        for f in known:
            while poly.rem(f)==0:poly=poly.exquo(f)
        return max(map(sum,poly))==0

    from itertools import combinations
    selected=[(i,r) for i,w,r in rows if w==2]
    mons=sorted(set(m for _,r in selected for m in r))
    M=[[K(row.get(m,R.zero)) for m in mons] for _,row in selected]
    piv=[]
    for rr in range(3):
     candidates=sorted((len(v.numer)+len(v.denom),i,j) for i in range(rr,len(M)) for j,v in enumerate(M[i]) if j not in [q[1] for q in piv] and v)
     _,i,j=next((n,i,j) for n,i,j in candidates if unit(M[i][j].numer) and unit(M[i][j].denom))
     M[rr],M[i]=M[i],M[rr];v=M[rr][j];piv.append((i,j,str(v)));M[rr]=[q/v for q in M[rr]]
     for i in range(rr+1,len(M)):
      v=M[i][j]
      if v:M[i]=[x-v*y for x,y in zip(M[i],M[rr])]
    remaining=[j for j in range(len(mons)) if j not in [q[1] for q in piv]]
    N=[]
    for row in M[3:]:
     den=R.one
     for j in remaining:den=den.lcm(row[j].denom)
     assert unit(den)
     N.append([row[j].numer*den.exquo(row[j].denom) for j in remaining])
    g=R.zero;record=[]
    for comb in combinations(range(len(remaining)),3):
     x,y,z=comb
     f=N[0][x]*(N[1][y]*N[2][z]-N[1][z]*N[2][y])-N[0][y]*(N[1][x]*N[2][z]-N[1][z]*N[2][x])+N[0][z]*(N[1][x]*N[2][y]-N[1][y]*N[2][x])
     if not f:continue
     for k in known:
      while f.rem(k)==0:f=f.exquo(k)
     old=g;g=g.gcd(f)
     if g!=old:
      record.append({'columns':[remaining[j] for j in comb],'polynomial':str(f.as_expr()),'gcd':str(g.as_expr())})
      print('GCD',len(record),'terms',len(g),'degree',max(map(sum,g)),'seconds',time.time()-start,flush=True)
     if max(map(sum,g))==0:break
    assert max(map(sum,g))==0
    paths=[Path(__file__),NOTE,root/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json',root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json',root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py']
    out={'id':'S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102','all_passed':True,
     'quadratic_rows':[i for i,_,_ in rows],'quadratic_row_count':6,'pure_light_variables':7,
     'quadratic_monomial_count':28,'scalar_row_contents_are_a_F_B_units':True,
     'unit_pivots':piv,'stripped_residual_three_minors':record,'final_gcd':str(g.as_expr()),
     'bad_locus_dimension_upper_on_three_core_open':1,
     'quadratic_rows_sparse':[[[list(mon),[[list(mm),str(c)] for mm,c in poly.items()]] for mon,poly in row.items()] for _,_,row in rows],
     'proof_note':str(NOTE.relative_to(root)),
     'limitations':['Only nuisance-zero on a*F*B*C nonzero.','No claim of rank6 at every point, nor cubic or quartic rank coverage.','No new ambient saturation removal.'],
     'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'seconds':time.time()-begin}
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert not args.output.exists();result=verify()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','bad_locus_dimension_upper_on_three_core_open','seconds']}))
