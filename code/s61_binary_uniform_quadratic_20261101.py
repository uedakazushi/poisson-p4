#!/usr/bin/env python3
"""Universal binary quadratic obstruction rank on a dominant nine-dimensional slice."""
import sys,json,time,hashlib,argparse
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring

from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

def verify():
    start=time.time();root=Path(__file__).resolve().parents[1];data_path=root/'reports/S61_BINARY_QUADRATIC_PROJECTION_DATA_20261101.json';data=json.loads(data_path.read_text());R,*gens=ring(','.join(data['variables']),s.QQ);FF=R.to_field();symbols=s.symbols(' '.join(data['variables']));locals=dict(zip(data['variables'],symbols))
    def pol(v):return R.from_expr(s.sympify(v,locals=locals))
    def frac(v):return FF.from_expr(s.sympify(v,locals=locals))
    project=[[frac(v) for v in row] for row in data['projection']];rows=data['jacobian_rows'];Tp=[[frac(v) for v in row] for row in data['tangent']]
    for j in range(10):
     den=R.one
     for i in range(45):den=den.lcm(Tp[i][j].denom)
     for i in range(45):Tp[i][j]=Tp[i][j].numer*den.exquo(Tp[i][j].denom)

    m=S61();u=s.symbols('u1:62');fixedbeta={s.Symbol(k):s.Rational(v) for k,v in data['fixed_beta'].items()};beta={i:u[i].subs(fixedbeta) for i,w in enumerate(m.reduced['parameter_weights']) if w==1};beta[11]=beta[12]=s.S.Zero;two,F=specialize(m,beta,True);fixed=data['chart_fixed_q'];cols=[i for i in range(45) if i!=fixed];q0=[s.S.Zero]*45;q0[fixed]=1
    J=[[R.zero]*45 for _ in range(135)];ht=[]
    for row,f in enumerate(F):
     terms=[]
     for mon,c in f.items():
      c=pol(c)
      if len(mon)==2:
       i,j=mon
       if i!=fixed and q0[j]:J[row][cols.index(i)]+=c
       if j!=fixed and q0[i]:J[row][cols.index(j)]+=c
      elif len(mon)==1 and q0[mon[0]]:J[row][-1]+=c
      if fixed not in mon:terms.append(tuple([cols.index(k) for k in mon]+[44]*(2-len(mon)))+(c,))
     ht.append(terms)
    project_num=[];project_den=[]
    for r in range(135):
     den=R.one
     for v in project[r]:den=den.lcm(v.denom)
     nums=[v.numer*den.exquo(v.denom) for v in project[r]];project_num.append(nums);project_den.append(den)
     for j in range(45):assert den*J[r][j]-sum((v*J[rows[k]][j] for k,v in enumerate(nums) if v),R.zero)==0,(r,j)
    for i,r in enumerate(rows):assert all(project[r][j]==FF(int(i==j)) for j in range(35))

    for i in range(135):
     for j in range(10):assert sum((a*Tp[k][j] for k,a in enumerate(J[i]) if a),R.zero)==0

    assert all(v == 0 for v in Tp[44])  # The tangent homogenizing coordinate is zero.
    mons=list(combinations_with_replacement(range(10),2));mi={m:i for i,m in enumerate(mons)};cache={}
    def prod(i,j):
     key=tuple(sorted((i,j)))
     if key not in cache:
      result=[R.zero]*55
      for a,v in enumerate(Tp[i]):
       if not v:continue
       for b,w in enumerate(Tp[j]):
        if w:result[mi[tuple(sorted((a,b)))]]+=v*w
      cache[key]=result
     return cache[key]
    Q=[]
    for row,terms in enumerate(ht):
     q=[R.zero]*55
     for i,j,c in terms:
      if i==44 or j==44:continue
      for n,v in enumerate(prod(i,j)):
       if v:q[n]+=c*v
     Q.append(q)

    A=[];nz=[]
    for r in range(135):
     row=[]
     for j in range(55):row.append(project_den[r]*Q[r][j]-sum((v*Q[rows[k]][j] for k,v in enumerate(project_num[r]) if v),R.zero))
     A.append(row)
     if any(row):nz.append(r)
    assert nz==[19,60,117,123,133]
    sample={u[i]:s.Rational(j+2) for j,i in enumerate([i for i,w in enumerate(m.reduced['parameter_weights']) if w==1])};sample.update(fixedbeta);vals=[s.QQ.from_sympy(sample[v]) for v in symbols]
    def ev(f):return sum((c*__import__("math").prod(v**e for v,e in zip(vals,mon)) for mon,c in f.items()),s.QQ.zero)
    Am=s.Matrix([[s.QQ.to_sympy(ev(f)) for f in A[i]] for i in nz]);pc=list(Am.to_DM().rref()[1]);assert len(pc)==5;det=Am[:,pc].det(method='domain-ge');assert det


    Js=s.Matrix([[s.QQ.to_sympy(ev(f)) for f in J[i]] for i in rows]);assert Js.rank()==35
    Ts=s.Matrix([[s.QQ.to_sympy(ev(f)) for f in row] for row in Tp]);assert Ts.rank()==10
    paths=[Path(__file__),data_path,root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py',root/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',root/'reports/S61_BINARY_UNIFORM_RANK_20261101.json']
    return dict(id='S61_BINARY_UNIFORM_QUADRATIC_20261101',all_passed=True,field='QQ(nine dominant slice coordinates)',projection_identities_checked=True,tangent_kernel_identities_checked=True,full_jacobian_rank=35,tangent_dimension=10,universal_nonzero_residual_rows=nz,universal_quadratic_obstruction_rank=5,sample_minor_columns=pc,sample_minor=str(det),sample={str(k):str(v) for k,v in sample.items()},proof_note='research/notes/S61_BINARY_UNIFORM_QUADRATIC_20261101.md',limitations=['Does not yet prove a universal cubic rank bound or transverse multiplicity65','No claim of t-nilpotence in the unsliced completed local ring'],source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();r=verify()
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({k:r[k] for k in ['all_passed','universal_quadratic_obstruction_rank','seconds']}))


if __name__=='__main__':main()
