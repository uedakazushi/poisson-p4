#!/usr/bin/env python3
"""Exact flat4 quotient at one free-M0 flag endpoint; no ambient removal claim."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_baseline_thickenings_20261101 import coords,enc
from weight1_structure import minor_witness
ROOT=Path(__file__).resolve().parents[1]
def run():
    st=time.time();m=S61();path=ROOT/'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json';r=json.loads(path.read_text())['strata']['free_M0']['local_data'];pt=s.Matrix(r['point61']).applyfunc(s.Rational);one=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==1];two,F=specialize(m,dict(zip(one,pt[one,:])),True);q0=pt[two,:];co=r['chart_columns'];fi=two.index(r['chart_fixed_u']-1);x,y,z,w=X[1:];G=-18*x*x*w-2*x*z*w;Q=s.zeros(4);Q[0,1]=x*x;Q[1,2]=-s.diff(G,w);Q[1,3]=s.diff(G,z);Q-=Q.T;dq=s.Matrix(coords(m,Q))[two,:];dq-=q0*dq[fi]/q0[fi];ss=s.Symbol('s');q=q0+ss*dq;M=s.zeros(135,45)
    for row,f in enumerate(F):
        val=0
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon;val+=c*q[i]*q[j]
                if i!=fi:M[row,co.index(i)]+=c*q[j]
                if j!=fi:M[row,co.index(j)]+=c*q[i]
            elif len(mon)==1:M[row,44]+=c*q[mon[0]]
        assert s.expand(val)==0
    N=M.to_DM().nullspace().to_Matrix().T.applyfunc(s.factor);assert N.cols==4;base=dq[co,:].col_join(s.zeros(1,1));C=base.row_join(N);piv=list(C.to_DM().rref()[1]);assert piv[0]==0 and len(piv)==4;N=C[:,piv]
    def valuation(v):
        if not v:return 1000
        a,b=s.fraction(s.cancel(v));return min(i[0]for i,c in s.Poly(a,ss).terms())-min(i[0]for i,c in s.Poly(b,ss).terms())
    for step in range(30):
        for j in range(1,4):
            v=min(valuation(x)for x in N[:,j]);N[:,j]=N[:,j].applyfunc(lambda x:s.cancel(x/ss**v))
        N0=N.subs(ss,0)
        if N0.rank()==4:break
        rel=N0.nullspace()[0];p=next(j for j in range(1,4)if rel[j]);N[:,p]=(N*rel/ss).applyfunc(s.cancel)
    else:raise AssertionError('lattice iteration bound')
    for j in range(1,4):
        content=s.factor(s.gcd_list(list(N[:,j])));assert content.subs(ss,0)!=0;N[:,j]=N[:,j].applyfunc(lambda v:s.cancel(v/content));row=next(i for i,v in enumerate(N[:,j])if v.subs(ss,0)!=0);scale=N[row,j].subs(ss,0);N[:,j]=N[:,j].applyfunc(lambda v:s.factor(v/scale))
    N0=N.subs(ss,0);assert(M.to_DM()*N.to_DM()).is_zero_matrix and N[:,0]==base;assert N0.rank()==4;witness=minor_witness(N0);assert witness['rank']==4 and all(valuation(v)>=0 for v in N);assert N[-1,0]==N[-1,1]==N[-1,2]==0 and valuation(N[-1,3])==3
    T=s.Matrix(r['tangent_basis']).applyfunc(s.Rational);rr=list(T.T.to_DM().rref()[1]);V=T[rr,:].inv()*N0[rr,:];assert T*V==N0
    paths=[Path(__file__),path,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_minimal_baseline_thickenings_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_MINIMAL_FREE_FLAG_FLAT4_PROBE_20261101',all_passed=True,field='QQ(s)',point61=list(map(str,pt)),q0=list(map(str,q0)),delta=list(map(str,dq)),chart_columns=co,jacobian=enc(M),generic_full_rank=41,regular_kernel=enc(N),endpoint_tangent=enc(V),embedding_minor=witness,actual_quotient='QQ[[s,e1,e2,e3]]/(ei*ej:1<=i<=j<=3)',homogenizer_normal_coefficients=list(map(str,N[-1,1:])),homogenizer_order_on_curve=3,scope=['One exact fixed-beta quotient of the full135 completed ring.','t is nonzero in this quotient and t^2=0 there.','No t-nilpotence or saturation-removal assertion for the ambient ring, and no relative H11 theorem.'],proof_note='research/notes/S61_MINIMAL_FREE_FLAG_FLAT4_PROBE_20261101.md',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.time()-st)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
