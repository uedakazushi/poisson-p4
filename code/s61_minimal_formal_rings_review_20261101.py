#!/usr/bin/env python3
"""Independent reconstruction of the minimal-infinity local and rank certificates."""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J
ROOT=Path(__file__).resolve().parents[1]
def mat(rows):return s.Matrix(rows).applyfunc(s.sympify)
def grouped(encoded,weights,beta,minimum=0):
    two=[i for i,w in enumerate(weights) if w==2];positions={a:i for i,a in enumerate(two)};cache={():s.S.One};result=[]
    for equation in encoded:
        row={}
        for indices,n,d in equation:
            q=tuple(positions[i] for i in indices if weights[i]==2)
            if len(q)<minimum:continue
            b=tuple(i for i in indices if weights[i]==1)
            if b not in cache:cache[b]=s.prod(beta[i] for i in b)
            row[q]=row.get(q,0)+s.Rational(n,d)*cache[b]
        result.append({key:s.cancel(value) for key,value in row.items() if value!=0})
    return two,result

def derivative(polys,q,columns):
    A=s.zeros(len(polys),45);t=s.zeros(len(polys),1)
    for row,f in enumerate(polys):
        top=0
        for mon,v in f.items():
            if len(mon)==2:
                i,j=mon;top+=v*q[i]*q[j];A[row,i]+=v*q[j];A[row,j]+=v*q[i]
            elif len(mon)==1:t[row,0]+=v*q[mon[0]]
        assert s.cancel(top)==0
    return A[:,columns].applyfunc(s.cancel),t.applyfunc(s.cancel)

def run():
    start=time.time();model=S61();weights=model.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1];two=[i for i,w in enumerate(weights) if w==2]
    extra_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';equations=model.reduced['equations']+[json.loads(extra_path.read_text())['universal_new_equation']]
    local_path=ROOT/'reports/S61_MINIMAL_INFINITY_LOCAL_DATA_20261101.json';universal_path=ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json';local=json.loads(local_path.read_text());universal=json.loads(universal_path.read_text());point_checks={}
    for kind in ['M0','M1']:
        point_checks[kind]={}
        for count in [134,135]:
            rec=local['cases'][kind][str(count)];point=list(map(s.Rational,rec['point61']));q=s.Matrix([point[i] for i in two]);columns=rec['chart_columns'];_,polys=grouped(equations[:count],weights,{i:point[i] for i in one});A,tc=derivative(polys,q,columns);full=A.row_join(tc)
            assert full==mat(rec['jacobian']);T=mat(rec['tangent_basis']);C=mat(rec['cokernel_rows']);assert full*T==s.zeros(count,T.cols) and C*full==s.zeros(C.rows,45)
            assert T.rank()==45-full.rank()==T.cols
            n=T.cols;mons=list(combinations_with_replacement(range(n),2));forms={qindex:list(T[i,:]) for i,qindex in enumerate(columns)};forms[next(i for i in range(45) if i not in columns)]=[s.S.Zero]*n;tf=list(T[-1,:]);H=s.zeros(count,len(mons));mi={a:i for i,a in enumerate(mons)}
            for row,f in enumerate(polys):
                for mon,v in f.items():
                    first,second=[forms[i] for i in mon]+[tf]*(2-len(mon))
                    for i in range(n):
                        for j in range(n):H[row,mi[tuple(sorted((i,j)))]]+=v*first[i]*second[j]
            obstruction=C*H;assert obstruction==mat(rec['quadratic_obstruction_basis']);r=n-1;assert obstruction.rank()==r*(r+1)//2
            Hess=s.Matrix.vstack(*[s.Matrix(n,n,lambda i,j:obstruction[row,mi[tuple(sorted((i,j)))]]*(2 if i==j else 1)) for row in range(obstruction.rows)]);assert Hess.rank()==r
            radical=mat(rec['quadratic_common_radical']);assert radical.rank()==1 and Hess*radical==s.zeros(Hess.rows,1) and T[-1,:]*radical==s.zeros(1,1)
            point_checks[kind][str(count)]={'jacobian_rank':full.rank(),'tangent_dimension':n,'transverse_length_upper_bound':n,'obstruction_rank':obstruction.rank(),'Hessian_radical_dimension':1}
    c=s.Symbol('c');bb=s.symbols('b0:9');universal_checks={};x,y,z,w=X[1:]
    for kind in ['M0','M1']:
        rec=universal['cases'][kind];Q=s.zeros(5);Q[1,2]=x*x;Q[2,3]=-c*x*z;Q[2,4]=c*x*w
        if kind=='M1':Q[3,4]=x*x;Q[2,3]-=x*z;Q[2,4]-=x*w
        Q-=Q.T;point=model.normalize(X[0]**2*J+Q)['parameters'];q=s.Matrix([point[i] for i in two]);assert q==mat(rec['canonical_q_coordinates'])
        BK=mat(rec['beta_kernel']);bv=BK*s.Matrix(bb);_,polys=grouped(equations,weights,dict(zip(one,bv,strict=True)),minimum=1);fixed=next(i for i,a in enumerate(q) if a);columns=[i for i in range(45) if i!=fixed];A,tc=derivative(polys,q,columns);A0=A[:70,:]
        K=mat(rec['first_kernel']);assert (A0*K).applyfunc(s.cancel)==s.zeros(70,K.cols);rr=rec['first_pivot_rows'];cc=rec['first_pivot_columns'];rank=rec['first_rank'];assert len(rr)==len(cc)==rank
        minor=A0.extract(rr,cc);det=minor.to_DM().det().as_expr();assert s.cancel(det-s.sympify(rec['first_pivot_determinant']))==0 and det!=0
        assert K.subs(c,s.Rational(2,7)).rank()==44-rank
        U=s.zeros(44,rank);iv=minor.inv()
        for j,i in enumerate(cc):U[i,:]=iv[j,:]
        assert (tc[:70,:]-A0*U*tc[rr,:]).applyfunc(s.cancel)==s.zeros(70,1)
        residual=(A[70:,:]*K).applyfunc(s.cancel);res_t=(tc[70:,:]-A[70:,:]*U*tc[rr,:]).applyfunc(s.cancel)
        assert residual==mat(rec['residual_matrix']) and res_t==mat(rec['residual_t_column'])
        universal_checks[kind]={}
        for count in [134,135]:
            info=rec['rank_certificates'][str(count)];F=residual[:count-70,:].row_join(res_t[:count-70,:]);N=mat(info['residual_kernel']);free=info['free_columns'];assert N[free,:]==s.eye(len(free));assert (F.to_DM()*N.to_DM()).is_zero_matrix
            sample={c:s.Rational(2,7),**dict(zip(bb,range(1,10),strict=True))};full_rank=rank+F.subs(sample).rank();leading_rank=rank+residual[:count-70,:].subs(sample).rank()
            assert full_rank==info['full_rank'] and leading_rank==info['leading_rank'];assert (N[-1,:]==s.zeros(1,N.cols))==info['homogenizer_zero_on_tangent']
            universal_checks[kind][str(count)]={'full_rank':full_rank,'leading_rank':leading_rank,'saved_kernel_independent_and_annihilated':True,'homogenizer_zero_on_tangent':info['homogenizer_zero_on_tangent']}
        print(json.dumps({'reviewed':kind,'universal_rank_identities':True}),flush=True)
    paths=[Path(__file__),ROOT/'code/s61_model.py',model.source/'symplectic_slice_reduced_equations.json',local_path,universal_path,extra_path,ROOT/'research/notes/S61_MINIMAL_INFINITY_FORMAL_RINGS_20261101.md']
    return {'id':'S61_MINIMAL_FORMAL_RINGS_REVIEW_20261101','all_passed':True,'point_checks':point_checks,'universal_checks':universal_checks,'proof_review':'Finite-module and generic length argument accepted; square-zero ideal and t assertions follow by torsion-freeness. Source-to-canonical M1 scaling convention explicitly checked in the corrected proof note.','proof_note':'research/notes/S61_MINIMAL_FORMAL_RINGS_REVIEW_20261101.md','seconds':time.time()-start,'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=run()
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
if __name__=='__main__':main()
