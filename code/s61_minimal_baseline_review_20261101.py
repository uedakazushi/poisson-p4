#!/usr/bin/env python3
"""Independent full-kernel and Hessian audit of the four baseline thickenings."""
import argparse,hashlib,json
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json'


def run():
    r=json.loads(INPUT.read_text());m=S61();matrix=lambda a:s.Matrix(a).applyfunc(s.sympify)
    one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2]
    index={i:j for j,i in enumerate(two)}
    source=m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    def equations(beta):
        result=[]
        for eq in source:
            terms={}
            for mon,n,d in eq:
                qmon=tuple(index[i] for i in mon if i in index)
                value=s.Rational(n,d)*s.prod(beta[one.index(i)] for i in mon if i in one)
                terms.setdefault(qmon,[]).append(value)
            result.append({mon:s.expand(s.Add(*values)) for mon,values in terms.items()})
        return result
    def jacobian(F,q,cols):
        A=s.zeros(135,45);tc=s.zeros(135,1)
        for row,f in enumerate(F):
            value=s.S.Zero
            for mon,c in f.items():
                if len(mon)==2:
                    i,j=mon;value+=c*q[i]*q[j];A[row,i]+=c*q[j];A[row,j]+=c*q[i]
                elif len(mon)==1:tc[row]+=c*q[mon[0]]
            assert s.expand(value)==0
        return A[:,cols].applyfunc(s.expand),tc.applyfunc(s.expand)
    results={}
    for kind,entry in r['cases'].items():
        un=entry['universal'];point=entry['point'];ld=point['local_data']
        K=matrix(un['beta_kernel']);b=s.symbols('b0:10');q=s.Matrix(list(map(s.Rational,un['q_point'])))
        beta_rows=[]
        for equation in source[:134]:
            row=[s.S.Zero]*16;constant=s.S.Zero
            for mon,num,den in equation:
                qmon=[index[i] for i in mon if i in index]
                if len(qmon)!=2:continue
                coefficient=s.Rational(num,den)*s.prod(q[i] for i in qmon)
                bmon=[one.index(i) for i in mon if i in one]
                assert len(bmon)<=1
                if bmon:row[bmon[0]]+=coefficient
                else:constant+=coefficient
            assert constant==0;beta_rows.append(row)
        beta_matrix=s.Matrix(beta_rows)
        assert beta_matrix.rank()==6 and beta_matrix*K==s.zeros(134,10)
        cols=un['chart_columns'];F=equations(K*s.Matrix(b));A,tc=jacobian(F,q,cols)
        assert K.rank()==10 and (A-matrix(un['leading_jacobian'])).applyfunc(s.expand)==s.zeros(135,44)
        assert tc==s.Matrix(list(map(s.sympify,un['t_column'])))
        C=matrix(un['first_kernel']);rr=un['first_pivot_rows'];cc=un['first_pivot_columns']
        assert A[:70,:]*C==s.zeros(70,C.cols) and C.rank()==C.cols
        U=s.zeros(44,len(cc));inverse=A.extract(rr,cc).inv()
        for i,j in enumerate(cc):U[j,:]=inverse[i,:]
        R=(A*C).applyfunc(s.cancel);ct=(tc-A*U*tc[rr,:]).applyfunc(s.cancel)
        assert R==matrix(un['residual_matrix'])
        assert ct==s.Matrix(list(map(s.sympify,un['residual_t_column'])))
        N=matrix(un['residual_kernel']);n=4 if kind=='M0' else 3
        assert N.cols==n and N.to_DM().rank()==n
        assert (R.row_join(ct).to_DM()*N.to_DM()).is_zero_matrix
        tangent=(C*N[:-1,:]-U*tc[rr,:]*N[-1,:]).col_join(N[-1,:]).applyfunc(s.cancel)
        assert tangent.to_DM().rank()==n
        assert (A.row_join(tc).to_DM()*tangent.to_DM()).is_zero_matrix
        assert (N[-1,:]==s.zeros(1,n))==(kind=='M1')
        assert (ct[:70,:]==s.zeros(70,1))==(kind=='M0')
        if kind=='M0':
            pivot=next(i for i,v in enumerate(N[-1,:]) if v!=0)
            leading_kernel=s.Matrix.hstack(*[tangent[:-1,j]-tangent[:-1,pivot]*N[-1,j]/N[-1,pivot]
                                             for j in range(n) if j!=pivot])
        else:leading_kernel=tangent[:-1,:]
        assert leading_kernel.to_DM().rank()==3
        assert (A.to_DM()*leading_kernel.to_DM()).is_zero_matrix
        sample={v:i+1 for i,v in enumerate(b)}
        assert A.subs(sample).rank()==41 and A.row_join(tc).subs(sample).rank()==45-n

        # Pointwise full135 Hessians, without the main local-data routine.
        p=list(map(s.Rational,ld['point61']));F=equations([p[i] for i in one])
        qp=s.Matrix([p[i] for i in two]);cols=ld['chart_columns'];fixed=next(i for i in range(45) if i not in cols)
        Ap,tp=jacobian(F,qp,cols);D=Ap.row_join(tp)
        assert D==matrix(ld['jacobian']) and D.rank()==45-n and Ap.rank()==41
        T=matrix(ld['tangent_basis']);assert D*T==s.zeros(135,n) and T.rank()==n
        monomials=list(combinations_with_replacement(range(n),2));H=s.zeros(135,len(monomials))
        full=s.zeros(46,n)
        for j,i in enumerate(cols):full[i,:]=T[j,:]
        full[45,:]=T[-1,:]
        for row,f in enumerate(F):
            for mon,c in f.items():
                if fixed in mon:continue
                i,j=[*mon,*([45]*(2-len(mon)))]
                for v,(a,bb) in enumerate(monomials):
                    H[row,v]+=c*(full[i,a]*full[j,bb]+(full[i,bb]*full[j,a] if a!=bb else 0))
        left=s.Matrix.hstack(*D.T.nullspace()).T
        projected=left*H
        assert projected.rank()==n*(n-1)//2
        independent=list(projected.T.rref()[1]);forms=projected[independent,:]
        hessians=[]
        for row in range(forms.rows):
            B=s.zeros(n)
            for i,(a,bb) in enumerate(monomials):
                B[a,bb]+=forms[row,i]*(2 if a==bb else 1)
                if a!=bb:B[bb,a]+=forms[row,i]
            hessians.append(B)
        stacked=s.Matrix.vstack(*hessians);assert stacked.rank()==n-1
        z=s.Symbol('z');curve=s.Matrix(list(map(s.sympify,point['tail_q_coordinates'])))
        assert curve.subs(z,1)==qp
        assert all(s.expand(sum(c*s.prod(curve[i] for i in mon) for mon,c in f.items() if len(mon)==2))==0 for f in F)
        derivative=curve.diff(z).subs(z,1)-qp*curve.diff(z).subs(z,1)[fixed]/qp[fixed]
        curve_vector=derivative[cols,:].col_join(s.zeros(1,1))
        v=s.Matrix(list(map(s.Rational,point['tail_tangent_in_kernel'])))
        assert v!=s.zeros(n,1) and T*v==curve_vector and stacked*v==s.zeros(stacked.rows,1)
        assert (T[-1,:]==s.zeros(1,n))==(kind=='M1')

        a,bb=s.symbols('a b');x,y,zz,w=X[1:];V=s.Matrix([x,y,zz,w]);c=2 if kind=='M0' else 1
        def tensor(a,bb,c):
            Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2
            Q[2,3]=0 if kind=='M0' else x*x/2
            Q[1,2]=-bb*x*x-c*x*zz-(x*zz if kind=='M1' else 0)
            Q[1,3]=a*x*x+c*x*w
            return Q-Q.T
        change=matrix(entry['normalization']['shear']);exchange=matrix(entry['normalization']['paired_exchange'])
        for S in (change,exchange):assert S*J[1:,1:]*S.T==J[1:,1:]
        transformed=change.inv()*tensor(a,bb,c).subs(dict(zip(V,change*V)),simultaneous=True)*change.inv().T
        assert (transformed-tensor(0,0,c)).applyfunc(s.expand)==s.zeros(4)
        paired=entry['normalization']['paired_source_c']
        transformed=exchange.inv()*tensor(0,0,c).subs(dict(zip(V,exchange*V)),simultaneous=True)*exchange.inv().T
        assert (transformed-tensor(0,0,paired)).applyfunc(s.expand)==s.zeros(4)
        results[kind]=dict(canonical_beta_equation_rank=6,universal_leading_rank=41,universal_full_rank=45-n,
            top70_t_compatibility_zero=kind=='M0',point_tangent_dimension=n,
            quadratic_obstruction_rank=n*(n-1)//2,common_Hessian_kernel_dimension=1,
            curve_is_common_radical=True,transverse_length=n,finite_normalizations_verified=True,
            accepted_t_behavior='nonzero with square zero' if kind=='M0' else 'zero')
    paths=[Path(__file__),INPUT,ROOT/'code/s61_minimal_baseline_thickenings_20261101.py',
           ROOT/'research/notes/S61_MINIMAL_BASELINE_THICKENINGS_20261101.md',
           ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',
           m.source/'symplectic_slice_reduced_equations.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_MINIMAL_BASELINE_REVIEW_20261101',all_passed=True,cases=results,
        scope='Independent135 Jacobians, full-kernel rank upper bounds, actual point Hessians/curve and finite symplectic changes. Formal finite-free argument reviewed in the note.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();result=run()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'output':str(args.output)}))
