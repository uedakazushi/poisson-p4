#!/usr/bin/env python3
"""Universal Jacobian ranks along both minimal nilpotent infinity families."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_infinity_flags import leading_equations_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]

def verify():
    start=time.time();c=s.Symbol('c');b=s.symbols('b0:9');m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];out={}
    x,y,z,w=X[1:]
    for kind in ['M0','M1']:
        Q=s.zeros(5);Q[1,2]=x*x;Q[3,4]=x*x if kind=='M1' else 0;Q[2,3]=-c*x*z-(x*z if kind=='M1' else 0);Q[2,4]=c*x*w-(x*w if kind=='M1' else 0);Q-=Q.T
        pt=m.normalize(X[0]**2*J+Q)['parameters'];beta_matrix=leading_equations_matrix(m,pt)[:,1:];BK=s.Matrix.hstack(*beta_matrix.nullspace()).applyfunc(s.cancel)
        assert BK.shape==(16,9) and (beta_matrix*BK).applyfunc(s.cancel)==s.zeros(134,9)
        bv=BK*s.Matrix(b);beta=dict(zip(one,list(bv)));two,polys=specialize(m,beta,True);q=s.Matrix([pt[i] for i in two]);fixed=next(i for i,v in enumerate(q) if v);cols=[i for i in range(45) if i!=fixed]
        A=s.zeros(135,45);col=s.zeros(135,1)
        for row,f in enumerate(polys):
         val=0
         for mon,coeff in f.items():
          if len(mon)==2:
           i,j=mon;val+=coeff*q[i]*q[j];A[row,i]+=coeff*q[j];A[row,j]+=coeff*q[i]
          elif len(mon)==1:col[row,0]+=coeff*q[mon[0]]
         assert s.cancel(val)==0
        A=A[:,cols].applyfunc(s.cancel);col=col.applyfunc(s.cancel)
        A0=A[:70,:];K=s.Matrix.hstack(*A0.nullspace()).applyfunc(s.cancel)
        cc=list(A0.rref()[1]);rr=list(A0[:,cc].T.rref()[1]);U=s.zeros(44,len(cc));inv=A0[rr,cc].inv()
        for j,i in enumerate(cc):U[i,:]=inv[j,:]
        compat=(col[:70,:]-A0*U*col[rr,:]).applyfunc(s.cancel);assert compat==s.zeros(70,1)
        res=(A[70:,:]*K).applyfunc(s.cancel);ct=(col[70:,:]-A[70:,:]*U*col[rr,:]).applyfunc(s.cancel)
        certificates={}
        for total in [134,135]:
            AA=res[:total-70,:];full=AA.row_join(ct[:total-70,:]);rf,piv=full.to_DM().rref();rf=rf.to_Matrix();free=[j for j in range(full.cols) if j not in piv];null=s.zeros(full.cols,len(free))
            for i,j in enumerate(free):
                null[j,i]=1
                for row,pc in enumerate(piv):null[pc,i]=-rf[row,j]
            assert (full.to_DM()*null.to_DM()).is_zero_matrix
            _,leadpiv=AA.to_DM().rref()
            actual=(len(cc)+len(leadpiv),len(cc)+len(piv),null[-1,:]==s.zeros(1,null.cols))
            expected={('M0',134):(41,41,False),('M0',135):(42,42,False),('M1',134):(42,42,False),('M1',135):(42,43,True)}
            assert actual==expected[kind,total]
            certificates[str(total)]=dict(leading_rank=actual[0],full_rank=actual[1],homogenizer_zero_on_tangent=actual[2],residual_kernel=[[str(v) for v in row] for row in null.tolist()],pivots=list(piv),free_columns=free)
        out[kind]=dict(canonical_q_coordinates=list(map(str,q)),beta_kernel=[[str(v) for v in row] for row in BK.tolist()],first_rank=len(cc),first_pivot_rows=rr,first_pivot_columns=cc,first_pivot_determinant=str(s.factor(A0[rr,cc].det())),first_kernel=[[str(v) for v in row] for row in K.tolist()],first_t_compatibility_zero=True,residual_matrix=[[str(v) for v in row] for row in res.tolist()],residual_t_column=list(map(str,ct)),rank_certificates=certificates)
        print(json.dumps(dict(component=kind,universal_ranks={n:v['full_rank'] for n,v in certificates.items()})),flush=True)
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_infinity_flags.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_MINIMAL_INFINITY_LOCAL_DATA_20261101.json',ROOT/'reports/S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101',all_passed=True,field='QQ(c,b0,...,b8)',cases=out,proof_note='research/notes/S61_MINIMAL_INFINITY_FORMAL_RINGS_20261101.md',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();result=verify()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
if __name__=='__main__':main()
