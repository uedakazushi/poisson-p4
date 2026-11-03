#!/usr/bin/env python3
"""Bounded paired-source odd order AFTER eight actual odd unit constraints.

No rawY=0 slice is substituted for the odd ideal. Runtime60s/address2.2GB.
The57by57 matrix is saved before any iterative jet computation.
"""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,vector,action
ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def frac(x):return Q(int(s.numer(x)),int(s.denom(x)))
def encode(obj):return json.dumps(obj,separators=(',',':'))+'\n'
class Limit(Exception):pass
def alarm(*_):raise Limit('60second bound')

def evaluator(point,basis,group_indices,degree):
    N=degree
    def zero():return[Q(0)]*(N+1)
    def const(v):return[Q(v)]+[Q(0)]*N
    def add(a,b):return[x+y for x,y in zip(a,b)]
    def neg(a):return[-x for x in a]
    def mul(a,b):
        c=zero()
        for i,x in enumerate(a):
            if x:
                for j in range(N+1-i):
                    if b[j]:c[i+j]+=x*b[j]
        return c
    def reciprocal(a):
        assert a[0];c=zero();c[0]=1/a[0]
        for i in range(1,N+1):c[i]=-sum((a[j]*c[i-j]for j in range(1,i+1)),Q(0))/a[0]
        return c
    def group(q,inverse=False):
        G=[[const(int(i==j))for j in range(5)]for i in range(5)]
        order=list(range(20));order=order[::-1]if inverse else order
        for k in order:
            i,j=divmod(group_indices[k],5);v=list(q[k][:N+1])
            if inverse:
                if i==j:v=add(reciprocal(add(const(1),v)),const(-1))
                else:v=neg(v)
            column=[mul(G[r][i],v)for r in range(5)]
            for r in range(5):G[r][j]=add(G[r][j],column[r])
        return G
    def run(q):
        G=group(q);Gi=group(q,True)
        polynomial=[const(v)for v in point]
        for k,b in enumerate(basis):
            for i,c in enumerate(b):
                if c:polynomial[i]=add(polynomial[i],[c*x for x in q[20+k][:N+1]])
        forms=[]
        for a,b in MONS:
            row=[]
            for r,t in MONS:
                v=mul(Gi[a][r],Gi[b][t])
                if r!=t:v=add(v,mul(Gi[a][t],Gi[b][r]))
                row.append(v)
            forms.append(row)
        middle=[]
        for k in range(10):
            row=[zero()for _ in range(15)]
            for m in range(15):
                coefficient=polynomial[15*k+m]
                if any(coefficient):
                    for n in range(15):
                        if any(forms[m][n]):row[n]=add(row[n],mul(coefficient,forms[m][n]))
            middle.append(row)
        output=[zero()for _ in range(150)]
        for k,(i,j)in enumerate(PAIRS):
            for l,(a,b)in enumerate(PAIRS):
                coefficient=add(mul(G[i][a],G[j][b]),neg(mul(G[i][b],G[j][a])))
                if any(coefficient):
                    for m in range(15):
                        if any(middle[l][m]):output[15*k+m]=add(output[15*k+m],mul(coefficient,middle[l][m]))
        return output
    return run

def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    signal.signal(signal.SIGALRM,alarm);signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
    status=dict(id='AMBIENT_CA_V1_ODD_ORDER_20261103_v63',all_passed=False,blocked=True)
    try:
        cpath=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';crit=load(cpath)
        C=s.Matrix(crit['fixed_fiber_linear_basis150_by20']);free=crit['affine_slice_free_indices']
        delta=s.zeros(20,1);values=[s.Rational(1,10),0,0,s.Rational(1,10),-s.Rational(1,200),0,0,-s.Rational(1,20),s.Rational(1,400),-s.Rational(1,24000)]
        for i,v in zip(free,values):delta[i]=v
        h,x,y,z,w=X;ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);Pzero=ey*D.T-D*ey.T;pv=vector(Pzero)+C*delta;P=s.zeros(5)
        for k,(i,j)in enumerate(PAIRS):
            P[i,j]=sum(pv[15*k+l]*X[a]*X[b]for l,(a,b)in enumerate(MONS));P[j,i]=-P[i,j]
        flag=crit['flag_positions'];leftrows=list(C.T.rref()[1]);left=C[leftrows,:].inv()
        whole=C.row_join(s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P))for i,j in flag]))
        operators=[s.eye(5)[:,i]*s.eye(5)[j,:]for i,j in itertools.product(range(5),repeat=2)]
        O=s.Matrix.hstack(*[vector(action(M,P))for M in operators]);gi=list(O.rref()[1]);assert len(gi)==20
        assert O.row_join(pv).rank()==20
        source=s.zeros(29,25)
        for k,(i,j)in enumerate(itertools.product(range(5),repeat=2)):
            if [i,j]in flag:source[20+flag.index([i,j]),k]=1
            else:source[:20,k]=left*O[leftrows,k]
        assert whole*source==O
        group_source=source[:,gi];extension=group_source.row_join(s.eye(29)[:,:20]);pc=list(extension.rref()[1]);extras=[j-20 for j in pc[20:]];assert len(extras)==9
        frame=group_source.row_join(s.eye(29)[:,extras]);assert frame.det()!=0
        J=O[:,gi].row_join(C[:,extras]);jc=list(J.rref()[1]);assert len(jc)==28
        kernel=J.nullspace()[0];freecol=next(i for i in range(29)if i not in jc);assert kernel[freecol]==1
        Jr=J[:,jc];S=s.eye(29)[:,jc]
        gs=[operators[i]for i in gi];Mk=sum((kernel[i]*gs[i]for i in range(20)),s.zeros(5));Pk=sum((kernel[20+k]*s.Matrix(5,5,lambda i,j:0)for k in []),s.zeros(5))
        def bivector(v):
            QP=s.zeros(5)
            for k,(i,j)in enumerate(PAIRS):
                QP[i,j]=sum(v[15*k+l]*X[a]*X[b]for l,(a,b)in enumerate(MONS));QP[j,i]=-QP[i,j]
            return QP
        bb=[bivector(C[:,i])for i in extras];Pk=sum((kernel[20+k]*bb[k]for k in range(9)),s.zeros(5));LMp=action(Mk,P)
        Hp=[]
        for j in range(29):
            Mj=gs[j]if j<20 else s.zeros(5);Pj=s.zeros(5)if j<20 else bb[j-20];N=s.zeros(5)
            if j<20:
                for i in range(20):
                    if i<j:N+=kernel[i]*gs[i]*gs[j]
                    if i>j:N+=kernel[i]*gs[j]*gs[i]
            K=N-(Mk*Mj+Mj*Mk)/2
            value=action(K,P)+(action(Mk,action(Mj,P))+action(Mj,LMp))/2+action(Mk,Pj)+action(Mj,Pk)
            Hp.append(vector(value))
        Horiginal=s.Matrix.hstack(*Hp);Hnew=Horiginal*kernel.row_join(S);W=Hnew[:,0]/2
        rg=list(O[:,gi].T.rref()[1]);assert len(rg)==20
        order=rg+[i for i in range(150)if i not in rg];yr=list(Jr[order,:].T.rref()[1]);Y=[order[i]for i in yr];assert Y[:20]==rg and len(Y)==28
        inv=Jr[Y,:].inv();wc=W-Jr*inv*W[Y,:];branch=next(i for i,v in enumerate(wc)if v);assert branch not in Y
        base=Y+[branch];Jn=s.zeros(150,1).row_join(Jr)
        even=Jn[rg,:].row_join(s.zeros(20,28));odd=Hnew.row_join(Jr)
        first=even.col_join(odd[base,:]);assert first.rank()==49
        rows=list(first.col_join(odd).T.rref()[1]);assert rows[:49]==list(range(49))and len(rows)==57
        selected=[i-49 for i in rows[49:]];block=first.col_join(odd[selected,:]);assert block.shape==(57,57)
        blockinv=block.inv();assert blockinv*block==s.eye(57)
        E=Hnew[:,1:]*inv;E=E-Jr*inv*E[Y,:];BY=E-wc/wc[branch]*E[branch,:]
        regular=s.eye(28)[:20,:].col_join(BY[selected,:]);assert regular.det()!=0
        point=list(map(frac,pv));fbase=[list(map(frac,C[:,i]))for i in extras]
        setup=dict(point150=list(map(str,pv)),group_GL_indices=gi,extra_fiber_indices=extras,source29_frame=mat(frame),source29_target_Jacobian=mat(J),kernel29=list(map(str,kernel)),source_r28_columns=jc,full150_by29_mixed_Hessian=mat(Hnew),orbit20_raw_rows=rg,actual_target_Y28=Y,actual_branch29th_row=branch,corrected_quadratic150=list(map(str,wc)),actual_odd_selected8_rows=selected,regular_odd8_plus_orbit20_minor=str(s.factor(regular.det())),full57_linear_block=mat(block),inverse57=mat(blockinv))
        (out/'setup.json.gz').write_bytes(gzip.compress(encode(setup).encode(),mtime=0))
        print('SETUP57/unit8 ready',time.monotonic()-start,flush=True)
        # Kernel coordinate has oddpart epsilon. Pairing +/-epsilon then
        # makes it exactly tracefree after all actual R outputs are even.
        q=[[Q(0)]*8 for _ in range(29)]
        for i in range(29):q[i][1]=frac(kernel[i])
        stages=[];previous=None;exact=None
        for degree in[3,5,7]:
            run=evaluator(point,fbase,gi,degree);baseline=run(q)
            if previous is not None:
                assert all(row[:len(old)]==old for row,old in zip(baseline,previous))
            rhs=s.Matrix([-s.Rational(baseline[i][degree-1])for i in rg]+[-s.Rational(baseline[i][degree])for i in base+selected]);solution=blockinv*rhs
            even_coeff=solution[:29,:];odd_coeff=solution[29:,:]
            delta_even=kernel*even_coeff[0]+S*even_coeff[1:,:];delta_odd=S*odd_coeff
            for i in range(29):q[i][degree-1]=frac(delta_even[i]);q[i][degree]=frac(delta_odd[i])
            corrected=[list(row)for row in baseline]
            output_even=Jn*even_coeff;output_odd=Hnew*even_coeff+Jr*odd_coeff
            for i in range(150):corrected[i][degree-1]+=frac(output_even[i]);corrected[i][degree]+=frac(output_odd[i])
            assert all(corrected[i][degree-1]==0 for i in rg)
            assert all(corrected[i][degree]==0 for i in base+selected)
            nonzero=[i for i,row in enumerate(corrected)if row[degree]]
            stages.append(dict(source_degree=degree,conductor_nu_degree=(degree-1)//2,all150_output_series=[[str(c)for c in row]for row in corrected],source_even_coeff=list(map(str,delta_even)),source_odd_coeff=list(map(str,delta_odd)),remaining_odd_nonzero_rows=nonzero))
            previous=corrected
            print('ACTUAL8ODD_ELIM order',degree,'remaining',len(nonzero),'sec',time.monotonic()-start,flush=True)
            if nonzero:exact=(degree-1)//2;break
        degree=stages[-1]['source_degree'];final=evaluator(point,fbase,gi,degree)(q);assert final==previous
        (out/'jets.json.gz').write_bytes(gzip.compress(encode(dict(stages=stages,final_source29_series=[[str(c)for c in row]for row in q],all150_final_direct_replay=True)).encode(),mtime=0))
        status.update(all_passed=True,blocked=False,CA_V1_actual_odd_regular_constraints=8,actual_source_pair_base_outputs_even=True,orbit20_target_parameters_fixed=True,
            source_kernel_odd_part_exact_epsilon=True,exact_normal_conductor_order=exact,lower_bound_if_not_detected=4 if exact is None else exact,
            all150_lower_odd_coefficients_checked=True,raw_Y_slice_not_substituted_for_odd_constraints=True,final_direct_source_series_replay_passed=True,
            scope='Bounded CA V1 odd order after8 actual odd constraints; whole genericV type requires the separately stated dense orbit and slice application proof.')
    except (Limit,MemoryError)as error:
        status.update(blocker=str(error),next_step='Reuse saved setup57 and evaluate paired-source jets; do not rebuild completed matrices.')
    finally:
        signal.alarm(0);status['seconds']=time.monotonic()-start;status['max_rss_kb']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        status['source_sha256']=sha(Path(__file__));status['evidence_sha256']={p.name:sha(p)for p in out.iterdir()if p.is_file()}
        (out/'report.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status,indent=2),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
