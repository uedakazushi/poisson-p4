#!/usr/bin/env python3
"""Bounded QQ(V) actual odd firstjets from full nonlinear CB source."""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from ambient_CA_V1_coefficient_hessian_20261103_v63 import operator
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(x)for x in row]for row in M.tolist()]
class Limit(Exception):pass
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    def alarm(*_):raise Limit('60second bound')
    signal.signal(signal.SIGALRM,alarm);signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
    result=dict(id='AMBIENT_CB_V_ACTUAL_ODD_LINEAR_20261103_v63',all_passed=False)
    framepath=ROOT/'reports/AMBIENT_CB_V_FAMILY_FRAME_20261103_v63/frame.json.gz'
    try:
        data=load(framepath);V=s.Symbol('V');K=s.QQ.frac_field(V)
        dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(K)
        p=s.Matrix(data['point150']);C=s.Matrix(data['fiber_jacobian150_by20']);G=s.Matrix(data['fiber_basis_G150_by3'])
        gi=data['ordered_group19_GL_indices'];extra=data['extra10_fiber_indices'];kap=s.Matrix(data['kernel29']);J=s.Matrix(data['fulltarget150_by29_J']);jc=data['complement28_columns']
        B=C[:,extra];ops=[operator(*divmod(i,5))for i in gi]
        def L(op,v):
            result=s.zeros(150,1)
            for row,col,c in op:result[row]+=c*v[col]
            return result.applyfunc(s.expand)
        LP=[L(op,p)for op in ops];LF=[[L(op,B[:,j])for j in range(10)]for op in ops]
        assert s.Matrix.hstack(*LP,B)==J
        LL={(i,j):L(ops[i],LP[j])for i in range(19)for j in range(i,19)}
        H=[]
        for j in range(19):
            col=sum((kap[i]*LL[min(i,j),max(i,j)]for i in range(19)),s.zeros(150,1))+sum((kap[19+f]*LF[j][f]for f in range(10)),s.zeros(150,1));H.append(col.applyfunc(s.expand))
        qcol=extra.index(3);fiberf={i:extra.index(i)for i in range(3)if i in extra}
        for f,index in enumerate(extra):
            col=sum((kap[i]*LF[i][f]for i in range(19)),s.zeros(150,1))
            if index==3:col-=sum((kap[19+j]*G[:,i]for i,j in fiberf.items()),s.zeros(150,1))
            if index in fiberf:col-=kap[19+qcol]*G[:,index]
            H.append(col.applyfunc(s.expand))
        H=s.Matrix.hstack(*H);W=(H*kap/2).applyfunc(s.expand)
        # Independent univariate source expansion includes the quadratic ratio*f.
        series=[p,(B*kap[19:,:]).applyfunc(s.expand),(-kap[19+qcol]*sum((kap[19+j]*G[:,i]for i,j in fiberf.items()),s.zeros(150,1))).applyfunc(s.expand),s.zeros(150,1)]
        for i in range(18,-1,-1):
            old=series;series=[row.copy()for row in old]
            for degree in range(4):
                v=old[degree]
                for power in range(1,4-degree):
                    v=L(ops[i],v);series[degree+power]+=kap[i]**power/s.factorial(power)*v
            series=[row.applyfunc(s.expand)for row in series]
        assert series[1]==s.zeros(150,1)and series[2]==W;Z=series[3]
        (out/'source_jets.json.gz').write_bytes(gzip.compress((json.dumps(dict(full150_by29_H_kernel=mat(H),pure_kernel_quadratic150=mat(W),pure_kernel_cubic150=mat(Z),all_ratio_f_mixed_terms_retained=True),separators=(',',':'))+'\n').encode(),mtime=0))
        D=dm(J[:,jc]);Hnew=dm(H)*dm(kap.row_join(s.eye(29)[:,jc]));Wr=dm(W);Zr=dm(Z)
        order=list(DomainMatrix.from_Matrix(J[:,:19].T).convert_to(K).rref()[1]);order+= [i for i in range(150)if i not in order]
        rr=D.extract(order,range(28)).transpose().rref()[1];Y=[order[i]for i in rr];assert len(Y)==28
        inv=D.extract(Y,range(28)).inv()
        rem=lambda M:M-D*inv*M.extract(Y,range(M.shape[1]))
        wc=rem(Wr);branch=next(i for i,x in enumerate(wc.to_Matrix())if x);weight=wc[branch,0].element
        mixed=Hnew.extract(range(150),range(1,29));E=rem(mixed*inv)
        zz=rem(Zr-mixed*inv*Wr.extract(Y,[0]))
        norm=lambda M:M-wc.scalarmul(K.one/weight)*M.extract([branch],range(M.shape[1]))
        BY=norm(E);BN=norm(zz).scalarmul(K.one/weight);Lrow=D.extract([branch],range(28))*inv
        A=(BY-BN*Lrow).hstack(BN);pc=list(A.rref()[1]);rank=len(pc)
        rows=list(A.extract(range(150),pc).transpose().rref()[1]);minor=s.factor(K.to_sympy(A.extract(rows,pc).det()))
        coeff=A.extract(rows,pc).inv()*A.extract(rows,range(29));assert A.extract(range(150),pc)*coeff==A
        critical=rem(Hnew);criticalrank=critical.rank()
        (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(dict(actual_base_Y28=Y,branch29th_row=branch,branch_weight=str(K.to_sympy(weight)),actual_all150_by29_odd_firstjet=mat(A.to_Matrix()),rank_rows=rows,rank_columns=pc,rank_minor=str(minor),whole_column_factorization=mat(coeff.to_Matrix()),critical_cotangent_rank=criticalrank),separators=(',',':'))+'\n').encode(),mtime=0))
        result.update(all_passed=True,field='QQ(V)',actual_full_odd_firstjet_rank=rank,source_critical_cotangent_rank=criticalrank,
          all_nonlinear_ratio_f_terms_in_Hessian_and_cubic_retained=True,source_pure_kernel_second_derivative_checked=True,
          rank_minor=str(minor),generic_V_modulus_retained=True,full_CB_completion_requires_separate_geometric_application=True)
    except(Limit,MemoryError)as error:result.update(blocked=str(error),next_step='Reuse source_jets.json.gz if present; finish QQ(V) rank2 trace projection from saved H/W/Z and frame.')
    finally:
        signal.alarm(0);result['seconds']=time.monotonic()-started
        result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),framepath,ROOT/'code/ambient_CA_V1_coefficient_hessian_20261103_v63.py']}
        result['evidence_sha256']={p.name:sha(p)for p in out.iterdir()if p.is_file()};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
