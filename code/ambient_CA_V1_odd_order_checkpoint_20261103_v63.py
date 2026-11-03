#!/usr/bin/env python3
"""Save only the cheap numeric frame after the bounded odd-order timeout."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,vector,action
ROOT=Path(__file__).resolve().parents[1]
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    old=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    crit=json.loads(gzip.decompress(old.read_bytes()));C=s.Matrix(crit['fixed_fiber_linear_basis150_by20']);free=crit['affine_slice_free_indices']
    delta=s.zeros(20,1);values=[s.Rational(1,10),0,0,s.Rational(1,10),-s.Rational(1,200),0,0,-s.Rational(1,20),s.Rational(1,400),-s.Rational(1,24000)]
    for i,v in zip(free,values):delta[i]=v
    h,x,y,z,w=X;ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T;pv=vector(P0)+C*delta
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2));P=s.zeros(5)
    for k,(i,j)in enumerate(pairs):P[i,j]=sum(pv[15*k+l]*X[a]*X[b]for l,(a,b)in enumerate(mons));P[j,i]=-P[i,j]
    flag=crit['flag_positions'];rr=list(C.T.rref()[1]);left=C[rr,:].inv();whole=C.row_join(s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P))for i,j in flag]))
    operators=[s.eye(5)[:,i]*s.eye(5)[j,:]for i,j in itertools.product(range(5),repeat=2)];O=s.Matrix.hstack(*[vector(action(M,P))for M in operators]);gi=list(O.rref()[1]);assert len(gi)==20
    source=s.zeros(29,25)
    for k,(i,j)in enumerate(itertools.product(range(5),repeat=2)):
        if [i,j]in flag:source[20+flag.index([i,j]),k]=1
        else:source[:20,k]=left*O[rr,k]
    assert whole*source==O
    group=source[:,gi];ext=group.row_join(s.eye(29)[:,:20]);piv=list(ext.rref()[1]);extra=[i-20 for i in piv[20:]];frame=group.row_join(s.eye(29)[:,extra]);assert len(extra)==9 and frame.det()!=0
    J=O[:,gi].row_join(C[:,extra]);cols=list(J.rref()[1]);assert len(cols)==28
    kernel=J.nullspace()[0];assert J*kernel==s.zeros(150,1)
    data=dict(id='AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63',all_passed=True,
        checkpoint_only_no_cusp_order_claim=True,point150=list(map(str,pv)),point_P=mat(P),fixedCA150_by20=mat(C),
        group20_GL_indices=gi,extra9_fiber_indices=extra,normalizer_source29_frame=mat(frame),
        full150_by29_target_Jacobian=mat(J),kernel29=list(map(str,kernel)),complement28_columns=cols,
        exact_GL_orbit_rank=20,source_frame_dimension=29,target_source_differential_rank=28,
        blocked_matrix_to_construct='H_kappa:150x29 mixed source Hessian; then57constraint block',
        finite_source_recipe='G=product_in_saved_order(I+g_i*E_i); F=G_*(PV+sum k_j*saved_fiber_basis_j)',
        actual_regular_odd8_rows_not_yet_identified=True,
        previous60second_failure_has_NO_saved_setup57=True,
        original_failure_source_sha256=hashlib.sha256((ROOT/'code/ambient_CA_V1_odd_order_20261103_v63.py').read_bytes()).hexdigest(),
        seconds=time.monotonic()-start)
    (out/'checkpoint.json.gz').write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    (out/'report.json').write_text(json.dumps({k:v for k,v in data.items()if k not in['point150','point_P','fixedCA150_by20','normalizer_source29_frame','full150_by29_target_Jacobian','kernel29']},indent=2)+'\n')
    print('CHEAP V1 FRAME SAVED; NO HESSIAN/57/ORDER RERUN',data['seconds'])
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
