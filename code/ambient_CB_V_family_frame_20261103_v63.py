#!/usr/bin/env python3
"""Cheap whole QQ(V) source frame; no cusp or full image-type inference."""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from ambient_CA_global_flag_model_20261103_v63 import X,canonical,vector
from ambient_CB_global_flag_model_20261103_v63 import raw
from ambient_CA_V1_coefficient_hessian_20261103_v63 import operator
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
    V=s.Symbol('V');K=s.QQ.frac_field(V)
    piv=lambda M:list(DomainMatrix.from_Matrix(M).convert_to(K).rref()[1])
    h,x,y,z,w=X;sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
    change=lambda P:(sw*P.subs({h:x,x:h},simultaneous=True)*sw.T).applyfunc(s.expand)
    cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
    F=[vector(change(canonical(raw(f,0,0))))for f in[h*h,h*x,x*x]]
    G=[vector(change(canonical(raw(0,f,0))))for f in[h*h,h*x,x*x]]
    HH=[vector(change(canonical(raw(0,0,H))))for H in cubics]
    C19=s.Matrix.hstack(*F,*HH);rows=piv(C19.T);left=C19[rows,:].inv()
    ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);point0=vector(ey*D.T-D*ey.T)
    center=left*point0[rows,:];assert C19*center==point0
    critical=ROOT/'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    old=json.loads(gzip.decompress(critical.read_bytes()));free=old['slice_free_indices']
    alpha=list(center[:3,:])+[s.S.Zero]+list(center[3:,:])
    values=[V,0,0,V/2,-V*V/4,0,0,-V/2,V*V/4,-V**3/24]
    variables=s.symbols('s0:10');substitutions=dict(zip(variables,values))
    assert all(s.expand(s.sympify(f).subs(substitutions))==0 for f in old['predicted_entire_critical_ideal'])
    for i,c in zip(free,values):alpha[i]+=c
    assert alpha[3]==0
    point=sum((alpha[i]*F[i]for i in range(3)),s.zeros(150,1))+sum((alpha[i+4]*HH[i]for i in range(16)),s.zeros(150,1))
    C=s.Matrix.hstack(*F,-sum((alpha[i]*G[i]for i in range(3)),s.zeros(150,1)),*HH)
    cr=piv(C.T);Cleft=C[cr,:].inv()
    def L(index,v):
        result=s.zeros(150,1)
        for row,col,c in operator(*divmod(index,5)):result[row]+=c*v[col]
        return result.applyfunc(s.expand)
    O=s.Matrix.hstack(*[L(i,point)for i in range(25)]);gi=piv(O)
    flag=old['flag_positions'];whole=C.row_join(s.Matrix.hstack(*[O[:,5*i+j]for i,j in flag]));source=s.zeros(29,25)
    for index in range(25):
        pos=list(divmod(index,5))
        if pos in flag:source[20+flag.index(pos),index]=1
        else:source[:20,index]=Cleft*O[cr,index]
    assert (whole*source-O).applyfunc(s.expand)==s.zeros(150,25)
    group=source[:,gi];pp=piv(group.row_join(s.eye(29)[:,:20]));extra=[i-len(gi)for i in pp[len(gi):]]
    frame=group.row_join(s.eye(29)[:,extra]);assert frame.shape==(29,29)
    frame_det=s.factor(K.to_sympy(DomainMatrix.from_Matrix(frame).convert_to(K).det()));assert frame_det
    J=O[:,gi].row_join(C[:,extra]);jc=piv(J)
    null=DomainMatrix.from_Matrix(J).convert_to(K).nullspace().to_Matrix();assert null.shape==(1,29)
    kernel=null.T;common=s.gcd_list([c for c in kernel if c])
    kernel=kernel.applyfunc(lambda c:s.cancel(c/common));assert (J*kernel).applyfunc(s.cancel)==s.zeros(150,1)
    assert len(gi)==19 and len(extra)==10 and len(jc)==28
    data=dict(parameter='V',field='QQ(V)',point150=list(map(str,point)),fiber_alpha20=list(map(str,alpha)),
      fiber_jacobian150_by20=mat(C),fiber_basis_F150_by3=mat(s.Matrix.hstack(*F)),fiber_basis_G150_by3=mat(s.Matrix.hstack(*G)),
      ordered_group19_GL_indices=gi,extra10_fiber_indices=extra,fullsource29_frame=mat(frame),source_frame_det=str(frame_det),
      fulltarget150_by29_J=mat(J),kernel29=list(map(str,kernel)),kernel_removed_common_factor=str(common),complement28_columns=jc,
      flag_positions=flag,source_ratio_times_f_mixed_term='-delta_q*sum(delta_f_i*G_i); retain in every higher-jet stage',
      finite_source_recipe='ordered GL19 exponentials applied to full CB ratio/f/H source restricted to saved ten fiber coordinates')
    (out/'frame.json.gz').write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    paths=[Path(__file__),critical,ROOT/'code/ambient_CA_V1_coefficient_hessian_20261103_v63.py',ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py']
    result=dict(id='AMBIENT_CB_V_FAMILY_FRAME_20261103_v63',all_passed=True,field='QQ(V)',GL_orbit_rank=19,
      all29_source_frame_verified=True,source_frame_det=str(frame_det),actual_target_J_rank=28,
      entire_polynomial_V_curve_lies_in_full_critical_scheme=True,
      modulus_retained_not_transported_from_V1=True,no_odd_ideal_or_cusp_type_claim=True,
      next_matrix='mixed kernel Hessian for ordered GL19+nonlinear CB fiber10, then actual odd linear matrix or57 paired constraints',
      input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'frame.json.gz':sha(out/'frame.json.gz')},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print('CB QQ(V) ACTUAL FRAME29 /GL19 /TARGET28 PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
