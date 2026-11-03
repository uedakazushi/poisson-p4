#!/usr/bin/env python3
"""Intrinsic ramification cotangents of the three finite normalizers at B18.

The CB mixed ratio Hessian is included. Exact rank7 excludes importing the
old Ac8/Ac7 types; it does not determine the nonlinear critical ideals.
"""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector,action
from ambient_CB_global_flag_model_20261103_v63 import raw as rawB
ROOT=Path(__file__).resolve().parents[1];h,x,y,z,w=X
def mj(M):return [[str(v)for v in row]for row in M.tolist()]
def rank_certificate(M):
 cols=list(M.rref()[1]);V=M[:,cols];rows=list(V.T.rref()[1]);D=V[rows,:];coeff=D.inv()*M[rows,:];assert V*coeff==M
 return dict(rank=len(cols),minor_rows=rows,minor_columns=cols,minor=str(s.factor(D.det())),basis= mj(V),all_column_coefficients=mj(coeff))
def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cu=[s.prod((h,x,z,w)[i]for i in m)for m in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in m)];sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
 def change(P):return(sw*P.subs({h:x,x:h},simultaneous=True)*sw.T).applyfunc(s.expand)
 A=[change(canonical(raw(f,0,0)))for f in[h*h,h*x,x*x]]+[change(canonical(raw(0,1,0)))]+[change(canonical(raw(0,0,H)))for H in cu]
 B=[change(canonical(rawB(f,0,0)))for f in[h*h,h*x,x*x]]+[change(canonical(rawB(0,h*h,0)))]+[change(canonical(rawB(0,0,H)))for H in cu];N=[A[0],A[1],*A[3:]]
 ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T;positions=[(i,j)for i,j in itertools.product(range(5),repeat=2)if(i==1 and j!=1)or(i==0 and j not in(0,1))or(j==2 and i!=2)];assert len(positions)==9;records={}
 for name,basis,m,dimension,oldc in [('CA',A,-2,29,8),('CB',B,s.Rational(1,2),29,8),('N',N,-2,28,7)]:
  C=s.Matrix.hstack(*map(vector,basis));rows=list(C.T.rref()[1]);left=C[rows,:].inv()
  def rem(M):return M-C*left*M[rows,:]
  Q0=rem(s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P0))for i,j in positions]));M=s.zeros(5);M[1,0]=1;M[0,3]=m;v=s.Matrix([M[i,j]for i,j in positions]);assert Q0*v==s.zeros(150,1)
  co=left*vector(action(M,P0))[rows,0];assert C*co==vector(action(M,P0));R=s.Matrix.hstack(*[vector(action(M,P))for P in basis]);hessian=s.zeros(150,len(basis))
  if name=='CB':
   # source g=-q*f. Kernel fiber velocity is -co. Its mixed Hessian
   # with deltaq=1 is +co_f, including the essential f1=5/2 term.
   assert co[0]==co[2]==co[3]==0 and co[1]==s.Rational(5,2)
   hessian[:,3]=vector(change(canonical(rawB(0,co[1]*h*x,0))));R+=hessian
  R=rem(R);qcert=rank_certificate(Q0);assert qcert['rank']==8;qc=qcert['minor_columns'];qr=qcert['minor_rows'];linear=(R-Q0[:,qc]*Q0[qr,qc].inv()*R[qr,:]).applyfunc(s.expand);cert=rank_certificate(linear);assert cert['rank']==7
  records[name]=dict(fiber_dimension=len(basis),cone_source_dimension=dimension,projective_source_dimension=dimension-1,flag_Lie_positions=positions,kernel_flag_direction=list(map(str,v)),kernel_source_fiber_velocity=list(map(str,-co)),Q0=mj(Q0),Q0_rank_certificate=qcert,fiber_derivative_before_flag_kernel_correction=mj(R),CB_mixed_Hessian_correction=mj(hessian),entire_ramification_linearization=mj(linear),rank7_certificate=cert,source_ramification_cotangent_codimension=7,cone_ramification_tangent_dimension=dimension-7,projective_ramification_tangent_dimension=dimension-8,old_Ac_pinched_type_expected_codimension=oldc+1,old_type_excluded_by_intrinsic_ramification=True)
 files=[Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py',ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json',ROOT/'research/notes/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63.md'];report=dict(id='AMBIENT_B18_RAMIFICATION_COTANGENT_20261103_v63',all_passed=True,field='QQ',families=records,exact_ramification_ideal_linear_rank=7,old_Qa_A8_A7_local_types_do_not_extend_to_generic_B18=True,nonlinear_odd_ideals_or_fullD_type_claimed=False,scope='Intrinsic Fitting ramification cotangents of finite normalizers; tangent dimensions do not assert actual critical-locus dimensions.',input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in files},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('B18 NORMALIZER RAMIFICATION RANK7 + OLD TYPES EXCLUDED PASS',report['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
