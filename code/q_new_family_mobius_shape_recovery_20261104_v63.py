#!/usr/bin/env python3
"""Freeze only the completed small exact B/q recovery; no failed cancellation retry."""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_MOBIUS_SHAPE_RECOVERY_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 r,s,q,B,R,Z0=S.symbols('r s q B R Z0');h=s*s-25*r*q*q;F=100*r*(r-1)*q*q-4*r*s+s*s+4*r
 base=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63';rev=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63_REVIEW/report.json';assert read(rev)['all_passed']
 num={};f={};paths=[]
 for name in['F4','p0','p1','p2','p3','p4']:
  p=base/(name+'.json.gz');paths.append(p);d=read(p);num[name]=sum(S.Rational(z)*r**a*s**b*q**c*B**e for z,a,b,c,e in d['terms'])/S.Integer(d['denominator_integer']);f[name]=num[name]/h**d['denominator_h_power']
 assert S.expand(f['F4']-S.Rational(2,405)*F)==0
 Rsource=S.cancel((f['p3']+4*f['p4'])/f['p0']);expected=(80*q*r*r-6*r*r*s+6*r*r-11*r*s+10*r+s)/(8*r*(r-1)*(5*q-s+1));assert S.cancel(Rsource-expected)==0 and not Rsource.has(B)
 T=(r-1)*R-2*r;L=8*R*r*(r-1)-6*r*r-11*r+1;Nq=(s-1)*L-(r-1);qback=Nq/(40*r*T);J=(10*r-1)*s-10*r
 assert S.cancel(Rsource.subs(q,qback)-R)==0 and S.cancel(qback.subs(R,Rsource)-q)==0
 assert S.cancel((5*q-s+1).subs(q,qback)-(r-1)*J/(8*r*T))==0
 D=2359296000*r**4*(r-1)*q*q*(5*q-s+1);Bback=Z0*h*F**2/D
 assert S.cancel((f['p0']/f['F4']**2).subs(B,Bback)-Z0)==0
 assert S.cancel(Bback.subs(Z0,f['p0']/f['F4']**2)-B)==0
 assert S.cancel(D.subs(q,qback)-184320*r*(r-1)**2*Nq**2*J/T**3)==0
 const=[num['p'+str(i)].subs(B,0)for i in range(1,5)]+[S.expand(f['F4']**2*h)];poly=[S.Poly(z,r,s,q)for z in const];mons=sorted(set().union(*(set(p.monoms())for p in poly)));M=S.Matrix([[p.coeff_monomial(m)for p in poly]for m in mons]);assert M.rank()==3
 kernel=[list(map(str,v))for v in M.nullspace()];assert kernel==[['0','0','1/4','1','0'],['0','0','400','0','1']]
 # The Mobius determinant is nonzero as a polynomial in the source shapes.
 moebius_det=-40*r*(r-1)**2*J
 a=80*r*r;b=-6*r*r*s+6*r*r-11*r*s+10*r+s;c=40*r*(r-1);d=8*r*(r-1)*(1-s)
 assert S.expand(a*d-b*c-moebius_det)==0
 cert=dict(normalized_target={'Z0':'p0/F4^2','R':'(p3+4p4)/p0'},R_source=str(expected),T=str(T),L=str(L),Nq=str(Nq),q_recovery=str(qback),B_recovery=str(Bback),
  J=str(J),minus_s_plus5q_plus1_after_recovery=str((r-1)*J/(8*r*T)),D_after_recovery=str(184320*r*(r-1)**2*Nq**2*J/T**3),
  Mobius_determinant=str(moebius_det),constant_B0_relation_columns=['p1','p2','p3','p4','F4^2'],constant_B0_rank=3,constant_B0_kernel=kernel,
  known_tau_fixed_field_reparametrization='QQ(xi,r,s,q,B)=QQ(xi,r,s,R,Z0)',
  displayed_dense_chart_units=['xi','r','s','Z0','T','r-1','Nq','J','h(qback)','F(qback)'],
  remaining_shape_variables=['r','s'],base_target_parameters=['Z0','R'],scope='Rational dense chart only; no global extension/properness/index1 or branch-degree claim')
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');signal.alarm(0)
 report=dict(id=STEM,all_passed=True,exact_source_Bfree_Mobius_ratio=True,exact_two_way_q_B_recovery=True,
  two_unknown_shapes_over_Z0_R=True,constant_B0_span_rank3_no_new_p1_p2_combination=True,
  no_normalized_p1_p2_N12_cancellation_elimination_source26_new_points_or_failed_job_retry=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),rev,*paths]},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],R=str(expected),q_recovery=str(qback)),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
