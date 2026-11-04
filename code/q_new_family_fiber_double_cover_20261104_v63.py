#!/usr/bin/env python3
"""The actual generic fiber as an explicit quadratic cover of the J-line."""
import argparse, hashlib, json, resource, signal, time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_FIBER_DOUBLE_COVER_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 old=ROOT/'reports/Q_NEW_FAMILY_EXPLICIT_FIBER_GENERATOR_20261104_v63';v0=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63_REVIEW/report.json'
 assert read(old/'report.json')['all_passed'] and read(v0)['all_passed']
 cert=read(old/'certificate.json');u,K,Iv,w=s.symbols('u K Iv w_fiber');r,ss,q,J=s.symbols('r0 s0 q0 J')
 source=s.sympify(cert['J_in_w'])
 sub={K:(8*ss+u*q)/(8*ss-u*q),Iv:u/(512*(8*ss+u*q))}
 h=ss*ss-25*r*q*q;A=200*r*q*(ss-1);B=ss-25*r*q*q;D=r*q*q*h;U=1600*r;C=s.Rational(256,25)
 expected=C*(A+B*w)**2/(D*(U-w*w))
 diff=s.cancel(source.subs(sub,simultaneous=True)-expected)
 n,d=s.fraction(diff)
 assert s.Poly(n,u).rem(s.Poly(u*u-U,u)).is_zero
 jcrit=256-256/(25*r*q*q)
 assert s.factor(C*(A*A-B*B*U)/(D*U)-jcrit)==0
 F=s.expand(J*D*(U-w*w)-C*(A+B*w)**2)
 disc=s.factor(s.discriminant(F,w));assert s.cancel(disc-(80*D)**2*r*J*(J-jcrit))==0
 p=s.Poly(F,w);assert p.degree()==2 and p.nth(2)!=0
 boundary=200*r*q/ss
 assert s.cancel(expected.subs(w,boundary)-256)==0
 wzero=-A/B;wcrit=-B*U/A
 assert s.cancel(expected.subs(w,wzero))==0
 assert s.cancel(expected.subs(w,wcrit)-jcrit)==0
 # The two branch points are rational on the actual rational curve.
 derivative=s.factor(s.diff(expected,w))
 assert s.cancel(s.diff(expected,w).subs(w,wzero))==0
 assert s.cancel(s.diff(expected,w).subs(w,wcrit))==0
 result=dict(new_v0_coordinates=['xi0','r0','s0','q0','B0'],
  WARNING='newv0 coordinates, not oldt2 same-named functions',
  scalar_curve_parameter='w_fiber from the explicit actual Luroth certificate',
  J_in_w=str(expected),quadratic_equation=str(F),discriminant=str(disc),
  branch_values=['0',str(jcrit)],branch_points_in_w=[str(wzero),str(wcrit)],
  actual_function_field='K_B(w_fiber)=K_B(J)(sqrt(r0*J*(J-Jcrit)))',
  degree=2,geometric_branch_count=2,actual_curve_genus=0,
  v0_actual_section_w=str(boundary),v0_actual_section_J=256,
  no_whole_special_fiber_or_morphism_classification=True,
  no_NEW_beta_equation_or_spectral_degree_claim=True)
 (out/'certificate.json').write_text(json.dumps(result,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),old/'report.json',old/'certificate.json',v0]
 report=dict(id=ID,all_passed=True,independent_review_pending=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
