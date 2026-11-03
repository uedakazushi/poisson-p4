#!/usr/bin/env python3
"""One bounded bivariate polynomial-module GB on the compact J0-safe source.

Distinct from PID Smith over QQ(e), no third inverse variable is introduced.
The annihilator support calculation is attempted only after the module GB.
Any result still requires explicit cofactor/lift certificates before promotion.
"""
import argparse,hashlib,json,resource,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def expression(s):
 z=s.split();n=int(z[0]);assert len(z)==1+4*n
 if not n:return '0'
 out=[]
 for k in range(n):
  c=z[1+4*k];a,e,d=map(int,z[2+4*k:5+4*k]);assert d==0
  t='('+c+')'
  if a:t+='*a^'+str(a)
  if e:t+='*e^'+str(e)
  out.append(t)
 return '('+'+'.join(out)+')'
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 src=ROOT/'work/nonsquare-J0-pid-v63/J0-matrix.txt';assert sha(src)=='4e8ed3e3c8eb89b614b15ab9a11c91ff51b90e5ef6cca826a27043eb6167e683'
 with src.open()as f:
  n,m,k=map(int,f.readline().split());assert[n,m,k]==[36,70,0];rows=[]
  for i in range(n):rows.append('{'+','.join(expression(f.readline())for _ in range(m))+'}')
 script='R=QQ[a,e,MonomialOrder=>GRevLex];F=500*a+90*e^2+3*e-24;Bs=300*a*e+54*e^3+9*e^2-24*e-4;C=100*a+18*e^2-9*e-2;U=a*F*Bs*C*(3*e+2);\nM=matrix{'+','.join(rows)+'};\n'
 script+='print "J0_POLYNOMIAL_MODULE_INPUT_READY";<<flush;G=gb M;N=gens G;print("J0_FULL_MODULE_GB_DONE",numColumns N,sum(flatten entries N,p->#terms p));<<flush;\n'
 script+='Q=coker M;I=ann Q;outI=openOut '+json.dumps(str(out/'annihilator.txt'))+';outI<<toExternalString gens gb I<<endl;close outI;print("J0_ANNIHILATOR_DONE",numgens I,dim(R/I));<<flush;\n'
 script+='S=saturate(I,ideal U);outS=openOut '+json.dumps(str(out/'saturated-annihilator.txt'))+';outS<<toExternalString gens gb S<<endl;close outS;print("J0_SATURATED_SUPPORT_DONE",numgens S,dim(R/S));<<flush;\n'
 (out/'run.m2').write_text(script)
 def limits():resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));resource.setrlimit(resource.RLIMIT_CPU,(120,121))
 t=time.monotonic()
 with(out/'M2.log').open('w')as f:
  try:p=subprocess.run(['M2','--script',str(out/'run.m2')],stdout=f,stderr=subprocess.STDOUT,preexec_fn=limits,timeout=125);status='completed'if p.returncode==0 else f'exit_{p.returncode}'
  except subprocess.TimeoutExpired:status='wall_timeout'
 log=(out/'M2.log').read_text();report={'id':'S61_NONSQUARE_J0_MODULE_GB_20261103_V63','status':status,'source_shape':[36,70],'source_sha256':sha(src),'coefficient_ring':'QQ[a,e]','no_inverse_variable_or_coefficient_function_field':True,'source_J_inversion_used':False,'only_after_gb_known_unit_saturation':'a*F*Bs*C*(3e+2)','module_gb_completed':'J0_FULL_MODULE_GB_DONE'in log,'annihilator_completed':'J0_ANNIHILATOR_DONE'in log,'saturated_support_completed':'J0_SATURATED_SUPPORT_DONE'in log,'CPU_limit_seconds':120,'AS_limit_bytes':2000000000,'solver_seconds':time.monotonic()-t,'seconds':time.monotonic()-start,'scope':'Diagnostic until explicit original-module cofactor identities and known-unit saturation witnesses are checked. Failure/cap is not a nonunit or positive-dimensional support proof.'}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
