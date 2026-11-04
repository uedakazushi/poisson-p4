#!/usr/bin/env python3
"""Exact six-variable Darboux matrix and the component normalization product."""
import argparse, hashlib, json, resource, signal, time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='S61_RATIONAL_COMPONENT_TRANSFER_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 out.mkdir(parents=True,exist_ok=False)
 inputs=[ROOT/'reports/S61_REBUILD_20261030.json',ROOT/'reports/S61_ADAPTER_20261030.json',
  ROOT/'reports/20261030-s61-weight1/finite_identity_certificate.json']
 for q in inputs:assert json.loads(q.read_text())['all_passed'],q
 a,b,c,d,e,f=s.symbols('a b c d e f');p=a*f-b*e+c*d
 C=s.Matrix([[0,a,b,c],[-a,0,d,e],[-b,-d,0,f],[-c,-e,-f,0]])
 R=s.Matrix([[1,0,0,0],[0,0,0,1/c],[-e/c,1,0,-a/c],[-f/p,0,c/p,-b/p]])
 J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
 assert (R*C*R.T-J).applyfunc(s.cancel)==s.zeros(4)
 assert s.cancel(R.det()-1/p)==0
 assert (R.inv()*J*R.inv().T-C).applyfunc(s.cancel)==s.zeros(4)
 assert (R.subs(c,1)*C.subs(c,1)*R.subs(c,1).T-J).applyfunc(s.cancel)==s.zeros(4)
 al=s.symbols('alpha0:4');be=s.symbols('beta0:4')
 G=s.eye(5);H=s.eye(5)
 for i in range(4):G[0,i+1]=al[i];H[0,i+1]=be[i]
 assert G*H==G.subs(dict(zip(al,[al[i]+be[i]for i in range(4)])),simultaneous=True)
 assert G.inv()==G.subs(dict(zip(al,[-v for v in al])),simultaneous=True)
 evidence=dict(skew_value_variables=list(map(str,[a,b,c,d,e,f])),pfaffian=str(p),
  Darboux_rows=[[str(v)for v in row]for row in R.tolist()],determinant='1/(a*f-b*e+c*d)',
  all16_Darboux_and_inverse_identities=True,
  affine_cone_product='U6 x A4 x Y',projective_product='U5 x A4 x Y',
  projective_constant_chart='c=1, a*f-b*e+d nonzero',
  normalization_inverse='record C and old-S61 perspective alpha; inverse perspective, then diag(1,R(C)^(-1)) pushforward',
  rationality_implication='Y rational dimension n => projective reduced component rational dimension n+9',
  new_component_Y16_rationality_assumed_not_proved_here=True)
 ep=out/'certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),ROOT/'code/s61_model.py',
  ROOT/'research/notes/S61_MODEL_AUDIT_20261030.md',*inputs]
 report=dict(id=ID,all_passed=True,exact_six_variable_Darboux_without_root=True,
  finite_S61_perspective_and_full_elimination_reused=True,
  rank4_component_projective_birational_product_n_plus9_proved=True,
  new_component16_or25_rationality_proved_here=False,
  no_new_point_source26_GB_or_model_rebuild=True,
  input_sha256={str(q.relative_to(ROOT)):sha(q)for q in paths},
  evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in ['all_passed','seconds']}))
if __name__=='__main__':
 q=argparse.ArgumentParser(description=__doc__);q.add_argument('--output',type=Path,required=True);main(q.parse_args().output)
