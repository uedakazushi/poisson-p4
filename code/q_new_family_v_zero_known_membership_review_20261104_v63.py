#!/usr/bin/env python3
"""Independent full-family binding and simple graph witnesses, without producer solver."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J,pushforward,quotient_coefficients
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_V_ZERO_KNOWN_MEMBERSHIP_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 primary=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_ELLIPTIC_GRAPH_MEMBERSHIP_20261104_v63';gp=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_GENERIC_KNOWN_EXCLUSIONS_20261104_v63/report.json'
 geo=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
 assert read(primary/'report.json')['all_passed'] and read(gp)['all_passed']
 saved=read(primary/'normalized-source.json');u=s.Symbol('u');n0,n1,n2=s.symbols('n0 n1 n2');x,y,z,w=s.symbols('x y z w')
 P4=s.Matrix(read(geo)['affine_P']).subs({n0:u,n1:1,n2:0});P=s.zeros(5)
 for i in range(4):
  for j in range(i+1,4):
   f=s.Poly(P4[i,j],x,y,z,w)
   for e,c in f.terms():P[i+1,j+1]+=c*X[0]**(2-sum(e))*s.prod(X[k+1]**e[k]for k in range(4))
   P[j+1,i+1]=-P[i+1,j+1]
 R=s.Matrix([[1,0,0,0],[0,0,0,s.Rational(9,8)],[0,1,0,-s.Rational(1,32)],[0,0,-s.Rational(9,2),0]])
 M=s.eye(5);M[1:,1:]=R
 value=P.subs(dict(zip(X,[1,0,0,0,0])))[1:,1:];assert R*value*R.T==J[1:,1:]
 Q=s.eye(5)
 for k,a in enumerate(saved['perspective4']):Q[0,k+1]=s.sympify(a)
 actual=quotient_coefficients(pushforward(P,Q*M))
 assert all(s.expand(a-s.sympify(b))==0 for a,b in zip(actual,saved['full150']))
 model=S61();values=[actual[k]for k in model.certificate['raw_free_indices']]
 assert all(s.expand(a-s.sympify(b))==0 for a,b in zip(values,saved['full61']))
 assert all(actual[k]==0 for k in model.certificate['perspective_orbit_pivots_raw'])
 # The exact point is bound to original P, so no normalization solver is needed.
 witnesses={};paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),primary/'report.json',primary/'normalized-source.json',primary/'graph-residuals.json',gp,geo,ROOT/'code/s61_model.py']
 variables=s.symbols('u1:62')
 for component,stem,ordinal,expected in [('T39','Q51',1,-u/16),('T40','Q52',0,s.Rational(9,64))]:
  p=ROOT/f'reports/{stem}_UNIVERSAL_SECTION_20261101.json';g=read(p);assert g['all_passed'];paths.append(p)
  one=g['weight1_u_indices_zero_based'];two=g['weight2_u_indices_zero_based'];sub={variables[i]:values[i]for i in one}
  residual=s.expand(values[two[ordinal]]-s.sympify(g['quadratic_section_q'][ordinal]).subs(sub,simultaneous=True))
  assert residual==expected
  witnesses[component]=dict(ordinal=ordinal,full61_index=two[ordinal],residual=str(residual),nonzero_on='all u!=0')
 result=dict(witnesses=witnesses,all_u_nonzero_misses_T01=True,T01_reason='accepted closed full T01 orbit dimension <=20, actual full boundary orbit23 for EVERY u!=0',
  all_u_nonzero_misses_T39_T40=True,actual_P_full150_bound_without_producer_solver=True,
  no_u_zero_type_or_other_component_exhaustion_claim=True)
 (out/'certificate.json').write_text(json.dumps(result,indent=2)+'\n')
 report=dict(id=ID,all_passed=True,independent_scope_and_actual_full_tensor_review=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
