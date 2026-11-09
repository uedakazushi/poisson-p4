"""Close the s=-1/3 rational exceptional NN slice with a singular plane.

No parameter is inverted. The remaining equations solve two coefficients
polynomially; the entire two-parameter family satisfies Jacobi, and its
nonzero quintic lies in (x2,x3)^2 identically.
"""
from pathlib import Path
import argparse,json,hashlib
from itertools import combinations
import sympy as s

def audit(model_path,jet_path,cert_path,out):
 model=json.loads(model_path.read_text());jet=json.loads(jet_path.read_text());cert=json.loads(cert_path.read_text())
 z=s.symbols('z0:40');x=s.symbols('x0:5');a,b=s.symbols('a b');t=s.Symbol('t')
 loc={str(v):v for v in z+x+(a,b,t)}
 checks={}
 remaining=[s.sympify(e,locals=loc) for e in cert['remaining_equations']]
 e1=z[11]/3-z[18]*z[20]
 e2=-z[1]/4+z[11]/12-z[18]*z[20]/4+3*z[20]**2/4
 checks['first_equation_present']=e1 in remaining
 checks['second_equation_present']=e2 in remaining
 checks['second_polynomial_elimination']=s.expand(-4*e2+e1-(z[1]-3*z[20]**2))==0
 add={z[11]:3*z[18]*z[20],z[1]:3*z[20]**2}
 checks['remaining_system_vanishes']=all(s.expand(e.xreplace(add))==0 for e in remaining)
 sub={s.sympify(k):s.sympify(v,locals=loc) for k,v in cert['final_substitution'].items()}
 sub.update(add)
 sub={k:s.expand(v.xreplace({z[18]:a,z[20]:b})) for k,v in sub.items()}
 sub.update({z[18]:a,z[20]:b})
 checks['all_40_parameters_accounted_for']=set(sub)==set(z)
 P=s.zeros(5)
 for ij,e in model['linear_model']['P_entries'].items():
  i,j=map(int,ij);P[i,j]=s.expand(s.sympify(e,locals=loc).xreplace(sub));P[j,i]=-P[i,j]
 pf=[]
 for k in range(5):
  i,j,l,m=[v for v in range(5) if v!=k]
  pf.append(s.expand((-1)**k*(P[i,j]*P[l,m]-P[i,l]*P[j,m]+P[i,m]*P[j,l])))
 F=s.expand(sum(x[i]*pf[i] for i in range(5))/5)
 oldF=s.expand(s.sympify(model['F'],locals=loc).xreplace(sub))
 checks['canonical_F_agrees']=s.expand(F-oldF)==0
 for i in range(5):
  checks[f'divergence_{i}']=s.expand(sum(s.diff(P[j,i],x[j]) for j in range(5)))==0
  checks[f'gradient_pfaffian_{i}']=s.expand(s.diff(F,x[i])-pf[i])==0
 for i,j,k in combinations(range(5),3):
  J=s.expand(sum(P[i,l]*s.diff(P[j,k],x[l])+P[j,l]*s.diff(P[k,i],x[l])+P[k,l]*s.diff(P[i,j],x[l]) for l in range(5)))
  checks[f'direct_Jacobi_{i}{j}{k}']=J==0
 checks['all_stored_jacobi_equations_zero']=all(s.expand(s.sympify(e,locals=loc).xreplace(sub))==0 for e in model['jacobi_equations'])
 checks['all_stored_e8_equations_zero']=all(s.expand(s.sympify(e,locals=loc).xreplace(sub))==0 for e in jet['all_e8_equations'])
 poly=s.Poly(F,*x)
 checks['F_in_singular_plane_square']=all(m[2]+m[3]>=2 for m,c in poly.terms())
 checks['F_uniformly_nonzero']=poly.coeff_monomial(x[1]**3*x[3]**2)==-s.Rational(1,120)
 for i in range(5):checks[f'partial_{i}_vanishes_on_plane']=s.expand(s.diff(F,x[i]).subs({x[2]:0,x[3]:0}))==0
 assert all(checks.values()),[k for k,v in checks.items() if not v]
 report={'name':'31_NN_k1_IB_rational','scope':'Exact complete two-parameter solution of the reduced exceptional slice; every member has a nonzero quintic singular along a plane',
 'parameter_order':['a','b'],'final_substitution':{str(k):str(v) for k,v in sub.items()},
 'elimination_equations':{'e1':str(e1),'e2':str(e2),'z11':str(add[z[11]]),'z1':str(add[z[1]])},
 'P_entries':{f'{i}{j}':str(P[i,j]) for i,j in combinations(range(5),2)},'F':str(F),
 'singular_plane_generators':['x2','x3'],'uniform_nonzero_coefficient':{'monomial':'x1**3*x3**2','coefficient':'-1/120'},
 'checks':checks,'check_count':len(checks),'all_passed':True,
 'input_sha256':{name:{'filename':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for name,p in [('model',model_path),('jets',jet_path),('certificate',cert_path)]}}
 out.write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({'check_count':len(checks),'all_passed':True,'parameters':['a','b'],'singular_plane':['x2','x3'],'F_nonzero':True}))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('model',type=Path);ap.add_argument('jet',type=Path);ap.add_argument('certificate',type=Path);ap.add_argument('output',type=Path)
 a=ap.parse_args();audit(a.model,a.jet,a.certificate,a.output)
