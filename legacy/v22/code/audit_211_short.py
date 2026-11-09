"""Independent exact audit of the 211-triple short exclusion certificate."""
from pathlib import Path
import argparse,json
import sympy as s

def audit(old_path,jet_path,certificate_path,out):
 old=json.loads(old_path.read_text()); jets=json.loads(jet_path.read_text()); cert=json.loads(certificate_path.read_text())
 z=s.symbols('z0:43');t=s.Symbol('t');loc={str(v):v for v in z+(t,)}
 initial={s.sympify(k):s.sympify(v) for k,v in cert['initial_substitution'].items()}
 checks={}
 conic_sub={z[37]:0,z[39]:0,z[41]:0,z[38]:1,z[40]:1,z[42]:1,z[31]:1,
 z[23]:0,z[26]:0,z[27]:0,z[24]:-z[33]**2/3,z[28]:z[32]**2/3,
 z[29]:z[32]*z[33]/3-2*z[36]/3,z[13]:z[34]**2,z[14]:z[34]*z[35],z[9]:-z[35]**2,
 z[11]:2*z[20]+2*z[32]*z[34],z[8]:2*z[18]-2*z[33]*z[35],z[12]:2*z[21]+z[32]*z[35]+z[33]*z[34]}
 def mat(which,k,sub):return s.Matrix([[s.expand(s.sympify(v,locals=loc).subs(sub,simultaneous=True)) for v in row] for row in jets['components'][k][which]])
 B=mat('B',0,conic_sub);v=s.Matrix([1,-t*t*z[34]-t*z[32],-t*t*z[35]-t*z[33]])
 for k,e in enumerate(B+v*v.T):checks[f'conic_rank_one_identity_{k}']=s.expand(e)==0
 B1=mat('B',1,conic_sub);B2=mat('B',2,conic_sub)
 expected1=s.Matrix([[0,t*z[32]**2/3,t*z[30]],[t*z[32]**2/3,-t*z[20]-z[32]/2,-t*z[22]],[t*z[30],-t*z[22],t*z[16]]])
 expected2=s.Matrix([[0,-t*z[33]**2/3,t*z[25]],[-t*z[33]**2/3,-t*z[18]+z[33]/2,-t*z[19]],[t*z[25],-t*z[19],t*z[16]]])
 for k,e in enumerate(B1-expected1):checks[f'first_leaf_matrix_{k}']=s.expand(e)==0
 for k,e in enumerate(B2-expected2):checks[f'second_leaf_matrix_{k}']=s.expand(e)==0
 post={**conic_sub,z[32]:0,z[33]:0,z[30]:0,z[25]:0}
 # Compose the zero substitutions into the prior nonlinear expressions.
 clear={z[32]:0,z[33]:0,z[30]:0,z[25]:0}
 post={k:s.sympify(v).subs(clear,simultaneous=True) for k,v in post.items()}
 for leaf,sign in [(1,-1),(2,1)]:
  H=mat('H',leaf,post); B=mat('B',leaf,post)
  H0=H.subs(t,0); target=s.zeros(3);target[1,2]=sign
  checks[f'leaf_{leaf}_constant_H']=H0==target
  center=(H+t*t*s.eye(3)/2)*B
  row=[s.expand(center[1,k]).coeff(t,1) for k in range(3)]
  expected=[0,-sign*(z[22] if leaf==1 else z[19]),sign*z[16]]
  checks[f'leaf_{leaf}_order_one_center']=all(s.expand(a-b)==0 for a,b in zip(row,expected))
 B1=mat('B',1,initial);checks['first_leaf_final_B']=B1==-t*z[20]*s.diag(0,1,0)
 q19=s.expand(s.sympify(old['jacobi_equations'][19],locals=loc).subs(initial,simultaneous=True))
 q143=s.expand(s.sympify(old['jacobi_equations'][143],locals=loc).subs(initial,simultaneous=True))
 checks['first_jacobi_restriction']=q19==z[29]/5
 checks['second_jacobi_restriction']=q143==z[5]-z[29]*z[34]
 H1=mat('H',1,initial);E=s.expand(((H1+t*t*s.eye(3)/2)*B1)[1,1]).coeff(t,3)
 checks['center_cubic_coefficient']=s.expand(E-z[20]*(z[5]-s.Rational(3,10)))==0
 checks['polynomial_certificate']=s.expand(E-z[20]*q143-5*z[20]*z[34]*q19+s.Rational(3,10)*z[20])==0
 assert all(checks.values()),[k for k,v in checks.items() if not v]
 report={'scope':'Exact identities for conic rank-one substitution, leaf matrix reductions and two-Jacobi/one-center certificate; closed-point and nonzero hypotheses are explained in proof',
 'checks':checks,'count':len(checks),'all_passed':True}
 out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'count':len(checks),'all_passed':True}))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('old_model',type=Path);ap.add_argument('jet_model',type=Path);ap.add_argument('full_certificate',type=Path);ap.add_argument('output',type=Path)
 a=ap.parse_args();audit(a.old_model,a.jet_model,a.full_certificate,a.output)
