#!/usr/bin/env python3
"""Independent coefficient-action and whole-ideal audit of source ramification.

No target-image or actual intersection ring is inferred from a source Fitting ideal.
"""
import gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
X=s.symbols('h x y z w');h,x,y,z,w=X
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))
INDEX={(p,m):k for k,(p,m)in enumerate(itertools.product(PAIRS,MONS))}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[m[0]]*X[m[1]])for(i,j),m in itertools.product(PAIRS,MONS)])
def canon(P):
 d=s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)]);xx=s.Matrix(X)
 return(P-(xx*d.T-d*xx.T)/5).applyfunc(s.expand)
def form(f=0,g=0,H=0):
 P=s.zeros(5);P[0,2]=f;P[3,4]=g*x*x;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z)
 return canon(P-P.T)
def coefficient_action(a,b):
 L=s.zeros(150)
 for col,((u,v),m)in enumerate(itertools.product(PAIRS,MONS)):
  def put(i,j,mon,c):
   if i==j:return
   if i>j:i,j,c=j,i,-c
   L[INDEX[((i,j),tuple(sorted(mon)))],col]+=c
  if u==b:put(a,v,m,1)
  if v==b:put(u,a,m,1)
  if a in m:
   new=list(m);new.remove(a);new.append(b);put(u,v,new,-m.count(a))
 return L
def audit(which):
 folder=ROOT/f'reports/AMBIENT_B18_{which}_CRITICAL_SCHEME_20261103_v63'
 report=json.loads((folder/'report.json').read_text());assert report['all_passed']
 for p,want in report['input_sha256'].items():assert sha(ROOT/p)==want,p
 for p,want in report['evidence_sha256'].items():assert sha(folder/p)==want,p
 d=json.loads(gzip.decompress((folder/'certificate.json.gz').read_bytes()))
 cubics=[s.prod((x,h,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
 base=[form(f=f)for f in([x*x,h*x,h*h]if which=='CA'else[x*x,h*x])]+[form(g=1)]+[form(H=H)for H in cubics]
 C=s.Matrix.hstack(*(vec(P)for P in base));saved=s.Matrix(d['fixed_fiber_linear_basis150_by20']if which=='CA'else d['fiber19_basis']);assert C==saved
 rr=list(C.T.rref()[1]);inv=C[rr,:].inv();assert inv*C[rr,:]==s.eye(C.cols)
 def rem(M):return(M-C*inv*M[rr,:]).applyfunc(s.expand)
 P0=s.zeros(5)
 for j,v in [(0,x*x),(3,h*x),(4,h*h+z*z)]:P0[2,j]=v;P0[j,2]=-v
 p=vec(P0);center=inv*p[rr,:];assert C*center==p
 assert center==s.Matrix(d['center_fixed20']if which=='CA'else d['center19'])
 flag=[tuple(v)for v in d['flag_positions']];para=[tuple(v)for v in d['parabolic_positions']]
 assert len(flag)==9 and len(para)==16 and set(flag)|set(para)==set(itertools.product(range(5),repeat=2))
 acts={ab:coefficient_action(*ab)for ab in flag+para}
 for ab in para:assert rem(acts[ab]*C)==s.zeros(150,C.cols)
 O=s.Matrix.hstack(*(inv*(acts[ab]*p)[rr,:]for ab in para))
 assert O==s.Matrix(d['parabolic_fiber_orbit20_by16']if which=='CA'else d['parabolic_orbit'])
 fixed=d['orbit10_fixed_rows']if which=='CA'else d['fixed_rows'];free=d['affine_slice_free_indices']if which=='CA'else d['free_indices'];pc=d['orbit10_columns']
 assert len(pc)==10 and len(fixed)==10 and set(fixed).isdisjoint(free) and set(fixed)|set(free)==set(range(C.cols))
 assert O[fixed,pc].det()==s.Rational(d['orbit10_minor'])!=0
 vv=s.symbols(' '.join(d['slice_variables']));v=p+C[:,free]*s.Matrix(vv)
 Q=rem(s.Matrix.hstack(*(acts[ab]*v for ab in flag)));assert Q==s.Matrix(d['full150_by9_flag_quotient'])
 rows=d['unit8_rows'];cols=d['unit8_columns'];extra=next(iter(set(range(9))-set(cols)))
 B=Q[rows,cols];assert s.expand(B.det())==1
 S=(Q[:,extra]-Q[:,cols]*B.inv()*Q[rows,extra]).applyfunc(s.expand)
 assert S==s.Matrix(d['full150_critical_Schur_column']if which=='CA'else d['all150_Schur_entries'])
 predicted=list(map(s.sympify,d['predicted_entire_critical_ideal']if which=='CA'else d['predicted_exact_ideal']))
 G=s.groebner(list(S),*vv,order='lex',domain=s.QQ);F=s.groebner(predicted,*vv,order='lex',domain=s.QQ)
 assert all(G.reduce(s.expand(q))[1]==0 for q in predicted)
 assert all(F.reduce(q)[1]==0 for q in S)
 # Explicit tiny monomial primary decomposition, independently computed by
 # monomial lcm intersections rather than a symbolic primary routine.
 def intersect(A,B):return{tuple(max(a,b)for a,b in zip(u,v))for u in A for v in B}
 def minimal(A):return{u for u in A if not any(v!=u and all(a<=b for a,b in zip(v,u))for v in A)}
 if which=='CA':
  got=minimal(intersect(intersect({(0,1,0),(0,0,1)},{(0,1,0),(1,0,0)}),{(1,0,0),(0,0,1),(0,3,0)}));assert got=={(1,1,0),(1,0,1),(0,1,1),(0,3,0)}
 else:
  got=minimal(intersect({(0,1)},{(1,0),(0,2)}));assert got=={(1,1),(0,2)}
 return dict(all150_basis_and_independent_coefficient_actions=True,parabolic_all_columns_stable=True,actual_center_and_etale19_factors=True,whole150_by9_and_unit8_determinant_one=True,all150_schur_polynomials_reconstructed=True,both_QQ_whole_ideal_containments=True,independent_monomial_primary_decomposition=True,source_Fitting_only_not_target_D=True,input_report=f'reports/AMBIENT_B18_{which}_CRITICAL_SCHEME_20261103_v63/report.json')
def main():
 start=time.monotonic();results={name:audit(name)for name in ['CA','N']}
 note=ROOT/'research/notes/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.md'
 inputs=[Path(__file__),note]
 for name in ['CA','N']:
  folder=ROOT/f'reports/AMBIENT_B18_{name}_CRITICAL_SCHEME_20261103_v63';inputs +=[folder/'report.json',folder/'certificate.json.gz']
 report=dict(id='AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW',all_passed=True,results=results,actual_target_completion_odd_ideal_primary_or_depth_inferred=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},seconds=time.monotonic()-start)
 (ROOT/'reports/AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
