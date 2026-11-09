"""Exact finite algebra checks for the line-line obstruction. SymPy required."""
from fractions import Fraction as Q
import itertools,json
from pathlib import Path
import argparse
ap=argparse.ArgumentParser();ap.add_argument("--output",type=Path,required=True);args=ap.parse_args()
import sympy as s
checks=[]
def check(label,truth):
 assert truth,label
 checks.append(label)
sets={(h,m):[Q(1),-Q(h,2),-Q(m,3),Q(h,2)+Q(m,3)-1] for h in [1,2] for m in [0,1,2]}
matched=[]
for A,sa in sets.items():
 check('tangent_multiplicity_'+str(A),sa.count(1)==1)
 for B,sb in sets.items():
  for ratio in set(sa)-{0}:
   if sorted(sa)==sorted(ratio*v for v in sb):
    check('spectral_match_'+str((A,B,ratio)),A==B and (ratio==1 or (ratio==-1 and A[0]==2)))
    matched.append([list(A),list(B),str(ratio)])
x,y,z,w=s.symbols('x y z w');v=(x,y,z,w)
a,b,c,d,e,f,g,h,i,j=s.symbols('a b c d e f g h i j')
params=(a,b,c,d,e,f,g,h,i,j)
def mons(n):
 return [s.prod(v[k]**ee[k] for k in range(4)) for ee in itertools.product(range(n+1),repeat=4) if sum(ee)==n]
def weight(mon,weights):
 ee=s.Poly(mon,*v).monoms()[0];return sum(ee[k]*weights[k] for k in range(4))
for m,weights in [(1,(3,-3,-1,1)),(2,(3,-3,-2,2))]:
 expected={(0,1):[x*y,z*w],(0,2):[x*z]+([w*w] if m==1 else []),(0,3):[x*w],(1,2):[y*z],(1,3):[y*w]+([z*z] if m==1 else []),(2,3):[x*y,z*w]}
 for I,J in itertools.combinations(range(4),2):
  actual=[mon for mon in mons(2) if weight(mon,weights)==weights[I]+weights[J]]
  check(f'quadratic_weight_template_{m}_{I}{J}',set(actual)==set(expected[I,J]))
 if m==2:
  for I in range(4):check(f'no_quadratic_vector_resonance_{I}',not [mon for mon in mons(2) if weight(mon,weights)==weights[I]])
  for I,J in itertools.combinations(range(4),2):
   actual=[mon for mon in mons(3) if weight(mon,weights)==weights[I]+weights[J]]
   expected_c={(0,2):[y*w*w],(1,3):[x*z*z]}.get((I,J),[])
   check(f'cubic_weight_template_{I}{J}',set(actual)==set(expected_c))
 def tensor(m):
  XZ=c*x*z+d*(w*w if m==1 else y*w*w)
  YW=g*y*w+h*(z*z if m==1 else x*z*z)
  return s.Matrix([[0,a*x*y+b*z*w,XZ,e*x*w],[-a*x*y-b*z*w,0,f*y*z,YW],[-XZ,-f*y*z,0,i*x*y+j*z*w],[-e*x*w,-YW,-i*x*y-j*z*w,0]])
 P=tensor(m)
 Jac={}
 for I,J,K in itertools.combinations(range(4),3):
  Jac[I,J,K]=s.Poly(s.expand(sum(P[I,L]*s.diff(P[J,K],v[L])+P[J,L]*s.diff(P[K,I],v[L])+P[K,L]*s.diff(P[I,J],v[L]) for L in range(4))),*v)
 if m==1:
  check('m1_bi',s.expand(Jac[0,1,2].coeff_monomial(x*y*z)-b*i)==0)
  check('m1_dh_obstruction',s.expand(Jac[0,1,2].coeff_monomial(w*z*z).subs(b,0)+2*d*h)==0)
  F=s.Poly(s.expand(P[0,1]*P[2,3]-P[0,2]*P[1,3]+P[0,3]*P[1,2]),*v)
  check('m1_cubic_x',F.coeff_monomial(x*z**3)==-c*h)
  check('m1_cubic_y',F.coeff_monomial(y*w**3)==-d*g)
 else:
  # The corrected m2 proof uses the radial cubic identity, not a condition on xyzw.
  P2=P.subs({d:0,h:0})
  quad=mons(2);u=s.symbols('u0:40')
  Z2=[sum(u[10*I+k]*mon for k,mon in enumerate(quad)) for I in range(4)]
  Rad=s.Matrix(4,4,lambda I,J:(v[I]*Z2[J]-v[J]*Z2[I])/5)
  PF5=s.Poly(s.expand(P2[0,1]*Rad[2,3]+Rad[0,1]*P2[2,3]-P2[0,2]*Rad[1,3]-Rad[0,2]*P2[1,3]+P2[0,3]*Rad[1,2]+Rad[0,3]*P2[1,2]),*v)
  check('m2_radial_no_x2z3',PF5.coeff_monomial(x*x*z**3)==0)
  check('m2_radial_no_y2w3',PF5.coeff_monomial(y*y*w**3)==0)
  check('m2_radial_no_any_weight0',all(weight(s.prod(v[I]**mm[I] for I in range(4)),weights)!=0 for mm in PF5.monoms()))

report={'all_passed':True,'checks':len(checks),'check_labels':checks,'spectral_matches':matched,'scope':'Exact algebra only; geometric line-flag proof in companion note'}
args.output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'all_passed':True,'checks':len(checks),'spectral_matches':matched},indent=2))
