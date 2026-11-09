"""Exact independent audits of the leading-center identity and line spectra.

The center identity is reconstructed from all basis-pair derivation equations
for the Heisenberg algebra. The line spectra are enumerated with Fractions.
"""
import json
from fractions import Fraction as Q
from pathlib import Path
import argparse
ap=argparse.ArgumentParser();ap.add_argument("--output",type=Path,required=True);args=ap.parse_args()
from itertools import product
import sympy as s

checks={}
d=s.symbols('d0:9'); D=s.Matrix(3,3,d); q,a=s.symbols('q a')
E=[s.eye(3)[:,i] for i in range(3)]
def bracket(x,y):return s.Matrix([0,0,x[0]*y[1]-x[1]*y[0]])
# Inputs minus output = q*[ , ]; this is the leading bracket weight k*alpha.
eq=[]
for i,j in product(range(3),repeat=2):
 eq.extend(bracket(D*E[i],E[j])+bracket(E[i],D*E[j])-D*bracket(E[i],E[j])-q*bracket(E[i],E[j]))
eq=[s.expand(x) for x in eq if x!=0]
G=s.groebner(eq,*d,q,domain=s.QQ)
checks['center_invariant_first_component']=G.reduce(D[0,2])[1]==0
checks['center_invariant_second_component']=G.reduce(D[1,2])[1]==0
checks['trace_minus_twice_center_equals_conformal_weight']=G.reduce(s.trace(D)-2*D[2,2]-q)[1]==0
GG=s.groebner(eq+[s.trace(D)+a],*d,q,a,domain=s.QQ)
checks['center_eigenvalue_equals_minus_half_slope_plus_weight']=GG.reduce(2*D[2,2]+a+q)[1]==0

def spectrum(h,m):return tuple(sorted([Q(1),-Q(h,2),-Q(m,3),Q(h,2)+Q(m,3)-1]))
types=list(product([1,2],[0,1,2]))
matches=[]
for A,B in product(types,repeat=2):
 SA=spectrum(*A);SB=spectrum(*B)
 # Since 1 belongs to SB, each possible scale must itself belong to SA.
 for ratio in set(SA)-{Q(0)}:
  if tuple(sorted(ratio*x for x in SB))==SA:
   matches.append({'left':A,'right':B,'ratio':str(ratio)})
checks['all_matches_classified']=all((m['ratio']=='1' and m['left']==m['right']) or
 (m['ratio']=='-1' and m['left']==m['right'] and m['left'][0]==2) for m in matches)
checks['exactly_nine_matches']=len(matches)==9
checks['tangent_eigenvalue_always_simple']=all(spectrum(*A).count(Q(1))==1 for A in types)
checks['zero_slope_not_hidden_in_nonzero_spectrum']=all(any(x!=0 for x in spectrum(*A)) for A in types)
checks['h2_spectra_exclude_other_tangent_ratio_two']=all(Q(2) not in spectrum(2,m) for m in range(3))
checks['line_conic_opposite_slopes_excluded']=all(spectrum(2,m)!=tuple(sorted([Q(1),Q(-1),Q(1,2),Q(-1,2)])) for m in range(3))

# Cubic projective obstruction: E wedge Z2 has the same weight as Z2.
weights=[3,-3,-2,2]
monomial2=[tuple(int(k==i)+int(k==j) for k in range(4)) for i in range(4) for j in range(i,4)]
resonant=[]
for powers in monomial2:
 for output in range(4):
  weight=sum(a*b for a,b in zip(powers,weights))-weights[output]
  if weight==0:resonant.append((powers,output))
checks['m2_no_quadratic_vectorfield_resonances']=not resonant
checks['m2_forced_cube_monomial_x2z3_weight_zero']=2*weights[0]+3*weights[2]==0
checks['m2_forced_cube_monomial_y2w3_weight_zero']=2*weights[1]+3*weights[3]==0

assert all(checks.values())
out={'checks':checks,'checks_count':len(checks),'all_passed':True,'spectral_matches':matches,
 'scope':'Exact local linear-algebra and finite rational spectrum checks; geometric extension arguments remain in the proof text'}
args.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
