"""Exact model audit for the fourfold modular-slope reduction.

This does not use a Groebner basis or claim the normalized full Jacobi ideal
is the unit ideal. It verifies the geometry-facing linear specialization.
"""
from pathlib import Path
import argparse,json
import sympy as s

ap=argparse.ArgumentParser();ap.add_argument('model',type=Path);ap.add_argument('output',type=Path);args=ap.parse_args()
d=json.loads(args.model.read_text())
m=d['linear_model']; x=s.symbols('x0:5'); z=s.symbols('z0:48')
local={str(a):a for a in x+z}
slope_sub={z[42]:0,z[44]:0,z[46]:0}
active_sub={z[41]:1,z[43]:1,z[45]:1,z[47]:1}
P={k:s.expand(s.sympify(v,locals=local).subs(slope_sub)) for k,v in m['P_entries'].items()}
checks={f'P{k}_independent_x0':s.diff(v,x[0])==0 for k,v in P.items()}
res=[s.expand(s.sympify(v,locals=local).subs(slope_sub)) for v in m['residue_expressions']]
checks['all_four_linear_slopes_zero']=all(res[i]==0 for i in (0,2,4,6))
checks['all_four_quadratic_coefficients_free']=res[1::2]==[z[41],z[43],z[45],z[47]]
eqs=[]
for e in d['jacobi_equations']:
 p=s.Poly(s.expand(s.sympify(e,locals=local).subs(slope_sub).subs(active_sub)),*[a for a in z if a not in slope_sub and a not in active_sub],domain=s.QQ)
 if p.as_expr()!=0:
  e=p.monic().as_expr()
  if e not in eqs:eqs.append(e)
normalized={k:str(s.expand(v.subs(active_sub))) for k,v in P.items()}
output={'scope':'Independent exact verification of slope-zero and active-residue normalization; not a full Jacobi decomposition',
 'checks':checks,'all_passed':all(checks.values()),'affine_slice_parameter_count':41,
 'remaining_distinct_equations':len(eqs),'P_entries':normalized,'jacobi_equations':[str(e) for e in eqs]}
args.output.write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({k:v for k,v in output.items() if k not in ('P_entries','jacobi_equations')},indent=2))
