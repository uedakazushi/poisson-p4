"""Finite exact local-ring elimination by variables having unit linear coefficient.
No formal truncation, standard basis, sample point, or nilpotence index search.
"""
import sympy as s

def eliminate(equations,implicit,retained,max_terms=100000):
 original=list(equations);current=list(equations);implicit=list(implicit);retained=list(retained);sub={};steps=[];denominators=set()
 def normalize(f,variables):
  num,den=s.fraction(s.cancel(f));assert den.subs(dict.fromkeys(variables,0))!=0
  if den!=1:denominators.add(str(den))
  return s.expand(num)
 current=[normalize(f,implicit+retained)for f in current]
 while implicit:
  variables=implicit+retained;candidates=[]
  for i,f in enumerate(current):
   if f==0:continue
   P=s.Poly(f,*variables)
   if len(P.terms())>max_terms:return dict(steps=steps,substitution=sub,remaining_implicit=implicit,retained=retained,residuals=current,denominators=denominators,stopped='termcap')
   for j,z in enumerate(implicit):
    if P.degree(z)!=1:continue
    coefficient=s.diff(f,z);value=coefficient.subs(dict.fromkeys(variables,0))
    if value==0:continue
    score=(int(P.total_degree()>1),int(coefficient.free_symbols!=set()),len(P.terms()),len(s.Poly(coefficient,*variables).terms()),j,i);candidates.append((score,i,z,coefficient,value))
  if not candidates:break
  _,i,z,coefficient,value=min(candidates,key=lambda c:c[0]);remainder=s.expand(current[i]-coefficient*z);assert not remainder.has(z);solution=s.cancel(-remainder/coefficient);assert solution.subs(dict.fromkeys(variables,0))==0;before=current[i];sub={t:s.cancel(f.subs(z,solution))for t,f in sub.items()};sub[z]=solution;implicit.remove(z);current=[normalize(f.subs(z,solution),implicit+retained)for f in current]
  steps.append(dict(equation_index=i,variable=str(z),equation_before=str(before),unit_coefficient=str(coefficient),unit_at_origin=str(value),remainder=str(remainder),rational_solution=str(solution),remaining_implicit=list(map(str,implicit)),all_residual_numerators_after=list(map(str,current))))
 return dict(steps=steps,substitution=sub,remaining_implicit=implicit,retained=retained,residuals=current,denominators=denominators,stopped='complete'if not implicit else'no_unitlinear_generator')
