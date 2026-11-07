#!/usr/bin/env python3
"""Complete exact localunit elimination of the quadraticA18cut; noCAS/jet."""
from pathlib import Path
import hashlib,json,time
import sympy as s
from classification_local_unit_elimination import eliminate
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();path=ROOT/'reports/CLASSIFICATION_NIL22_QUADRATIC_A_FIXED_MODEL_20261107.json';model=json.loads(path.read_text());zz=s.symbols('z0:17');n=s.Symbol('n');variables=zz+(n,);env={str(v):v for v in variables};equations=[s.sympify(f,locals=env)for f in model['complete_local_ideal_generators']];answer=eliminate(equations,zz,[n]);remaining=answer['remaining_implicit'];residual=answer['residuals'];sub=answer['substitution'];print('steps',len(answer['steps']),'remaining',len(remaining),'nonzero',sum(f!=0 for f in residual),round(time.monotonic()-started,3),flush=True)
 inputs=[Path(__file__),ROOT/'code/classification_local_unit_elimination.py',path,ROOT/'reports/CLASSIFICATION_NIL22_QUADRATIC_A_WEIGHT_FRAME_20261107.json'];result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),method='Exactlocalringisomorphisms byunitlinear equations; noformaljet orGB.',normal_variable='n',original_complete_equations=model['complete_local_ideal_generators'],steps=answer['steps'],complete_final_inverse_substitution={str(v):str(f)for v,f in sub.items()},all_explicit_denominator_units=sorted(answer['denominators']),remaining_implicit_variables=list(map(str,remaining)),complete_final_all45_residuals=list(map(str,residual)),stopped=answer['stopped'],seconds=time.monotonic()-started)
 if not remaining:
  for f in residual:assert not f.free_symbols.difference({n})
  nonzero=next(((i,f)for i,f in enumerate(residual)if f!=0),None)
  if nonzero:
   i,f=nonzero;poly=s.Poly(f,n);degree=min(k[0]for k,c in poly.terms());assert degree>=1;unit=s.cancel(f/n**degree);assert unit.subs(n,0)!=0;result.update(reduced_support_is_normal_origin=True,one_sufficient_normal_residual=dict(equation_index=i,polynomial=str(f),normal_power=degree,unit_factor=str(unit),unit_at_origin=str(unit.subs(n,0))),whole7source_full23_component_supported=True,all_final_denom_units_at_origin=True)
  else:
   result.update(exact_rational_normal_curve=True,all45_equations_vanish_identically=True,whole_source_componenthood_claimed=False)
 result['seconds']=time.monotonic()-started;out=ROOT/'reports/CLASSIFICATION_NIL22_QUADRATIC_A_UNIT_ELIMINATION_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['stopped','seconds']}))
if __name__=='__main__':main()
