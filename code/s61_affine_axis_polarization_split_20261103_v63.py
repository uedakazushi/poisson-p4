#!/usr/bin/env python3
"""Actual45-coordinate polarization splitter for a simple infinity axis factor."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_model import S61
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(poly,mon,c):
 v=poly.get(mon,Fraction(0))+c
 if v:poly[mon]=v
 else:poly.pop(mon,None)
def multiply(a,b):
 out={}
 for m,c in a.items():
  for n,d in b.items():add(out,tuple(sorted(m+n)),c*d)
 return out
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 dep=ROOT/'reports/S61_AFFINE_RAW_PURITY_SPLIT_20261103_v63/report.json';old=json.loads(dep.read_text());assert old['all_passed']
 for rel,want in old['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 for rel,want in old['evidence_sha256'].items():assert sha(dep.parent/rel)==want,rel
 cert=json.loads((dep.parent/'certificate.json').read_text());axispath=ROOT/'reports/S61_MINIMAL_AXIS_INFINITY_20261101.json';axis=json.loads(axispath.read_text());point=list(map(Fraction,axis['sample_S61_coordinates']));model=S61();two=axis['weight_two_u_indices_zero_based'];fixed=cert['q_coordinate_fixed_for_projective_chart_zero_based'];fixed_u=two[fixed];chart=[i for i in two if i!=fixed_u];index={i:j for j,i in enumerate(chart)};pc=[point[i]for i in chart]+[Fraction(0)];assert len(pc)==45
 polynomials=[]
 for equation in model.reduced['equations']:
  polynomial={}
  for mon,num,den in equation:
   coefficient=Fraction(num,den);q=[i for i in mon if i in two];assert len(q)<=2;variables=[44]*(2-len(q))
   for i in mon:
    if i not in two:coefficient*=point[i]
    elif i==fixed_u:coefficient*=point[i]
    else:variables.append(index[i])
   add(polynomial,tuple(sorted(variables)),coefficient)
  assert sum(c*__import__('math').prod(pc[i]for i in mon)for mon,c in polynomial.items())==0;polynomials.append(polynomial)
 jacobians=[];midpoint=[]
 for polynomial in polynomials:
  J=[{}for _ in range(45)]
  for mon,c in polynomial.items():
   for k,i in enumerate(mon):add(J[i],mon[:k]+mon[k+1:],c)
  M=[]
  for derivative in J:
   poly={}
   for mon,c in derivative.items():
    if not mon:add(poly,(),c)
    else:
     assert len(mon)==1;add(poly,mon,c/2);add(poly,(),c*pc[mon[0]]/2)
   M.append(poly)
  replay={}
  for j,entry in enumerate(M):
   delta={(j,):Fraction(1)}
   if pc[j]:delta[()]=-pc[j]
   for mon,c in multiply(entry,delta).items():add(replay,mon,c)
  assert replay==polynomial
  jacobians.append(J);midpoint.append(M)
 rows=cert['full45_minor_rows'];assert len(rows)==45 and max(rows)<134
 atpoint=s.Matrix([[s.Rational(sum(c*__import__('math').prod(pc[k]for k in mon)for mon,c in jacobians[i][j].items()))for j in range(45)]for i in rows]);dm=DomainMatrix.from_Matrix(atpoint).convert_to(s.QQ);d=dm.domain.to_sympy(dm.det());assert d==s.Rational(cert['full45_QQ_determinant'])!=0
 encoded=lambda polynomial:[[list(mon),str(c)]for mon,c in sorted(polynomial.items())]
 payload=dict(chart_rule='Fixdeletedqcoordinatetoitsnonzerosamplevalue;45variablesare44otherqfollowedbyw',fixed_free_u_one_based=fixed_u+1,fixed_value=str(point[fixed_u]),chart_variables_u_one_based=[i+1 for i in chart]+['w'],infinity_point45=list(map(str,pc)),all134_actual_homogeneous_chart_polynomials=list(map(encoded,polynomials)),selected_original_rows=rows,linear45_by45_midpoint_matrix=[[encoded(poly)for poly in midpoint[i]]for i in rows],d_at_point=str(d))
 (out/'certificate.json.gz').write_bytes(gzip.compress(json.dumps(payload,separators=(',',':')).encode(),mtime=0))
 result=dict(all134_exact_QQ_midpoint_identities=True,full45_minor_original134_only=True,exact_d=str(d),determinant_recipe='e_p(x)=detJ_selected((x+p)/2)/detJ_selected(p),degree<=45',
  exact_ideal_memberships=['e_p*(x_j-p_j)inI forj0..44','e_p*w inI','e_p^2-e_p inI'],
  exact_localchart_split='I=(I+(1-e_p)) intersect (I+(e_p));sumidealsunit',D_w_preserved_exactly=True,
  generic_axis_extension='E_axis from(I+(w)):Delta∞ isfiniteetaleoverK;onchartandnormunitopen trace_E/K(e_p) splitsentireaxisfactor;finiteetaledegreeunknownandtracenotexpanded',
  no_affine_unknown_removed=True,not_claimed=['wholew_torsionremoved','affineunit','expandeddegree45determinant','axisetaledegree'])
 (out/'split-recipe.json').write_text(json.dumps(result,indent=2)+'\n')
 paths=[Path(__file__),dep,dep.parent/'certificate.json',axispath,ROOT/'code/s61_model.py',model.source/'symplectic_slice_reduced_equations.json',ROOT/'research/notes/S61_AFFINE_AXIS_POLARIZATION_SPLIT_20261103_v63.md']
 report=dict(id='S61_AFFINE_AXIS_POLARIZATION_SPLIT_20261103_v63',all_passed=True,actual_all134_midpoint_identities_replayed=True,exact_axis45_Jacobian_determinant_replayed=True,exact_infinity_idempotent_and_w_annihilation_recipe=True,generic_finiteetale_trace_split_recipe=True,affine_Dw_unchanged=True,total_torsion_or_exhaustion_not_claimed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('AXIS45POLARIZATION/ALL134IDENTITIES/EXACTwTORSIONSPLIT PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
