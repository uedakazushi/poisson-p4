#!/usr/bin/env python3
"""Raw-purity obstruction and exact sheaf-colon split of isolated infinity."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_model import S61
from s61_known_fiber_audit import equation_jacobian
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 oldpath=ROOT/'reports/S61_MINIMAL_AXIS_INFINITY_20261101.json';old=json.loads(oldpath.read_text());assert old['all_passed']and old['all_135_identities']
 for rel,want in old['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 model=S61();two=old['weight_two_u_indices_zero_based'];assert two==[i for i,w in enumerate(model.reduced['parameter_weights'])if w==2];point=list(map(s.Rational,old['sample_S61_coordinates']));source=model.reduced['equations'];assert len(source)==134
 top=[[term for term in equation if sum(i in two for i in term[0])==2]for equation in source]
 assert all(model.evaluate(eq,point)==0 for eq in top)
 leading=equation_jacobian(top,point)[:,two];fixed=next(j for j,i in enumerate(two)if point[i]!=0);chart=[j for j in range(45)if j!=fixed]
 tcolumn=s.Matrix([model.evaluate([term for term in eq if sum(i in two for i in term[0])==1],point)for eq in source]);full=leading[:,chart].row_join(tcolumn)
 rows=old['independent_full_homogenized_minor']['rows'];assert len(rows)==45 and max(rows)<134;minor=DomainMatrix.from_Matrix(full[rows,:]).convert_to(s.QQ);det=minor.domain.to_sympy(minor.det());assert det==s.Rational(old['independent_full_homogenized_minor']['determinant'])!=0
 dominance=s.Matrix([[s.Rational(c)for c in row]for row in old['dominance_minor']['matrix']]);domdet=DomainMatrix.from_Matrix(dominance).convert_to(s.QQ).det();assert s.Rational(domdet)==s.Rational(old['dominance_minor']['determinant'])==569
 dimpath=ROOT/'reports/S61_AFFINE_BOUNDARY_DIMENSION_20261102.json';dim=json.loads(dimpath.read_text());assert dim['all_passed']and dim['raw_generic_infinity_exact_dimension']==3
 certificate=dict(original134_rows_only=True,all134_homogenized_equations_zero_at_axis_sample=True,
  q_coordinate_fixed_for_projective_chart_zero_based=fixed,deleted_full_coordinate_u_one_based=two[fixed]+1,
  full45_minor_rows=rows,full45_minor_columns_q_zero_based_then_homogenizer=chart+['w'],full45_QQ_determinant=str(det),beta_dominance16_QQ_determinant=str(domdet),
  global_delta_recipe='Homogenize each original134 to qdegree2 usingw, take45Jacobian rows selected and44qcolumns excludingfixed plusw;Delta homogeneousdegree45.',
  raw_Y134_has_generic_reduced_isolated_infinity_component=True,raw_Y135_contains_generic_dimension3_boundary=True,raw_purity3_and_pure3_CM_strategy_false=True,
  projective_sheaf_ideal_recipe=dict(I='actual divided-normal compactincidence Y, all134 retained',w='original homogenizer',Delta='saved original134 globalJacobian45minor pulledback to Y',isolated_infinity_ideal='I_inf=(I+(w)):Delta^infinity',complement_ideal='I_rest=I:I_inf',scheme_equalities=['I=I_inf intersect I_rest','I_inf+I_rest=O_ambient as sheafideals','Y=Y_inf disjointunion Y_rest','Y_rest intersectD(w)=Y intersectD(w)']),
  finite_etale_isolated_infinity_removed_only=True,colon_Groebner_basis_or_degree_of_etale_factor_not_computed=True,affine_unknown_or_nonreduced_support_preserved=True,
  not_claimed=['rawYpure3','remainingYpure3','fullw_torsionremoved','affineresidualempty'])
 (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 paths=[Path(__file__),oldpath,dimpath,ROOT/'code/s61_model.py',ROOT/'code/s61_known_fiber_audit.py',model.source/'symplectic_slice_reduced_equations.json',ROOT/'research/notes/S61_AFFINE_RAW_PURITY_SPLIT_20261103_v63.md']
 report=dict(id='S61_AFFINE_RAW_PURITY_SPLIT_20261103_v63',all_passed=True,fresh_original134_full45_QQ_minor_replayed=True,fresh_beta_dominance16_minor_replayed=True,raw_purity3_strategy_refuted=True,exact_finite_sheaf_colon_split_recipe=True,affine_Dw_unchanged=True,colon_basis_or_total_w_saturation_not_claimed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('RAW PURITY3 FALSE; ORIGINAL134 AXISMINOR; EXACTSHEAFCOLON SPLIT PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
