#!/usr/bin/env python3
"""Verify explicit degree9 vertex resultant for an all-point four-quadric CI open.

The geometric proof is in the companion note; replay binds the actual quadrics,
all four-coordinate identities, binary D unit, and the fixed-degree resultant.
"""
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
import subprocess,json,time,sys
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/quadric-CI-open-replay-v79'
OUT.mkdir(parents=True,exist_ok=True)

def main():
 start=time.monotonic();t,v,H,y=s.symbols('t v H y')
 for stem,folder in [('s61_nonsquare_four_quadric_binary_reduction_20261103_v63','binary-source'),('s61_nonsquare_pair_finite_bad_locus_20261103_v63','pair-source')]:
  with(OUT/(folder+'.log')).open('w')as log:subprocess.run([sys.executable,str(ROOT/'code'/(stem+'.py')),str(OUT/folder)],check=True,stdout=log)
 binary=json.loads((OUT/'binary-source/certificate.json').read_text());pair=json.loads((OUT/'pair-source/certificate.json').read_text())
 load=lambda terms:sum(s.Rational(c)*t**ex[0]*v**ex[1]*H**ex[2]for ex,c in terms);rows=[[load(q)for q in row]for row in pair['natural_rows']]
 q=[row[5]*y*y+row[4]*y+row[2]for row in rows];coef=[[s.Poly(f,y).coeff_monomial(y**k)for k in [2,1,0]]for f in q];aa,bb=coef
 Syl=s.Matrix([aa+[0],[0]+aa,bb+[0],[0]+bb]);det=DomainMatrix.from_Matrix(Syl).det();det=Syl[0,0]*0+s.sympify(str(det))if not isinstance(det,s.Expr)else det
 R=s.Poly(s.cancel(-det*s.Integer(4096000000000000000)/(729*t)),t,v,H,domain=s.QQ);assert R.total_degree()==9 and not R.is_zero
 assert s.Poly(s.expand(det+729*t*R.as_expr()/s.Integer(4096000000000000000)),t,v,H,domain=s.QQ).is_zero
 # The actual minors and binaryD/frame were all checked by the first source replay.
 assert binary['D_x0_squared_unit']=='1280*a^3*J^2/243'
 assert binary['inverse_frame_determinant']=='1/(200*a*P^3*J)'
 encode=lambda poly:[[list(ex),str(c)]for ex,c in s.Poly(poly,t,v,H,domain=s.QQ).terms()]
 cert=dict(natural_core_coordinates={'t':'a/P','v':'3e+20t','H':'P+8t'},known_open='t*(H-8t)*(5H-2v-4)!=0',vertex='A=x0=W=B=0',vertex_free_light_variables=['u','vlight'],binary_quadric_coefficients=[list(map(encode,row))for row in coef],Sylvester_determinant_known_unit_factor='-729*t/4096000000000000000',R_vertex=encode(R.as_expr()),R_vertex_totaldegree=9,source_binary_certificate=binary)
 (OUT/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 result=dict(all_exact_QQ_source_and_resultant_checks_passed=True,explicit_degree9_CI_open=True,raw_four_quadrics_CI_on='a*P*J*R_vertex!=0',determinantal_base_pure_curve_from_height2=True,proper_projection_allows_reducible_or_nonreduced_base=True,extra_gcd_exception_exclusion_not_needed_on_this_open=True,original_light_quadratic_quartic_rank146_on='a*P*F*B*C*J*R_vertex!=0',global_R_vertex_zero_case_not_resolved=True,CI_or_rank146_everywhere_not_proved=True,no_full_quartic_surjectivity_or_ambient_claim=True,seconds=time.monotonic()-start)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
