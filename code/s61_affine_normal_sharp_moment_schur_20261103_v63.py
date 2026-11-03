#!/usr/bin/env python3
"""Sharp532 controlled Schur determinant and135-generator scalar alternative."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def determinant(M):
 D=DomainMatrix.from_Matrix(M);return s.expand(D.domain.to_sympy(D.det()))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 deps=[ROOT/'reports/S61_AFFINE_NORMAL_MOMENT_SCHUR_20261103_v63/report.json',ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/report.json']
 for dep in deps:
  data=json.loads(dep.read_text());assert data['all_passed']
  for rel,want in data['input_sha256'].items():assert sha(ROOT/rel)==want,rel
  for rel,want in data.get('evidence_sha256',{}).items():assert sha(dep.parent/rel)==want,rel
 h,a=s.symbols('h a');fixtures=[]
 for name,G,h0,sources in [('simple',s.prod(h-i for i in range(1,6)),1,{0:h,1:h-1}),('multiple',(h-2)**3*(h-3)*(h+1),2,{0:h,1:s.S.One}),('sharp532',h**5,0,{0:h,133:s.S.One})]:
  G=s.expand(G);g4=s.div(G,h-h0,h)[0];S=sum(a**i*f for i,f in sources.items());r=s.expand((h-h0)*S);columns=[]
  for j in range(5):
   rem=s.Poly(s.rem(r*h**j,G,h),h);columns.append(s.Matrix([rem.coeff_monomial(h**k)for k in range(5)]))
  M=s.Matrix.hstack(*columns);delta=determinant(M[1:,:4]);expected=s.expand(s.resultant(g4,S,h));assert s.expand(delta-expected)==0 and delta!=0
  degree=s.Poly(delta,a).degree();assert degree<=532
  if name=='sharp532':assert delta==a**532 and degree==532
  fixtures.append(dict(name=name,G=str(G),h0=h0,quotient_S=str(S),multiplication5_by5=[[str(v)for v in row]for row in M.tolist()],delta0=str(delta),norm_D=str(expected),degree_a=degree))
 # A coefficient-ring fixture for the135-generator algebra reduction.
 x,b,c=s.symbols('x b c');G=h**5+c*h+b;r=h-x;cols=[]
 for j in range(5):
  rem=s.Poly(s.rem(r*h**j,G,h),h);cols.append(s.Matrix([rem.coeff_monomial(h**k)for k in range(5)]))
 M=s.Matrix.hstack(*cols);P4=M[1:,:4];delta=determinant(P4);assert delta==1
 numerator=s.expand(-(M[0,:4]*P4.adjugate())[0]);norm=determinant(M);assert numerator==x and s.expand(norm+G.subs(h,x))==0
 fixtures.append(dict(name='coefficient_ring_cyclic_algebra',G=str(G),r=str(r),delta0=str(delta),Cramer_h_numerator=str(numerator),Norm_C_r=str(norm),exact_cyclic_quotient='QQ[x,b,c,h]/(G,h-x)=QQ[x,b,c]/(x^5+c*x+b)'))
 (out/'fixtures.json').write_text(json.dumps(fixtures,indent=2)+'\n')
 recipe=dict(moment='r(a)=sum_i0..133 a^i*r_i',pivot_columns='b_j=sum_i a^i*Phi_column(5*i+j),j0..3',fixed_rows=[1,2,3,4],delta0='detbottom4bcolumns',
  exact_fiber_identity='delta0=Norm_(k[h]/(G/(h-h0)))(sum_i a^i*(r_i/(h-h0)))',annihilator_kernel='spanMONICdegree4G/(h-h0);poly<=3alwayscomplement',degree_a_bound=532,
  primitive_and_regular='coefficientsdelta0generateunitidealinB;delta0regularonB[a],B(a);B->B(a)_delta0injective',single_generic_unit_equivalent_to_B_empty=True,
  full_cover_options=dict(three_independent_generic_constants=True,requires_dimB_atmost2=True,fixed_QQ_constants=list(range(533)),fixed_cover_no_dimension_bound=True),
  degree10_option=dict(scalar_generators=670,formula='det[b0,b1,b2,b3,Phi_column_l]',formal_degree_atmost=10),
  degree34_option=dict(scalar_generators=135,norm='N=detM_r,all5columnsT^j*r(a)',Cramer='tau=U/delta0,U=-(toprow_b*adj(bottom4_b))_0',equations=['N','sum_k0..4 r_i,k*U^k*delta0^(4-k),i0..133'],exact_algebra='A(a)_delta0=Rprime_K(a)[delta0^-1]/(N,cleared134evaluations)',formal_degree_atmost=34,norm_formal_degree_atmost=10,U_and_delta_formal_degree_atmost=8),
  old536_construction_remains_valid=True,no_actual_scalar_unit_or_exhaustion=True)
 (out/'criterion.json').write_text(json.dumps(recipe,indent=2)+'\n')
 paths=[Path(__file__),*deps,ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/presentation.json.gz',ROOT/'research/notes/S61_AFFINE_NORMAL_SHARP_MOMENT_SCHUR_20261103_v63.md'];report=dict(id='S61_AFFINE_NORMAL_SHARP_MOMENT_SCHUR_20261103_v63',all_passed=True,sharp_degree532_norm_fixtures_verified=True,one_generic_regular_delta_exact_empty_criterion=True,exact135_scalar_algebra_recipe=True,cover_three_generic_or533_fixed=True,actual_units_not_computed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('SHARP532REGULARSCHUR/EXACT135SCALARALGEBRA PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
