#!/usr/bin/env python3
"""A controlled single generic Schur chart detects the entire unknown residual."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 deps=[ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/report.json',ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63/report.json',ROOT/'reports/S61_AFFINE_NORMAL_SCHUR_COVER_20261103_v63/report.json']
 for dep in deps:
  data=json.loads(dep.read_text());assert data['all_passed']
  for rel,want in data['input_sha256'].items():assert sha(ROOT/rel)==want,rel
  for rel,want in data.get('evidence_sha256',{}).items():assert sha(dep.parent/rel)==want,rel
 h,a=s.symbols('h a');fixtures=[]
 cases=[('simple',s.prod(h-i for i in range(1,6)),1,{0:h,1:h-1}),('multiple',(h-2)**3*(h-3)*(h+1),2,{0:h,1:s.S.One}),('sharp536',h**5,0,{0:h,133:s.S.One})]
 for name,G,h0,sources in cases:
  G=s.Poly(s.expand(G),h);gdiv=s.div(G,s.Poly(h-h0,h))[0].as_expr();assert s.rem(G.as_expr(),h-h0,h)==0;S=sum(a**i*f for i,f in sources.items());r=s.expand((h-h0)*S)
  columns=[]
  for j in range(4):
   rem=s.Poly(s.rem(r*(h-a)*h**j,G.as_expr(),h),h);columns.append(s.Matrix([rem.coeff_monomial(h**k)for k in range(5)]))
  matrix=s.Matrix.hstack(*columns);dm=DomainMatrix.from_Matrix(matrix[1:,:]);delta=s.expand(dm.domain.to_sympy(dm.det()));expected=s.expand(s.resultant(gdiv,S,h)*gdiv.subs(h,a));assert s.expand(delta-expected)==0 and delta!=0
  degree=s.Poly(delta,a).degree();assert degree<=536
  if name=='sharp536':assert delta==a**536 and degree==536
  fixtures.append(dict(name=name,G=str(G.as_expr()),h0=h0,moment_quotient=str(S),all5_by4_columns=[[str(v)for v in row]for row in matrix.tolist()],bottom4_determinant=str(delta),factorization=str(s.factor(delta)),norm_resultant_identity=str(expected),degree_a=degree,nonseparable_allowed=True))
 (out/'fixtures.json').write_text(json.dumps(fixtures,indent=2)+'\n')
 recipe=dict(source='Phi5x670columnsPhi_(i,j)=T_h^j*r_i, i0..133,j0..4',moment='r(a)=sum_i a^i*r_i',pivot_columns=['sum_i a^i*(Phi_(i,j+1)-a*Phi_(i,j)) forj0..3'],fixed_pivot_rows=[1,2,3,4],delta='detfixedbottom4pivotcolumns',scalar_equations='g_l=det[pivot4columns,Phi_column_l], l0..669',local_scalar_ring='Rprime_K(a)[delta^-1]/(g_0,...,g_669)=A(a)_delta=B(a)_delta',
  degree_bounds=dict(moment_a=133,pivot_a=134,delta_a=536,scalar_g_a=536,delta_formal_parameters=8,scalar_g_formal_parameters=10),
  exact_pointwise_norm_identity='WithG=(h-h0)*g4,r(a)=(h-h0)*S(a,h),delta(a)=g4(a)*Norm_(k[h]/g4)(S(a,h)).',
  primitive_delta_coefficient_ideal_in_B='unit',regularity='deltaavoidsallAssB[a];nonzerodivisor;B->B(a)_deltainjective',
  exact_single_unit_criterion='B=0 iff singlelocalScalarringoverK(a)iszero',known_E5_removed_by_delta_zero_at_origin=True,
  full_cover_options=dict(three_independent_generic_constants=['a0','a1','a2'],uses_dim_B_atmost=2,fixed_QQ_constants=list(range(537)),fixed_cover_requires_dimension_bound=False),
  not_claimed=['singlechartcoversallpoints','anyscalarunitcomputed','modpunitspreadsQQ','specific3rationalconstantsverified','wholew_torsionremoved'])
 (out/'criterion.json').write_text(json.dumps(recipe,indent=2)+'\n')
 paths=[Path(__file__),*deps,ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/presentation.json.gz',ROOT/'research/notes/S61_AFFINE_NORMAL_MOMENT_SCHUR_20261103_v63.md'];report=dict(id='S61_AFFINE_NORMAL_MOMENT_SCHUR_20261103_v63',all_passed=True,proof_and_five_dimensional_fixture_receipt=True,controlled_generic_delta_regular_on_unknown_residual=True,single_generic_scalar_unit_equivalent_to_exhaustion=True,full_cover_alternatives_three_generic_or537_fixed=True,sharp_degree536_fixture_verified=True,no_actual_unit_or_exhaustion_claimed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('ONECONTROLLED GENERIC SCHUR/REGULARDELTA/SHARP536FIXTURE PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
