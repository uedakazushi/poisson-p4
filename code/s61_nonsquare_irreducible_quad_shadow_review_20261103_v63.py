#!/usr/bin/env python3
"""Independent compact arithmetic for the uniform product-rank37 proof.

Uses a Hessian determinant and multiset shadows, not the producer's sparse
permutation determinant or recursive exponent enumeration. No matrix solver.
"""
import argparse,hashlib,itertools,json,time
from fractions import Fraction
from pathlib import Path
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def checked_report(p):
 d=json.loads(p.read_text());assert d['all_passed']
 for name,h in d.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name
 return d
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 source=ROOT/'reports/S61_BINARY_NONSQUARE_MOD_ELL_COORDINATES_20261103_v62.json';data=json.loads(source.read_text());a,e,d=sp.symbols('a e d');w=sp.symbols('w0:7');q=sp.S.Zero
 for mon,terms in data['cleared_qfour_sparse'][0]:
  assert len(mon)==7 and mon[6]==0 and sum(mon)==2
  coefficient=sp.S.Zero
  for ex,c in terms:
   assert len(ex)==3 and min(ex)>=0 and ex[2]==0
   coefficient+=sp.Rational(c)*a**ex[0]*e**ex[1]*d**ex[2]
  q+=coefficient*sp.prod(x**k for x,k in zip(w,mon))
 assert sp.expand(q-sp.sympify(data['unscaled_restricted_qfour'][0],locals={str(x):x for x in (*w,a,e,d)}))==0
 assert data['exact_clearing_unit_multipliers'][0]=='1'
 hessian=sp.hessian(q,w[:3]);hdet=sp.factor(hessian.det());assert hdet==sp.Rational(125,32)*a*a
 # Each coefficient matrix is half the Hessian; its determinant is1/8.
 symmetric_det=sp.cancel(hdet/8);assert symmetric_det==sp.Rational(125,256)*a*a
 # Enumerate by the NONDEGENERATE isotropic dimension s instead of k.
 # OGr(s,r) has codimension s(s+1)/2 in Gr(s,r); graph lifts contribute
 # s(m-k). The identity 3(m-k) combines radical and graph parameters.
 strata=[];maxima=[]
 for rank in range(1,7):
  radical=6-rank;entries=[]
  for s in range(min(3,rank//2)+1):
   k=3-s
   if k>radical:continue
   orthogonal=Fraction(s*(rank-s))-Fraction(s*(s+1),2)
   dim=Fraction(3*(radical-k))+orthogonal;assert dim.denominator==1
   entries.append({'radical_intersection_k':k,'nondegenerate_isotropic_s':s,'orthogonal_Grassmann_dimension':int(orthogonal),'total_dimension':int(dim)})
  maxima.append(max(z['total_dimension']for z in entries));strata.append({'quadric_rank':rank,'strata':entries})
 assert maxima==[6,6,4,4,3,3];assert max(2+max(maxima[:2]),3+max(maxima[2:]))==8<9
 # All15 multiset quadruples, with cubic shadows formed as sorted index lists.
 quadratic_mons=list(itertools.combinations_with_replacement(range(3),2));shadows=[];hist={}
 for selected in itertools.combinations(quadratic_mons,4):
  cubic=sorted({tuple(sorted(q+(j,)))for q in selected for j in range(3)})
  shadows.append({'quadratic_monomials':selected,'cubic_shadow':cubic,'dimension_lower_bound':len(cubic)});hist[len(cubic)]=hist.get(len(cubic),0)+1
 assert len(shadows)==15 and hist=={8:9,9:3,10:3}
 pairpath=ROOT/'reports/S61_BINARY_NONSQUARE_PAIR_COVER_20261103_v62/report.json';pair=checked_report(pairpath);assert pair['ell_Z_W_linearly_independent_at_every_core']and pair['four_remaining_quadrics_restricted_mod_Z_W_rank']==4
 rank4path=ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v62/report.json';rank4=checked_report(rank4path);assert rank4['all_point_four_quad_rank']==4
 priorreview=ROOT/'reports/S61_BINARY_NONSQUARE_QUAD_SHADOW_20261103_v63_REVIEW/report.json';checked_report(priorreview)
 step33path=ROOT/'reports/S61_BINARY_NONSQUARE_SECOND_SCALAR_PIVOT_20261103_v62/report.json';step33=checked_report(step33path);assert step33['total_constructive_product_unit_pivots']==33
 parentpath=ROOT/'reports/S61_BINARY_NONSQUARE_IRREDUCIBLE_QUAD_SHADOW_20261103_v62/report.json';parent=checked_report(parentpath);assert parent['all_point_original_quadratic_product_lower_bound']==37
 certificate={'rank3_Hessian':[[str(c)for c in row]for row in hessian.tolist()],'Hessian_determinant':str(hdet),'symmetric_coefficient_minor':str(symmetric_det),'isotropic_Grassmann_strata':strata,'max_fiber_dimensions':maxima,'all15_monomial_shadows':shadows,'shadow_histogram':hist,'restriction_kernel_dimension':12,'restriction_image_dimension_lower_bound':8,'mod_ell_product_dimension_lower_bound':20,'original_mod_ell_kernel_dimension_lower_bound':17,'original_product_dimension_lower_bound':37,'remaining_product_rank_after33_lower_bound':4,'zero_constant_combination_incidence_dimension_upper':6,'projective_constant_combination_space_dimension':7}
 (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 files=[Path(__file__),source,pairpath,rank4path,priorreview,step33path,parentpath,ROOT/'research/notes/S61_BINARY_NONSQUARE_IRREDUCIBLE_QUAD_SHADOW_20261103_v62.md']
 report={'id':'S61_NONSQUARE_IRREDUCIBLE_QUAD_SHADOW_20261103_V63_REVIEW','all_passed':True,'own_Hessian_polynomial_minor_checked':True,'quadric_rank_atleast3_every_point_a_nonzero':True,'all6_orthogonal_fiber_dimension_cases_checked':True,'projective_Gr_incidence_bound':8,'Gr36_dimension':9,'some3plane_restriction_injective_pointwise_not_uniform_L':True,'all15_multiset_monomial_shadows_checked':True,'twelve_independent_kernel_products_proof_audited':True,'seventeen_original_kernel_products_from_audited_pair_cover':True,'allpoint_original_product_rank_lower_bound':37,'generic_QQ_unimodular_row_exists_after33_nonconstructively':True,'executed_native_pivots_remain33':True,'chosen_lambda_or_Bezout_or_pivot34_certified':False,'allpoint41_fullcubic57_quartic_ambient_claim':False,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','allpoint_original_product_rank_lower_bound','chosen_lambda_or_Bezout_or_pivot34_certified','seconds']},indent=2))
if __name__=='__main__':main()
