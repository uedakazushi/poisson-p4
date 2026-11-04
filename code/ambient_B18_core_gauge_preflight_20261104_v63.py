#!/usr/bin/env python3
"""Exact tiny-core cut, linear-pivot and algebraic-gauge recipe preflight."""
import argparse,collections,gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_CORE_GAUGE_PREFLIGHT_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(M):return [[str(x)for x in row]for row in M.tolist()]
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 rp=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json';r=read(rp);assert r['all_passed']
 wp=ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63/report.json';wr=read(wp);wc=read(wp.parent/'certificate.json');assert wr['all_frozen_source_weight_and_sparsity_checks_passed']
 for rel,h in wr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 assert sha(wp.parent/'certificate.json')==wr['certificate_sha256']
 gp=ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json';gc=read(gp);weights=gc['physical150_torus_weights']
 op=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz';oc=read(op);O=s.Matrix(oc['physical150_by19_scalar_then_GL']);P0=s.Matrix(oc['physical_point150'])
 cuts=wc['common_physical_group_cuts19'];base=wc['remaining_actual_even_base10_raw_coordinates'];R9=[j for j in base if weights[j]==0];assert R9==[27,108,109,112,113,117,118,119,61]
 zcuts=[j for j in cuts if weights[j]==0];gzero=[j for j in range(19)if all(not O[i,j]or weights[i]==0 for i in range(150))];assert len(zcuts)==len(gzero)==15
 Q=O[zcuts,gzero];Qinv=Q.inv();assert Q*Qinv==s.eye(15);qdet=s.factor(Q.det());assert qdet
 params=s.symbols('p k0:9');p=params[0];Ks=list(params[1:]);K0=[k for k in Ks if k!=Ks[2]];origin=dict.fromkeys(params,0)
 records={}
 for name in['CA','CB']:
  fam=r['families'][name];U=s.Matrix([s.expand(s.sympify(v))for v in fam['core_full150_coordinates']]);assert U.subs(origin)==P0
  cutpolys=[s.expand(U[j]-P0[j])for j in cuts];active=[j for j,f in zip(cuts,cutpolys)if f]
  assert active==[15,16,22,23,50,65,132]
  U0=U.subs(Ks[2],0);rawbase=U0[R9,:];rawK=rawbase.jacobian(K0);assert all(not x.free_symbols.intersection(set(K0))for x in rawK)
  rawrank0=rawK.subs(p,0).rank();assert rawrank0==7
  minors=[s.factor(rawK[[j for j in range(9)if j!=drop],:].det())for drop in range(9)]
  nonzero=[f for f in minors if f];gcd=s.Poly(nonzero[0],p)
  for f in nonzero[1:]:gcd=s.gcd(gcd,s.Poly(f,p))
  gcd=gcd.monic().as_expr();genericrank=8 if nonzero else rawrank0
  L=U0.jacobian(K0).subs(origin);assert L.shape==(150,8)
  lin_gauge=-Qinv*L[zcuts,:]
  corrected=L[R9,:]+O[R9,gzero]*lin_gauge;assert corrected.rank()==8
  piv=list(corrected.T.rref()[1]);assert len(piv)==8
  fold=next(j for j in range(9)if j not in piv);B=corrected[piv,:];Binv=B.inv();assert B*Binv==s.eye(8)
  W=s.Matrix([s.Poly(v,*params).coeff_monomial(p*p)for v in U0]);assert W==s.Matrix(fam['raw_second_coefficient'])
  Wg=W[R9,:]-O[R9,gzero]*Qinv*W[zcuts,:]
  foldc=s.factor(Wg[fold]-(corrected[fold,:]*Binv*Wg[piv,:])[0]);assert foldc
  # Independent23x23 block verifies the sequential constant Schur recipe.
  block=O[zcuts,gzero].row_join(L[zcuts,:]).col_join(O[[R9[j]for j in piv],gzero].row_join(L[[R9[j]for j in piv],:]))
  blockdet=s.factor(qdet*B.det());assert blockdet and block.det()==blockdet
  records[name]={'exact19_group_cut_deviations':[[j,str(f)]for j,f in zip(cuts,cutpolys)],'nonzero_group_cuts_on_core':active,
   'exact_ungauged_zero_R9_functions':[[j,str(U0[j])]for j in R9],'ungauged_zero_R9_K8_coefficient_matrix':enc(rawK),
   'ungauged_zero_R9_K8_rank_at_p0':rawrank0,'ungauged_rank_over_QQ_p':genericrank,'all9_ungauged8minors':list(map(str,minors)),'monic_gcd_of_nonzero_ungauged8minors':str(gcd),
   'p_is_fold_parameter_NOT_a_unit':True,'zero_group15_linear_gauge_K8_coefficients':enc(lin_gauge),
   'corrected_gauged_R9_K8_matrix':enc(corrected),'selected8_base_raw_rows':[R9[j]for j in piv],'selected8_K_minor':str(s.factor(B.det())),'remaining_fold_raw_row':R9[fold],
   'full23_implicit_Jacobian_determinant':str(blockdet),'fold_quadratic_coefficient_after23_unit_eliminations':str(foldc),
   'fold_linear_coefficient_at_origin':0,'all_actual_FiniteRank2_coefficients_unexpanded':True}
 # Exact circuit for all19 or only15 zero-weight gauge equations.
 generators=[('tensor_scalar',None)]+[('GL_matrix_unit',[j//5,j%5])for j in oc['GL18_generator_indices']]
 assert len(generators)==19
 circuit={'source':'frozen full150 core polynomial vector U_i(p,K9)',
 'G':'I5+sum_{j=1..18} gamma_j E_(GL18_generator_indices[j-1])',
 'tensor_multiplier':'1+gamma0','inverse_G':'adj(G)/det(G)',
 'cleared_pushforward':'(1+gamma0)*G*U_i(adj(G)*X)*G^T',
 'denominator':'det(G)^2','gauge_equation_for_each_selected_raw_row_c':'coefficient_c(cleared_pushforward) - P0_c*det(G)^2 = 0',
 'all19_cut_rows':cuts,'all19_generator_order':generators,'zero_sector_positive_source_parameter':'k2=0',
 'zero_sector_positive_group_parameters_identically0':[j for j in range(19)if j not in gzero],
 'zero_sector15_group_parameter_indices':gzero,'zero_sector15_cut_rows':zcuts,
 'unit_factors_at_origin':['1+gamma0','det(G)','det(partial15zero_gauge_equations/partial15gamma)'],
 'zero_sector_group_Jacobian_determinant_at_origin':str(qdet),'individual_equation_group_degree_upper_bound':11,
 'scheme_or_formal_branch':'take the unique pointed algebraic implicit branch gamma=0 at sourceorigin; full polynomial system not expanded',
 'not_claimed':'rational or triangular closed gamma formulas; a finite/proper ordinary-Zariski base model; allzero17 coefficients'}
 evidence={'source_order':list(map(str,params)),'zero_source_order':['p',*map(str,K0)],'common_group_cuts19':cuts,'common_zero_base9_raw_rows':R9,'records':records,'exact_finite_polynomial_circuit_recipe':circuit,
 'algebraic_rank2_recipe':'For each family solve the15zero gauge equations plus8 selected base equations by the23x23 unit Jacobian in(gamma15,K8), over targetx9 andp. Substitute into remaining base equationF(x,p)=0. Nonzero recorded quadraticp coefficient implies Weierstrass degree2. Keep the unique monic quadratic pointed Hensel factor; its coefficients and allnormalizer trace coefficients are algebraic power series, not assumed rational.',
 'completed_normalizer':'QQ[[x9]][tau]/(tau^2-A(x)*tau-B(x)), finitefree2; A,B algebraic power series, explicit coefficient expansion UNCOMPUTED',
 'algebraicity_justification':'The unit-Jacobian solution is an algebraic implicit branch of a finite polynomial system; the last fold is a locallyquasifinite algebraic germ. Coordinates are algebraic overQQ(x9); the two localformal conjugates give algebraic trace/norm. Completion rank2 is certified, not a proper uncompleted model.',
 'no_GB_field_newpoint_fullgauge_expansion_or_primary_solver':True}
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),rp,wp,wp.parent/'certificate.json',gp,op]
 report={'id':ID,'all_tiny_core_cut_linearity_and_fold_checks_passed':True,'19cuts_identically_zero_on_source_cores':False,'ungauged_zero_base_differential_rank':7,'gauged_zero_base_differential_rank':8,'common_zero_group_unit_dimension':15,'unit_implicit_block_dimension':23,'pointed_Weierstrass_degree':2,'finite_polynomial_circuit_recipe_saved':True,'rank2_normalizer_algebraic_power_series_recipe_proved':True,'explicit_rational_gauge_or_expanded_trace_coefficients_computed':False,'finite_proper_uncompleted_base_model_claimed':False,'no_GB_elimination_newpoint_or_generic_solver':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())
