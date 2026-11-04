#!/usr/bin/env python3
"""Source-pinned combined finiteBbad corollary; no new solver/determinant/points."""
import hashlib,json,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_FINITE_BAD_LOCUS_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def pin(stem,file='certificate.json'):
 base=ROOT/'reports'/stem;p=base/file;assert sha(p)==load(base/'report.json')['evidence_sha256'][file];return base,p,load(p)
def degree(p):return max([int(m)for m,c in p if Q(c)],default=-1)
def main():
 t=time.monotonic();rb,rp,rc=pin('S61_NONSQUARE_RAW_H_B_PROJECTION_CURVE_RESIDUAL_INPUT_20261104_v63');assert rc['independent_exact_gcd_h1_k1_is1'];oldbound=rc['isolated_B_bad_points_bound_after_H4'];assert oldbound==1600
 lb,lp,lc=pin('S61_NONSQUARE_RAW_H_B_COMPONENT_COMMON_ROOT_L_20261104_v63');assert lc['generic_common_roots_all_U_boundary'];QL,QQ=None,None
 qb,qp,qc=pin('S61_NONSQUARE_RAW_H_B_COMPONENT_EXACT_QUOTIENT_Q_20261104_v63');assert len(qc['final']['common_gcd'])==2 and qc['final']['common_gcd'][0]['d_power']==1
 ab,ap,ac=pin('S61_NONSQUARE_RAW_H_B_Q_ROOT_ALL_MINORS_20261104_v63','native-certificate.json');nz=ac['nonzero_remainder_indices'];assert nz and nz[0]==0;norm=ac['all22_evaluations'][0]['remainder_quadratic_norm_e'];assert norm and all(m[0]==m[2]==0 for m,c in norm);norm_e=[[m[1],c]for m,c in norm];assert degree(norm_e)==38
 lpolys=lc['finite_e_exception_polynomials'];qpolys=qc['finite_e_exception_polynomials'];D=[[3,'486'],[2,'243'],[1,'-297'],[0,'-54']];ldeg=sum(map(degree,lpolys));qdeg=sum(map(degree,qpolys))+degree(D)+degree(norm_e);assert ldeg==65 and qdeg==361
 lbound=4*ldeg;qbound=8*qdeg;total=oldbound+lbound+qbound;assert total==4748
 fb,fp,fc=pin('S61_NONSQUARE_RAW_H_B_COMMON_CURVE_FACTOR_20261104_v63');hb,hp,hc=pin('S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63');assert max(m[2]for m,c in hc['monic_H_laurent'])==4
 cert={'theorem':'The actual7x3 B rank<3 locus onFULLU=D(adFBCJ) is finite or empty overQQbar/C.','does_NOT_solve':['emptiness/globalB3','rawH7 onp=0','fullcubic57','fullquartic239'],'all_B_bad_contained_in_three_covers':[
  {'name':'oldresidual','equations':['H=0','h1=0','k1=0'],'inputs':str(rp.relative_to(ROOT)),'reason_finite':'h1/k1coprime,totaldegrees16/25; Hmonic4wherea!=0','geometric_point_bound':oldbound},
  {'name':'L_component_exceptions','equations':['H=0','L=0','product(E_L)=0'],'e_polynomial_factors':lpolys,'total_e_degree_upper':ldeg,'reason_finite':'OutsideE_L,3sourcecommonGCDrootsareallU-boundary; everybadpointinsideFULLUmustbeexception. Llinearinagivesatmost1a,pere,andHmonic4.','geometric_point_bound':lbound},
  {'name':'Q_component_exceptions','equations':['H=0','Q=0','product(E_Q)*D*Norm_R0=0'],'e_polynomial_factors':qpolys+[D,norm_e],'total_e_degree_upper':qdeg,'reason_finite':'OutsideE_QandD,3sourcecommonrootisunique d=-N/D; actualminor0forcesR0=0,soitsnonzeroquadraticnormNorm_R0(e)=0. Qmonica2andHmonic4giveatmost8pointspere.','geometric_point_bound':qbound}],
  'total_geometric_point_bound':total,'bound_scope':'numberofgeometricpoints,NOTschemelength/notflatcoverdegree; overlapsandmultiplicitiesareovercounted','H_monic_d4_on_aunit':hc['monic_H_laurent'],'factor_G_Q_L4_source':str(fp.relative_to(ROOT)),'all_denominator_zero_cases_covered_not_dropped':True,'additional_FULLU_localizer_inverse_finiteness_not_needed_for_finite_point_statement':True,'next_small_input_gate':'TakeunivariatesquarefreepartofproductE_L andofproductE_Q*D*NormR0 (degrees65/361) beforefinitealgebra construction; no suchgcd/resultant computedhere.','no_new_GB_resultant_point_prime_or_field_job':True}
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md')]
 for base,p in[(rb,rp),(lb,lp),(qb,qp),(ab,ap),(fb,fp),(hb,hp)]:paths.extend([base/'report.json',p])
 report={'id':ID,'actual_B_rankbad_finite_or_empty':True,'three_geometric_point_upper_bounds':[oldbound,lbound,qbound],'geometric_point_upper_bound_total':total,'rawH7_p0_and_full57_239_NOT_inferred':True,'exception_univariate_product_degrees_upper':[ldeg,qdeg],'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{'certificate.json':sha(out/'certificate.json')},'seconds':time.monotonic()-t,'frozen':True,'independent_review_pending':True};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
