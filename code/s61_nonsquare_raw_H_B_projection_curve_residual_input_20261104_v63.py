#!/usr/bin/env python3
"""Savednorm GCD decomposition and residual e-resultant COST ONLY."""
import hashlib,json,math,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_PROJECTION_CURVE_RESIDUAL_INPUT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def poly(p):return{tuple(m):Q(c)for m,c in p if Q(c)}
def ser(p):return[[list(m),str(c)]for m,c in sorted(p.items(),reverse=True)]
def scale(p,c):return{m:c*v for m,v in p.items()if c*v}
def profile(p):
 den=math.lcm(*(x.denominator for x in p.values()));g=math.gcd(*(int(x*den)for x in p.values()));z={m:int(x*den)//g for m,x in p.items()};n=max(m[1]for m in p);lc={m:c for m,c in p.items()if m[1]==n};assert len(lc)==1 and (0,n,0)in lc
 return {'terms':len(p),'degrees_a_e':[max(m[i]for m in p)for i in range(2)],'totaldegree':max(sum(m)for m in p),'primitive_integer_coefficient_bits':max(abs(v).bit_length()for v in z.values()),'max_e_coefficient_support':max(sum(m[1]==j for m in p)for j in range(n+1)),'QQ_scale_to_primitive':str(Q(den,g)),'constant_e_leading_coefficient':str(lc[(0,n,0)])}
def main():
 t=time.monotonic();base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_ADDITIONAL_MINOR_NORM_20261104_v63';np=base/'native-certificate.json';assert sha(np)==load(base/'report.json')['evidence_sha256'][np.name];n=load(np);G=poly(n['gcd_newK_oldh']);h1=poly(n['oldh_quotient_by_gcd']);k1=poly(n['newK_quotient_by_gcd']);p0,p1=map(profile,[h1,k1]);pg=profile(G)
 RR,a,e,d=ring('a,e,d',s.QQ);rp=lambda p:RR.from_dict({m:s.QQ(str(v))for m,v in p.items()});g=rp(h1).gcd(rp(k1));assert g/g.LC==RR.one
 m=p0['degrees_a_e'][1];nn=p1['degrees_a_e'][1];da=nn*p0['degrees_a_e'][0]+m*p1['degrees_a_e'][0];N=m+nn;height=math.factorial(N)*((1<<p0['primitive_integer_coefficient_bits'])*p0['max_e_coefficient_support'])**nn*((1<<p1['primitive_integer_coefficient_bits'])*p1['max_e_coefficient_support'])**m;bounds={'e_degrees':[m,nn],'Sylvester_size_if_used':N,'univariate_a_resultant_degree_upper':da,'univariate_support_upper':da+1,'primitive_integer_resultant_coefficient_bits_upper':height.bit_length(),'conservative_univariate_output_bytes_upper':(da+1)*(100+(height.bit_length()*30103+99999)//100000),'exact_resultant_NOT_computed':True}
 mb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63';mp=mb/'certificate.json';assert sha(mp)==load(mb/'report.json')['evidence_sha256'][mp.name];H=load(mp)['monic_H_laurent'];Emonic=scale(G,1/Q(pg['constant_e_leading_coefficient']));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 cert={'common_G':ser(G),'common_G_monic_e12':ser(Emonic),'oldh_quotient_h1':ser(h1),'newK_quotient_k1':ser(k1),'source_profiles_G_h1_k1':[pg,p0,p1],'independent_exact_gcd_h1_k1_is1':True,'set_decomposition':'V(oldh,K)=V(G) union V(h1,k1); retain bothpieces','curve_cover_monic_H_d4':H,'curve_cover_free_basis_e_i_d_j':[[i,j]for j in range(4)for i in range(12)],'curve_cover_free_rank48':True,'isolated_projection_points_bound_Bezout':p0['totaldegree']*p1['totaldegree'],'isolated_B_bad_points_bound_after_H4':4*p0['totaldegree']*p1['totaldegree'],'residual_resultant_cost_only':bounds,'no_resultant_points_GB_newprime_or_badexistence':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',np,mb/'report.json',mp];report={'id':ID,'common_curve_G_terms':len(G),'common_curve_e_degree':12,'finite_free_curve_cover_rank':48,'residual_h1_k1_terms':[len(h1),len(k1)],'residual_degrees_e':[m,nn],'residual_intersections_finite_or_empty':True,'isolated_B_bad_bound':4*p0['totaldegree']*p1['totaldegree'],'resultant_cost_bounds':bounds,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{'certificate.json':sha(out/'certificate.json')},'seconds':time.monotonic()-t,'frozen':True,'independent_review_pending':True};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
