#!/usr/bin/env python3
"""ONE newminor monic4 norm input only; no norm or nested resultant calculation."""
import hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_ADDITIONAL_MINOR_NORM_INPUT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_input_20261104_v63.py';spec=importlib.util.spec_from_file_location('fixednormhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';bp=base/'certificate.json';br=load(base/'report.json');assert sha(bp)==br['evidence_sha256'][bp.name];c=load(bp);A=c['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling'];profiles=[]
 for i,p in enumerate(A):
  if i in[10,15]:continue
  profiles.append({'index':i,'terms':len(p),'degrees_a_e_d':[max(m[j]for m,v in p)for j in range(3)],'totaldegree':max(sum(m)for m,v in p)})
 profiles.sort(key=lambda p:(p['degrees_a_e_d'][2],p['terms'],p['totaldegree'],p['index']));ki=profiles[0]['index'];assert ki==20;k=ar.poly(A[ki])
 mb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63';mp=mb/'certificate.json';assert sha(mp)==load(mb/'report.json')['evidence_sha256'][mp.name];H=ar.poly(load(mp)['monic_H_laurent'])
 M=[[{}for j in range(4)]for i in range(4)];Qs=[];tr=[]
 for j in range(4):
  q,r,t=ar.divide_monic(ar.shift(k,(0,0,j)),H);Qs.append(q);tr.append(t)
  for i in range(4):M[i][j]={(m[0],m[1],0):v for m,v in r.items()if m[2]==i}
  assert ar.add(ar.mul(q,H),*(ar.shift(M[i][j],(0,0,i))for i in range(4)))==ar.shift(k,(0,0,j))
 allp=[p for row in M for p in row];terms=sum(map(len,allp));assert terms<150000
 mina=min(m[0]for p in allp for m in p);maxa=max(m[0]for p in allp for m in p);maxe=max(m[1]for p in allp for m in p);clear=max(0,-mina);den=math.lcm(*(x.denominator for p in allp for x in p.values()));mx=max(abs(int(x*den))for p in allp for x in p.values());st=max(map(len,allp));height=24*(st*mx)**4
 bounds={'matrix_uniform_a_clear_power':clear,'matrix_uniform_QQ_denominator':str(den),'cleared_matrix_degree_a':maxa+clear,'cleared_matrix_degree_e':maxe,'cleared_matrix_integer_coefficient_bits':mx.bit_length(),'cleared_det_degree_a_upper':4*(maxa+clear),'det_degree_e_upper':4*maxe,'cleared_det_support_box_upper':(4*(maxa+clear)+1)*(4*maxe+1),'cleared_det_integer_coefficient_bits_upper':height.bit_length(),'det_permutation_products':24,'det_laurent_a_exponent_lower':4*mina,'det_laurent_a_exponent_upper':4*maxa}
 bounds['conservative_cleared_det_serialized_bytes_upper']=bounds['cleared_det_support_box_upper']*(100+(height.bit_length()*30103+99999)//100000+len(str(den**4)));assert bounds['conservative_cleared_det_serialized_bytes_upper']<32*1024*1024
 # A different source-bound GCD only certifies nonzero proposed norm, not coprime projection equations.
 RR,a,e,d=ring('a,e,d',s.QQ);rp=lambda p:RR.from_dict({m:s.QQ(str(v))for m,v in p.items()});ha=-min(m[0]for m in H);HH=ar.shift(H,(ha,0,0));gcd=rp(k).gcd(rp(HH));gcd=gcd/gcd.LC;dd=max(m[2]for m in gcd);assert dd==0
 fc=ROOT/'reports/S61_NONSQUARE_RAW_H_B_FINITE_CURVE_COVER_20261104_v63';fp=fc/'certificate.json';assert sha(fp)==load(fc/'report.json')['evidence_sha256'][fp.name];E=ar.poly(load(fp)['monic_E_e28']);lc={m:v for m,v in k.items()if m[2]==4}
 cert={'selection_rule':'Excludeindices10,15; order by ddegree,inputterms,totaldegree,index; selectfirst only','ordered_remaining_profiles':profiles,'selected_generator_index':ki,'source_minor35_indices':[i for i,x in enumerate(c['minor35_to_distinct_generator_index'])if x==ki],'source_k':A[ki],'source_H':ar.ser(H),'actual_monic_d_remainder_k':ar.ser(ar.add(*(ar.shift(M[i][0],(0,0,i))for i in range(4)))),'matrix_rows':[[ar.ser(p)for p in row]for row in M],'column_quotients_by_H':[ar.ser(p)for p in Qs],'all4_full_monic_division_identities_verified':True,'division_steps_by_column':tr,'matrix_entry_stats':[[ar.stats(p)for p in row]for row in M],'cost_bounds_only':bounds,'source_k_Hclear_gcd':ar.ser({m:Q(str(v))for m,v in gcd.items()}),'gcd_ddegree_is0_so_proposed_norm_nonzero':True,'existing_E':ar.ser(E),'normK_h_coprimality_NOT_tested_or_inferred':True,'future_coprimality_goal':'Compute only this secondnormK, stripknownU-supportedfactors with exactreconstruction, then gcd(K,h) inQQ[a,e]. gcd1 would make their commonprojection0D, hence B-badfinite viafiniteH. Sharedcurvewouldbe onlynecessary, notbadexistence.','unknown_localizer_inverses':[],'norm_nestedresultant_GB_point_prime_jobs':0}
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');paths=[Path(__file__),helper,ROOT/'research/notes'/(ID+'.md'),base/'report.json',bp,mb/'report.json',mp,fc/'report.json',fp];report={'id':ID,'selected_k_index':ki,'source_k_terms':len(k),'source_k_ddegree':4,'monic_remainder_k_terms':sum(len(M[i][0])for i in range(4)),'matrix4x4_total_terms':terms,'entry_term_counts':[[len(p)for p in row]for row in M],'cost_bounds':bounds,'proposed_norm_nonzero_from_exact_source_gcd':True,'norm_h_coprimality_NOT_computed':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{'certificate.json':sha(out/'certificate.json')},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
