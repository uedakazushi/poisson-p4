#!/usr/bin/env python3
"""ONE deterministic fixed22 pair-GCD search; resultant output COST ONLY."""
import argparse,hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_COPRIME_PAIR_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(p):return {tuple(m):Q(c)for m,c in p if Q(c)}
def mul(a,b):
 r={}
 for m,c in a.items():
  for n,d in b.items():
   k=tuple(x+y for x,y in zip(m,n));v=r.get(k,Q())+c*d
   if v:r[k]=v
   else:r.pop(k,None)
 return r
def serial(p):return [[list(m),str(c)]for m,c in sorted(p.items(),reverse=True)]
def intprim(p):
 den=math.lcm(*(x.denominator for x in p.values()));z={m:int(c*den)for m,c in p.items()};g=math.gcd(*z.values());return {m:c//g for m,c in z.items()},Q(g,den)
def profile(i,p):
 z,scale=intprim(p);dd=max(m[2]for m in z);slices=[sum(m[2]==k for m in z)for k in range(dd+1)]
 return {'generator_index':i,'terms':len(p),'totaldegree':max(sum(m)for m in p),'degrees_a_e_d':[max(m[j]for m in p)for j in range(3)],'integer_primitive_scale_QQ':str(scale),'integer_coefficient_bits':max(abs(c).bit_length()for c in z.values()),'max_d_coefficient_support':max(slices),'d_coefficient_supports':slices,'d_leading_coefficient':serial({m:c for m,c in p.items()if m[2]==dd})}
def bounds(p,q):
 m=p['degrees_a_e_d'][2];n=q['degrees_a_e_d'][2];N=m+n;ba=n*p['degrees_a_e_d'][0]+m*q['degrees_a_e_d'][0];be=n*p['degrees_a_e_d'][1]+m*q['degrees_a_e_d'][1];btotal=n*p['totaldegree']+m*q['totaldegree']-m*n
 # N! * (2^bF * max_slice_supportF)^n * (2^bG * max_slice_supportG)^m
 upper=math.factorial(N)*((1<<p['integer_coefficient_bits'])*p['max_d_coefficient_support'])**n*((1<<q['integer_coefficient_bits'])*q['max_d_coefficient_support'])**m
 return {'Sylvester_size':N,'d_degrees':[m,n],'resultant_degree_a_upper':ba,'resultant_degree_e_upper':be,'resultant_totaldegree_a_e_upper':btotal,'resultant_support_box_upper':(ba+1)*(be+1),'integer_resultant_abs_coefficient_bits_upper':upper.bit_length(),'integer_resultant_abs_coefficient_upper':str(upper)}
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';ip=base/'certificate.json';br=read(base/'report.json');assert sha(ip)==br['evidence_sha256'][ip.name];c=read(ip);A=c['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling'];assert len(A)==22
 P=[poly(p)for p in A];profiles=[profile(i,p)for i,p in enumerate(P)];pairs=[(i,j)for i in range(22)for j in range(i+1,22)]
 def key(ij):
  i,j=ij;b=bounds(profiles[i],profiles[j]);return(b['Sylvester_size'],b['resultant_support_box_upper'],profiles[i]['terms']+profiles[j]['terms'],i,j)
 pairs.sort(key=key);native=out/'native-input.txt'
 with native.open('w')as f:
  f.write('22\n')
  for p in A:f.write(str(len(p))+''.join(' '+str(c)+' '+' '.join(map(str,m))for m,c in p)+'\n')
  f.write(str(len(pairs))+'\n');f.writelines(str(i)+' '+str(j)+'\n'for i,j in pairs)
 (out/'profiles-and-order.json').write_text(json.dumps({'variable_order':['a','e','d'],'profiles':profiles,'ordered_pairs':pairs,'sort_key':'Sylvester_size, support_box_upper, sum_input_terms, first_index, second_index'},indent=2)+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_coprime_pair_preflight_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-coprime-pair-preflight'
 with(out/'compile.log').open('w')as f:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as f:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','B-pair-native',str(exe),str(native),str(out/'native-certificate.json')],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=read(out/'native-certificate.json');tried=n['tried_ordered_pairs'];assert tried and [x['pair']for x in tried]==[list(x)for x in pairs[:len(tried)]]
 for trial in tried:
  i,j=trial['pair'];g=poly(trial['gcd']);assert mul(g,poly(trial['literal_quotients'][0]))==P[i]and mul(g,poly(trial['literal_quotients'][1]))==P[j]
 selected=tried[-1]['pair']if n['found_coprime_pair']else None;cert={'source22_profile_count':22,'all_tried_literal_quotient_reconstructions_own_Fraction':True,'selected_pair':selected,'resultant_NOT_computed':True,'GB_factorization_or_point_jobs':0}
 if selected:
  i,j=selected;assert poly(tried[-1]['gcd'])=={(0,0,0):Q(1)};R,a,e,d=ring('a,e,d',s.QQ);rp=lambda p:R.from_dict({m:s.QQ(str(c))for m,c in p.items()});g=rp(P[j]).gcd(rp(P[i]));assert g/g.LC==R.one
  for k in selected:z,scale=intprim(P[k]);assert {m:scale*c for m,c in z.items()}==P[k]
  b=bounds(profiles[i],profiles[j]);cert.update({'selected_pair_source':{'indices':selected,'polynomials':[A[i],A[j]],'profiles':[profiles[i],profiles[j]]},'independent_reverse_pair_SymPy_exact_gcd_is1':True,'resultant_cost_bounds_only':b,'field_resultant_nonzero_by_Gauss':True,'all_point_projection_necessary_condition':'For any actual B-rankbad point the selectedpair vanishes; the formal d-resultant vanishes on its(a,e) projection. No finiteprojection asserted; LC-zero cases retained.','d_leading_zero_cases_explicit':[profiles[i]['d_leading_coefficient'],profiles[j]['d_leading_coefficient']]})
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');exe.unlink(missing_ok=True)
 paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),base/'report.json',ip]
 report={'id':ID,'selected_pair':selected,'tried_pairs':len(tried),'native_seconds':n['native_seconds'],'selected_cost':cert.get('resultant_cost_bounds_only'),'all_tried_reconstructions_and_independent_selected_coprimality':True,'resultant_GB_factorization_point_jobs':0,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
