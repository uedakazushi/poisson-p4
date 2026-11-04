#!/usr/bin/env python3
"""ONE R46 exactfactor+wholeproduct, then CLOSEDfield17x24input+chartunitpreflightONLY."""
import hashlib,importlib.util,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_THETA_P80_FACTOR_INPUT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(exist_ok=False,parents=True)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integerhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 base=ROOT/'reports/S61_NONSQUARE_THETA_P80_CUBIC_MINOR_UNITS_20261104_v63';ip=base/'exact-known-unit-strip.json';assert sha(ip)==load(base/'report.json')['evidence_sha256'][ip.name];src=load(ip);R0=src['remaining_polynomial'];R={(int(m),0,0):Q(v)for m,v in R0};assert len(R)==47 and max(m[0]for m in R)==46
 mb=ROOT/'reports/S61_NONSQUARE_THETA_P80_CUBIC_GENERIC_20261104_v63';mp=mb/'substituted17x24-before-rank.json';mr=load(mb/'report.json');assert sha(mp)==mr['evidence_sha256'][mp.name];M=load(mp);assert M['variable']=='zeta'and M['ordinary_polynomial_terms']==5182;matrix=M['integer_matrix17x24'];assert len(matrix)==17 and all(len(row)==24 for row in matrix)
 cb=ROOT/'reports/S61_NONSQUARE_THETA_RANK_DROP_CURVES_20261104_v63';cp=cb/'certificate.json';assert sha(cp)==load(cb/'report.json')['certificate_sha256'];curve=next(x for x in load(cp)['curve_profiles']if x['m']==80);openfactors=[[[int(m[0]),v]for m,v in p]for p in curve['exact_univariate_open_factors']];assert len(openfactors)==6
 def mul(a,b):
  aa,ad=ar.intp(a);bb,bd=ar.intp(b);r={}
  for m,c in aa.items():
   for z,d in bb.items():
    k=tuple(x+y for x,y in zip(m,z));r[k]=r.get(k,0)+c*d
  return{k:Q(v,ad*bd)for k,v in r.items()if v}
 up=lambda p:{(int(m),0,0):Q(v)for m,v in p if Q(v)};unit={(0,0,0):Q(1)}
 for p in openfactors:unit=mul(unit,up(p))
 pre={'R_terms':47,'R_degree':46,'R_primitive_integer_bits':404,'source_integer_matrix_shape':[17,24],'source_terms':5182,'max_source_degree':34,'closed_factor_total_input_coeff_bound':408*46,'rank_or_rightinverse_job_authorized_here':False};(out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 factorinput=out/'factor-input.txt';factorinput.write_text(ar.line(R));factorcpp=ROOT/'code/s61_nonsquare_raw_H_B_common_curve_factor_20261104_v63.cpp';factorexe=ROOT/'work/P80-R46-factor';reducpp=ROOT/'code/s61_nonsquare_Theta_P80_factor_input_20261104_v63.cpp';reduexe=ROOT/'work/P80-closedfactor-input'
 with(out/'compile.log').open('w')as log:
  for cpp,exe in[(factorcpp,factorexe),(reducpp,reduexe)]:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 hashes=[sha(factorexe),sha(reduexe)]
 with(out/'factor.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','P80factor-native',str(factorexe),str(factorinput),str(out/'factor-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 f=load(out/'factor-certificate.json');prod={(0,0,0):Q(f['factor_constant'])};factors=[]
 for rec in f['factors']:
  p=ar.poly(rec['factor']);assert all(m[1]==m[2]==0 for m in p)
  for j in range(rec['exponent']):prod=mul(prod,p)
  factors.append([[m[0],str(v)]for m,v in sorted(p.items(),reverse=True)])
 assert prod==R;sumdegree=sum(max(m for m,v in p)for p in factors);assert sumdegree<=46
 def line(p):return str(len(p))+''.join(' '+str(m)+' '+str(v)for m,v in p)+'\n'
 inp=out/'reduction-input.txt'
 with inp.open('w')as stream:
  stream.write(str(len(factors))+'\n');stream.writelines(map(line,factors));stream.write(line([[m[0],str(v)]for m,v in sorted(unit.items(),reverse=True)]));stream.write('17 24\n')
  for row in matrix:stream.writelines(map(line,row))
 with(out/'reduction.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','P80input-native',str(reduexe),str(inp),str(out/'component-inputs.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'component-inputs.json');profiles=[]
 for i,rec in enumerate(n['components']):
  F=up(rec['factor']);assert F==up(factors[i]);ar.check([(up(rec['chart_unit_Bezout1'][0]),F),(up(rec['chart_unit_Bezout1'][1]),unit)],{(0,0,0):Q(1)})
  for r in range(17):
   for c in range(24):ar.check([(up(rec['matrix_quotients17x24'][r][c]),F),(up(rec['matrix_remainders17x24'][r][c]),{(0,0,0):Q(1)})],up(matrix[r][c]))
  deg=rec['field_degree'];profiles.append({'factor_index':i,'degree':deg,'exponent':f['factors'][i]['exponent'],'factor_terms':len(factors[i]),'input_field_coefficients':rec['reduced_rational_coefficient_count'],'nonzero_matrix_entries':rec['nonzero_matrix_entries'],'input_max_num_plus_den_bits':rec['max_numerator_plus_denominator_bits'],'rank_problem':'K^24->K^17','rightinverse_shape_if_rank17':[24,17],'rightinverse_Kbasis_coeff_upper':24*17*deg,'QQ_instead_large_matrix_shape':[17*deg,24*deg],'all_chart_FULLU_factors_units_certified':True})
 selected=min(profiles,key=lambda p:(p['degree'],p['factor_terms'],p['input_max_num_plus_den_bits'],p['factor_index']))
 cert={'variable':'zeta','native_factor_embedding_variable0':'zeta','original_R46':R0,'factor_constant':f['factor_constant'],'irreducible_factor_records':f['factors'],'whole_R46_product_independent_integer_checked':True,'actual_input_matrix17x24_source':str(mp.relative_to(ROOT)),'all_closed_factor_matrix_divisions_independent_integer_checked':True,'all_6_chart_open_factors':openfactors,'unit_product':[[m[0],str(v)]for m,v in sorted(unit.items(),reverse=True)],'all_closed_factor_chartunit_Bezout1_independent_integer_checked':True,'factor_profiles':profiles,'selected_smallest_factor_index':selected['factor_index'],'scheme_clopen':'R46coprimeeachcertifiedchartlinearunit, soeveryclosedfactorliesinFULLUchart; no factorcomponentdiscarded.','source_rank_preservation':'Commonzeta24*a*B androwLCM/contentclearsareallQQ/FULLUunits; sourceboundactual17x24rankpersistsunderfieldreduction.','future_exact_rightinverse_protocol':'Onlyifauthorized: rowGaussian17x24overK=QQ[zeta]/fi, certifyeachpivotinversebyQQ[zeta]Bezout1 modulo fi, save24x17rightinverse B andverifyALL289entries C*B=I17. No sourcepointtest.','no_rank_rightinverse_GB_point_prime_or_new57expansion_jobs':True};(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');factorexe.unlink(missing_ok=True);reduexe.unlink(missing_ok=True);paths=[Path(__file__),factorcpp,reducpp,helper,ROOT/'research/notes'/(ID+'.md'),base/'report.json',ip,mb/'report.json',mp,mb/'source-denominator-and-operator-binding.json',cb/'report.json',cp]
 report={'id':ID,'R46_factor_degrees_exponents':[[p['degree'],p['exponent']]for p in profiles],'selected_smallest_factor_index':selected['factor_index'],'selected_profile':selected,'total_closed_matrix_field_coefficients':sum(p['input_field_coefficients']for p in profiles),'all_independent_product_source_reductions_chartunits_checked':True,'factor_native_seconds':f['native_seconds'],'reduction_native_seconds':n['native_seconds'],'rank_rightinverse_GB_point_prime_jobs':0,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHAs_not_committed':hashes,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
