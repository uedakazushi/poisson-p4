#!/usr/bin/env python3
"""Independent fixed modular product/normality and source/Gauss proof audit."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1];AUTHOR='S61_NONSQUARE_PREPOWER_R434_SAUX_DIVISION_20261104_v63';STEM=AUTHOR+'_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));hashes={}
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 work=ROOT/'work/prePOWER-R434-Saux-fixed-divrem-independent-review';work.mkdir(parents=True,exist_ok=True)
 def load(p):hashes[p]=sha(p);return read(p)
 def get(st,n='report.json'):return load('reports/'+st+'/'+n)
 r=get(AUTHOR);assert r['all_passed'] and hashes['reports/'+AUTHOR+'/report.json']=='05c3fd88a40071a6b97506ea52659fd7f04d713f7fc2fa96c6ae7c56d1476489'
 for p,h in r['input_sha256'].items():assert sha(p)==h;hashes[p]=h
 for n,h in r['evidence_sha256'].items():p='reports/'+AUTHOR+'/'+n;assert sha(p)==h;hashes[p]=h
 P=get('S61_NONSQUARE_PREPOWER_U4_T15_MODP_POLYNOMIAL_20261104_v63_REVIEW');assert P['all_passed']
 assert hashes['reports/S61_NONSQUARE_PREPOWER_U4_T15_MODP_POLYNOMIAL_20261104_v63_REVIEW/report.json']=='fe28fa127d5134991fd22fc7a60347ad69e2990dc925e1fabfe90a9f68f7702b'
 I=get('S61_NONSQUARE_PREPOWER_R434_MODP_IRREDUCIBILITY_20261104_v63_REVIEW');assert I['QQ_R434_irreducible_independently_proved']
 assert hashes['reports/S61_NONSQUARE_PREPOWER_R434_MODP_IRREDUCIBILITY_20261104_v63_REVIEW/report.json']=='9daa160d5e2f4f89a6b11f5ee1b1be318a016cb390187c14844306d2f029bb38'
 ch=get('S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63_REVIEW');assert ch['all_passed']
 strip=get('S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63','native.json');entry=get('S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_ENTRY_20261104_v63','native.json')
 rank=get('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63','certificate.json');rawpeer=get('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63_REVIEW');assert rawpeer['all_passed']
 native=get(AUTHOR,'native.json');R=strip['remaining_polynomial'];S=entry['exact_Saux_Laurent_polynomial']
 assert len(R)==434 and len(S)==7022 and native['R_leading_monomial']==[2,15,3]
 assert native['R_leading_coefficient']==27306 and len(native['quotient'])==3212 and len(native['remainder'])==6729
 def modpoly(pol):
  d={}
  for m,c in pol:
   q=Q(c);assert q.denominator%32003
   a=q.numerator*pow(q.denominator,-1,32003)%32003
   if a:d[tuple(m)]=a
  return d
 Rbar=modpoly(R);Sbar=modpoly(S)
 assert max(map(sum,Rbar))==max(sum(m) for m,c in R)==20 and max(map(sum,Sbar))==max(sum(m) for m,c in S)==61
 leading=max(Rbar,key=lambda m:(sum(m),-m[2],-m[1],-m[0]));assert leading==(2,15,3) and Rbar[leading]==27306
 # Actual source polynomiality/p-integrality, independent of trusting its
 # full public QQ expansion: propagate conservative a-valuations on the
 # already source-bound768 operation input. No coefficient multiplication.
 oldinput='work/prePOWER-direct-one-entry-independent-review/input.txt'
 assert sha(oldinput)=='b8b29b2e9239548811f64e45d828b49d0e3a6e15969d369422e56e1ecd5a3e74';hashes[oldinput]=sha(oldinput)
 token=iter((ROOT/oldinput).read_text().split());reserve=int(next(token));count=int(next(token));assert count==768
 low=[];denchecks=0;INF=10**8
 for i in range(768):
  op=next(token)
  if op=='p':
   n=int(next(token));ex=[]
   for _ in range(n):
    a,e,d,c=int(next(token)),int(next(token)),int(next(token)),Q(next(token));assert c.denominator%32003;denchecks+=1
    if c:ex.append(a)
   low.append(min(ex) if ex else INF)
  else:
   a,b=int(next(token)),int(next(token));assert a<i and b<i and op in ('a','s','m')
   low.append(INF if op=='m' and (low[a]==INF or low[b]==INF) else low[a]+low[b] if op=='m' else min(low[a],low[b]))
 assert low[767]>=0
 degree_a=lambda pol:max(m[0] for m,c in pol if Q(c))
 deg={'R434':degree_a(R),'p':degree_a(rank['p']),'Delta407':degree_a(rank['nonzero_selected_minor_Delta']),'h22':degree_a(strip['actual_h22']),'B':1}
 assert deg=={'R434':11,'p':2,'Delta407':8,'h22':4,'B':1}
 assert rank['pairdet_formula']=='-(6/5)*B*p^2'
 # All named FULLU factors have a-degree <=1 (or zero for d), so the same
 # irreducible R cannot disappear under the FULLU localization.
 text=[]
 def emit(pol):
  text.append(str(len(pol)))
  for m,c in pol:text.append(' '.join(map(str,m))+' '+str(c))
 emit(R);emit(S);emit(native['quotient']);emit(native['remainder']);text.append('2 15 3 27306')
 ip=work/'input.txt';ip.write_text('\n'.join(text)+'\n');assert ip.stat().st_size<32<<20
 pre={'no_division_or_factorization_call':True,'source768_conservative_a_valuation':low[767],'source768_QQ_denominator_unit_checks':denchecks,'source_actual_Saux_in_Zp_polynomial_ring':True,'source_R_and_Saux_degree_preserved':[20,61],'normality_terms':6729,'quotient_terms':3212,'a_degrees':deg,'source_fp_witness':hashes['reports/S61_NONSQUARE_PREPOWER_U4_T15_MODP_POLYNOMIAL_20261104_v63_REVIEW/report.json']}
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 cpp=Path(__file__).with_suffix('.cpp');exe=work/'checker';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,53-(time.monotonic()-started)))
 with (work/'native.log').open('w') as log:subprocess.run([str(exe),str(ip),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,57-(time.monotonic()-started)))
 n=read('reports/'+STEM+'/native.json');assert n['all_passed'] and n['division_calls']==n['factor_calls']==0
 proof={'normal_remainder':'All6729 terms avoid divisibility by LM(R)=(2,15,3), under DEGREVLEX. In a field polynomial ring the leading term of any nonzero multiple of R is divisible by LM(R), so this nonzero normal remainder proves Rbar does not divide the actual symbolic source Sauxbar.',
  'Gauss_transfer':'Rprimitive is p-primitive over Z_(32003), and actual source Saux is a p-integral polynomial by the768 atom denominator/valuation proof. If R divides actual Saux overQQ, Gauss gives a p-integral quotient, whose reduction contradicts the nonzero normal remainder.',
  'QQ_gcd':'Independent QQ irreducibility of R plus nondivisibility gives gcd_QQ(R,actualSaux)=1.',
  'complex_base_extension':'R and Saux form a regular sequence overQQ polynomial ring (R irreducible, Saux nonzero modR). Flat field extension toC preserves the second multiplication injection, so they have no common complex divisor even though R itself is not proved absolutely irreducible.',
  'additional_coprimality':'a-degreeR=11; nonzero p,Delta407,h22,B have a-degrees2,8,4,1, each strictly smaller. Irreducibility forbids R dividing any. Delta5=-(6/5)B*p²*Delta407, so gcd_QQ(R,Delta5*h22)=1. No absolute irreducibility premise.',
  'actual_direct_P5':'OnD(Delta5) first eliminate H5 only. Its u4 row is pair-only, so the corrected pure row equals rho_u/(delta2*rhoLCM)*Bhat. The actual5border determinant is rho_u/(delta2*rhoLCM)*Saux, a unit multiple of Saux. NoD4/h22 inverse is used.',
  'generic_old_R_hypersurface_coverage':'Thus actualP5+H5 is invertible onFULLU D(Delta5*Saux), yielding the direct10pivot chart. No commonCdivisor of R andDelta5*Saux implies this chart is dense on each geometricR component; remaining intersections have codimension at least2 if nonempty.',
  'not_claimed':['Unitideal or disjointness of R and Saux','Absolute/C irreducibility of R','Residual rank5 or global full57','All chart complements removed','Independent exactQQ7022 expansion']}
 (out/'independent-normality-Gauss-and-chart-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 for p in [str(Path(__file__).relative_to(ROOT)),str(cpp.relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'fixed_qR_plus_nonzero_normal_remainder_verified':True,'actual_source_Saux_QQ_R_nondivisibility_independently_proved':True,'QQ_gcd_R_actualSaux_equals1':True,'no_common_complex_divisor_by_flat_regular_sequence_extension':True,'QQ_gcd_R_Delta5_h22_equals1':True,'actual_old_R_hypersurface_generic_directP5_coverage':True,'unitideal_disjoint_absolute_irreducibility_or_full57_claimed':False,'no_new_division_factor_point_prime_or_QQ_retry':True,'input_sha256':hashes,'ignored_reproducible_input':{'path':str(ip.relative_to(ROOT)),'sha256':sha(str(ip.relative_to(ROOT)))},'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')}),flush=True)
if __name__=='__main__':main()
