#!/usr/bin/env python3
"""Independent primitive normalization, source signs, and typed13th pivot audit."""
import hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from functools import reduce
from pathlib import Path
import s61_nonsquare_prePOWER_kappa11_schur_review_20261104_v63 as U
import s61_nonsquare_prePOWER_full57_source_lift_review_20261104_v63 as P
ROOT=Path(__file__).resolve().parents[1];BASE='S61_NONSQUARE_PREPOWER_K120_PIVOT13_SCHUR_20261104_v63';ID=BASE+'_REVIEW';OUT=ROOT/'reports'/ID
ENTRY='S61_NONSQUARE_PREPOWER_K120_T14_NESTED_ENTRY_20261104_v63';OLD='S61_NONSQUARE_PREPOWER_G115_PIVOT12_SCHUR_20261104_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(n,v):(OUT/n).write_text(json.dumps(v,indent=2)+'\n')
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={};cache={}
 def load(base,name='report.json'):
  p=ROOT/'reports'/base/name;rel=str(p.relative_to(ROOT));pins[rel]=sha(p)
  if rel not in cache:cache[rel]=json.loads(p.read_text())
  return cache[rel]
 r=load(BASE);assert pins[f'reports/{BASE}/report.json']=='63a29774840ba840347f0e9003096c37e540f328331dbfdba6c5f610713f5fb0'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/BASE/n)==h
 assert load(ENTRY+'_REVIEW')['actual_N7_NONZERO_char0']and pins[f'reports/{ENTRY}_REVIEW/report.json']=='0571cdcd3f71e2a2e853373b85a72eb25a2158f6d5bcbf3b7ba2b774cd70a763'
 assert load(OLD+'_REVIEW')['all_passed'] and load('S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63_REVIEW')['all_passed']
 bar=load(BASE,'primitive-N7bar.json');src=load(ENTRY,'native.json')['normalized_N7'];oldpol=src['integer_coefficients'];newpol=bar['polynomial'];D=int(src['global_denominator']);content=reduce(math.gcd,(abs(int(c))for m,c in oldpol));amin=min(m[0]for m,c in oldpol)
 assert len(oldpol)==len(newpol)==14986 and content>0 and amin==9 and F(content,D)==F(16,949680497646331787109375)
 for (a,x),(b,y) in zip(oldpol,newpol):assert int(x)==content*int(y)and a==[b[0]+9,b[1],b[2]]
 assert reduce(math.gcd,(abs(int(c))for m,c in newpol))==1 and min(m[0]for m,c in newpol)==0
 assert bar['content_positive']==str(content)and bar['original_global_denominator']==str(D)and F(bar['rational_scalar'])==F(content,D)
 # Source-specific K120 small row genuinely annihilates ALL middle3 columns.
 hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json');b=hc['raw_compact14x5'];g=P.dec(b[6][1]);hp=P.add(P.dec(b[11][1]),P.dec(b[9][1]),-1)
 for j in [1,2,3]:assert not P.add(P.mult(g,P.add(P.dec(b[11][j]),P.dec(b[9][j]),-1)),P.mult(hp,P.dec(b[6][j])),-1)
 c11=F(hc['known_row_units_rho14'][11]['constant']);assert hc['known_row_units_rho14'][11]==hc['known_row_units_rho14'][9]
 op=load(BASE,'operator.json');prev=load(OLD,'operator.json');N=op['nodes'];assert len(N)==78
 for i,n in enumerate(N):
  assert n['id']==i and all(j<i for j in n['args']);sh=n['shape'];args=n['args'];ar=[N[j]['shape']for j in args];typ=n['op']
  if typ in ['external_frozen_reference','external_integral_polynomial']:
   ref=n['reference'];p=ROOT/ref['path'];assert sha(p)==ref['sha256'];pins[ref['path']]=sha(p);z=cache.setdefault(ref['path'],json.loads(p.read_text()))
   for k in ref['json_pointer'].strip('/').split('/'):z=z[int(k)]if isinstance(z,list)else z[k]
   if ref['json_pointer'].startswith('/exports/'):assert prev['nodes'][z]['shape']==sh
   else:assert sh==[1,1]
  elif typ=='matrix_product':assert ar[0][1]==ar[1][0]and sh==[ar[0][0],ar[1][1]]
  elif typ in ['difference','sum']:assert ar==[sh,sh]
  elif typ=='central_scalar_multiple':assert ar[0]==[1,1]and ar[1]==sh
  elif typ in ['negative','transpose']:assert ar==[sh if typ=='negative'else sh[::-1]]
  elif typ=='submatrix':assert sh==[len(n['rows']),len(n['columns'])]and all(j in range(ar[0][0])for j in n['rows'])and all(j in range(ar[0][1])for j in n['columns'])
  elif typ in ['zero','identity']:assert not ar and(typ=='zero'or sh[0]==sh[1])
  elif typ in ['row_permutation','column_permutation']:assert not ar and sh[0]==sh[1]and sorted(n['order'])==list(range(sh[0]))
  elif typ in ['certified_chart_inverse','known_FULLU_parameter_unit_inverse']:assert ar==[[1,1]]and sh==[1,1]
  elif typ=='literal_known_FULLU_parameter_unit':assert not ar and sh==[1,1]and n['polynomial_monomial']==[9,0,0]and n['variable_order']==['a','e','d']and F(n['QQ_scalar'])==c11*F(content,D)
  elif typ=='block_matrix':
   rr=n['blocks'];assert args==[j for row in rr for j in row];S=[[N[j]['shape']for j in row]for row in rr];assert all(len({s[0]for s in row})==1 for row in S)and all(len({row[j][1]for row in S})==1 for j in range(len(S[0])))and sh==[sum(row[0][0]for row in S),sum(s[1]for s in S[0])]
  else:raise AssertionError(typ)
 def flat(i):return sum((flat(j)for j in N[i]['args']),[])if N[i]['op']=='matrix_product'else[i]
 def sub(i,p,rr,cc):assert N[i]['op']=='submatrix'and N[i]['args']==[p]and N[i]['rows']==rr and N[i]['columns']==cc
 for i,pointer in enumerate(['/exports/residual3x10','/exports/fullsource15','/exports/fullsource15_inverse','/exports/fulltarget22','/exports/fulltarget22_inverse']):assert N[i]['reference']['json_pointer']==pointer
 assert N[5]['reference']['json_pointer']=='/exact_normalized_D6'and N[6]['reference']['json_pointer']=='/polynomial'
 assert N[7]['args']==[6]and N[7]['chart']=='D(N7bar)'and N[9]['args']==[8]and flat(11)==[5,9,7]
 assert [i for i,n in enumerate(N)if n['op'] in ['certified_chart_inverse','known_FULLU_parameter_unit_inverse']]==[7,9]
 assert N[12]['order']==[1,0,2]and N[14]['order']==[7,0,1,2,3,4,5,6,8,9]and flat(17)==[12,0,14]
 sub(18,17,[0],[0]);sub(19,17,[0],list(range(1,10)));sub(20,17,[1,2],[0]);sub(21,17,[1,2],list(range(1,10)))
 assert N[22]['args']==[11,20]and N[23]['args']==[11,19]and N[25]['args']==[21,24]and flat(24)==[22,19]
 assert N[29]['blocks']==[[11,26],[27,28]]and N[27]['args']==[22]and N[32]['blocks']==[[18,30],[20,31]]
 assert N[37]['blocks']==[[33,34],[35,36]]and N[34]['args']==[23]and N[41]['blocks']==[[38,23],[39,40]]
 sub(43,42,[1,2],[0,1,2]);assert flat(42)==[29,12]
 assert N[48]['blocks']==[[44,46],[47,45]]and flat(45)==[29,12]and flat(49)==[48,1]
 assert N[54]['blocks']==[[50,52],[53,51]]and flat(51)==[13,32]and flat(55)==[2,54]
 assert N[60]['blocks']==[[56,58],[59,57]]and flat(57)==[14,37]and flat(61)==[3,60]
 assert N[66]['blocks']==[[62,64],[65,63]]and flat(63)==[41,15]and flat(67)==[66,4]
 for i,pointer in zip(range(68,73),['surviving3_source_original25','surviving3_corrected_full57','surviving3_corrected_original63','surviving3_actual_pure84','surviving3_raw7']):assert N[i]['reference']['json_pointer']=='/exports/'+pointer
 for i,old in zip(range(73,78),range(68,73)):assert flat(i)==[43,old]
 # Block identity and both left/right inverse products are universal in the
 # free noncommuting block algebra; no parameter sample or scalar expansion.
 U.nc_block_proof()
 assert op['source_remaining_current15']==[9,14]and op['target_remaining_original32']==[1,3,4,6,7,8,9,17,18]
 assert op['source_pivot_order_current15']==prev['source_pivot_order_current15']+[12]and op['target_pivot_original32']==prev['target_pivot_original32']+[14]
 assert op['chart']=='FULLU D(Delta407*Saux*T*N6*N7bar)'and op['all7_remaining_extra_tail_effects_retained']and op['new_expanded_parameter_coefficients']==0
 proof=dict(all14986_primitive_normalizations_literal=True,q=str(F(content,D)),a_power=9,all3_K120_middle_heavy_identities=True,all78_types_hashes_and_wiring=True,both_Schur_inverse_products=True,source_actual_current15_not_bare135=True,all_corrected57_63_raw7_and_pure84_carried_by_same_W=True,seven_remaining_tails_and_prior_tail2_effects_kept=True,source_ratio='s*N6=c11*N7',ratio_proof='Middle3 killed by K120; known raw4x7 proportionality permits the literal G20 pair replacement. Cramer numerator D7/g=a²N7 and old D6=a²N6; a² cancels. Dense-field identity extends to the original localized domain.',inverse='N6/(c11*q*a^9*N7bar)',no_new_g_D5_kappa_p_R_h22_inverse=True,N7_peer_scope='symbolic source polynomial NONZERO, not full QQ coefficient identity',samecokernel=[2,9],remainingrank2_NOT_proved=True)
 save('independent-normalization-source-and-Schur-proof.json',proof)
 for p in[Path(__file__),Path(U.__file__),Path(P.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 rep=dict(id=ID,all_passed=True,frozen=True,typed_nodes=78,samecokernel_shape=[2,9],chart=op['chart'],all_source_and_scalar_normalization_checks=True,no_new_parameter_coefficient_or_rank_search=True,full57_or239_NOT_proved=True,input_sha256=pins,evidence_sha256={'independent-normalization-source-and-Schur-proof.json':sha(OUT/'independent-normalization-source-and-Schur-proof.json')},seconds=time.monotonic()-st)
 save('report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
