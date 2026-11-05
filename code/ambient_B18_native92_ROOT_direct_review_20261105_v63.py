#!/usr/bin/env python3
"""Independent direct universal-64 trace convolution for the fixed native92.

Do not form producer Gamma*L or intrinsic Res4(S) coefficient matrices.
Multiply N-mask, G8-mask and original input-mask simultaneously, then route
the unique monic collision. Split its location only to verify all322 pieces.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NATIVE92_ROOT_DIRECT_REVIEW_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v)for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,str):
  try:F(x);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(x,(int,float,F)))
def main():
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={}
 def read(p):q=ROOT/p;pins[p]=sha(q);return json.loads(q.read_text())
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 gate='reports/AMBIENT_B18_ELL2_NATIVE_PIVOT23_SECOND_SOURCE_20261105_v63_V2/'
 r=read(gate+'report.json');assert pins[gate+'report.json']=='ab2e8d901a32da8e2f0e410f86e5bf46512b74fcc3342ecc12edbb52ccb316a0'
 d=read(gate+'factored-native23-Gamma-pivot-plane2-source-and-full-L2S0.json')
 primary='reports/AMBIENT_B18_NATIVE_PIVOT23_SECOND_VALUES_20261105_v63/'
 pr=read(primary+'report.json');assert pins[primary+'report.json']=='8eea9b0748d2d230dc6495c910197db61f15e95b6d92721075bbc97191b4048d'
 actual=read(primary+'actual92-Gamma-nativepivot23-plane2-and322-full-product-pieces.json')
 assert r['all_passed']and pr['all_passed']
 mono=[(0,0),(1,0),(0,1),(2,0),(1,1)];index={m:i for i,m in enumerate(mono)}
 rridx={(0,3):0,(1,1):1,(3,0):2};reidx={(0,4):0,(1,2):1,(2,1):2,(4,0):3}
 # Numeric source vectors plus direct92/322 prefix and failure duplication.
 # No recursive global circuit or transformed L/S matrices are retained.
 source_count=count([r,d,pr,actual,pins])
 owner_prior=source_count+2880+160+10+112+2*414+4096
 assert owner_prior<21034,(source_count,owner_prior)
 gd=math.prod(int(x['D'])for x in d['source_positive_classes_D_U'].values())
 gu=math.prod(int(x['U'])for x in d['source_positive_classes_D_U'].values())*d['source_path_count_upper']
 assert max(gd.bit_length(),gu.bit_length())==280
 dump('pre-arithmetic-SOURCE.json',dict(input_sha256=pins,whole_numeric_prior=owner_prior,whole_numeric_cap=21034,source_bits=280,temporary_bits=841,hard_bits=[384,1153],RAW_cap=4_882_000,soft_seconds=45,hard_seconds=60,RAM=2_000_000_000,direct64_mask_rule='XOR output and floor per trace-bit; doublecollision starts degree4; full monic carries kept',no_physical_or_honest_or_upper_query=True))
 N=[[[[F(q)for q in mask]for mask in a]for a in row]for row in d['cached_N_12_by12_plane5_free4_UNEVALUATED_IN_PIVOT_RESPONSE']]
 G=[[[F(q)for q in mask]for mask in a]for a in d['canonical_G8_plane5_free4_UNEVALUATED_IN_PIVOT_RESPONSE']]
 nu=[[F(q)for q in a]for a in d['canonical_nu_two_plane5']]
 Gamma=[{i:F(q)for i,q in row if F(q)}for row in d['Gamma_two']]
 assert all(row[:3]==[F(0)]*3 for row in nu)
 results=[];pieces=[];ops=0;maximum=0;cursor={}
 def ck(q):
  nonlocal maximum
  maximum=max(maximum,abs(q.numerator).bit_length(),q.denominator.bit_length())
  if maximum>280 or time.monotonic()-start>45:
   dump('guard-frontier.json',dict(cursor=cursor,uncommitted_fraction=str(q),results=results,pieces=pieces,current_target=target,current_pieces=pp,input_sha256=pins,no_retry=True));raise RuntimeError('source280 or45second STOP')
  return q
 # Verify all64 routes from the integer trace-bit exponents themselves.
 routes={}
 for a in range(4):
  for c in range(4):
   for b in range(4):
    powers=[((a>>j)&1)+((c>>j)&1)+((b>>j)&1)for j in range(2)]
    om=sum((v%2)<<j for j,v in enumerate(powers));collision=sum((v//2)<<j for j,v in enumerate(powers))
    assert om==a^b^c
    s_carry=c&b;l_carry=a&(c^b)
    assert collision==s_carry|l_carry and not(s_carry&l_carry)
    routes[a,c,b]=(om,collision,bool(s_carry))
 for ri,recipe in enumerate(d['native23_pivot_recipes']):
  target=[[F(0),F(0)]for _ in range(2)];pp=[[[F(0)]*3,[F(0)]*4]for _ in range(2)];b=recipe['input_trace_mask']
  for entry in recipe['S_source_entries']:
   if entry[0]=='G8':_,source,g,sign=entry;field=G[g];sign=F(sign)
   else:
    _,source,constant=entry;field=[[F(constant),F(0),F(0),F(0)]]+[[F(0)]*4 for _ in range(4)];sign=F(1)
   for physical in range(12):
    for ni,nm in enumerate(mono):
     for a,nv in enumerate(N[source][physical][ni]):
      if not nv:continue
      for gi,gm in enumerate(mono):
       for c,gv in enumerate(field[gi]):
        if not gv:continue
        om,collision,in_S=routes[a,c,b]
        if collision==3:continue
        carry=[(0,F(1))]if not collision else[(k,q)for k,q in enumerate(nu[collision-1])if q]
        for k,q in carry:
         km=mono[k];combined=(nm[0]+gm[0]+km[0],nm[1]+gm[1]+km[1])
         if combined not in [(2,0),(1,1)]:continue
         channel=0 if combined==(2,0)else 1
         lm=(nm[0]+(0 if in_S else km[0]),nm[1]+(0 if in_S else km[1]))
         sm=(gm[0]+(km[0]if in_S else 0),gm[1]+(km[1]if in_S else 0))
         pi=(rridx if channel==0 else reidx)[index[lm],index[sm]]
         for t in range(2):
          gamma=Gamma[t].get(4*physical+om,F(0))
          if not gamma:continue
          cursor=dict(native=ri,source12=source,physical12=physical,Nbase=ni,Nmask=a,Gbase=gi,Gmask=c,inputmask=b,carrybase=k,channel=t,plane=channel,piece=pi)
          v=ck(ck(ck(ck(nv*gv)*sign)*q)*gamma);ops+=5
          target[t][channel]=ck(target[t][channel]+v);pp[t][channel][pi]=ck(pp[t][channel][pi]+v);ops+=2
  tv=[[str(q)for q in row]for row in target];pv=[[[str(q)for q in row]for row in t]for t in pp]
  assert tv==actual['Gamma_Xp2_native23_twochannels_rr_rell'][ri],ri
  assert pv==actual['all23_two_by2_full_product_piece_values'][ri],ri
  for t in range(2):
   for ch in range(2):
    q=target[t][ch];assert q==sum(pp[t][ch],F(0))and gd%q.denominator==0 and abs(q.numerator)*gd<=gu*q.denominator
  results.append(tv);pieces.append(pv)
 assert len(results)==23 and ops<500_000
 ep='direct-universal64-native92-and322-pieces.json';dump(ep,dict(native_sources=actual['original768_native23_sources'],Gamma_Xp2=results,product_pieces=pieces,all92_and322_match=True))
 report=dict(id=STEM,all_passed=True,independent_method='direct triplemask floor-bitcarry convolution before Gamma, no factorGammaL or intrinsicRes4S matrices',all92_and322_exact=True,native648_both_channels_zero=(results[2]==[['0','0'],['0','0']]),source_positive_bits280_temporary841=True,actual_max_bits=maximum,whole_numeric_prior=owner_prior,numeric_cap=21034,operations=ops,universal_routes_checked=64,no_new_honest_physical_gauge_N3_or_upper_values=True,input_sha256=pins,code_sha256=sha(Path(__file__)),evidence_sha256={ep:sha(out/ep),'pre-arithmetic-SOURCE.json':sha(out/'pre-arithmetic-SOURCE.json')},seconds=time.monotonic()-start)
 dump('report.json',report);print(json.dumps(report))
if __name__=='__main__':main()
