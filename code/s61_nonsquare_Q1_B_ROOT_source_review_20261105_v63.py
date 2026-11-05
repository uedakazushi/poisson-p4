#!/usr/bin/env python3
"""ROOT source/code audit for ordinary named-B valuation, no root or Hasse values."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B_HASSE_WORKER_CUSTODY_20261105_v63'
OUT=ROOT/'reports/S61_NONSQUARE_PREPOWER_Q1_B_HASSE_ROOT_SOURCE_REVIEW_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(exist_ok=False);p=BASE/'report.json';assert sha(p)=='1d90c857126e5f209a5f4f495b8b50aaed903de78624c8fadd01ed39c434c146'
 r=json.loads(p.read_text());assert r['all_passed'] and r['value_worker_NOT_run'];pins={str(p.relative_to(ROOT)):sha(p)}
 for group in ['input_sha256','evidence_sha256']:
  for name,want in r[group].items():
   p=ROOT/name if group=='input_sha256'else BASE/name;assert sha(p)==want;pins[str(p.relative_to(ROOT))]=want
 native=(BASE/'native-input.txt').read_bytes().splitlines(keepends=True)
 old=(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt').read_bytes().splitlines(keepends=True)
 assert native[0]==b'Q1BH001 4811 18 0\n' and native[1]==old[1] and native[2:]==old[3:]
 D=int(native[1]);degrees=[];l1=0
 for line in native[2:]:
  key,c=map(int,line.split());a,rem=divmod(key,1048576);e,d=divmod(rem,1024)
  assert 0<=a<=16 and 0<=e<=33 and 0<=d<=18 and a+e+d<=50 and c
  l1+=abs(c);degrees.append(d)
 assert len(degrees)==4811 and l1<2**301 and D.bit_length()==197
 # Exact small polynomial proof of20b1+Rw=1000a; no b1 inverse/chart.
 b1={(1,0):50,(0,2):18,(0,1):-9,(0,0):-2};RR={(0,0):40,(0,1):-60};w={(0,1):6,(0,0):1}
 det={k:20*v for k,v in b1.items()}
 for x,c in RR.items():
  for y,v in w.items():k=(x[0]+y[0],x[1]+y[1]);det[k]=det.get(k,0)+c*v
 assert{k:v for k,v in det.items()if v}=={(1,0):1000}
 assert -1*(-9)+4*(-2)==1
 cache=[19*(19-h)for h in range(19)];assert sum(cache)==3610
 layer=[sum(1 for a in range(35)for e in range(70)if a+e<=86-k)for k in range(19)];assert max(layer)==2297
 assert 2*(4811+3610+2297)+16000==37436
 assert 2*(4811+3610+190+361+324)+16000==34592
 assert 2**301*179**18<2**436 and 100**18<2**120
 inner=[sum(cache[d-k]for d in degrees if d>=k)for k in range(19)]
 eligible=[sum(d>=k for d in degrees)for k in range(19)]
 construction=342+sum((t+1)**2*(19-t)for t in range(19))+4*sum((t+1)**2 for t in range(18));assert construction==22078
 schedules=[construction+sum(inner[:k+1])+2*inner[k]+sum(eligible[:k+1])+2*eligible[k]for k in range(19)]
 assert max(schedules)==12768337
 cpp=(ROOT/'code/s61_nonsquare_Q1_B_hasse_20261105_v63.cpp').read_text()
 driver=(ROOT/'code/s61_nonsquare_prePOWER_Q1_B_hasse_entry_20261105_v63.py').read_text()
 for text in ['in.reserve(n)','magic!="Q1BH001"||n!=4811||dm!=18||shift!=0','int rr[2]={40,-60}',
              'bcoef[4]={50,18,-9,-2}','hc=18-t','update(S[hc],pkey+r,v,ns,0)',
              'R.clear();R.shrink_to_fit();nr=0;P.clear();np=0;N.clear();nn=0',
              'hc=d-k','choose[d][k]','update(A,offset+key,out,na,2)','a+e<=86-kc',
              'convolution(k,6)','convolution(k,7)','base_ready','Q1BHP001',
              'symbolic_b1_power18','b1_NOT_assumed_FULLU_unit',
              'std::stoull(argv[4])+5000000+200000>8000000','rlimit ram{2000000000,2000000000},cpu{60,60}']:
  assert text in cpp,text
 assert 'timeout=58' in driver and r['RAW_total_upper']==6169921<8000000
 facts=dict(input4811_and_D197_byte_identical=True,
  ordinary_valuation_lemma='B=b1*d-R primitive in QQ[a,e][d]; any common factor of b1,R divides1000a, while a does not divide R. Thus coefficientgcd1, prime B by Gauss. Local fraction root detects exact ordinary B order. Auxiliary binary matrix determinant1000a is verified, b1 is never a FULLUunit assumption.',
  primitive_ZZ_identity=True,exact_identity='P_k=sum I_ijl binomial(l,k)R^(l-k)b1^(18-l+k)a^i e^j; H_k=P_k/(D*b1^18)',
  symbolic_b1_denominator_retained=True,cache3610=3610,accumulator2297=2297,whole_owned37436=37436,
  construction_owned34592=34592,source_integer_height436=436,hard518=518,all_products_upper=12768337,
  RAW_upper=6169921,code_audit='one Rcache/current-next b1 workspace; all19 S cached then release R/P/N; firstnonzero forward/reverse/reproduction in one accumulator; full15field prefix incl activebase and all maps; no double full I/O output',
  no_root_cache_Hasse_Bvalue_quotient_or_GCD_evaluated=True)
 (OUT/'source-proof.json').write_text(json.dumps(facts,indent=2)+'\n');pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 result=dict(id=OUT.name,all_passed=True,frozen=True,SOURCE_only=True,ONE_Q1_B_valuation_ready=True,
  input_sha256=pins,evidence_sha256={'source-proof.json':sha(OUT/'source-proof.json')})
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(facts,indent=2));print('report_sha256',sha(OUT/'report.json'))
if __name__=='__main__':main()
