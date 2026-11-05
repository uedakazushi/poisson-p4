#!/usr/bin/env python3
"""SOURCE/code/custody preparation ONLY for one Q1 named-J case; compile/help, no J values."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_Q1_J_HASSE_WORKER_CUSTODY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 source=load('S61_NONSQUARE_PREPOWER_TWO_TAIL_CANONICAL_SOURCE_20261105_v63');assert source['all_passed']and source['source_only']
 pre=load('S61_NONSQUARE_PREPOWER_TWO_TAIL_CANONICAL_SOURCE_20261105_v63','ONE-known-J-case-PREFLIGHT-ONLY.json')
 assert pre['whole_owned29336']==29336 and pre['source_integer_bits372']==372
 peer=load('S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63')
 assert peer['all_passed']and peer['full_remainder_zero']and peer['all4811_quotient_coefficients_match']
 inp=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt')
 assert sha(inp)==peer['evidence_sha256'][inp.name]
 L=0;mx=0;count=0;last=None;hist=[0]*19;body=hashlib.sha256()
 with inp.open()as f,(out/'native-input.txt').open('w')as g:
  assert f.readline().strip()=='T1C4QN001';Dline=f.readline();D=int(Dline);assert D>0 and D.bit_length()==197
  assert int(f.readline())==4811;g.write('Q1JH001 4811 18 0\n');g.write(Dline)
  for line in f:
   key,c=map(int,line.split());assert c and(last is None or last<key);last=key
   a,x=divmod(key,1048576);e,d=divmod(x,1024)
   assert 0<=a<=16 and 0<=e<=33 and 0<=d<=18 and a+e+d<=50
   g.write(line);body.update(line.encode());L+=abs(c);mx=max(mx,abs(c).bit_length());count+=1;hist[d]+=1
 assert count==4811 and L.bit_length()<=301
 # Count exact source index multiplicities; no root polynomial coefficients or Hasse values.
 H=[sum(n*(l-k+1)for l,n in enumerate(hist)if l>=k)for k in range(19)]
 N=[sum(n for l,n in enumerate(hist)if l>=k)for k in range(19)]
 all_schedule=[sum(H[:k+1])+2*H[k]+2*(sum(N[:k+1])+2*N[k])+342+18+1 for k in range(19)]
 assert max(all_schedule)==441936 and max(all_schedule)<=1200000
 counts=[sum(a+e<=68-k for a in range(35)for e in range(52-k))for k in range(19)]
 assert max(counts)==1667 and 2*(4811+190+1667)+16000==29336
 cpp=pin(ROOT/'code/s61_nonsquare_Q1_J_hasse_20261105_v63.cpp')
 drv=pin(ROOT/'code/s61_nonsquare_prePOWER_Q1_J_hasse_entry_20261105_v63.py');pin(Path(__file__))
 note=pin(ROOT/f'research/notes/{ID}.md')
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu -lflint','Q1J-build',str(cpp),str(out/'native-worker')]
 with(out/'compile.log').open('w')as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=30)
 for name,cmd in [('worker-help',[str(out/'native-worker'),'--help']),('driver-help',[str(ROOT/'.venv/bin/python'),str(drv),'--help'])]:
  r=subprocess.run(cmd,text=True,capture_output=True,check=True,timeout=5);(out/(name+'.log')).write_text(r.stdout+r.stderr)
 proof=dict(ONLY_canonical_Q1_named_J=True,value_worker_NOT_run=True,
  Q1='T1/(a3*C4)=canonical4811 numeratorI/D197; not the original S61 quadric q1',
  unchanged4811body_sha256=body.hexdigest(),header_only_changed=True,no_a_C_or_otherfactor_redivision=True,
  input_actual_maxbits=mx,input_actual_L1_bits=L.bit_length(),input_D197=D,
  J='5*a*d-6*e-4',root='d0=(6e+4)/(5a)',known_denominator_unit='5a, named FULLU a only',
  primitive_irreducible_proof='coefficientgcd(5a,6e+4)=1 because a does not divide nonzero6e+4; primitive linear in d irreducible over QQ(a,e), hence irreducible in QQ[a,e,d] byGauss. This is gcd1, NOT a claim the two coefficients generate the unitideal.',
  valuation_proof='in QQ(a,e)[d],J=5a(d-d0), so first nonzero Hasseorder is exactvJ; primitiveness makes this the ordinary QQ[a,e,d] multiplicity as well. All0..18 cannot vanish for nonzero degree18input.',
  fixed_formula='P_k=sum_{l>=k}I_ijl*choose(l,k)*(4+6e)^(l-k)*5^(18-l+k)*a^(i+18-l+k)*e^j',
  exact_denominator='D197*5^18*a^18; numericclearD*5^18 plus SYMBOLIC externala exponent18, never a numerical evaluation',
  exact_identity='[z^k]Q1(a,e,d0+z)=P_k/(D197*5^18*a^18)',
  orders='0..18 stop FIRST nonzero only',
  output_support='a<=34,e<=51-k,a+e<=68-k; a and e nonnegative after fixedmonomialclear',all19_support_bounds=counts,
  input_terms4811=4811,root_cache190=190,single_accumulator_max1667=1667,
  whole_owned29336=29336,
  owned_equation='2*(4811+190+1667)+16000; coefficient+packedkey countedseparately; input reservedonce, no expectedpolynomial or secondaccumulator',
  ALL_numeric_metadata16000=16000,
  metadata='19x19choose,19powers5,sourceexponentdecode/key/cursors/counters,temporaryZs,L1/D/denominator records,prefix/output formatting scratch',
  all_source_integer_height372=372,hard_all_integer518=518,
  height='L1I<2^301. For fixedk, coefficient/binomial/rootweight<=5^(18-l)*15^l<=15^18. Thus allpartial sums/sourceproducts<2^372; reverse subtraction is the remaining subset, no 3x acc copy. Numericdenominator<=239bits.',
  integer_rootcoeff_product_scope='root(4+6e)^h constants only; NO sourceQ1 powers',
  source_d_histogram=hist,onepass_inner_products_by_k=H,eligible_input_counts_by_k=N,
  all_schedule_firstnonzero_k0to18_upper=all_schedule,
  ALL_products_scope_upper441936=max(all_schedule),existing_ALL_products_bound1200000_UNCHANGED=True,
  products_scope='all orders through FIRSTnonzero plus two repeat passes at thatorder; inner source*root, two scale products pereligibleinput,342rootcache,18fivepowers,1numericDclear',
  verification='forwardfirstnonzero, reversesource/reverseroot subtractionZERO, forwardreproduction; oneaccumulator and allcoeffreadback. Independentdifferentarithmetic peer remains required afteractual.',
  prefix_version='Q1JHP001',allstate='input/root arrays/accumulator/D/L1; phase,k,inputi,rootr,rootcacheh/j,products/scalings/updates/peak; completedzeroorders externalfile retained',
  fixed_restore='choose19x19 and5powers0..18 regeneratedfromfixedliteral metadata only; no oldunknown source replay',
  guard_order='owned/bit/L1 checkedbeforecommit; phase5 afterallthreepasses clockguard beforeoutput; mathematicalSTOP savesone fullprefix',
  IO='successoutput ONLY afterreverse+reproduction; IO failure never triggerssecond fullprefix; noatomicdoublebuffer',
  soft45_hard60_RAM2GB=True,no_Jvalue_quotient_GCD_otherfactor_tail_or_badlocus=True)
 dump('fixed-J-formula-height-lifetime-and-cursors.json',proof)
 files={str(p.relative_to(ROOT)):p for p in out.iterdir()};files.update({n:ROOT/n for n in pins})
 raw=sum(p.stat().st_size for p in files.values())+20000;total=raw+4000000+200000;assert total<=8000000
 dump('RAW-closure.json',dict(actual_prepared_all_read_sources_codes_metadata=raw,fullprefix_OR_output_upper4000000=4000000,
  all_future_receipts_logs_reports200000=200000,total_upper=total,hard8000000=8000000,
  no_extra_body_copy_or_IO_doubleprefix=True))
 r=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,value_worker_NOT_run=True,
  canonical4811_literal_BODY_unchanged=True,ordinary_ring_Jvaluation_lemma=True,
  source372_hard518=True,whole_owned29336=29336,ALL_products_upper441936=441936,
  existing1200000_product_bound_UNCHANGED=True,RAW_total_upper=total,
  compiled_help_ONLY=True,fresh_ROOT_source_code_audit_and_ONE_GO_REQUIRED=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',r);signal.alarm(0)
 print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
