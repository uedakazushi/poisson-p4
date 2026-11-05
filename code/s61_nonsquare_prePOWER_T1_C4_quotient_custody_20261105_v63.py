#!/usr/bin/env python3
"""Prepare immutable T1/C4 quotient input/code; no quotient coefficient evaluated."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_WORKER_CUSTODY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):
  p=pin(ROOT/'reports'/st/name);return json.loads(p.read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 peer=load('S61_NONSQUARE_PREPOWER_T1_C_HASSE_HORNER_REVIEW_20261105_v63');assert peer['all_passed']
 pre=load('S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_SOURCE_20261105_v63','finite-operator-and-preflight.json');assert pre['owned_MAX']==81776
 add=load('S61_NONSQUARE_PREPOWER_T1_C4_QUOTIENT_SOURCE_ADDENDUM_20261105_v63','proof.json');assert add['C4_valuation_independent_peer_fulfilled']and add['ALL_scaled_root_and_polynomial_products_upper']==2712113
 old=load('S61_NONSQUARE_PREPOWER_T1_C_HASSE_WORKER_CUSTODY_20261105_v63');assert old['all_passed']and old['value_worker_NOT_run']
 inp=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C_HASSE_WORKER_CUSTODY_20261105_v63/native-input.txt');assert sha(inp)==old['evidence_sha256'][inp.name]
 count=0;L=0;previous=None
 with inp.open()as f,(out/'native-input.txt').open('w')as g:
  assert f.readline().strip()=='T1CH001 7595 20 3';Dline=f.readline();D=int(Dline);assert D>0 and D.bit_length()==197
  g.write('T1C4Q001 7595 20 3 16\n');g.write(Dline)
  for line in f:
   key,c=map(int,line.split());assert c and(previous is None or previous<key);previous=key
   a,x=divmod(key,1048576);e,d=divmod(x,1024);assert 3<=a<=23 and 0<=e<=41 and 0<=d<=18 and a+e+d<=61
   g.write(line);count+=1;L+=abs(c)
 assert count==7595 and L.bit_length()<=363
 cpp=pin(ROOT/'code/s61_nonsquare_T1_C4_quotient_20261105_v63.cpp');drv=pin(ROOT/'code/s61_nonsquare_prePOWER_T1_C4_quotient_entry_20261105_v63.py');me=pin(Path(__file__))
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu -lflint','T1C4-build',str(cpp),str(out/'native-worker')]
 with (out/'compile.log').open('w')as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=30)
 for name,cmd in [('worker-help',[str(out/'native-worker'),'--help']),('driver-help',[str(ROOT/'.venv/bin/python'),str(drv),'--help'])]:
  r=subprocess.run(cmd,text=True,capture_output=True,check=True,timeout=5);(out/(name+'.log')).write_text(r.stdout+r.stderr)
 proof=dict(ONLY_T1=True,only_known_units_removed='a^3*C^4',no_quotient_values_evaluated=True,
  input_header_only_changed_body7595lines_identical=True,source_C4_independent_Horner_premise=sha(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C_HASSE_HORNER_REVIEW_20261105_v63/report.json'),
  formula='Qint_i=sum_{j>=i+4} I_j*binom(j-i-1,3)*R^(j-i-4)*100^(20-j+i); Dq=D1*100^20',
  negative_a_remainder='NOT silently dropped: the operator is Laurent polynomial-part; ALL coefficient identity C^4*Qint=100^20*U1int verifies the entire negative-a rational remainder is zero',
  source_slabs='do not truncate before whole slab; preserve all e<=73,d<=18 contributions; fail if completed nonzero outside provedQ a16/e33/d18/total50',
  final_Q_support_upper=10013,partial_slab_support_upper=1406,remainder_support_upper=14991,
  whole_owned_MAX=81776,formula_owned_MAX=54606,coefficient_key_cells_per_entry=2,ALL_metadata_scalar_reserve16000=16000,
  lifecycle='release ALL289 Rpower slots and empty slab BEFORE allocating Cpower cache/remainder',
  Cpower_actual_support_upper55=55,constant_cache_reserve289=289,
  input_l1_bits363=363,Qint_l1_bits472=472,full_multiplyback_partial_l1_bits502=502,
  proof='L1Q<4*100^16*L1I; ||C||1=129,||C^4||1<2^29; ||100^20 I||1+||C^4 Qint||1<2^502',
  hard_all_integer_temp_bits518=518,
  ALL_scalar_root_and_polynomial_product_upper=2712113,
  mathematical_guard_cursor='inputi BEFOREcommit; sourcei/inputj/rootr BEFOREsum; completedslab phase3 then transferkey BEFOREcopy; Cpowerh/key/t; inputremainderi; pendingQkey/Ckey before multiplyback; phase7 allidentitycompletebeforeoutput',
  prefix_version='T1C4QP001',all_active_arrays_in_prefix=['input','Rpower slots','Cpower maps','Q','whole slab','remainder','D1,inputL1','all phase/i/j/r/key cursors'],
  fixed_metadata_reconstruction='100powers0..20 andchoose(h+3,3)0..16 from literal constants; never unknown source reevaluation',
  no_second_product_or_accumulator_no_IO_doubleprefix=True,
  soft45_hard60_RAM2GB=True,no_GCD_otherfactor_T0_or_tail=True)
 dump('source-height-lifetime-and-guards.json',proof)
 raw=sum(p.stat().st_size for p in out.iterdir())+sum(p.stat().st_size for p in [cpp,drv,me])+20000
 total=raw+8410960+200000;assert total<10610960
 dump('RAW-closure.json',dict(actual_prepared_and_codes_plus_metadata=raw,full_state_upper8410960=8410960,all_future_receipts_upper200000=200000,total_upper=total,strict_hard_limit10610960=10610960,
  success_output_only_AFTER_whole_multiplyback=True,IO_error_no_second_full_prefix=True,one_state_no_atomic_doublebuffer=True))
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 r=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,value_worker_NOT_run=True,
  input_literal7595=True,independent_C4_premise_pinned=True,negative_a_remainder_bound=True,
  whole_owned81776=81776,hard_integer_bits518=518,RAW_total_upper=total,ALL_products_upper2712113=2712113,
  compiled_help_only=True,fresh_ROOT_source_code_audit_REQUIRED=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',r);signal.alarm(0);print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
