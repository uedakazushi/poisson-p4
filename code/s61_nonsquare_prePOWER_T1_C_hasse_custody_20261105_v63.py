#!/usr/bin/env python3
"""SOURCE ONLY: T1's copied7595 integers and prepared known-C Hasse worker."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T1_C_HASSE_WORKER_CUSTODY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):
  p=pin(ROOT/'reports'/st/name);return json.loads(p.read_text())
 def dump(n,d):(out/n).write_text(json.dumps(d,indent=2)+'\n')
 old=load('S61_NONSQUARE_PREPOWER_TAIL1_CACHED_KAPPA_ENTRY_20261105_v63');assert old['all_passed']
 peer=load('S61_NONSQUARE_PREPOWER_TAIL1_AXIS_REVIEW_20261105_v63');assert peer['all_passed']and peer['T1_polynomial_nonzero_independently_proved']
 pre=load('S61_NONSQUARE_PREPOWER_TAIL1_FULL57_UNION_SOURCE_20261105_v63','NEXT-known-C-valuation-PREFLIGHT-ONLY.json')
 assert pre['integer_coefficient_bits_upper']==518 and pre['owned_coefficient_key_metadata_upper']==35188
 source=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_TAIL1_CACHED_KAPPA_ENTRY_20261105_v63/integer-column4-and-result.txt')
 assert sha(source)==old['evidence_sha256'][source.name]
 nn=out/'native-input.txt';L=0;bits=0;last=None;count=0
 with source.open()as f,nn.open('w')as g:
  assert f.readline().strip()=='J1KV001'
  for j in range(5):
   node,D=map(int,f.readline().split());n=int(f.readline())
   if j<4:
    for _ in range(n):f.readline()
    continue
   assert node==155 and n==7595 and D.bit_length()==197
   g.write(f'T1CH001 7595 20 3\n{D}\n')
   for _ in range(n):
    line=f.readline();k,c=map(int,line.split());assert c and(last is None or last<k);last=k
    a,r=divmod(k,1048576);e,d=divmod(r,1024);assert 3<=a<=23 and 0<=e<=41 and 0<=d<=18 and a+e+d<=61
    # Exact stored text bytes, no new coefficient rescaling or substitution.
    g.write(line);L+=abs(c);bits=max(bits,abs(c).bit_length());count+=1
  assert not f.read().strip()
 assert count==7595 and L.bit_length()<=359
 cpp=pin(ROOT/'code/s61_nonsquare_T1_C_hasse_20261105_v63.cpp')
 drv=pin(ROOT/'code/s61_nonsquare_prePOWER_T1_C_hasse_entry_20261105_v63.py');me=pin(Path(__file__))
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu -lflint','T1-C-build',str(cpp),str(out/'native-worker')]
 with (out/'compile.log').open('w')as f:subprocess.run(cmd,check=True,stdout=f,stderr=subprocess.STDOUT,timeout=30)
 for name,cmd in [('worker-help',[str(out/'native-worker'),'--help']),('driver-help',[str(ROOT/'.venv/bin/python'),str(drv),'--help'])]:
  r=subprocess.run(cmd,text=True,capture_output=True,check=True,timeout=5);(out/(name+'.log')).write_text(r.stdout+r.stderr)
 formula=dict(ONLY_T1=True,input='literal7595 common-D integer coefficients of stored T1',a_shift3='U1=T1/a^3; C does not divide a, so v_C(T1)=v_C(U1)',
  C='100*a+18*e^2-9*e-2',a0='(-18*e^2+9*e+2)/100',root_integer_coefficients_low_to_high=[2,9,-18],
  scalar_formula='D1*100^20*[z^k]U1(a0+z,e,d)=sum I_ijl*choose(i,k)*100^(20-i+k)*R(e)^(i-k)*e^j*d^l',
  each_Hasse_order='0..20 sequential, stop FIRST nonzero; no C quotient materialized',
  exact_valuation_lemma='C=100(a-a0) and100 is a QQ unit; the first nonzero z coefficient has exact C order',
  proof_Hasse_coefficient_bits='L1(input)<2^363, choose<=2^20, ||R||1=29<100, hence all integer partial sums <=L1*2^20*100^20 <2^518',
  proof_bivariate_support='e<=41+2*20=81,d<=18 =>1558 slots',
  proof_root_cache='sum_h(2h+1),h0..20=441 slots; no root powers evaluated during preparation',
  strict_owned=2*(7595+441+1558)+16000,
  strict_owned_detail='inputcoefficient+packedkey15190; conservative rootcoefficient+index882; acccoefficient+key3116; allmetadata/scratch16000',
  coefficient_bits518=518,denominator_bits_upper=331,
  first_pass_products_upper=7595*441,
  verification_first_nonzero_only_products_upper=2*7595*41,
  combined_products_upper=7595*441+2*7595*41,
  prefix='T1CHP001: full input vector, all root arrays inclzero slots, current accumulator, D/inputL1, phase/k/i/r and root-build cursors; fixed choose/100-powers reconstructed from literal metadata',
  no_second_accumulator_or_full_product=True,
  prefix_pending='convolution i/r points to the UNCOMMITTED next scalar term; phase5 full layer complete, pre-output guard',
  completed_zero_orders='separate small text file; prior k all complete before next layer',
  mathematical_guards_before_commit=True,output_IO_failure='bounded partial output only, no second full prefix',
  no_GCD_factor_other_tail_or_parameter_point=True)
 assert formula['strict_owned']==35188 and formula['combined_products_upper']<6528123
 dump('formula-ownership-and-cursors.json',formula)
 dump('literal-custody.json',dict(input_source_sha256=sha(source),source_node155=True,
  copied_count7595=7595,coefficient_and_key_lines_BYTE_identical=True,
  maximum_integer_coefficient_bits=bits,input_integer_l1_bits=L.bit_length(),common_denominator_bits197=197,
  coefficient_substitution_NOT_run=True,root_power_polynomials_NOT_run=True))
 raw=sum(p.stat().st_size for p in out.iterdir())+sum(p.stat().st_size for p in [cpp,drv,me])+20000
 total=raw+2200000+100000;assert total<19442896
 dump('RAW-closure.json',dict(actual_prepared_and_codes_plus_metadata_reserve=raw,prefix_OR_success_upper=2200000,
  reports_logs_reserve=100000,total_upper=total,hard_limit19442896=19442896,
  alloutputs_common_denominator_once_no_QQ_duplicate=True,one_state_file_no_atomic_doublebuffer=True))
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 report=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,ONLY_T1=True,
  copied_literal_coefficients7595=7595,exact_C_root_normalization_pinned=True,
  no_root_power_or_Hasse_layer_evaluated=True,value_worker_NOT_run=True,
  strict_owned35188=35188,coefficient_bits518=518,integer_products_upper=formula['combined_products_upper'],
  hard_seconds60=60,soft_seconds45=45,RAM2GB=2000000000,RAW_total_upper=total,
  worker_compiled_help_only=True,independent_source_code_audit_REQUIRED_before_any_valuation=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
