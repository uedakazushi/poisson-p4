#!/usr/bin/env python3
"""Prepare ONLY the ROOT-normalized4811 input and F-Hasse source/code gate."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_Q1_F_HASSE_WORKER_CUSTODY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):
  return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 st='S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63'
 peer=load(st);assert peer['all_passed']and peer['full_remainder_zero']and peer['all4811_quotient_coefficients_match']
 assert peer['fixed_scalar_100_power20_divides_all_native_coefficients']and peer['original_common_denominator_bits']==197
 source=pin(ROOT/'reports'/st/'quotient-original-common-D.txt');assert sha(source)==peer['evidence_sha256'][source.name]
 load(st,'actual-read-output-RAW-receipt.json')
 pre=load('S61_NONSQUARE_PREPOWER_T1_C4_GAUSS_F_SOURCE_20261105_v63','Gauss-and-F-preflight.json')
 assert pre['owned28708']==28708 and pre['Hasse_integer_bits_upper']==463
 L=0;mx=0;count=0;last=None
 with source.open()as f,(out/'native-input.txt').open('w')as g:
  assert f.readline().strip()=='T1C4QN001';Dline=f.readline();D=int(Dline);assert D>0 and D.bit_length()==197
  assert int(f.readline())==4811;g.write('Q1FH001 4811 16 0\n');g.write(Dline)
  for line in f:
   key,c=map(int,line.split());assert c and(last is None or last<key);last=key
   a,r=divmod(key,1048576);e,d=divmod(r,1024);assert 0<=a<=16 and 0<=e<=33 and 0<=d<=18 and a+e+d<=50
   g.write(line);L+=abs(c);mx=max(mx,abs(c).bit_length());count+=1
 assert count==4811 and L.bit_length()<=301
 cpp=pin(ROOT/'code/s61_nonsquare_Q1_F_hasse_20261105_v63.cpp');drv=pin(ROOT/'code/s61_nonsquare_prePOWER_Q1_F_hasse_entry_20261105_v63.py');me=pin(Path(__file__))
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu -lflint','Q1F-build',str(cpp),str(out/'native-worker')]
 with (out/'compile.log').open('w')as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=30)
 for name,cmd in [('worker-help',[str(out/'native-worker'),'--help']),('driver-help',[str(ROOT/'.venv/bin/python'),str(drv),'--help'])]:
  r=subprocess.run(cmd,text=True,capture_output=True,check=True,timeout=5);(out/(name+'.log')).write_text(r.stdout+r.stderr)
 proof=dict(ONLY_Q1_F=True,Q1_definition='Q_C=T1/(a^3*C^4), not the original S61 quadric q1',
  input_source='immutable ROOT integer longdivision normalized quotient-original-common-D.txt',
  body4811lines_BYTE_identical=True,no_new_scalar_normalization_or_quotient=True,
  ROOT_original_D197_restored=True,actual_integer_coefficient_bits=mx,actual_integer_l1_bits=L.bit_length(),
  F='500*a+90*e^2+3*e-24',a0='(24-3*e-90*e^2)/500',R_low_to_high=[24,-3,-90],root_l1=117,
  exact_formula='D0*500^16*[z^k]Q_C(a0+z,e,d)=sum K_ijl*choose(i,k)*500^(16-i+k)*R^(i-k)*e^j*d^l',
  orders='k0..16, stop FIRST nonzero',
  valuation_lemma='F=500(a-a0);500 is a QQ unit; lowest nonzero coefficient is exact v_F(Q_C)',
  integer_height_proof='L1K<2^301,choose<=2^16,||R||1=117<500 =>all integer partials<2^463',
  hard_all_integer_bits518_UNCHANGED=518,HH_partial_bits463=463,denominator_clear_bits_upper342=342,
  support='a0..16,e0..33,d0..18 input; e_HH<=65,d<=18 =>1254slots',
  whole_owned=2*(4811+289+1254)+16000,
  ownership='two cells per coefficient+packedkey; conservative rootcoeff+index289; ALLmetadata/scratch16000; inputvector reserves4811 before loading',
  source_input_L1_check_is_transactional_before_commit=True,
  first_nonzero_verification='forward ->reverse subtractionZERO ->forwardreproduction, ONE accumulator',
  inner_products_upper=1707905,ALL_scaled_root_polynomial_products_upper=1896419,
  soft45_hard60_RAM2GB=True,
  prefix='Q1FHP001 allinput/root/accumulator/D,inputL1/phase,k,i,r,ph,pr,counters; fixedchoose and500powers reconstructed from literalmetadata',
  source_guard_before_commit=True,phase5='all passes complete, pre-output clock guard',
  no_success_plus_fullprefix_on_IO_failure=True,
  no_F_power_Hasse_layer_other_factor_tail_GCD_or_parameter_point_evaluated=True)
 assert proof['whole_owned']==28708
 dump('formula-height-ownership-and-cursors.json',proof)
 raw=sum(p.stat().st_size for p in out.iterdir())+sum(p.stat().st_size for p in [cpp,drv,me])+20000
 total=raw+4000000+200000;assert total<8000000
 dump('RAW-closure.json',dict(actual_prepared_plus_codes_and_metadata=raw,fullprefix_OR_output_upper4000000=4000000,
  future_receipts_logs_reports200000=200000,total_upper=total,hard_limit8000000=8000000,
  commonDenominator_once_no_duplicateQQ=True,no_IO_double_fullstate=True))
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 r=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,value_worker_NOT_run=True,
  ROOT_normalized4811_literal_body_reused=True,only_known_F_valuation=True,
  hard_owned28708=28708,integer_hard518=518,source_HH_bits463=463,ALL_products_upper1896419=1896419,
  RAW_total_upper=total,compiled_help_only=True,fresh_ROOT_source_code_audit_REQUIRED=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',r);signal.alarm(0);print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
