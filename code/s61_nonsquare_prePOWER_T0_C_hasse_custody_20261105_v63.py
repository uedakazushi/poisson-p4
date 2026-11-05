#!/usr/bin/env python3
"""Prepare fixed saved-T0 C-Hasse input/code only; no valuation execution."""
import argparse,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T0_C_HASSE_WORKER_CUSTODY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):
  return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 sourcegate=load('S61_NONSQUARE_PREPOWER_T0_C_HASSE_SOURCE_20261105_v63');assert sourcegate['all_passed']and sourcegate['C_valuation_NOT_computed']
 plan=load('S61_NONSQUARE_PREPOWER_T0_C_HASSE_SOURCE_20261105_v63','T0-literal-frame-and-preflight.json');D=int(plan['common_D197']);assert D.bit_length()==197
 old=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63');assert old['all_passed']and old['output_terms']==7208
 peer=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63_REVIEW');assert peer['all_passed']
 p=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63/root239.poly');assert sha(p)==old['evidence_sha256'][p.name]
 jp=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63/root239.json');assert sha(jp)==old['evidence_sha256'][jp.name]
 n=L=mx=0;seen=set();bodyhash=hashlib.sha256()
 with p.open()as f,(out/'native-input.txt').open('w')as g:
  offset=int(f.readline());assert offset==4;assert int(f.readline())==7208
  g.write(f'T0CH001 7208 20 0\n{D}\n')
  for line in f:
   a,e,d,c=line.split();a,e,d=map(int,[a,e,d]);cc=Q(c);assert cc and D%cc.denominator==0
   scale=D//cc.denominator;assert scale*cc.denominator==D
   I=cc.numerator*scale;assert I
   assert 0<=a<=20 and 0<=e<=40 and 0<=d<=18 and a+e+d<=57
   key=a*1048576+e*1024+d;assert key not in seen;seen.add(key)
   # BODY IS ALREADY U0=T0/a^4. No extra a subtraction is permitted.
   text=f'{key} {I}\n';g.write(text);bodyhash.update(text.encode())
   n+=1;L+=abs(I);mx=max(mx,abs(I).bit_length())
 assert n==7208 and L.bit_length()==326 and mx==320
 cpp=pin(ROOT/'code/s61_nonsquare_T0_C_hasse_20261105_v63.cpp');drv=pin(ROOT/'code/s61_nonsquare_prePOWER_T0_C_hasse_entry_20261105_v63.py');me=pin(Path(__file__))
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu -lflint','T0C-build',str(cpp),str(out/'native-worker')]
 with (out/'compile.log').open('w')as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=30)
 for name,cmd in [('worker-help',[str(out/'native-worker'),'--help']),('driver-help',[str(ROOT/'.venv/bin/python'),str(drv),'--help'])]:
  r=subprocess.run(cmd,text=True,capture_output=True,check=True,timeout=5);(out/(name+'.log')).write_text(r.stdout+r.stderr)
 proof=dict(ONLY_OLD_T0=True,T0='oldT/root239=B*kappa*tail2',old_source_full_QQ_peer_pinned=True,
  literal_poly_header=dict(offset=4,terms=7208),
  body_already_U0='bodyexponents mean a0..20 of U0=T0/a^4; variable a itself is unchanged',
  future_worker_header_shift0=True,forbidden='subtract4 again or treat bodya0 as originalT0a0',
  literal7208_monomial_order_unchanged=True,normalized_body_sha256=bodyhash.hexdigest(),
  fixed_QQ_conversion='d divides D; scale=D/d; scale*d=D; I=n*scale; commonD197. No polynomial/contentGCD or rational cross-products.',
  input_exact_integer_L1_bits326=326,input_max_integer_bits320=320,
  C='100a+18e^2-9e-2',a0='(2+9e-18e^2)/100',root_low_to_high=[2,9,-18],root_l1=29,
  exact_formula='D*100^20*[z^k]U0(a0+z,e,d)=sum I_ijl*choose(i,k)*100^(20-i+k)*R^(i-k)*e^j*d^l',
  exact_valuation='firstnonzero k0..20 equalsv_C(U0)=v_C(T0); C doesnot dividea',
  signed_height='L1input<2^326,choose<=2^20,29<100 =>all integerpartials<2^481',
  HH_bits481=481,all_integer_hard518=518,denominator_bits_upper331=331,
  input7208=7208,root_cache441=441,one_accumulator1539=1539,
  eHH_max80=80,dHH_max18=18,
  whole_owned34376=2*(7208+441+1539)+16000,
  ownership='two cells per coefficient+packedkey; conservative rootcoefficient+index; allmetadata/scratch16000; reserve7208 inputvector',
  preparation_owned_upper='one streamed QQrecord/scalarden+keyset7208+allmetadata16000 <25000; no fullJSON/second coefficientpayload',
  native_inner_products_upper3769784=3769784,
  ALL_conversion_scaled_root_products_upper4117089=4117089,
  firstnonzero_verification='forward,reverse subtractionZERO,forwardreproduce inONE accumulator',
  source_input_L1_transaction_before_commit=True,
  prefix='T0CHP001 fullinput/root/accumulator,D,L1,all phase/k/i/r/root/counter cursors; fixedchoose/100powers literalmetadata reconstruction',
  before_commit_guards=True,phase5='allfirstlayerpassescomplete;clockguardbeforeoutput',
  no_IO_doubleprefix=True,soft45_hard60_RAM2GB=True,
  C_powers_layers_valuation_NOT_run=True,no_otherfactor_GCD_tail_or_badlocus_inference=True)
 assert proof['whole_owned34376']==34376
 dump('literal-frame-height-ownership-and-cursors.json',proof)
 raw=sum(p.stat().st_size for p in out.iterdir())+sum(p.stat().st_size for p in [cpp,drv,me])+20000
 total=raw+4400000+200000;assert total<8000000
 dump('RAW-closure.json',dict(actual_prepared_and_codes_plus_metadata=raw,full_state_OR_output_upper4400000=4400000,
  allfuture_receipts_logs_reports200000=200000,total_upper=total,hard_limit8000000=8000000,
  one_state_no_atomic_doublebuffer=True,no_success_before_allfirstlayerpasses=True,no_IO_doublefullprefix=True))
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 r=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,value_worker_NOT_run=True,
  ONLY_T0_root239=True,all7208_QQvalues_preserved_fixedD=True,old_offset4_removed_once_in_literalbody=True,
  future_shift0_no_double_shift=True,integer_hard518=518,HH_bits481=481,hard_owned34376=34376,
  ALL_products_upper4117089=4117089,RAW_total_upper=total,
  compiled_help_only=True,fresh_ROOT_source_code_audit_REQUIRED=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',r);signal.alarm(0);print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
