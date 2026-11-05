#!/usr/bin/env python3
"""Prepare fixed T0/C4 canonical-D integer operator; compile/help ONLY, no quotient values."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_T0_C4_QUOTIENT_WORKER_CUSTODY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):
  return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 source=load('S61_NONSQUARE_PREPOWER_T0_C4_NEXT_SOURCE_20261105_v63');assert source['all_passed']and source['source_only']
 pre=load('S61_NONSQUARE_PREPOWER_T0_C4_NEXT_SOURCE_20261105_v63','fixed-C4-integer-quotient-SOURCE.json')
 assert pre['verification_owned79558']==79558 and pre['ALL_scaled_root_polynomial_products_upper2581132']==2581132
 peer=load('S61_NONSQUARE_PREPOWER_T0_C_ROOT_JET_HORNER_REVIEW_20261105_v63');assert peer['all_passed']
 publication=load('S61_NONSQUARE_PREPOWER_T0_C_ROOT_JET_HORNER_REVIEW_PUBLICATION_20261105_v63');assert publication['all_passed']
 primary=load('S61_NONSQUARE_PREPOWER_T0_C_HASSE_ENTRY_20261105_v63');assert primary['result']['C_valuation']==4
 old=load('S61_NONSQUARE_PREPOWER_T0_C_HASSE_WORKER_CUSTODY_20261105_v63');assert old['all_passed']and old['value_worker_NOT_run']
 inp=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T0_C_HASSE_WORKER_CUSTODY_20261105_v63/native-input.txt')
 assert sha(inp)==old['evidence_sha256'][inp.name]
 count=0;L=0;keys=set();body=hashlib.sha256()
 with inp.open()as f,(out/'native-input.txt').open('w')as g:
  assert f.readline().strip()=='T0CH001 7208 20 0';Dline=f.readline();D=int(Dline)
  assert D==pre['common_D0_D1_integer']and D>0 and D.bit_length()==197
  g.write('T0C4Q001 7208 20 0 16\n');g.write(Dline)
  for line in f:
   key,c=map(int,line.split());assert c and key not in keys;keys.add(key)
   a,x=divmod(key,1048576);e,d=divmod(x,1024)
   assert 0<=a<=20 and 0<=e<=40 and 0<=d<=18 and a+e+d<=57
   assert abs(c).bit_length()<=320
   g.write(line);body.update(line.encode());count+=1;L+=abs(c)
 assert count==7208 and L.bit_length()<=326
 cpp=pin(ROOT/'code/s61_nonsquare_T0_C4_quotient_20261105_v63.cpp')
 drv=pin(ROOT/'code/s61_nonsquare_prePOWER_T0_C4_quotient_entry_20261105_v63.py');me=pin(Path(__file__))
 note=pin(ROOT/f'research/notes/{ID}.md')
 # No --run call is present in this preparation program.
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu -lflint','T0C4-build',str(cpp),str(out/'native-worker')]
 with(out/'compile.log').open('w')as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=30)
 for name,cmd in [('worker-help',[str(out/'native-worker'),'--help']),('driver-help',[str(ROOT/'.venv/bin/python'),str(drv),'--help'])]:
  r=subprocess.run(cmd,text=True,capture_output=True,check=True,timeout=5);(out/(name+'.log')).write_text(r.stdout+r.stderr)
 proof=dict(ONLY_T0=True,prepared_only=True,independent_Cvaluation4_peer_gate_fulfilled=True,
  source_body_SHA256=body.hexdigest(),input_header_only_changed_body7208_lines_identical=True,
  input_body='U0=T0/a4 already; original file offset4 is provenance, new native shift0; NO second a4 removal',
  C='100*a+18*e^2-9*e-2=100*a-R(e)',R='2+9*e-18*e^2',
  input_D_integer=D,input_D_bits197=197,input_integer_L1_bits=L.bit_length(),
  quotient='J0=I0/C4 in ZZ[a,e,d]; Q0=J0/D0=T0/(a4*C4)',
  formula='S_i=sum_{j>=i+4}I_j*binom(j-i-1,3)*R^(j-i-4)*100^(20-j+i), i0..16; S_i=10020*J0_i',
  primitive_Gauss='-1*(-9)+4*(-2)=1; Cprimitive; accepted Cvaluation4 implies rational division, Gauss implies J0 integral',
  whole_slab_before_any_final_support_filter=True,partial_slab_e_max72=72,partial_slab_d_max18=18,partial_slab_terms_upper1387=1387,
  slab_transfer='only after fullsum: assert e<=32,totaldegree<=49; exact fixed10020 divide EVERY coefficient; reverse multiply; store canonical integer',
  no_independent_contents_or_GCD_search=True,
  final_support_upper9690=9690,final_a16_e32_d18_total49=True,
  negative_a_remainder='full C4*J0=I0 coefficient identity, using originalI0 unscaled; NO assumed negative Laurent truncation',
  remainder_a20_e40_d18_total57_support14592=True,
  formula_owned53148=53148,verification_owned79558=79558,
  coefficient_and_key_cells_per_entry2=2,ALL_numeric_metadata_reserve16000=16000,
  metadata_includes=['inputD/L','100powers21','binomial17','temporaryZs','all size/counter/cursor/exponent records','publication buffers and formatting scratch'],
  lifetime='input7208+Rcache289+J0<=9690+ONEslab1387; release ALLRcache and emptyslab before Cpower55/reserve289+remainder14592; no product polynomial',
  temporary_L1_bits=dict(input326=326,scaledslab435=435,canonicalJ303=303,scalar_reverse436=436,whole_multiplyback333=333,hard518=518),
  ALL_scalar_root_polynomial_product_upper2581132=2581132,exact_scalar_division_checks_upper9690=9690,
  scalar_divisor_fixed10020_not_searched=True,success_denominator='originalD197, no10020 factor',
  prefix_version='T0C4QP001',prefix_allstate=['inputD/L and literal7208','Rcache incl zero slots','Cpower maps','Q','whole slab','remainder','phase/sourcei/rootr/transferkey/Ccursor/Qcursor','product/scalar/update/owned counters'],
  pending_cursor='beforecommit at source/root multiplication; completed slab before scalar exactdivide/transferkey; constantC h/key/t; remainder inputi; Qkey/Ckey; phase7 fullmultiplybackdone',
  fixed_metadata_reconstruction='literal100powers0..20 andbinomial(h+3,3)0..16 only; no unknown source replay',
  no_IO_doubleprefix='success output only after wholeidentity; IO failure no secondfull prefix; mathSTOP writes ONEstate, noatomicduplicate',
  soft45_hard60_RAM2GB=True,no_actual_quotient_GCD_T1_otherfactor_tail_or_badlocus=True)
 dump('literal-frame-height-lifetime-and-guards.json',proof)
 # This is RAW union rather than a claim that all source files are live in memory.
 base_files={str(p.relative_to(ROOT)):p for p in out.iterdir()}
 base_files.update({n:ROOT/n for n in pins})
 actual=sum(p.stat().st_size for p in base_files.values())+20000
 total=actual+8300000+200000;assert total<=11000000
 dump('RAW-closure.json',dict(actual_prepared_codes_external_inputs_plus_metadata=actual,full_state_or_output_upper8300000=8300000,
  state_bound='(31779polycoeff+16000metadata)*170bytes+100000=8222430<8300000; each coeff and numeric key separately owned',
  receipts_upper200000=200000,total_upper=total,hard11000000=11000000,
  no_atomic_duplicate_state=True,success_only_after_fullmultiplyback=True,IO_failure_no_secondfullprefix=True))
 r=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,value_worker_NOT_run=True,
  independent_T0_Cvaluation4_gate_fulfilled=True,canonical_original_D197=True,
  body7208_unchanged_shift0=True,whole_owned79558=79558,hard_all_integer518=518,
  ALL_products_upper2581132=2581132,RAW_total_upper=total,compiled_help_ONLY=True,
  fresh_ROOT_source_code_audit_and_ONE_GO_REQUIRED=True,input_sha256=pins,
  evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',r);signal.alarm(0)
 print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
