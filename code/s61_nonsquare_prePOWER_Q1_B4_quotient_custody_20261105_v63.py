#!/usr/bin/env python3
"""Prepare FOUR fixedB lex50ad divisions; compile/help ONLY, no quotient value."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_WORKER_CUSTODY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def load(st,name='report.json'):return json.loads(pin(ROOT/'reports'/st/name).read_text())
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 peer=load('S61_NONSQUARE_PREPOWER_Q1_B_ROOT_JET_HORNER_REVIEW_20261105_v63');assert peer['all_passed']
 load('S61_NONSQUARE_PREPOWER_Q1_B_ROOT_JET_HORNER_REVIEW_PUBLICATION_20261105_v63')
 source=load('S61_NONSQUARE_PREPOWER_Q1_B_NEXT_SOURCE_20261105_v63');assert source['all_passed']
 pre=load('S61_NONSQUARE_PREPOWER_Q1_B_NEXT_SOURCE_20261105_v63','fixed-primitive-B4-lex-quotient-SOURCE-ONLY.json')
 assert pre['whole_owned52826']==52826 and pre['final_full_B4multiplyback_partial_bits386']==386
 assert pre['all_product_upper1249351']==1249351
 old=load('S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63');assert old['all_passed']
 inp=pin(ROOT/'reports/S61_NONSQUARE_PREPOWER_T1_C4_LONGDIVISION_REVIEW_20261105_v63/quotient-original-common-D.txt');assert sha(inp)==old['evidence_sha256'][inp.name]
 count=0;L=0;prev=None;body=hashlib.sha256()
 with inp.open()as f,(out/'native-input.txt').open('w')as g:
  assert f.readline().strip()=='T1C4QN001';Dline=f.readline();D=int(Dline);assert D>0 and D.bit_length()==197
  assert int(f.readline())==4811;g.write('Q1B4Q001 4811 4 0\n');g.write(Dline)
  for line in f:
   key,c=map(int,line.split());assert c and(prev is None or prev<key);prev=key
   a,x=divmod(key,1048576);e,d=divmod(x,1024);assert 0<=a<=16 and 0<=e<=33 and 0<=d<=18 and a+e+d<=50
   g.write(line);body.update(line.encode());L+=abs(c);count+=1
 assert count==4811 and L.bit_length()<=301
 cpp=pin(ROOT/'code/s61_nonsquare_Q1_B4_quotient_20261105_v63.cpp');drv=pin(ROOT/'code/s61_nonsquare_prePOWER_Q1_B4_quotient_entry_20261105_v63.py');pin(Path(__file__))
 note=pin(ROOT/f'research/notes/{ID}.md');history=[]
 for folder in sorted((ROOT/'reports').glob(ID+'_ATTEMPT*')):
  for p in sorted(folder.rglob('*')):
   if p.is_file():pin(p);history.append(str(p.relative_to(ROOT)))
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu -lflint','Q1B4-build',str(cpp),str(out/'native-worker')]
 with(out/'compile.log').open('w')as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=30)
 for name,cmd in [('worker-help',[str(out/'native-worker'),'--help']),('driver-help',[str(ROOT/'.venv/bin/python'),str(drv),'--help'])]:
  r=subprocess.run(cmd,text=True,capture_output=True,check=True,timeout=5);(out/(name+'.log')).write_text(r.stdout+r.stderr)
 proof=dict(prepared_ONLY_no_B4quotient=True,independent_vB4_gate_fulfilled=True,
  input_canonical4811_BODY_byteidentical=True,body_sha256=body.hexdigest(),original_D197=D,input_L1_bits=L.bit_length(),
  Q2='Q1/B4=T1/(a3*C4*B4), integernumerator/originalD197; nottheoriginalS61quadricq2',
  B='50ad+18e2d-9ed-2d+60e-40',primitive_ZZ='-1*(-9)+4*(-2)=1; Gauss + exactvB4 impliesintegerquotient',
  fixed_order='lex a>e>d, packedkey a*1048576+e*1024+d injectiveoverboxes; LM50ad',
  FOURpasses='leadinga,d>=1; integerLC50 division +50backcheck; addnewquotientonce, subtractALL6terms; FULLremainderZERO eachpass; moveQintoR withoutcopy fornextpass',
  input_lifetime='original4811 are streamed directlyintoR, NEVER retained as a vector alongsideR/Q. Atfinalverification immutableinput is stream-reread intoR ONLY afterfourpasses andfactorcachebuild.',
  stage_support_bounds=pre['stage_support_bounds'],stage_boxes=pre['true_stage_degree_boxes_a_e_d_total'],
  exact_L1_bounds=pre['exact_L1_bounds'],stage_R363_fullback386_hard518=True,
  main_owned52826=52826,finalverify_owned45706=45706,
  owned='2*(Rold<=10013+Qnext<=8400)+16000; factorbuild holdsQfinal4615+twofixedpowers<=450, thenonlyFcache225+Q4615+Rverify10013',
  numeric_metadata16000='D/L/Lr/activequotient,quotientkey/monomial/source/read/cofactorcursors,pass/6termposition,counters,stagecounts/FNVtelemetry andstrings/IOscratch; allcoeff+packedkey two cells',
  active_state='0 none;1 checkedquotientscalarbeforeQinsert;2 Qinserted/Bterm nextpending. Prefixkeepsactivecoefficient aswellasR/Q.',
  final_B4_source='buildfixedB^4 from6literalterms,onecurrent/nextfactoronly; releaseoldfactorbeforeloadingoriginalinput. ALLcoeff B4*finalQ=originalI checked, no negative exponent orremainder discard.',
  source_ALL_coefficient_product_upper1249351=1249351,proposed_productguard1300000=1300000,
  hash_telemetry='FNV1a64canonicalkey/coefficientline records in descendingquotientproductionorder, unsignedmod2^64; informationalONLY, notcryptographicorproof. No intermediate polynomialfiles. sourceSHA256 checkedbefore/afterbydriver.',
  metadata_hash_char_upper4368490=170*sum(pre['stage_support_bounds'][1:]),
  prefix_version='Q1B4QP001',all_prefix_arrays='R/Q/fixedcurrentF/nextF; D/L/Lr/activequotient; stagecounts/digests,all16cursors/counters/sourcefingerprints',
  cursor='phase,pass,Bterm,factorh,sourcei,active_ready,activeQkey,factorkey,verifyQkey,verifyFkey,products,scalars,updates,peak,inputFNV,rereadFNV',
  reread='sameimmutable4811/nativeinput,D/header/order/L1/FNVchecked; driver allSHA256pinsverifiedagain afterjob. Not anoracle orold572 replay.',
  lifecycle_raw='onefullprefix OR finalsuccessoutput; noatomicdoublebuffer/noIOsecondfullprefix; outputonlyafterALLB4multiplybackzero',
  prefix_bound5950210=5950210,reserve6000000=6000000,output_max4615=True,
  history=history,soft45_hard60_RAM2GB=True,no_b1inverse_polyGCD_otherfactor_tail_or_badlocus=True)
 assert proof['metadata_hash_char_upper4368490']==4368490
 dump('fixedB-fourpass-lifetime-height-cursors-and-guards.json',proof)
 files={str(p.relative_to(ROOT)):p for p in out.iterdir()};files.update({n:ROOT/n for n in pins})
 raw=sum(p.stat().st_size for p in files.values())+20000;assert raw<=1800000
 total=raw+6000000+200000;assert total<=8000000
 dump('RAW-closure.json',dict(all_prepared_inputs_history_codes_andmetadata=raw,fullprefix_OR_output_reserve6000000=6000000,
  complete_prefix_source_bound5950210=5950210,future_receipts200000=200000,total_upper=total,hard8000000=8000000,noIOdoubleprefix=True))
 r=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,value_worker_NOT_run=True,
  independent_Bvaluation4_gate_fulfilled=True,originalD197_and_BODY4811_unchanged=True,no_original_vector_retained=True,
  whole_owned52826=52826,source363_386_hard518=True,ALL_coefficient_products_upper1249351=1249351,
  RAW_total_upper=total,compiled_help_ONLY=True,fresh_ROOT_source_code_audit_and_ONE_GO_REQUIRED=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',r);signal.alarm(0);print(json.dumps({k:v for k,v in r.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
