#!/usr/bin/env python3
"""Prepare fixed40 CRT/exactfactor worker; compile/help/synthetic IO fixture ONLY, no target arithmetic."""
from pathlib import Path
import argparse,hashlib,json,resource,signal,subprocess,time,re
ROOT=Path(__file__).resolve().parents[1]
ID='S61_TWO_TAIL_CRT_LIFT_WORKER_CUSTODY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(p):pins[str(p.relative_to(ROOT))]=sha(p);return p
 def dump(n,x):(out/n).write_text(json.dumps(x,indent=2)+'\n')
 srcpath=pin(ROOT/'reports/S61_TWO_TAIL_CRT_LIFT_SOURCE_20261105_v63/report.json');src=json.loads(srcpath.read_text());assert src['all_passed']and src['source_only']
 for n,h in src['input_sha256'].items():assert sha(pin(ROOT/n))==h,n
 for n,h in src['evidence_sha256'].items():assert sha(pin(srcpath.parent/n))==h,n
 sched=json.loads((srcpath.parent/'fixed-primes-degreefilter-and-CRT-normalization-SOURCE.json').read_text())
 budget=json.loads((srcpath.parent/'whole-ownership-binary-prefix-and-RAW-budget.json').read_text())
 assert budget['integer_owned54090']==54090 and budget['RAW_union_upper26MiB']==26*2**20
 math_paths=[x['path']for x in sched['inputs']]+['reports/S61_TWO_TAIL_FIXED_KRONECKER_ENTRY_20261105_v63/monic-gcd-mod32003.txt']
 for n in math_paths:pin(ROOT/n)
 cpp=pin(ROOT/'code/s61_two_tail_CRT_lift_20261105_v63.cpp');hpp=pin(ROOT/'code/s61_two_tail_CRT_state_20261105_v63.hpp')
 drv=pin(ROOT/'code/s61_two_tail_CRT_lift_entry_20261105_v63.py');pin(Path(__file__))
 note=pin(ROOT/f'research/notes/{ID}.md')
 actual_primes=list(map(int,re.search(r'P\[40\]=\{([^}]+)\}',hpp.read_text()).group(1).split(',')))
 assert actual_primes==[r['p']for r in sched['fixed40_primes']]
 history=[]
 for folder in sorted((ROOT/'reports').glob(ID+'_ATTEMPT*')):
  for p in sorted(folder.rglob('*')):
   if p.is_file():pin(p);history.append(str(p.relative_to(ROOT)))
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -Wl,-rpath,/workspace/cloud-setup/prefix/usr/lib/x86_64-linux-gnu -lflint -lgmp','two-tail-CRT-build',str(cpp),str(out/'native-worker')]
 with(out/'compile.log').open('w')as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=30)
 for name,cmd in [('worker-help',[str(out/'native-worker'),'--help']),('driver-help',[str(ROOT/'.venv/bin/python'),str(drv),'--help'])]:
  r=subprocess.run(cmd,text=True,capture_output=True,check=True,timeout=5);(out/(name+'.log')).write_text(r.stdout+r.stderr)
 # A schema fixture has no source paths, no target input opened, and invokes no gcd/CRT.
 with(out/'fixture.log').open('w')as f:subprocess.run([str(out/'native-worker'),'--prefix-fixture',str(out/'prefix-schema-fixture.bin')],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=5)
 r=subprocess.run([str(ROOT/'.venv/bin/python'),str(drv),'--inspect-prefix',str(out/'prefix-schema-fixture.bin')],text=True,capture_output=True,check=True,timeout=5)
 fixture=json.loads(r.stdout);assert fixture['mode']==99 and fixture['max_map_integer_bits']==8192 and fixture['max_scalar_bits']==8193
 assert fixture['source_paths']==['','','']and fixture['current_owned']<60000 and fixture['binary_roundtrip_SHA256']==sha(out/'prefix-schema-fixture.bin')
 dump('prefix-fixture-inspection.json',fixture)
 proof=dict(preparation_ONLY_NO_TARGET_ARITHMETIC=True,
  math_inputs=math_paths,input_arrays='4500+2317 immutable sparse key/integer records, reserved capacities exactlycharged; <=301bitcoefficients andD197 unchanged',
  fixed_prime_list=actual_primes,seed32003_reused=True,newprime_Euclid_calls_atmost39=True,
  source_gamma226_bit_scalar='computedfromtwoinputLCs onlyduringfutureworker; allfixedprimes/LCs/gamma/sourcepreflightpinned',
  prime_control='no reset/noextra prime; higherdegree3899skip/lowerSTOP; fieldf/g buffers freedBEFOREeachprimeallocation, so capacitytotal18521 notold+new overlap',
  CRT_control='3900denseentriesINPLACE; currenti oldM/newM andCRSUB0..4 stored; MchangesONLYafterallentriescommit; anyuncommitted scalar isState z[PENDING]/cursors',
  content='positiveintegergcd ofALLbalanced3900, everydivisionexactchecked/backmultiplied, leadinggamma exactwhenM>2gamma; nocandidateintegerdivisionbeforecontent',
  selected_candidate='shape<=207 witha6/e9/d4,total11,uniqueW13at(6,1,4); supportmayrestorezero-mod32003coeffs; otherwiseSTOPno alternativecandidate',
  height='exactsourceformula computedfromactualcandidateLC/L1tail. BOTH BR<=8192 requiredbeforeintegerwork. WeightedtriangularL1 controlsfalsecandidatepartialR/Qandfreshmultiplybacks.',
  partial_work_boxes='N0R a16/d18/W42 andQ a10/d14/W29;N2R a12/d14/W35andQ a6/d10/W22. No assumptionsuccessfule-degreedrop.',
  new_coefficient_cap8192=True,uncommitted_add_scratch8193='mapcoeffs NEVER exceed8192; an uncommittedadd mayhave8193bits andisStatePENDING, binaryscalarcap8193. 1030serializedbytes stillwithin1040recordreserve. Old518jobsuntouched.',
  key='W*2^26+legacy(a*2^20+e*2^10+d); integerpartial e42 fitslegacy, NOT usingKron646/19/1forintegerproducts',
  all_exact_ids='integerLCdivisibility/constantreversecheck,thenALLmultivariableH*Q0=N0andH*Q2=N2 againstoriginalinputarrays; notonlyunivariatedivisions',
  field_owned52055=52055,integer_owned54090=54090,hard_owned60000=60000,
  lifetime='no copied localsparseR/Q/candidate outsideState. H/R/Q allState maps. Fieldvectors+denseCRT freedbeforeintegerR/Q; verifiedfirstcofactor savedthenQclearedbeforesecond.',
  modular_product_cap6000000000=6000000000,bigint_product_cap2100000=2100000,
  prefix='TTCRT001 LE typedbinary; all64uintcursors/14bigscalars/40x8ledger/3sourcepaths/inputswithcapacities/f+gwithcapacities/CRTwithcapacity/H/R/Q. Pendingdeltaandnewsum bothState scalars.',
  schema_length32='binarybigint useslen32 (paddingstill<=1040/48/96/80reservedrecords); signs/magnitude canonical; NO libinternallayout',
  prefix_fixture='syntheticinput-3,600bitCRT,8192bitmaprecord,8193bituncommittedscalar; originalsourcefilesnotopened; no finitefieldgcd/CRT/targetcoefficient work',
  fixture_stream_roundtrip_verified=True,no_double_fullprefix_roundtrip_file=True,
  STOP='onefullbinaryprefixattempt; sourcewriterhasNOclockrethrowwhensavingtimeSTOP. PrefixIO failure nevertriggerssecondprefix; postdriverhandlesmissingreportasfailure,neverfalseproof.',
  soft45_hard58driver60CPU_RAM2GB=True,history=history,
  fresh_ROOT_whole_code_source_audit_and_ONE_GO_REQUIRED=True)
 dump('source-code-lifetimes-cursors-and-guards.json',proof)
 maximal=dict(CONDITIONAL_SOURCE_ONLY_no_maximality_claim_yet=True,
  conditions='bothchar0identitiesexact + H_Krondegree3899 + Hnormalizedmod32003equalverifiedseed andLCunit (driverchecksreadback)',
  saved_seed_quotient_peer='S61_TWO_TAIL_FIXED_KRONECKER_QUOTIENT_REVIEW_V2_20261105_v63: quotientdegrees6645/4076coprime',
  deduction='cofactorKroneckerimagesdifferonlybyunitfromsavedcoprimequotients;theirLCs survive32003. Anychar0commoncofactorfactorwouldmapnonconstantanddividemodgcd1. ThusactualcofactorscoprimeandHmaximalgcd uptoQQscalar,WITHOUT claimingHirreducible.',
  note='finaltruecofactorordinarye/ddegreeboundsfollowfromexactpolynomialproduct; doNOTusebroadpartialworke42asKroneckerinjectionbound',
  adoption='onlyafteractualidentitiesANDindependentvaluepeer; producerreportdoesnotpromotethisdeduction')
 dump('conditional-maximality-lemma.json',maximal)
 files={str(p.relative_to(ROOT)):p for p in out.iterdir()};files.update({n:ROOT/n for n in pins})
 base=sum(p.stat().st_size for p in files.values())+20000;assert base<=3*2**20
 total=base+23*2**20;assert total<32*2**20
 dump('RAW-closure.json',dict(allprepared_history_mathinputs_codes_andmetadata=base,
  oneprefix16MiB=16*2**20,cofactoroutputs5MiB=5*2**20,CRTcandidate1MiB=2**20,reportslogs1MiB=2**20,
  RAW_actualbase_plus_reserves=total,SOURCE_upper26MiB=26*2**20,hard32MiB=32*2**20,
  noatomic_or_secondfullprefix=True,fixture_and_allpreparationhistory_in_BASE=True))
 result=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,arithmetic_worker_NOT_run=True,
  compiled_help_and_synthetic_schema_fixture_ONLY=True,fixed40_primes=True,math_input_paths=math_paths,
  hard_owned60000=60000,conditional_integer8192=True,RAW_total_upper=total,
  fresh_ROOT_code_source_audit_and_ONE_GO_REQUIRED=True,input_sha256=pins,
  evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump('report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
