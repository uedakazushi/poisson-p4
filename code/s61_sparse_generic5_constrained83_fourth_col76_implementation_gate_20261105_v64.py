#!/usr/bin/env python3
"""Freeze allcustody/parse/partials/implementation/storage and safe --help."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd,prod
from itertools import product
import argparse,json,hashlib,ast,subprocess,sys,time,signal,resource
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_COL76_IMPLEMENTATION_GATE_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0

def main(a):
 assert a.freeze
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);worker=ROOT/'code/s61_sparse_generic5_constrained83_fourth_col76_coefficients_20261105_v64.py';ast.parse(worker.read_text());ast.parse(Path(__file__).read_text());target=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_COL76_COEFFICIENTS_20261105_v64';assert not target.exists()
 for name,code in [('worker-help.txt',worker),('helper-help.txt',Path(__file__))]:
  with(out/name).open('wb')as f:subprocess.run([sys.executable,str(code),'--help'],stdout=f,stderr=subprocess.STDOUT,timeout=10,check=True)
 assert not target.exists()
 sp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_UNIT_COL76_SOURCE_20261105_v64/report.json';assert sha(sp)=='349d1606feb479950eb690241bfdb95c49586731465b9f456e6656d20ea467cc';s=read(sp);assert s['all_passed'];bp=sp.parent/'col76-fixed-source-operands.json';pp=sp.parent/'all20-cleared-col76-support-height-plans.json'
 for p in(bp,pp):assert sha(p)==s['evidence_sha256'][p.name]
 bundle=read(bp)['source_operands'];plans=read(pp)['all20_source_plans'];assert len(plans)==20
 auditp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_SOURCE_OWNERSHIP_AUDIT_20261105_v64/report.json';assert sha(auditp)=='ff6eb57226988b1a3a4c8dda9c5e5a8f1621ff8bd1ef372ae40b37ed508c8614';audit=read(auditp);assert audit['all_passed']and audit['whole_source_numeric_upper']==87029
 inputs={ROOT/p:h for p,h in s['input_sha256'].items()};inputs.update({sp:sha(sp),bp:sha(bp),pp:sha(pp),auditp:sha(auditp),worker:sha(worker),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_COL76_COEFFICIENTS_20261105_v64.md':sha(ROOT/'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_COL76_COEFFICIENTS_20261105_v64.md')})
 for n,h in s['evidence_sha256'].items():inputs[sp.parent/n]=h
 for n,h in audit['evidence_sha256'].items():inputs[auditp.parent/n]=h
 for p,h in inputs.items():assert sha(p)==h
 # Bind all63 leaves to their ORIGINALcached sources, not just suppliednumbers.
 from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
 planp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64/actual-fixed24-entry-source-and-exponent-plan.json';v=read(planp);E={(z['row'],z['column']):z for z in v['actual_entry_sources576']};del v
 cachep=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR3_COL48_COEFFICIENTS_20261105_v64/actual-residual-col48-and-structural-frontier.json';cache={z['row']:z for z in read(cachep)['actual21_residual_col48_entries']};profiles={};custody=[];scalarprep=0
 for label,z in bundle.items():
  p=z['original_literal'];ref=z['source']
  if label=='Q':assert p=={'content':'1','shift':[0]*5,'terms':[['5',0,0,0,1,0],['1',0,0,0,0,0]]}
  elif label.startswith('U'):
   row=int(label[1:]);c=cache[row];path=cachep.parent/c['file'];assert ref=={'path':str(path.relative_to(ROOT)),'sha256':c['sha256']}and sha(path)==c['sha256'];blob=path.read_bytes();actual=decode(blob[9:])if blob[8]else{'content':'0','shift':[0]*5,'terms':[]};assert p==actual;assert Q(z['QQ_transform_scalar'])==Q(6075,256)and z['transform_monomial_exponent5']==[-2,0,0,0,1]
  else:
   row,col=(31,50)if label=='P'else((83,76)if label=='R'else(int(label[1:]),76 if label.startswith('W')else 50));e=E[row,col]
   if e['terms']:
    assert ref==e['source'];path=ROOT/ref['file']
    with path.open('rb')as f:f.seek(ref['offset']);blob=f.read(ref['bytes'])
    assert hashlib.sha256(blob).hexdigest()==ref['payload_sha256']and p==decode(blob)
   else:assert ref=={'kind':'sourcezero','row':row,'column':col}and p=={'content':'0','shift':[0]*5,'terms':[]}
  if not label.startswith('U'):assert Q(z['QQ_transform_scalar'])==1 and z['transform_monomial_exponent5']==[0]*5
  rawq=Q(p['content']);transform=Q(z['QQ_transform_scalar']);scalarprep=max(scalarprep,abs(rawq.numerator).bit_length()+abs(transform.numerator).bit_length()+4,rawq.denominator.bit_length()+transform.denominator.bit_length()+4);q=rawq*transform;shift=[p['shift'][j]+z['transform_monomial_exponent5'][j]for j in range(5)];terms=[(int(t[0]),tuple(shift[j]+t[j+1]for j in range(5)))for t in p['terms']];assert len({e for c,e in terms})==len(terms);profiles[label]={'q':q,'terms':terms,'max':max([abs(c)for c,e in terms]+[0])};custody.append({'label':label,'actual_original_source_equal':True,'terms':len(terms)})
 alltuples=0;verified=[];partialbits=0
 for plan in plans:
  D=1;pre=[];row=plan['row'];source_union=set();M=0
  for p in plan['products']:
   labels=p['factor_labels'];assert labels in [['Q','W'+str(row)],['U'+str(row),'P','R'],['Q','V'+str(row),'R']];q=Q(1)if len(labels)==2 else(Q(405,8)if labels[0].startswith('U')else Q(-405,8));assert p['prefactor_monomial_exponent5']==([0]*5 if len(labels)==2 else[0,0,0,0,-2]);num_b=q.numerator.bit_length();den_b=q.denominator.bit_length()
   for k in labels:q*=profiles[k]['q'];num_b+=abs(profiles[k]['q'].numerator).bit_length();den_b+=profiles[k]['q'].denominator.bit_length()
   scalarprep=max(scalarprep,num_b+4,den_b+4);assert q==Q(p['source_product_content']);D=D*q.denominator//gcd(D,q.denominator);pre.append((p,labels,q))
  assert D==int(plan['integer_common_D'])
  for p,labels,q in pre:
   ns=[len(profiles[k]['terms'])for k in labels];n=prod(ns);assert n==p['literal_tuple_count'];alltuples+=n;assign=n//max(ns);h=assign*prod(profiles[k]['max']for k in labels);weight=q.numerator*(D//q.denominator);assert h==int(p['coefficient_product_integer_abs_upper'])and str(weight)==p['exact_scalar_integer_weight_after_D'];M+=abs(weight)*h
   support=set()
   for combo in product(*(profiles[k]['terms']for k in labels)):
    e=tuple(p['prefactor_monomial_exponent5'][j]+sum(z[1][j]for z in combo)for j in range(5));support.add(e)
   assert len(support)==p['exact_product_support_cardinality'];source_union|=support
  assert M==int(plan['integer_numerator_abs_upper'])and source_union=={tuple(e)for e in plan['exact_source_support_sorted']};partialbits=max(partialbits,M.bit_length()+2);verified.append({'row':row,'products':len(pre),'source_support_exact':True,'integer_D_M_recomputed':True,'absolute_all_partial_integer_coefficients_bound':'Everypartial coefficient sums a subset of primitive tuples; each fixedoutput has atmostassign tuples. SourceM sums absoluteweights*assign*productmax, so allproducts and accumulatedcoefficients stay<=M.','future_owned_source_limit':plan['future_owned_extra_including25kmeta']})
 assert alltuples==1068
 parsecount=numeric(bundle)+numeric(plans)+numeric(custody)+numeric(verified)+sum(7*len(z['terms'])for z in profiles.values())+2000;assert parsecount<25000
 # Fixed perrow working maps charge allbefore/after/block/back and encode,
 # 7numericfields permonomial plus128scratch;sourcelimits keepfurthermargin.
 for p in plans:assert 25000+28*p['exact_union_support_cardinality']+256<=p['future_owned_extra_including25kmeta']<=28968
 base=sum(p.stat().st_size for p in inputs);caps=sum(p['uncompressed_binary_RAW_upper']for p in plans);future=base+caps+300000+131072;assert future<=16658923 and future<=32*1024*1024
 proof={'all_passed':True,'source_only':True,'ONE20_cleared_column_job_safe':True,'worker_sha256':sha(worker),'operand_bundle':str(bp.relative_to(ROOT)),'source_plans':str(pp.relative_to(ROOT)),'actual_READ_input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'actual_READ_input_RAW_bytes':base,'all63_leaf_custody_checks':custody,'all20_support_D_M_and_partial_checks':verified,'all_source_tuple_count1068':alltuples,'coefficient_fixed_metadata_limit':25000,'actual_parsed_metadata_and_shared_leaf_numeric_bound':parsecount,'coefficient_owned_limit':28968,'source_producer_postGC_tracked38348_NOT_whole_parse':True,'source_full_parse_numeric_upper_reference':87029,'source_full_parse_audit_sha256':sha(auditp),'scalar_content_preparation_unreduced_bit_upper':scalarprep,'family_integer_partial_bit_upper':partialbits,'all_new_values_representation':'LaurentnumeratorF/(5q0+1), no newcolumnoperation; commonQcomplementOPEN. P290onlyonce. Other19columns/RHS/full173retained.','helper_and_worker_help_safe':True,'binary_output_upper':caps,'other_artifact_upper':300000,'prepared_gate_HELP_and_addendum_input_upper':131072,'fresh_RAW_upper':future,'hard_seconds':60,'hard_address_space_bytes':2000000000,'no_fieldGCD_points_prime_rank_or_other19updates':True}
 (out/'implementation-custody-storage-bit-addendum.json').write_bytes(enc(proof));report={'id':ID,'all_passed':True,'source_only':True,'ONE20_cleared_column_job_safe':True,'owned_upper':28968,'RAW_upper':future,'all63_leaf_custody':True,'all20_support_partial_height':True,'safe_worker_helper_help':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));assert sum(p.stat().st_size for p in out.iterdir())<=131072;signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true',help='freeze only source/custody/storage proof, never coefficientvalues');main(p.parse_args())
