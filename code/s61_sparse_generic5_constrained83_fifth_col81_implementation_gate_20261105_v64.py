#!/usr/bin/env python3
"""Freeze80-source custody/reuse, VALUEparse/height/storage for ONE19bundle."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd,prod
from itertools import product
import argparse,json,hashlib,ast,subprocess,sys,time,signal,resource
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIFTH_COL81_IMPLEMENTATION_GATE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def poly(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']}
def plus(A,B,c=Q(1)):
 R=dict(A)
 for e,q in B.items():
  R[e]=R.get(e,Q(0))+c*q
  if not R[e]:del R[e]
 return R
def timesQ(A):
 R=dict(A)
 for e,q in A.items():
  f=list(e);f[3]+=1;f=tuple(f);R[f]=R.get(f,Q(0))+5*q
 return{e:q for e,q in R.items()if q}

def main(a):
 assert a.freeze
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);worker=ROOT/'code/s61_sparse_generic5_constrained83_fifth_col81_coefficients_20261105_v64.py';ast.parse(worker.read_text());ast.parse(Path(__file__).read_text());target=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIFTH_COL81_COEFFICIENTS_20261105_v64';assert not target.exists()
 for name,code in [('worker-help.txt',worker),('helper-help.txt',Path(__file__))]:
  with(out/name).open('wb')as f:subprocess.run([sys.executable,str(code),'--help'],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=10)
 assert not target.exists()
 sp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIFTH_UNIT_COL81_SOURCE_20261105_v64/report.json';assert sha(sp)=='d668d98214a0f755afe96990eb143f355c90ecaf8df7a9c9a988eeb4509b11fe';s=read(sp);assert s['all_passed']and s['whole_source_numeric_upper_including_parse']==96956;bp=sp.parent/'next-col81-fixed-source-operands.json';pp=sp.parent/'all19-cleared-col81-support-height-plans.json';fp=sp.parent/'existing-commonQ-factor-and-remainder-witnesses.json'
 for p in(bp,pp,fp):assert sha(p)==s['evidence_sha256'][p.name]
 bundle=read(bp)['source_operands'];plans=read(pp)['all19_source_plans'];factorproof=read(fp);assert len(bundle)==80 and len(plans)==19
 csp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIFTH_SOURCE_CUSTODY_PLAN_20261105_v64/report.json';assert sha(csp)=='40cf0a104292997bb86016e152fa61dd53af38109c12c8bef14bf6b305783b31';cs=read(csp);assert cs['all_passed']and cs['accepted_exact_old_reuse']==39
 inputs={ROOT/p:h for p,h in s['input_sha256'].items()};inputs.update({ROOT/p:h for p,h in cs['input_sha256'].items()});inputs.update({sp:sha(sp),bp:sha(bp),pp:sha(pp),fp:sha(fp),csp:sha(csp),worker:sha(worker),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_FIFTH_COL81_COEFFICIENTS_20261105_v64.md':sha(ROOT/'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_FIFTH_COL81_COEFFICIENTS_20261105_v64.md')})
 for root,r in [(sp,s),(csp,cs)]:
  for name,h in r['evidence_sha256'].items():inputs[root.parent/name]=h
 assert not any(str(p).endswith('current-matrix.bin')for p in inputs),'old14MBmatrix NOTin this actualREAD set'
 for p,h in inputs.items():assert sha(p)==h
 oldp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_UNIT_COL76_SOURCE_20261105_v64/col76-fixed-source-operands.json';old=read(oldp)['source_operands'];acceptedp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_COL76_IMPLEMENTATION_GATE_20261105_v64/report.json';assert sha(acceptedp)=='1c5f79ea5b373e2c8b7f95aff5ccefe9a058c7d4d7bd20df6bc7429bd6c14b23';accepted=read(acceptedp);acp=acceptedp.parent/'implementation-custody-storage-bit-addendum.json';assert sha(acp)==accepted['evidence_sha256'][acp.name];oldchecks={z['label']:z for z in read(acp)['all63_leaf_custody_checks']}
 # Only small alreadyacceptedcachebinaries are read below. No recursion into
 # acceptedgate.input_sha256, which historicallycontainsoldmatrix data.
 from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
 n2p=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR_COL81_COEFFICIENTS_20261105_v64/actual-residual-col81-and-structural-frontier.json';n2={z['row']:z for z in read(n2p)['actual22_residual_col81_entries']}
 c4p=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_COL76_COEFFICIENTS_20261105_v64/actual-cleared-col76-and-commonQ-representation.json';c4={z['row']:z for z in read(c4p)['actual20_col76_numerator_values']}
 def cached(path,h):
  p=ROOT/path;assert p in inputs and inputs[p]==h and sha(p)==h;b=p.read_bytes();assert b[:8]==b'Q5CELL1\n';return decode(b[9:])if b[8]else{'content':'0','shift':[0]*5,'terms':[]}
 custody=[];profiles={};reused=0;small=0;quot=0;aux=0;scalarprep=0
 for label,z in bundle.items():
  p=z['original_literal'];src=z['source']
  if label=='Q'or label.startswith('U')or label.startswith('V'):
   assert label in oldchecks and oldchecks[label]['actual_original_source_equal']and z==old[label];reused+=1;method='EXACT39objectreuse_accepted1c5f79physicalproof_NOoldmatrixread'
  elif label in('R','T')or label.startswith('W'):
   row=83 if label=='R'else(31 if label=='T'else int(label[1:]));orig=n2[row];path=str((n2p.parent/orig['file']).relative_to(ROOT));assert src=={'kind':'actual_stage2_col81','row':row,'path':path,'sha256':orig['sha256']};assert p==cached(path,orig['sha256']);small+=1;method='literalacceptedstage2binary'
  elif label.startswith('C'):
   row=int(label[1:]);orig=c4[row];path=str((c4p.parent/orig['file']).relative_to(ROOT));assert src['kind']=='known_Q_quotient_existing_col76'and src['row']==row and src['original_numerator_path']==path and src['original_numerator_sha256']==orig['sha256'];assert timesQ(poly(p))==poly(cached(path,orig['sha256']));quot+=1;method='knownlinearQremultiply_existingvalue'
  else:
   assert label=='G'and src['kind']=='fixed_source_row_sum'and src['rows']==[14,35]and src['QQ_coefficients']==['1','6'];ps=[]
   for row,term in zip([14,35],src['source_terms']):
    orig=n2[row];path=str((n2p.parent/orig['file']).relative_to(ROOT));assert term=={'path':path,'sha256':orig['sha256']};ps.append(poly(cached(path,orig['sha256'])))
   assert poly(p)==plus(ps[0],ps[1],Q(6));aux+=1;method='literalactualrow14plus6row35'
  if not(label.startswith('U')):assert Q(z['QQ_transform_scalar'])==1 and z['transform_monomial_exponent5']==[0]*5
  raw=Q(p['content']);tr=Q(z['QQ_transform_scalar']);scalarprep=max(scalarprep,abs(raw.numerator).bit_length()+abs(tr.numerator).bit_length()+4,raw.denominator.bit_length()+tr.denominator.bit_length()+4);q=raw*tr;shift=[p['shift'][j]+z['transform_monomial_exponent5'][j]for j in range(5)];terms=[(int(t[0]),tuple(shift[j]+t[j+1]for j in range(5)))for t in p['terms']];assert len({e for c,e in terms})==len(terms);profiles[label]={'q':q,'terms':terms,'M':max([abs(c)for c,e in terms]+[0])};custody.append({'label':label,'source_method':method,'terms':len(terms),'literal_identity_or_acceptedobject_equal':True})
 assert(reused,small,quot,aux)==(39,21,19,1)
 Gproof=factorproof['G81_Q_not_divisible'];assert poly(Gproof['G_exact'])==poly(bundle['G']['original_literal']);rem=poly(Gproof['known_Q_remainder_NONZERO']);assert rem and all(e[3]==0 for e in rem);assert plus(timesQ(poly(Gproof['known_Q_quotient_part'])),rem)==poly(bundle['G']['original_literal'])
 tuples=0;verified=[];familybits=0
 for plan in plans:
  i=plan['row'];expected=[(['Q','W'+str(i)],Q(1),[0]*5),(['U'+str(i),'T'],Q(-1),[0]*5),(['Q','U'+str(i),'R'],Q(-1,40),[1,0,0,0,-2]),(['Q','V'+str(i),'R'],Q(-405,8),[0,0,0,0,-2]),(['C'+str(i),'G'],Q(-6075,2),[-1,0,0,0,0])];lookup={tuple(l):(q,sh)for l,q,sh in expected};D=1;pre=[]
  for p in plan['products']:
   labels=p['factor_labels'];q,sh=lookup[tuple(labels)];assert sh==p['prefactor_monomial_exponent5'];nb=abs(q.numerator).bit_length();db=q.denominator.bit_length()
   for k in labels:q*=profiles[k]['q'];nb+=abs(profiles[k]['q'].numerator).bit_length();db+=profiles[k]['q'].denominator.bit_length()
   scalarprep=max(scalarprep,nb+4,db+4);assert q==Q(p['source_product_content']);D=D*q.denominator//gcd(D,q.denominator);pre.append((p,labels,q))
  assert D==int(plan['integer_common_D']);union=set();M=0
  for p,labels,q in pre:
   ns=[len(profiles[k]['terms'])for k in labels];n=prod(ns);assign=n//max(ns);h=assign*prod(profiles[k]['M']for k in labels);weight=q.numerator*(D//q.denominator);assert n==p['literal_tuple_count']and h==int(p['product_primitive_integer_abs_upper'])and weight==int(p['exact_scalar_integer_weight_after_D']);tuples+=n;M+=abs(weight)*h;S=set()
   for combo in product(*(profiles[k]['terms']for k in labels)):S.add(tuple(p['prefactor_monomial_exponent5'][j]+sum(z[1][j]for z in combo)for j in range(5)))
   assert len(S)==p['exact_product_support_cardinality'];union|=S
  assert union=={tuple(e)for e in plan['exact_source_support_sorted']}and M==int(plan['integer_numerator_abs_upper']);assert 25000+28*len(union)+256<=plan['future_owned_extra_including25kmeta']<=30408;familybits=max(familybits,M.bit_length()+2);verified.append({'row':i,'products':len(pre),'exact_support_D_M_and_all_partial_bound_verified':True,'source_owned_limit':plan['future_owned_extra_including25kmeta']})
 assert tuples==2928
 parse=numeric(bundle)+numeric(plans)+numeric(custody)+numeric(verified)+numeric(factorproof)+sum(7*len(z['terms'])for z in profiles.values())+2000;assert parse<25000
 base=sum(p.stat().st_size for p in inputs);caps=sum(p['uncompressed_binary_RAW_upper']for p in plans);future=base+caps+300000+131072;assert future<=2901958 and future<=32*1024*1024
 proof={'all_passed':True,'source_only':True,'ONE19_cleared_column_job_safe':True,'worker_sha256':sha(worker),'operand_bundle':str(bp.relative_to(ROOT)),'source_plans':str(pp.relative_to(ROOT)),'actual_READ_input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'actual_READ_input_RAW_bytes':base,'all80_source_custody_checks':custody,'accepted39_original_physical_proof_reuse':{'report':str(acceptedp.relative_to(ROOT)),'sha256':sha(acceptedp),'exact_object_count':39,'oldmatrix_not_rehashed_or_read':True},'all19_support_D_M_partial_checks':verified,'all_source_tuples2928':tuples,'G81_knownQ_nondivisibility_preserved':True,'VALUEparse_and_shared_source_metadata_numeric_bound':parse,'coefficient_fixed_metadata_limit':25000,'coefficient_owned_limit':30408,'prior_SOURCE_whole_parse_numeric_upper96956_separate':True,'scalar_content_unreduced_temporary_bit_upper':scalarprep,'all_family_integer_partial_bit_upper':familybits,'source_partial_height_proof':'Each outputpartial sumusesasubset ofthefixedtuples. SourceM=sumabsoluteweight*maxassignments*productprimitive_max, so everyterm,productblock,partialaccumulator andbackadd stays<=M.','commonQ_representation_NOT_columnoperation':True,'all_other18_columns_full172_87plus66_RHS_and_oneP290_kept':True,'all_actual_READ_paths_explicit_NOold14MBmatrix':True,'worker_and_helper_help_safe':True,'binary_output_upper':caps,'other_artifact_upper':300000,'prepared_GATE_HELP_input_upper':131072,'fresh_RAW_upper':future,'hard_seconds':60,'hard_address_space_bytes':2000000000,'no_fieldGCD_rank_points_prime':True}
 (out/'implementation-custody-storage-bit-addendum.json').write_bytes(enc(proof));report={'id':ID,'all_passed':True,'source_only':True,'ONE19_cleared_column_job_safe':True,'owned_upper':30408,'RAW_upper':future,'all80_custody_and_G_nondivisibility':True,'VALUEparse_numeric_bound':parse,'no_old14MBmatrix_read':True,'worker_helper_help_safe':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));assert sum(p.stat().st_size for p in out.iterdir())<=131072;signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())
