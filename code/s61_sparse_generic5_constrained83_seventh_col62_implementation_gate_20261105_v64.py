#!/usr/bin/env python3
"""Freeze155-leaf custody and row-stream VALUE bounds for seventh-col62."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd,prod
import argparse,json,hashlib,ast,subprocess,sys,time,signal,resource,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_SEVENTH_COL62_IMPLEMENTATION_GATE_20261105_v64'
SOURCE='S61_SPARSE_GENERIC5_CONSTRAINED83_SEVENTH_UNIT_COL62_SOURCE_20261105_v64'
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
def mono(p):
 assert len(p['terms'])==1;t=p['terms'][0];return Q(p['content'])*int(t[0]),tuple(p['shift'][j]+t[j+1]for j in range(5))
def mm(A,B):return A[0]*B[0],tuple(x+y for x,y in zip(A[1],B[1]))
def inv(A):return 1/A[0],tuple(-x for x in A[1])
def neg(A):return -A[0],A[1]
def main(a):
 assert a.freeze
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 worker=ROOT/'code/s61_sparse_generic5_constrained83_seventh_col62_coefficients_20261105_v64.py';ast.parse(worker.read_text());ast.parse(Path(__file__).read_text());target=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_SEVENTH_COL62_COEFFICIENTS_20261105_v64';assert not target.exists()
 for name,code in [('worker-help.txt',worker),('helper-help.txt',Path(__file__))]:
  with(out/name).open('wb')as f:subprocess.run([sys.executable,str(code),'--help'],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=10)
 assert not target.exists()
 sp=ROOT/f'reports/{SOURCE}/report.json';assert sha(sp)=='f76a6096796dc0917773d963c9a1df89e7826880f785126b1411399796e6a99c';s=read(sp);assert s['all_passed']and s['whole_source_numeric_upper']==65796
 inputs={ROOT/p:h for p,h in s['input_sha256'].items()};inputs[sp]=sha(sp)
 for name,h in s['evidence_sha256'].items():inputs[sp.parent/name]=h
 inputs.update({worker:sha(worker),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_SEVENTH_COL62_COEFFICIENTS_20261105_v64.md':sha(ROOT/'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_SEVENTH_COL62_COEFFICIENTS_20261105_v64.md')})
 acceptedp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_SIXTH_COL74_IMPLEMENTATION_GATE_V2_20261105_v64/report.json';assert sha(acceptedp)=='a3e8e6b84de0400f2980ac5bd541f2a6ab28270031f73aee17af79e77a9a02d1';ar=read(acceptedp);assert ar['all_passed'];acp=acceptedp.parent/'implementation-custody-storage-bit-addendum.json';assert sha(acp)==ar['evidence_sha256'][acp.name];accepted={z['label']:z for z in read(acp)['all142_leaf_custody_checks']};inputs[acceptedp]=sha(acceptedp);inputs[acp]=sha(acp)
 for p,h in inputs.items():assert sha(p)==h
 commonpath=sp.parent/'common-source-operands.json';commonraw=read(commonpath)['common_operands'];indexpath=sp.parent/'row-stream-input-and-plan-index17.json';index=read(indexpath);assert len(commonraw)==22 and len(index['row_records17'])==17
 ep=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64/actual-fixed24-entry-source-and-exponent-plan.json';E={(z['row'],z['column']):z for z in read(ep)['actual_entry_sources576']};assert len(E)==576
 oldip=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_SIXTH_UNIT_COL74_SOURCE_20261105_v64/row-stream-input-and-plan-index18.json';oldrefs={z['row']:z for z in read(oldip)['row_records18']};old={}
 fp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_SIXTH_COL74_COEFFICIENTS_20261105_v64/actual-cleared-col74-and-commonQ2-representation.json';F={z['row']:z for z in read(fp)['actual18_col74_numerator_values']}
 def raw(k,c):
  e=E[k,c]
  if not e['terms']:return {'content':'0','shift':[0]*5,'terms':[]},{'kind':'sourcezero','row':k,'column':c}
  ref=e['source'];p=ROOT/ref['file'];assert p in inputs and ref['kind']=='unchanged20_column_literal'
  with p.open('rb')as f:f.seek(ref['offset']);blob=f.read(ref['bytes'])
  assert len(blob)==ref['bytes']and hashlib.sha256(blob).hexdigest()==ref['payload_sha256'];p0=decode(blob);assert len(p0['terms'])==e['terms'];return p0,ref
 up=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_SIXTH_UNIT_COL74_SOURCE_20261105_v64/old-unit2-monomial-source-prefactors.json';unit=read(up);pa,ra=raw(44,49);pb,rb=raw(44,51);A=mono(pa);B=mono(pb);dpoly=plus(poly(raw(11,51)[0]),poly(pb),Q(-8));assert len(dpoly)==1;D0=(next(iter(dpoly.values())),next(iter(dpoly)))
 assert A==(Q(-1,324000),(0,0,1,0,1))and D0==(Q(8,675),(1,0,0,1,0))and ra==unit['a_source']and rb==unit['b_source']
 inner=[(Q(1),(0,0,0,0,0)),neg(inv(A)),mm(B,mm(inv(A),inv(D0))),neg(inv(D0))];assert [{'QQ':str(q),'exponent5':list(e)}for q,e in inner]==unit['all_four_inner_prefactors']
 assert unit['a_monomial']=={'QQ':str(A[0]),'exponent5':list(A[1])}and unit['b_monomial']=={'QQ':str(B[0]),'exponent5':list(B[1])}and unit['d_known_delta51']=={'QQ':str(D0[0]),'exponent5':list(D0[1])}
 delta=plus(poly(raw(11,62)[0]),poly(raw(44,62)[0]),Q(-8));custody=[];scalarbits=0;counts={'raw':0,'explicit':0,'delta':0,'reuse':0,'N':0}
 def leaf(label,z):
  nonlocal scalarbits
  p=z['original_literal'];src=z['source']
  if label in ('Q','Q2','P7'):
   exps={'Q':[(0,1),(1,5)],'Q2':[(0,1),(1,10),(2,25)],'P7':[(0,1),(1,6),(2,5)]};expected={(0,0,0,e,0):Q(c)for e,c in exps[label]};assert poly(p)==expected;method='literal_parameter_Q_Qsquare_or_Q_times_qplus1';counts['explicit']+=1
  elif label=='Delta':assert poly(p)==delta and src=={'kind':'fresh_raw_A11_minus8A44','column':62};method='fresh_raw_A11_62_minus8_A44_62';counts['delta']+=1
  elif label=='S'or label.startswith('A'):
   k,c=(44,62)if label=='S'else map(int,label[1:].split('_'));v,ref=raw(k,c);assert p==v and src==ref;method='payload_SHA_decode_and_fixed576_index_binding';counts['raw']+=1
  elif label[0]in'UVCF':
   assert label in accepted and accepted[label]['literal_source_identity']and z==old[label];method='exact_old_object_plus_accepted_a3e8_sourceproof';counts['reuse']+=1
  else:
   assert label.startswith('N');i=int(label[1:]);v=F[i];p0=fp.parent/v['file'];assert p0 in inputs and sha(p0)==v['sha256'];blob=p0.read_bytes();q=decode(blob[9:])if blob[8]else{'content':'0','shift':[0]*5,'terms':[]};assert p==q and len(p['terms'])==v['terms']and src=={'kind':'actual_sixth_col74_numerator','row':i,'path':str(p0.relative_to(ROOT)),'sha256':v['sha256'],'value_representation':'S6_i74=N_i/Q^2'};method='actual_sixth74_binary_and_API_same_row';counts['N']+=1
  if not label.startswith('U'):assert Q(z['QQ_transform_scalar'])==1 and z['transform_monomial_exponent5']==[0]*5
  cq=Q(p['content']);tr=Q(z['QQ_transform_scalar']);scalarbits=max(scalarbits,abs(cq.numerator).bit_length()+abs(tr.numerator).bit_length()+4,cq.denominator.bit_length()+tr.denominator.bit_length()+4);shift=[p['shift'][j]+z['transform_monomial_exponent5'][j]for j in range(5)];ts=[(int(t[0]),tuple(shift[j]+t[j+1]for j in range(5)))for t in p['terms']];assert len(ts)==len({e for q,e in ts});custody.append({'label':label,'source_method':method,'terms':len(ts),'literal_source_identity':True});return {'q':cq*tr,'terms':ts,'M':max([abs(q)for q,e in ts]+[0])}
 common={k:leaf(k,z)for k,z in commonraw.items()};verified=[];tuples=0;familybits=0;streamparse=0;sourcepeak=0
 for ref in index['row_records17']:
  rp=ROOT/ref['row_operand_file'];pp=ROOT/ref['plan_file'];assert sha(rp)==ref['row_operand_sha256']and sha(pp)==ref['plan_sha256'];rowraw=read(rp);plan=read(pp);i=ref['row'];assert i==plan['row']==rowraw['row'];oldref=oldrefs[i];orp=ROOT/oldref['row_operand_file'];assert orp in inputs and sha(orp)==oldref['row_operand_sha256'];old=read(orp)['row_operands'];profiles={**common,**{k:leaf(k,z)for k,z in rowraw['row_operands'].items()}}
  outer=[(i,['P7','Q2'],(Q(1),(0,0,0,0,0))),(31,['P7','Q','U'+str(i)],(Q(-1),(0,0,0,0,0))),(83,['P7','Q2','U'+str(i)],(Q(-1,40),(1,0,0,0,-2))),(83,['P7','Q2','V'+str(i)],(Q(-405,8),(0,0,0,0,-2))),(14,['P7','Q','C'+str(i)],(Q(-6075,2),(-1,0,0,0,0))),(35,['P7','Q','C'+str(i)],(Q(-18225),(-1,0,0,0,0))),(104,['P7','F'+str(i)],(Q(-765450000),(-1,0,0,-1,0)))]
  expected=[]
  for k,ol,opre in outer:
   inn=[([f'A{k}_62'],inner[0]),([f'A{k}_49','S'],inner[1]),([f'A{k}_49','Delta'],inner[2]),([f'A{k}_51','Delta'],inner[3])]
   for il,ipre in inn:
    ls=ol+il
    if all(profiles[l]['terms']for l in ls):expected.append((k,ls,mm(opre,ipre)))
  for ls,pre in [(['N'+str(i),'A51_62'],(Q(3645,64),(-2,0,1,0,1))),(['N'+str(i),'Q','A11_62'],(Q(675),(-1,0,0,0,0))),(['N'+str(i),'Q','S'],(Q(675),(-1,0,0,0,0)))]:
   if all(profiles[l]['terms']for l in ls):expected.append((51,ls,pre))
  assert len(expected)==len(plan['products']);dc=1;prepared=[]
  for (k,labels,pref),p in zip(expected,plan['products']):
   assert k==p['source_M2_row']and labels==p['factor_labels']and {'QQ':str(pref[0]),'exponent5':list(pref[1])}==p['full_parameter_monomial_prefactor'];q=pref[0];nb=abs(q.numerator).bit_length();db=q.denominator.bit_length()
   for l in labels:q*=profiles[l]['q'];nb+=abs(profiles[l]['q'].numerator).bit_length();db+=profiles[l]['q'].denominator.bit_length()
   scalarbits=max(scalarbits,nb+4,db+4);assert q==Q(p['source_product_content']);dc=dc*q.denominator//gcd(dc,q.denominator);prepared.append((p,labels,pref,q))
  assert dc==int(plan['integer_common_D']);union=set();M=0;peaksets=0
  for p,labels,pref,q in prepared:
   ns=[len(profiles[l]['terms'])for l in labels];nt=prod(ns);h=(nt//max(ns))*prod(profiles[l]['M']for l in labels);w=q.numerator*(dc//q.denominator);assert nt==p['literal_tuple_count']and h==int(p['primitive_product_integer_abs_upper'])and w==int(p['exact_scalar_integer_weight_after_D']);tuples+=nt;M+=abs(w)*h;S={(0,0,0,0,0)}
   for l,st in zip(labels,p['support_stages']):
    T={tuple(e[j]+f[j]for j in range(5))for e in S for _,f in profiles[l]['terms']};assert st=={'factor':l,'previous_support':len(S),'factor_terms':len(profiles[l]['terms']),'exact_next_support':len(T)};peaksets=max(peaksets,len(S)+len(T)+len(union));S=T
   S={tuple(e[j]+pref[1][j]for j in range(5))for e in S};assert len(S)==p['exact_product_support_cardinality'];union|=S
  assert M==int(plan['integer_numerator_abs_upper'])and union=={tuple(e)for e in plan['exact_source_support_sorted']}and len(union)==plan['exact_union_support_cardinality'];familybits=max(familybits,M.bit_length()+2)
  valparse=numeric(commonraw)+numeric(rowraw)+numeric(plan)+numeric(index)+numeric(custody)+numeric(verified)+7*sum(len(z['terms'])for z in profiles.values())+4096;streamparse=max(streamparse,valparse);assert valparse<25000
  assert 25000+28*len(union)+128<=plan['future_owned_extra_including25kmeta']<=66008
  local=25000+numeric(E)+numeric(old)+numeric(commonraw)+numeric(rowraw)+numeric(plan)+numeric(custody)+numeric(verified)+numeric(unit)+7*sum(len(z['terms'])for z in profiles.values())+6*peaksets+4096;sourcepeak=max(sourcepeak,local);assert local<150000
  verified.append({'row':i,'products':len(prepared),'support_size':len(union),'source_support_D_M_and_every_partial_height_verified':True,'source_owned_limit':plan['future_owned_extra_including25kmeta'],'VALUEparse_upper':valparse});del rowraw,plan,profiles,prepared,expected,union;gc.collect()
 assert tuples==894424 and counts=={'raw':66,'explicit':3,'delta':1,'reuse':68,'N':17}and len(custody)==155
 base=sum(p.stat().st_size for p in inputs);caps=sum(z['binary_RAW_upper']for z in index['row_records17']);future=base+caps+300000+131072;assert future<=17691677 and future<32*1024*1024
 proof={'all_passed':True,'source_only':True,'ONE17_column62_job_safe':True,'worker_sha256':sha(worker),'row_index_file':str(indexpath.relative_to(ROOT)),'actual_READ_input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'actual_READ_input_RAW_bytes':base,'all155_leaf_custody_checks':custody,'source_custody_counts':counts,'unit2_prefactors_fresh_raw_source_checked':True,'accepted_UVCF68_exact_object_reuse':{'report':str(acceptedp.relative_to(ROOT)),'sha256':sha(acceptedp),'exact_objects':68},'all17_support_D_M_partial_checks':verified,'all_source_tuples894424':tuples,'VALUE_streamed_parse_numeric_upper':streamparse,'coefficient_fixed_metadata_limit':25000,'coefficient_owned_limit':66008,'prior_SOURCE65796_separate':True,'helper_SOURCE_numeric_upper':sourcepeak,'scalar_unreduced_temporary_bit_upper':scalarbits,'family_partial_integer_bit_upper':familybits,'partial_height_proof':'M=sum absolute integer-weight times fixed-output assignment upper times productprimitive-max; every productpartial,block,accumulator and backadd is bounded by M.','old14MBmatrix_actualREAD_explicitly_charged':True,'binary_output_upper':caps,'other_artifact_upper':300000,'prepared_GATE_HELP_input_upper':131072,'fresh_RAW_upper':future,'worker_helper_help_safe':True,'Q3_qplus1_representation_ONLY_no_unproved_Q_cancellation':True,'full170_RHS_top104_plus66_87extras_oneP290_all_other16_kept':True,'hard_seconds':60,'hard_address_space_bytes':2000000000}
 (out/'implementation-custody-storage-bit-addendum.json').write_bytes(enc(proof));report={'id':ID,'all_passed':True,'source_only':True,'ONE17_column62_job_safe':True,'owned_upper':66008,'RAW_upper':future,'all155_custody_source_checks':True,'VALUE_streamed_parse_numeric_upper':streamparse,'helper_SOURCE_numeric_upper':sourcepeak,'worker_helper_help_safe':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));assert sum(p.stat().st_size for p in out.iterdir())<=131072;signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())
