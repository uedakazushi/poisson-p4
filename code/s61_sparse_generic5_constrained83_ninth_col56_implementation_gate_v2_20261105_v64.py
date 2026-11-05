#!/usr/bin/env python3
"""SOURCE-only full custody/code/lifetime gate for the prepared staged15 job."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,json,hashlib,ast,re,time,signal,resource,subprocess,sys,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_IMPLEMENTATION_GATE_V2_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def digest(v):return hashlib.sha256(enc(v)).hexdigest()
def scalars(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',v)else 0
 if isinstance(v,list):return sum(map(scalars,v))
 if isinstance(v,dict):return sum(map(scalars,v.values()))
 return 0
def poly(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']}
def main(a):
 assert a.freeze
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);worker=ROOT/'code/s61_sparse_generic5_constrained83_ninth_col56_staged_coefficients_v2_20261105_v64.py';target=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64';assert not target.exists()
 for code in(worker,Path(__file__)):
  tree=ast.parse(code.read_text());fn=next(z for z in tree.body if isinstance(z,ast.FunctionDef)and z.name=='main');defs={z.name for z in fn.body if isinstance(z,ast.FunctionDef)};names=set()
  def scan(z):
   if isinstance(z,(ast.FunctionDef,ast.ClassDef,ast.Lambda)):return
   if isinstance(z,ast.Name)and isinstance(z.ctx,ast.Store):names.add(z.id)
   for v in ast.iter_child_nodes(z):scan(v)
  for z in fn.body:scan(z)
  assert not defs&names,(code,defs&names)
 for name,code in [('worker-help.txt',worker),('helper-help.txt',Path(__file__))]:
  with(out/name).open('wb')as f:subprocess.run([sys.executable,str(code),'--help'],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=10)
 assert not target.exists()
 tp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_THIN_CONTROL_SOURCE_20261105_v64/report.json';assert sha(tp)=='201a0b9f1c6141c15681c96f9f63d3219586dbc46f4b074a6b715fd75f12b403';tr=read(tp);assert tr['all_passed'];inputs={ROOT/p:h for p,h in tr['input_sha256'].items()};inputs[tp]=sha(tp)
 for name,h in tr['evidence_sha256'].items():inputs[tp.parent/name]=h
 inputs.update({worker:sha(worker),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64.md':sha(ROOT/'research/notes/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64.md')})
 codec=ROOT/'code/s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64.py';inputs[codec]=sha(codec);cp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STREAM_CODEC_SOURCE_GATE_20261105_v64/report.json';assert sha(cp)=='c9921d8f8db77f5a4cf2739b26d88a4728cb871c4e58b4ed9aa748a6d09e2d6f';cr=read(cp);assert cr['all_passed'];inputs[cp]=sha(cp)
 for name,h in cr['evidence_sha256'].items():inputs[cp.parent/name]=h
 ap=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_IMPLEMENTATION_GATE_20261105_v64/report.json';assert sha(ap)=='00170ebc42ae16908f3fe3178386d23b15b5adc5697d7f5704a9d37a89d5260b';ar=read(ap);ac=ap.parent/'implementation-custody-storage-bit-addendum.json';assert sha(ac)==ar['evidence_sha256'][ac.name];accepted={z['label']:z for z in read(ac)['all165_leaf_custody_checks']};inputs[ap]=sha(ap);inputs[ac]=sha(ac)
 for p,h in inputs.items():assert sha(p)==h
 program_path=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64/staged-literal-col56-program-and-positive-profiles.json';program=read(program_path);ops=program['literal_ops'];profiles=program['positive_support_and_height_profiles'];control_path=tp.parent/'thin-controller.json';control=read(control_path);thinproof=read(tp.parent/'thin-custody-lifetime-and-RAW-proof.json');assert thinproof['future_fixed_metadata_upper']==20000;fragments=[];records=[]
 for ref in control['phase_index16']:
  p=ROOT/ref['file'];assert sha(p)==ref['sha256'];f=read(p);fragments.append(f)
  for z in f['literal_records']:records.append(z)
 assert len(records)==len(ops)==700
 # This all-fragment parse belongs ONLY to SOURCE, never to VALUE.
 sourcepeak=25000+scalars(program)+scalars(fragments)+scalars(control)+scalars(accepted);assert sourcepeak<150000
 Epath=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64/actual-fixed24-entry-source-and-exponent-plan.json';E={(z['row'],z['column']):z for z in read(Epath)['actual_entry_sources576']};oldidx=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_UNIT_COL82_SOURCE_20261105_v64/row-stream-input-and-plan-index16.json';oldrefs={z['row']:z for z in read(oldidx)['row_records16']};api_path=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64/actual-cleared-col82-and-commonQ3RE-representation.json';api={z['row']:z for z in read(api_path)['actual16_col82_numerator_values']};scope_path=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL82_NINTH_UNIT_SOURCE_20261105_v64/actual-row35-ninth-unit-and-full168x15-operator.json';scope=read(scope_path)
 for p in(Epath,oldidx,api_path,scope_path):assert p in inputs
 def original_raw(row,col):
  z=E[row,col]
  if not z['terms']:return {'content':'0','shift':[0]*5,'terms':[]},{'kind':'sourcezero','row':row,'column':col}
  src=z['source'];p=ROOT/src['file'];assert p in inputs
  with p.open('rb')as f:f.seek(src['offset']);b=f.read(src['bytes'])
  assert hashlib.sha256(b).hexdigest()==src['payload_sha256'];return decode(b),src
 delta=poly(original_raw(11,56)[0]);rr=poly(original_raw(44,56)[0])
 for e,q in rr.items():delta[e]=delta.get(e,Q(0))-8*q
 delta={e:q for e,q in delta.items()if q};counts={'raw':0,'explicit':0,'delta':0,'F9':0,'cached':0,'Z':0};leafchecks=[];scalarbits=0;familybits=0
 for i,(rec,op)in enumerate(zip(records,ops)):
  assert rec['index']==i and rec['original_op']==op and rec['original_op_sha256']==digest(op);name=op['output'];pout=profiles[name];assert rec['output_profile']==pout and rec['output_profile_sha256']==digest(pout);n=pout['terms_upper'];M=int(pout['primitive_integer_abs_upper']);q=Q(pout['known_QQ_content']);familybits=max(familybits,M.bit_length())
  if op['op']=='input':
   ref=rec['input_leaf'];p=ROOT/ref['file'];assert sha(p)==ref['sha256'];leaf=read(p);z=leaf['operand'];binding=op['binding'];parent=ROOT/leaf['parent_path'];assert parent in inputs and sha(parent)==leaf['parent_sha256'];data=read(parent);old=data if binding.get('entire_operand')else(data['common_operands']if'common_operands'in data else data['row_operands'])[binding['label']];assert z==old and leaf['original_binding']==binding and leaf['original_input_op_sha256']==digest(op)
   lit=z['original_literal'];src=z['source'];method=''
   if name in('Q','Q2','P7','E'):
    expected={'Q':{(0,0,0,0,0):Q(1),(0,0,0,1,0):Q(5)},'Q2':{(0,0,0,0,0):Q(1),(0,0,0,1,0):Q(10),(0,0,0,2,0):Q(25)},'P7':{(0,0,0,0,0):Q(1),(0,0,0,1,0):Q(6),(0,0,0,2,0):Q(5)},'E':{(0,0,0,0,0):Q(35),(0,0,1,0,0):Q(32),(0,0,0,1,0):Q(15)}}[name];assert poly(lit)==expected;counts['explicit']+=1;method='literal_existing_parameter_unit'
   elif name=='F9':
    assert lit=={'content':'1','shift':[0]*5,'terms':[[c,*e]for e,c in scope['primitive_F9_integer_terms75']]};counts['F9']+=1;method='accepted_primitive_F9_75'
   elif name=='Delta':assert poly(lit)==delta;counts['delta']+=1;method='fresh_old_raw_A11_minus8A44'
   elif name=='S'or name.startswith('A'):
    row,col=(44,56)if name=='S'else map(int,name[1:].split('_'));v,ref0=original_raw(row,col);assert lit==v and src==ref0;counts['raw']+=1;method='fresh_E576_and_authenticated_payload'
   elif name.startswith('Z'):
    row=int(name[1:]);v=api[row];path=api_path.parent/v['file'];assert src=={'kind':'actual_existing_N8_i82','row':row,'path':str(path.relative_to(ROOT)),'sha256':v['sha256']};assert sha(path)==v['sha256'];b=path.read_bytes();assert lit==(decode(b[9:])if b[8]else{'content':'0','shift':[0]*5,'terms':[]});counts['Z']+=1;method='actual_primary_col82_binary'
   else:
    assert name[0]in'UVCFNM'and name in accepted and accepted[name]['literal_source_identity'];row=int(name[1:]);ref0=oldrefs[row];oldpath=ROOT/ref0['row_operand_file'];assert sha(oldpath)==ref0['row_operand_sha256'];assert z==read(oldpath)['row_operands'][name];counts['cached']+=1;method='exact_old_object_plus_00170_actual_sourceproof'
   if name[0]not in'UVCFNM'or name=='F9':assert Q(z['QQ_transform_scalar'])==1 and z['transform_monomial_exponent5']==[0]*5
   qq=Q(lit['content'])*Q(z['QQ_transform_scalar']);S={tuple(lit['shift'][j]+t[j+1]+z['transform_monomial_exponent5'][j]for j in range(5))for t in lit['terms']};assert qq==q and len(S)==len(lit['terms'])==n and max([abs(int(t[0]))for t in lit['terms']]+[0])==M;assert digest([list(e)for e in sorted(S)])==pout['support_sha256'];box=[[min(e[j]for e in S),max(e[j]for e in S)]for j in range(5)]if S else[[0,0]]*5;assert box==pout['exponent_box5'];leafchecks.append({'op':i,'node':name,'terms':n,'method':method,'source_and_positive_profile_exact':True});sourcepeak=max(sourcepeak,25000+scalars(program)+scalars(fragments)+scalars(E)+scalars(data)+scalars(leaf)+5*n+scalars(leafchecks)+4096);assert sourcepeak<150000;del leaf,data,z,lit,S
  elif op['op']=='copy':
   A=profiles[op['A']];assert(q,M,n,pout['support_sha256'],pout['exponent_box5'])==(Q(A['known_QQ_content']),int(A['primitive_integer_abs_upper']),A['terms_upper'],A['support_sha256'],A['exponent_box5'])
  elif op['op']=='multiply':
   A=profiles[op['A']];B=profiles[op['B']];qa=Q(A['known_QQ_content']);qb=Q(B['known_QQ_content']);c=Q(op['scalar']);assign=min(A['terms_upper'],B['terms_upper']);assert q==qa*qb*c and M==assign*int(A['primitive_integer_abs_upper'])*int(B['primitive_integer_abs_upper']);assert op['positive_fixed_output_assignment_upper']==assign and op['primitive_pair_upper']==A['terms_upper']*B['terms_upper'];box=[[A['exponent_box5'][j][0]+B['exponent_box5'][j][0]+op['shift5'][j],A['exponent_box5'][j][1]+B['exponent_box5'][j][1]+op['shift5'][j]]for j in range(5)]if assign else[[0,0]]*5;assert box==pout['exponent_box5'];scalarbits=max(scalarbits,abs(qa.numerator).bit_length()+abs(qb.numerator).bit_length()+abs(c.numerator).bit_length()+4,qa.denominator.bit_length()+qb.denominator.bit_length()+c.denominator.bit_length()+4)
  elif op['op']=='add':
   A=profiles[op['A']];B=profiles[op['B']];qa=Q(A['known_QQ_content']);qb=Q(B['known_QQ_content'])*Q(op['scalarB']);d=qa.denominator*qb.denominator//gcd(qa.denominator,qb.denominator);ai=(qa*d).numerator;bi=(qb*d).numerator;g0=gcd(abs(ai),abs(bi));outq=Q(g0,d)if g0 else Q(0);wa=ai//g0 if g0 else 0;wb=bi//g0 if g0 else 0;assert q==outq and str(wa)==op['integer_weight_A']and str(wb)==op['integer_weight_B']and M==abs(wa)*int(A['primitive_integer_abs_upper'])+abs(wb)*int(B['primitive_integer_abs_upper']);assert max(A['terms_upper'],B['terms_upper'])<=n<=A['terms_upper']+B['terms_upper'];active=[z for z in(A,B)if z['terms_upper']];box=[[min(z['exponent_box5'][j][0]for z in active),max(z['exponent_box5'][j][1]for z in active)]for j in range(5)]if active else[[0,0]]*5;assert box==pout['exponent_box5'];scalarbits=max(scalarbits,qa.denominator.bit_length()+qb.denominator.bit_length()+abs(ai).bit_length()+abs(bi).bit_length()+8)
  else:raise AssertionError(op['op'])
 assert len(leafchecks)==181 and counts=={'raw':64,'explicit':4,'delta':1,'F9':1,'cached':96,'Z':15}
 # Independent complete liveness and all per-op buffer forecasts, with20k META.
 uses={k:0 for k in profiles};outputs={z['node']:z for z in program['row15_output_nodes']}
 for op in ops:
  for key in('A','B'):
   if key in op:uses[op[key]]+=1
 for z in outputs:uses[z]+=1
 live={};maximum=20000;outputcaps={};phasechecks=[];encoders={z['row']:z for z in program['streaming_binary_phase_forecasts']};currentPhase=0;phaseMax=20000
 def rawcap(p):
  q=Q(p['known_QQ_content']);m=int(p['primitive_integer_abs_upper'])*abs(q.numerator);box=p['exponent_box5'];eb=sum(max(1,((b-a).bit_length()+6)//7)for a,b in box);return 2048+(q.denominator.bit_length()+7)//8+p['terms_upper']*((m.bit_length()+7)//8+eb+16)
 control=read(tp.parent/'thin-controller.json');phaseEnds={x['phase']['end']:i for i,x in enumerate(control['phase_index16'])}
 for i,(rec,op)in enumerate(zip(records,ops)):
  name=op['output'];n=profiles[name]['terms_upper'];reuse=op['op']=='add'and uses[op['A']]==1;assert reuse==rec['last_use']['reuse_left_last_use_buffer'];growth=n-profiles[op['A']]['terms_upper']if reuse else n;assert growth>=0;peak=20000+7*(sum(live.values())+growth)+1024;maximum=max(maximum,peak);phaseMax=max(phaseMax,peak);live[name]=n;release=[]
  if uses[name]==0:release.append(name);live.pop(name,None)
  for key in('A','B'):
   if key in op:
    uses[op[key]]-=1
    if uses[op[key]]==0:release.append(op[key]);live.pop(op[key],None)
  o=outputs.get(name)
  if o:
   expected=[k for k in release if k!=name];assert rec['last_use']['release_after_arithmetic_before_encoding']==expected and rec['last_use']['release_output_after_encoding']==name;ep=20000+7*(sum(live.values())-n)+10*n+2048;maximum=max(maximum,ep);phaseMax=max(phaseMax,ep);outputcaps[str(o['row'])]=rawcap(profiles[name]);uses[name]-=1;live.pop(name,None)
  else:assert rec['last_use']['release_after_arithmetic_before_encoding']==release and rec['last_use']['release_output_after_encoding']is None
  if i in phaseEnds:
   if phaseEnds[i]==0:
    aux=control['shared_N8_35_56_node'];naux=profiles[aux]['terms_upper'];ep=20000+7*(sum(live.values())-naux)+10*naux+2048;maximum=max(maximum,ep);phaseMax=max(phaseMax,ep);outputcaps['shared35']=rawcap(profiles[aux])
   phasechecks.append({'phase':phaseEnds[i],'last_op':i,'prearithmetic_whole_owned_upper':phaseMax});phaseMax=20000
 assert all(v==0 for v in uses.values())and not live and maximum<=130791
 # Compact proof, never duplicate all700 positive profiles into VALUE metadata.
 proof={'all_passed':True,'prepared_one15_staged_job_safe':True,'worker_sha256':sha(worker),'stream_codec_sha256':sha(codec),'thin_controller_file':str((tp.parent/'thin-controller.json').relative_to(ROOT)),'actual_READ_inputs_source':'GATE report.input_sha256','actual_READ_inputs_map_sha256':digest({str(p.relative_to(ROOT)):h for p,h in inputs.items()}),'actual_READ_input_RAW_bytes':sum(p.stat().st_size for p in inputs),'fixed_metadata_limit':20000,'coefficient_owned_upper':130791,'all181_input_custody':leafchecks,'source_custody_counts':counts,'all700_contents_positive_integer_partial_bounds_and_lastuse_checked':True,'scalar_content_unreduced_temporary_bit_upper':scalarbits,'all_integer_coefficient_positive_bound_bits':familybits,'per_phase_prearithmetic_owned_upper':phasechecks,'output_caps':outputcaps,'one_shared_auxiliary_encoding_in_budget':True,'all16_stream_encoder_roundtrips_required':True,'VALUE_NO_full700_program_or_rowbundle_or_oldE_JSON_parse':True,'ROOT_independent_SOURCE_code_review_REQUIRED_before_execute':True,'ROOT_review_raw_limit':65536,'ROOT_review_numeric_limit':2000,'prepared_GATE_HELP_upper':131072,'coefficient_report_backadd_prefix_artifact_upper':300000,'coefficient_seconds':60,'coefficient_address_space_bytes':2000000000,'all168x15_87plus66_RHS_other14_oneP290_complements_kept':True}
 # The thin proof reserved10k for these actual proof/report/runtime structures.
 extra=scalars(proof)+2000+4*700+20*16+2048+100;assert extra<=10000,extra
 future=proof['actual_READ_input_RAW_bytes']+sum(outputcaps.values())+131072+65536+300000;assert future<=25200267 and future<32*1024*1024;proof['future_extra_metadata_reserved10k_actual_upper']=extra;proof['fresh_RAW_upper']=future
 (out/'implementation-custody-streaming-lifetime-proof.json').write_bytes(enc(proof));report={'id':ID,'all_passed':True,'source_only':True,'prepared_one15_staged_job_safe':True,'worker_SHA':sha(worker),'all181_source_custody_and_all700_positive_buffer_checks':True,'owned_upper':130791,'VALUE_fixed_metadata_upper':20000,'future_extra_metadata_upper':extra,'fresh_RAW_upper':future,'SOURCE_owned_upper':sourcepeak,'ROOT_review_REQUIRED_before_coefficients':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};assert sum(p.stat().st_size for p in out.iterdir())+len(enc(report))<=131072;assert not target.exists();(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())
