#!/usr/bin/env python3
"""SOURCE-only staged col56 contraction: supports, scalar heights and last use."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,json,hashlib,time,signal,resource,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
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
def profile(z):
 p=z['original_literal'];s=[p['shift'][j]+z['transform_monomial_exponent5'][j]for j in range(5)];return {'S':{tuple(s[j]+t[j+1]for j in range(5))for t in p['terms']},'q':Q(p['content'])*Q(z['QQ_transform_scalar']),'M':max([abs(int(t[0]))for t in p['terms']]+[0])}
def supporthash(S):return hashlib.sha256(enc([list(e)for e in sorted(S)])).hexdigest()
def main(a):
 assert a.build
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);(out/'rows').mkdir();inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py')}
 def gate(stem,h):
  p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];inputs[p]=h;return p,r
 def ev(p,r,name):
  f=p.parent/name;assert sha(f)==r['evidence_sha256'][name];inputs[f]=sha(f);return read(f)
 np,nr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_COL82_NINTH_UNIT_SOURCE_20261105_v64','7dbaf504004a61b9c2554b6c9c31bfffaf426bb324d15dca31221d365c009f95');scope=ev(np,nr,'actual-row35-ninth-unit-and-full168x15-operator.json');cost=ev(np,nr,'remaining15-literal-source-cost-comparison.json');assert cost['minimum_in_this_exact_source_recipe']==56
 ep,er=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64','a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a');E={(z['row'],z['column']):z for z in ev(ep,er,'actual-fixed24-entry-source-and-exponent-plan.json')['actual_entry_sources576']}
 sp,sr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_UNIT_COL82_SOURCE_20261105_v64','6925f8f7072565b98804d3c87262a5a20b5bdc05cc5f3171884f6f73d97bfdfb');oldcommon=ev(sp,sr,'common-source-operands.json')['common_operands'];idx=ev(sp,sr,'row-stream-input-and-plan-index16.json');refs={z['row']:z for z in idx['row_records16']}
 cp,cr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64','f7f9bc898a54b9a457a59ba4146f38b0e5a8ce9f2cbda906379a71db20e5f9cb');api=ev(cp,cr,'actual-cleared-col82-and-commonQ3RE-representation.json');values={z['row']:z for z in api['actual16_col82_numerator_values']}
 checked=set()
 def raw(k,c):
  e=E[k,c]
  if not e['terms']:return {'content':'0','shift':[0]*5,'terms':[]},{'kind':'sourcezero','row':k,'column':c}
  ref=e['source'];p=ROOT/ref['file'];assert ref['kind']=='unchanged20_column_literal'
  if p not in checked:assert sha(p)==er['input_sha256'][str(p.relative_to(ROOT))];inputs[p]=sha(p);checked.add(p)
  with p.open('rb')as f:f.seek(ref['offset']);b=f.read(ref['bytes'])
  assert hashlib.sha256(b).hexdigest()==ref['payload_sha256'];z=decode(b);assert len(z['terms'])==e['terms'];return z,ref
 def operand(p,src):return {'original_literal':p,'source':src,'QQ_transform_scalar':'1','transform_monomial_exponent5':[0]*5}
 common={k:oldcommon[k]for k in('Q','Q2','P7','E')};j=56
 z,r=raw(44,j);common['S']=operand(z,r)
 # Only actual old input coefficients are added for the literal Delta leaf.
 z11,_=raw(11,j);A={tuple(z['shift'][h]+t[h+1]for h in range(5)):Q(z['content'])*int(t[0])for t in z['terms']}
 B={tuple(z11['shift'][h]+t[h+1]for h in range(5)):Q(z11['content'])*int(t[0])for t in z11['terms']}
 for e,q in A.items():B[e]=B.get(e,Q(0))-8*q
 B={e:q for e,q in B.items()if q};dd=1
 for q in B.values():dd=dd*q.denominator//gcd(dd,q.denominator)
 zz={e:(q*dd).numerator for e,q in B.items()};gg=0
 for v in zz.values():gg=gcd(gg,abs(v))
 lo=[min(e[h]for e in zz)for h in range(5)]if zz else[0]*5;p={'content':str(Q(gg,dd))if zz else'0','shift':lo,'terms':[[str(v//gg),*[e[h]-lo[h]for h in range(5)]]for e,v in sorted(zz.items(),reverse=True)]};common['Delta']=operand(p,{'kind':'literal_old_raw_A11_minus8A44','column':j});del A,B,zz,z,z11
 for k in(31,83,14,35,104):
  for c in(j,49,51):z,r=raw(k,c);common[f'A{k}_{c}']=operand(z,r)
 for k in(51,11,52):z,r=raw(k,j);common[f'A{k}_{j}']=operand(z,r)
 f9=scope['primitive_F9_integer_terms75'];common['F9']=operand({'content':'1','shift':[0]*5,'terms':[[c,*e]for e,c in f9]},{'kind':'accepted_fixed_row35_primitive_F9_75','report':str((np/'report.json').relative_to(ROOT))if np.is_dir()else str(np.relative_to(ROOT))})
 (out/'common-source-operands.json').write_bytes(enc({'all_passed':True,'common_operands':common,'parameter_order':['H','L','s0','q0','t0']}));P={};nodes={};ops=[];inputrefs=[];phases=[];sourcepeak=25000;tupleops=0;resident_n=0;meta_numeric=0;parse_extra=0;context={}
 base=25000+numeric(E)+numeric(scope)+numeric(cost)+numeric(api)+numeric(idx)+numeric(common)+4096;assert base<150000
 def own(extra=0):
  nonlocal sourcepeak
  n=base+parse_extra+5*resident_n+meta_numeric+extra;sourcepeak=max(sourcepeak,n)
  if n>=150000:
   (out/'failure.json').write_bytes(enc({'all_passed':False,'source_only':True,'stop':'prospective SOURCE owned guard before new support insertion','forecast':n,'limit':150000,'context':context,'completed_rows':[z['row']for z in phases if 'row'in z],'new_matrix_coefficients':0,'no_automatic_retry':True}))
   with(out/'partial-source-program.json').open('w')as f:json.dump({'ops':ops,'profiles':nodes,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()}},f,separators=(',',':'))
   raise RuntimeError('SOURCE prospective owned guard '+str(n))
  return n
 def save(name,z,op):
  nonlocal resident_n,meta_numeric
  assert name not in nodes;P[name]=z;n=len(z['S']);resident_n+=n;box=[[min(e[h]for e in z['S']),max(e[h]for e in z['S'])]for h in range(5)]if n else[[0,0]]*5;nodes[name]={'terms_upper':n,'known_QQ_content':str(z['q']),'primitive_integer_abs_upper':str(z['M']),'support_sha256':supporthash(z['S']),'exponent_box5':box};ops.append({'output':name,**op});meta_numeric+=numeric(nodes[name])+numeric(ops[-1]);own()
 def leaf(name,z,ref):
  nonlocal context
  context={'op':'input_support','output':name};own(5*len(z['original_literal']['terms'])+64);save(name,profile(z),{'op':'input','binding':ref})
 def drop(*names):
  nonlocal resident_n
  for name in names:
   z=P.pop(name,None)
   if z is not None:resident_n-=len(z['S'])
 def copy_node(dest,A):
  nonlocal context
  context={'op':'copy_support','output':dest,'A':A};own(5*len(P[A]['S'])+64);save(dest,{'S':set(P[A]['S']),'q':P[A]['q'],'M':P[A]['M']},{'op':'copy','A':A})
 def mul(dest,A,B,q=Q(1),shift=(0,0,0,0,0)):
  nonlocal tupleops,context
  context={'op':'multiply_support','output':dest,'A':A,'B':B};x=P[A];y=P[B];S=set()
  for e in x['S']:
   for f in y['S']:
    g=tuple(e[h]+f[h]+shift[h]for h in range(5))
    if g not in S:own(5*(len(S)+1)+64);S.add(g)
  n=min(len(x['S']),len(y['S']));M=n*x['M']*y['M'];qq=x['q']*y['q']*q
  if not S:qq=Q(0);M=0
  tuples=len(x['S'])*len(y['S']);tupleops+=tuples;save(dest,{'S':S,'q':qq,'M':M},{'op':'multiply','A':A,'B':B,'scalar':str(q),'shift5':list(shift),'primitive_pair_upper':tuples,'positive_fixed_output_assignment_upper':n})
 def add(dest,A,B,q=Q(1),consume=False):
  nonlocal context
  context={'op':'add_support','output':dest,'A':A,'B':B};x=P[A];y=P[B];a=x['q'];b=y['q']*q;d=a.denominator*b.denominator//gcd(a.denominator,b.denominator);ai=(a*d).numerator;bi=(b*d).numerator;g=gcd(abs(ai),abs(bi));qq=Q(g,d)if g else Q(0);wa=ai//g if g else 0;wb=bi//g if g else 0;M=abs(wa)*x['M']+abs(wb)*y['M']
  if consume:
   assert A.endswith('_ninth_A')and B.endswith('_ninth_B')and dest.endswith('_ninth_output');S=x['S'];T=y['S'];nB=len(T);drop(A,B);own(5*(len(S)+len(T))+64);processed=0
   while T:
    own(5*(len(S)+len(T))+64);e=T.pop();processed+=1
    if e not in S:S.add(e)
   assert processed==nB
  else:
   S=set(x['S']);own(5*len(S)+64)
   for e in y['S']:
    if e not in S:own(5*(len(S)+1)+64);S.add(e)
  save(dest,{'S':S,'q':qq,'M':M},{'op':'add','A':A,'B':B,'scalarB':str(q),'integer_weight_A':str(wa),'integer_weight_B':str(wb),'lastuse_left_buffer_can_be_stolen_if_A_not_used_again':True,'SOURCE_consumed_both_lastuse_support_buffers':consume})
 def madd(dest,A,B,q,shift,tag):
  block=dest+'_block_'+tag;new=dest+'_sum_'+tag;mul(block,A,B,q,shift);add(new,dest,block);drop(dest,block);return new
 def m2(dest,k):
  copy_node(dest,f'A{k}_{j}');name=dest
  for a,b,q,s,t in [(f'A{k}_49','S',Q(324000),(0,0,-1,0,-1),'a'),(f'A{k}_49','Delta',Q(36000),(0,0,-1,0,-1),'b'),(f'A{k}_51','Delta',Q(-675,8),(-1,0,0,-1,0),'c')]:name=madd(name,a,b,q,s,t)
  return name
 for name,z in common.items():leaf(name,z,{'scope':'common-source-operands.json','label':name})
 del common
 shared={k:m2('shared_M2_'+str(k),k)for k in(31,83,14,35,104)}
 copy_node('shared_G0',shared[14]);add('shared_G','shared_G0',shared[35],Q(6));drop('shared_G0')
 copy_node('shared_B51_0',f'A51_{j}');add('shared_rows11plus44',f'A11_{j}','S');mul('shared_B51_block','Q','shared_rows11plus44',Q(320,27),(1,0,-1,0,-1));add('shared_B51','shared_B51_0','shared_B51_block');drop('shared_B51_0','shared_rows11plus44','shared_B51_block');mul('shared_P7F9','P7','F9')
 sharedkeys=set(P);rowbindings=[]
 def stage8(row):
  nonlocal parse_extra
  prefix=f'row{row}_';r=refs[row];f=ROOT/r['row_operand_file'];assert sha(f)==r['row_operand_sha256'];inputs[f]=sha(f);old=read(f)['row_operands'];bundle={}
  for c in(j,49,51):
   label=f'A{row}_{c}'
   if label not in P:z,src=raw(row,c);bundle[label]=operand(z,src)
  for k in('U','V','C','F','N','M'):bundle[k+str(row)]=old[k+str(row)]
  f=out/'rows'/f'row-{row}-source-operands.json';f.write_bytes(enc({'all_passed':True,'row':row,'row_operands':bundle}));rowbindings.append({'row':row,'file':str(f.relative_to(ROOT)),'sha256':sha(f)})
  parse_extra=numeric(old)+numeric(bundle)+4096;own()
  for name,z in bundle.items():leaf(name,z,{'scope':str(f.relative_to(out)),'label':name})
  del old,bundle;parse_extra=0;B=m2(prefix+'M2',row);H=prefix+'H4a';mul(H,'Q',B);drop(B);H=madd(H,'U'+str(row),shared[31],Q(-1),(0,0,0,0,0),'u31');mul(prefix+'QU','Q','U'+str(row));H=madd(H,prefix+'QU',shared[83],Q(-1,40),(1,0,0,0,-2),'u83');drop(prefix+'QU');mul(prefix+'QV','Q','V'+str(row));H=madd(H,prefix+'QV',shared[83],Q(-405,8),(0,0,0,0,-2),'v83');drop(prefix+'QV');H=madd(H,'C'+str(row),'shared_G',Q(-6075,2),(-1,0,0,0,0),'c');mul(prefix+'F6','Q',H);drop(H);F=madd(prefix+'F6','F'+str(row),shared[104],Q(-765450000),(-1,0,0,-1,0),'f');mul(prefix+'N7','P7',F);drop(F);N=madd(prefix+'N7','N'+str(row),'shared_B51',Q(3645,64),(-2,0,1,0,1),'n');mul(prefix+'N8','E',N);drop(N);N8=madd(prefix+'N8','M'+str(row),f'A52_{j}',Q(-273375,32768),(-3,0,0,0,2),'m');return N8
 pivot=stage8(35);assert pivot in P
 for name in list(P):
  if name not in sharedkeys and name!=pivot:drop(name)
 futureglobals={'Q','P7','E','S','Delta',f'A52_{j}',shared[31],shared[83],shared[104],'shared_G','shared_B51','shared_P7F9',pivot}
 for name in list(P):
  if name not in futureglobals:drop(name)
 sharedkeys=set(P);phases.append({'phase':'shared_pivot35','N8_35_56_literal_node':pivot,'computed_once_SOURCE_only':True,'scalar_function_actual_coefficients_NOT_evaluated':True})
 outputs=[]
 for row in scope['selected15_row_ids']:
  beginops=len(ops);N8=stage8(row)
  for name in list(P):
   if name not in sharedkeys and name!=N8:drop(name)
  v=values[row];f=cp.parent/v['file'];assert sha(f)==v['sha256'];inputs[f]=sha(f);b=f.read_bytes();p=decode(b[9:])if b[8]else{'content':'0','shift':[0]*5,'terms':[]};z=operand(p,{'kind':'actual_existing_N8_i82','row':row,'path':str(f.relative_to(ROOT)),'sha256':v['sha256']});zf=out/'rows'/f'row-{row}-existing-col82.json';zf.write_bytes(enc(z));leaf('Z'+str(row),z,{'scope':str(zf.relative_to(out)),'entire_operand':True});del p,z,b
  A=f'row{row}_ninth_A';B=f'row{row}_ninth_B';C=f'row{row}_ninth_output';mul(A,'shared_P7F9',N8);drop(N8);mul(B,'Z'+str(row),pivot,Q(-62208000000),(0,0,1,1,2));drop('Z'+str(row));add(C,A,B,consume=True);drop(A,B);outputs.append({'row':row,'node':C,'phase_begin_op':beginops,'phase_end_op':len(ops)-1});phases.append({'phase':'row_output','row':row,'output':C,'source_support':len(P[C]['S']),'lastuse_all_row_leaves_then_output_persist':True})
  for name in list(P):
   if name not in sharedkeys:drop(name)
  (out/'SOURCE-prefix.json').write_bytes(enc({'all_passed':True,'completed_rows':[z['row']for z in outputs],'no_matrix_coefficient_values':True}));gc.collect();assert time.monotonic()-start<55
 # Exact abstract last-use schedule from the frozen expression program.
 uses={name:0 for name in nodes}
 for op in ops:
  for k in('A','B'):
   if k in op:uses[op[k]]+=1
 for o in outputs:uses[o['node']]+=1
 live={};valuepeak=50000;schedule=[];rawcaps=0;encodepeaks=[];maxbits=0;outset={z['node']:z for z in outputs}
 for i,op in enumerate(ops):
  dest=op['output'];n=nodes[dest]['terms_upper'];reuse=op['op']=='add'and uses[op['A']]==1
  if op.get('SOURCE_consumed_both_lastuse_support_buffers'):assert uses[op['A']]==uses[op['B']]==1
  current=sum(live.values());growth=n-nodes[op['A']]['terms_upper']if reuse else n;assert growth>=0;peak=50000+7*(current+growth)+1024;valuepeak=max(valuepeak,peak);live[dest]=n
  if uses[dest]==0:live.pop(dest,None)
  for k in('A','B'):
   if k in op:
    uses[op[k]]-=1
    if uses[op[k]]==0:live.pop(op[k],None)
  schedule.append({'op':i,'output':dest,'reuse_left_lastuse_buffer':reuse,'whole_coefficient_owned_upper':peak,'resident_terms_after':sum(live.values())})
  if dest in outset:
   q=Q(nodes[dest]['known_QQ_content']);M=int(nodes[dest]['primitive_integer_abs_upper'])*abs(q.numerator);maxbits=max(maxbits,M.bit_length(),q.denominator.bit_length());box=nodes[dest]['exponent_box5'];eb=sum(max(1,((hi-lo).bit_length()+6)//7)for lo,hi in box);rawcap=2048+(q.denominator.bit_length()+7)//8+n*((M.bit_length()+7)//8+eb+16);rawcaps+=rawcap
   # Planned SAME primitive bytes, with in-place heapsort of alias key refs
   # and incremental canonical varint write/readback; no value/decode array.
   other=sum(live.values())-n;ep=50000+7*other+10*n+2048;encodepeaks.append({'row':outset[dest]['row'],'output_terms_upper':n,'planned_streaming_binary_owned_upper':ep,'binary_RAW_upper':rawcap});valuepeak=max(valuepeak,ep);uses[dest]-=1;live.pop(dest,None)
 program={'all_passed':True,'source_only':True,'parameter_order':['H','L','s0','q0','t0'],'common_source_operand_file':'common-source-operands.json','row_bindings16':rowbindings,'literal_ops':ops,'positive_support_and_height_profiles':nodes,'literal_shared_N8_35_56_node':pivot,'row15_output_nodes':outputs,'phases':phases,'source_no_coefficient_arithmetic':True,'lastuse_schedule':schedule,'streaming_binary_phase_forecasts':encodepeaks,'known_scalar_content_gcd_only_integer_rational_NO_polynomial_GCD':True,'full168x15_87plus66_RHS_other14_oneP290_complements_kept':True}
 (out/'staged-literal-col56-program-and-positive-profiles.json').write_bytes(enc(program));baseRAW=sum(p.stat().st_size for p in inputs);sourceRAW=sum(p.stat().st_size for p in out.rglob('*')if p.is_file());futureRAW=baseRAW+sourceRAW+rawcaps+2*1024*1024;safe=valuepeak<150000 and futureRAW<32*1024*1024
 report={'id':ID,'all_passed':True,'source_only':True,'all16_staged_SOURCE_rows':16,'actual_new_matrix_coefficients':0,'literal_ops':len(ops),'staged_pair_majorant':tupleops,'direct_unroll_majorant_NOT_run':995410878,'SOURCE_owned_peak':sourcepeak,'future_coefficient_whole_owned_upper':valuepeak,'future_RAW_upper':futureRAW,'future_fixed_metadata_allowance':50000,'future_streaming_encoder_source_safe_under150k_32MiB':safe,'height_integer_bits_upper':maxbits,'new_source_independent_gate_and_custody_REQUIRED':True,'no_coefficient_worker_GO':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',action='store_true');main(p.parse_args())
