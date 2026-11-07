#!/usr/bin/env python3
"""ONE future two-dot job on previously validated37 t; no root recomposition."""
import argparse,hashlib,json,signal,time,traceback
from pathlib import Path
from fractions import Fraction as Q
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_GAMMA12_CACHED_VALUES_20261105_v63'
GATE='reports/AMBIENT_B18_C_KBAR_GAMMA12_VALUE_GATE_20261105_v63/two-gamma-positive-source-and-validated37-input-only-gate.json'
GATE_SHA='18a9e953eaaf8f7966b2ae4455a60644214996e15fef3ba2a9f538318f83abfc'
MAXH=0;MAXPRE=0;OUT=None;READ=0;WRITTEN=0
SCI={'phase':'initialized','row':None,'term':None,'pending':None,'parsed_t37':[],'values':[],'terms':[],'accumulator':'0'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def h(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):
  ss=x.split('/');return int(len(ss)<=2 and all(s.lstrip('-').isdigit()for s in ss))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.keys()))+sum(map(count,x.values()))
 return 0
def mark(op,**kw):SCI['pending']={'operation':op,'completed':False,**kw}
def checked(q):
 global MAXH
 SCI['pending']['computed_candidate']=str(q);MAXH=max(MAXH,h(q));assert h(q)<=136
 SCI['pending'].update(completed=True,result=str(q));return q
def pre(n):
 global MAXPRE
 MAXPRE=max(MAXPRE,n);assert n<=273
def parse(s):
 mark('parse',token=str(s));ss=str(s).split('/');assert len(ss)<=2
 n=int(ss[0]);d=int(ss[1])if len(ss)==2 else 1;assert d
 assert max(abs(n).bit_length(),abs(d).bit_length())<=136
 return checked(Q(n,d))
def add(a,b):mark('add',operands=[str(a),str(b)]);pre(2*max(h(a),h(b))+1);return checked(a+b)
def mul(a,b):mark('multiply',operands=[str(a),str(b)]);pre(h(a)+h(b));return checked(a*b)
def timeout_stop(signum,frame):raise TimeoutError('wall60 STOP')
def preserve(exc):
 global OUT
 signal.alarm(0)
 if OUT is None:OUT=ROOT/'reports'/(STEM+'_PREOUTPUT_FAILURE');OUT.mkdir(exist_ok=False)
 assert count(SCI)<=1000
 data={'status':'STOP_COMPLETE_SCIENTIFIC_PREFIX','scientific_state':SCI,'exception_type':type(exc).__name__,'exception':str(exc),'traceback':traceback.format_exc(),'max_stored':MAXH,'max_pre':MAXPRE,'gate_sha256':GATE_SHA,'code_sha256':sha(Path(__file__)),'no_automatic_replay_or_capraise':True}
 txt=json.dumps(data,indent=2,default=str)+'\n';assert READ+WRITTEN+len(txt.encode())+20000<400000
 (OUT/'STOP-scientific-prefix.json').write_text(txt);(OUT/'phase.json').write_text(json.dumps({'phase':'STOP','scientific_prefix_saved':True})+'\n')
def main(peerpath):
 global OUT,READ,WRITTEN
 start=time.monotonic();signal.signal(signal.SIGALRM,timeout_stop);signal.alarm(60);SCI['phase']='validation'
 p=ROOT/GATE;assert sha(p)==GATE_SHA;gate=json.loads(p.read_text())
 peerpath=peerpath.resolve();assert peerpath.stat().st_size<=60000;peer=json.loads(peerpath.read_text());assert peer['all_passed']and count(peer)<=1000
 assert peer['input_sha256'][GATE]==GATE_SHA
 src=gate['validated_t_foreign_reference'];tp=ROOT/src['path'];assert sha(tp)==src['sha256']=='6dce213c1da5bf7cf2a957d6dd7d7e7d5e19b3f090fdd574fb466b18865b755e'
 ip=ROOT/src['ROOT_integer_postpeer']['path'];assert sha(ip)==src['ROOT_integer_postpeer']['sha256']=='57878946b9c677f4ed12ef9ae00e61fc9979f70f4c57db6acda45d8101163c72'
 integer=json.loads(ip.read_text());assert integer['all_passed']and integer['all37_rr_values_match']
 assert integer['input_sha256'][src['path']]==src['sha256']
 # Read existing37 values, without recomputing ANY earlier C0/root dot.
 original=json.loads(tp.read_text());given=original[src['field']];assert len(given)==37;del original
 READ=sum(x.stat().st_size for x in[p,peerpath,tp,ip]);OUT=ROOT/'reports'/STEM;OUT.mkdir(exist_ok=False)
 saved=0
 def save(name,data):
  nonlocal saved
  global WRITTEN
  owned=count(gate)+count(peer)+count(integer)+count(given)+count(SCI)+saved+count(data)+2000;assert owned<=8000
  txt=json.dumps(data,indent=2)+'\n';WRITTEN+=len(txt.encode());assert READ+WRITTEN+40000<400000
  (OUT/name).write_text(txt);saved+=count(data)
  (OUT/'phase.json').write_text(json.dumps({'phase':name,'owned':owned,'seconds':time.monotonic()-start,'max_stored':MAXH,'max_pre':MAXPRE})+'\n')
 SCI['phase']='parse_validated_t37'
 for j,s in enumerate(given):SCI['term']=j;SCI['parsed_t37'].append(parse(s))
 save('phase0-validated-existing37-input-reference.json',{'validated37path':src['path'],'sha256':src['sha256'],'independent_integer_peer_sha256':src['ROOT_integer_postpeer']['sha256'],'all37_existing_inputs_parsed':True,'earlier_roots_or37_dots_recomputed':False})
 SCI['phase']='only_two_gamma_dots'
 for row,w in enumerate(gate['gamma_weights2_by37']):
  assert len(w)==37;SCI['row']=row;terms=[];SCI['terms'].append(terms);z=Q(0);SCI['accumulator']='0'
  for j,c in enumerate(w):
   SCI['term']=j;a=parse(c);v=SCI['parsed_t37'][j];value=mul(a,v);terms.append([j,str(a),str(v),str(value)]);z=add(z,value);SCI['accumulator']=str(z)
  SCI['values'].append(z)
 save('two-actual-filtered-gamma-values-and-all74-products.json',{'gamma1_gamma2':list(map(str,SCI['values'])),'all74_termwise_products':SCI['terms'],'same37_inputs_reused_without_root_or_tdot_replay':True,'C_F_Delta_H_not_computed':True,'scope':'Twooldfilteredq2coefficientsONLY; no H0/nil/typechoice.'})
 ep=OUT/'two-actual-filtered-gamma-values-and-all74-products.json'
 report={'id':STEM,'all_passed':True,'actual_gamma1_gamma2':list(map(str,SCI['values'])),'only_two_gamma_dots_evaluated':True,'input_sha256':{GATE:GATE_SHA,str(peerpath.relative_to(ROOT)):sha(peerpath),src['path']:src['sha256'],src['ROOT_integer_postpeer']['path']:src['ROOT_integer_postpeer']['sha256']},'code_sha256':sha(Path(__file__)),'evidence_sha256':{ep.name:sha(ep)},'max_stored_H':MAXH,'max_pre_H':MAXPRE,'static_owned':8000,'seconds':time.monotonic()-start,'C_F_Delta_H_or_H0_nil_claimed':False}
 save('report.json',report);print(json.dumps(report));signal.alarm(0)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-peer-report',required=True,type=Path);a=p.parse_args()
 try:main(a.source_peer_report)
 except BaseException as exc:preserve(exc);raise
