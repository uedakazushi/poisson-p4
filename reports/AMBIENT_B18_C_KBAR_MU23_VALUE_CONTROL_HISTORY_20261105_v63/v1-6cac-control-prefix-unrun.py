#!/usr/bin/env python3
"""ONE cache-only37+14 dot job. Requires fresh source peer and parent GO."""
import argparse,hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_MU23_CACHED_VALUES_20261105_v63'
SOURCE='reports/AMBIENT_B18_C_KBAR_MU23_CACHE_SOURCE_20261105_v63'
COMPILED=SOURCE+'/complete-old-Pg-second-cache-source-and-fourteen-action-positive-gate.json'
COMPILED_SHA='317b2c5e491e7f8369e73a7578982731415ac45d77444219308dcc5cc34b0127'
STORED=0;PRE=0;COMMITTED=136;TEMPORARY=273
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def h(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def check(q):
 global STORED
 STORED=max(STORED,h(q));assert h(q)<=COMMITTED;return q
def parse(s):
 ss=str(s).split('/');assert len(ss)<=2
 n=int(ss[0]);d=int(ss[1])if len(ss)==2 else 1;assert d
 assert max(abs(n).bit_length(),abs(d).bit_length())<=COMMITTED
 return check(Q(n,d))
def pre(z):
 global PRE
 PRE=max(PRE,z);assert z<=TEMPORARY and z<=1153
def add(a,b):pre(2*max(h(a),h(b))+1);return check(a+b)
def mul(a,b):pre(h(a)+h(b));return check(a*b)
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):
  a=x.split('/');return int(1<=len(a)<=2 and all(c.lstrip('-').isdigit()for c in a))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.keys()))+sum(map(count,x.values()))
 return 0
def dot(weights,values,terms):
 z=Q(0)
 for root,c in weights:
  w=parse(c);v=values[int(root)];p=mul(w,v);z=add(z,p);terms.append([int(root),str(w),str(v),str(p)])
 return z
def main(peerpath):
 global COMMITTED,TEMPORARY
 started=time.monotonic();signal.alarm(60);assert (ROOT/COMPILED).exists()
 assert sha(ROOT/COMPILED)==COMPILED_SHA
 sourcepath=ROOT/SOURCE/'report.json';sr=json.loads(sourcepath.read_text());assert sr['all_passed']
 assert sr['evidence_sha256'][Path(COMPILED).name]==COMPILED_SHA
 peerpath=peerpath.resolve();peer=json.loads(peerpath.read_text());assert peer['all_passed']
 pins=peer['input_sha256'];assert pins[COMPILED]==COMPILED_SHA
 D=json.loads((ROOT/COMPILED).read_text());assert D['all21_regular_second_opaque_gradients_empty']
 COMMITTED=D['positive_whole_future_height']['committed'];TEMPORARY=D['positive_whole_future_height']['temporary']
 assert (COMMITTED,TEMPORARY)==(136,273)
 profile=D['future_numeric_profile'];static=profile['static_owned_upper'];assert static==14076 and static<60000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);saved=0;readbytes=(ROOT/COMPILED).stat().st_size+sourcepath.stat().st_size+peerpath.stat().st_size;written=0
 def save(name,data,live=()):
  nonlocal saved,written
  owned=count(D)+count(sr)+count(peer)+saved+count(data)+count(live)+2500;assert owned<=static and owned<60000 and owned<150000
  text=json.dumps(data,indent=2)+'\n';written+=len(text.encode());assert readbytes+written+20000<=800000 and readbytes+written<33554432
  (out/name).write_text(text);saved+=count(data)
  (out/'phase.json').write_text(json.dumps({'phase':name,'seconds':time.monotonic()-started,'owned':owned,'stored_bits':STORED,'pre_bits':PRE,'bytes_withreserve':readbytes+written+20000})+'\n')
 C={i:{int(k):parse(c)for k,c in row}for i,row in D['all21_given_scalar_cache_rows_not_selected_rr']};assert len(C)==21
 mon=D['actual_cache_monomials_original8'];assert mon[19]==[7,7]
 ri=mon.index([7]);assert ri==5
 # Mandatory complete origin/first readback BEFORE selecting requested rr.
 lower=[]
 for j,weights in enumerate(D['all37_Pg_source_weights']):
  terms0=[];termsr=[]
  z0=dot(weights,{i:p.get(0,Q(0))for i,p in C.items()},terms0);zr=dot(weights,{i:p.get(ri,Q(0))for i,p in C.items()},termsr)
  assert not z0 and not zr;lower.append([j,'0','0'])
 rr={i:p.get(19,Q(0))for i,p in C.items()};tvals=[];terms37=[]
 for weights in D['all37_Pg_source_weights']:
  terms=[];tvals.append(dot(weights,rr,terms));terms37.append(terms)
 save('phase1-all37-actual-Pg-second-values-and-lower-readbacks.json',{'all37_origin_and_r_readback_zero':lower,'all37_actual_t_head4r2':list(map(str,tvals)),'all37_termwise_cache_products':terms37,'no_source_derivative_solve':True},live=[C,rr,tvals])
 actual=[];terms14=[];direct_terms=[]
 for weights,direct in zip(D['all14_same37_weights'],D['all14_direct_Pg_source_weights']):
  terms=[];a=Q(0)
  for j,c in enumerate(weights):
   if parse(c):p=mul(parse(c),tvals[j]);a=add(a,p);terms.append([j,c,str(tvals[j]),str(p)])
  directterms=[];b=dot(direct,rr,directterms);assert a==b
  actual.append(a);terms14.append(terms);direct_terms.append(directterms)
 witness=next((i for i,x in enumerate(actual)if x),None);mu=2 if witness is not None else 3
 lam=parse(D['known_lambdaS']);minor=mul(mul(lam,lam),actual[witness])if witness is not None else Q(0)
 save('phase2-all14-actual-actions-two-compositions-and-generator-count.json',{'all14_z_action_labels':D['all14_z_action_labels'],'all14_actual_q2_mod_rKbar':list(map(str,actual)),'all14_termwise_M2_compositions':terms14,'all14_termwise_direct21root_compositions':direct_terms,'all14_two_compositions_agree':True,'actual_mu_R0_K':mu,'fixed_nonzero_action_witness':witness,'witness_action_label':D['all14_z_action_labels'][witness]if witness is not None else None,'witness_prolonged_minor_lambdaS_squared_times_action':str(minor),'scope':'Minimal BASE R0-module generatorcount for ORIGINALK; DmuK7unchanged. ActualfullKliftsEXIST,notconstructed. NooriginalH0/nil inference.'},live=[C,tvals,actual])
 ep={p.name:sha(p)for p in out.glob('phase[12]-*.json')};report={'id':STEM,'all_passed':True,'actual_mu_R0_K':mu,'full_D_module_generator_count_unchanged':7,'all37_rawM2_and14_action_values_complete':True,'witness_index':witness,'full_K_lifts_constructed':False,'original_H0_or_nil_changed':False,'input_sha256':{COMPILED:COMPILED_SHA,str(sourcepath.relative_to(ROOT)):sha(sourcepath),str(peerpath.relative_to(ROOT)):sha(peerpath)},'code_sha256':sha(Path(__file__)),'evidence_sha256':ep,'static_numeric_bound':static,'maximum_stored_bits':STORED,'maximum_preoperation_bits':PRE,'seconds':time.monotonic()-started}
 save('report.json',report);print(json.dumps(report));signal.alarm(0)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-peer-report',required=True,type=Path);a=p.parse_args();main(a.source_peer_report)
