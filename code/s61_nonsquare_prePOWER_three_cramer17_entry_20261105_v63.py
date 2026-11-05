#!/usr/bin/env python3
"""Immutable preparation for exactly THREE Cramer17 numerators; separate execute."""
import argparse,gc,hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def cap(b):
 z=b['box'];lo,hi=b['degree'];n=sum(max(0,min(z[2][1],hi-a-e)-max(z[2][0],lo-a-e)+1)for a in range(z[0][0],z[0][1]+1)for e in range(z[1][0],z[1][1]+1));rows=(z[0][1]-z[0][0]+1)*(z[1][1]-z[1][0]+1)+1
 return n+rows+17
def prepare(out):
 out.mkdir(parents=True,exist_ok=False);pins={};phasecells=[]
 def load(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 def subtree(rel,key):
  p=ROOT/rel;pins[rel]=sha(p);text=p.read_text();tag=json.dumps(key)+':';assert text.count(tag)==1;i=text.index(tag)+len(tag)
  while text[i].isspace():i+=1
  return json.JSONDecoder().raw_decode(text,i)[0]
 src='reports/S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_SOURCE_20261105_v63/'
 r=load(src+'report.json');assert r['all_passed'];schedule=load(src+'three-cramer-schedule.json')
 hr='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json';raw=subtree(hr,'raw_compact14x5');rho=subtree(hr,'known_row_units_rho14')
 # Literal arrays are retained only for the six selected entries.
 A=[];sourceproof=[]
 for i,rr in enumerate([1,3,4]):
  assert rho[rr]['a_d_F_B_C_J_powers']==[1,0,1,0,2,0];c=Q(rho[rr]['constant']);A.append([[(m,str(Q(v)*c))for m,v in raw[rr][j]]for j in[0,4]])
  sourceproof.append(dict(current_row=[2,4,5][i],compact_row=rr,rho=rho[rr],H27H31_QQ_scale=str(c),known_unit_carry='rho=c*aFC²; column27 factor1/(aFC²), column31 factor1/(aFC² gamma); gamma is exactly original H31 column unit'))
 del raw,rho;gc.collect()
 br='reports/S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63/native.json';bot=subtree(br,'bottom4x6')
 for i in range(3):A[i].append(bot[i][0])
 del bot;gc.collect();V=[]
 for i in range(3):
  x=load(f'reports/S61_NONSQUARE_PREPOWER_KAPPA17_ENTRY_20261105_v63/source-root{i}.json');V.append([[m,str(Q(int(v),int(x['denominator'])))]for m,v in x['integer_coefficients']])
 del x;gc.collect();assert sum(len(p)for row in A for p in row)+sum(map(len,V))==1647
 N=[]
 def atom(poly,bound):
  D=int(bound['D']);terms=[]
  for m,c in poly:
   z=Q(c)*D;assert z.denominator==1 and Q(int(z),D)==Q(c);terms.append((m,str(int(z))))
  i=len(N);N.append(dict(id=i,op='i',args=[-1,-1],bound=bound,terms=terms));return i
 ar=[[atom(A[i][j],schedule['A3_bounds'][i][j])for j in range(3)]for i in range(3)];vr=[atom(V[i],schedule['v17_bounds'][i])for i in range(3)]
 def node(op,a,b,bound):i=len(N);N.append(dict(id=i,op=op,args=[a,b],bound=bound));return i
 roots=[];by={(x['output'],x['term']):x for x in schedule['steps']}
 for j in range(3):
  acc=None
  for i in range(3):
   s=by[j,i];rr=s['minor_rows'];cc=s['minor_columns'];a=node('m',ar[rr[0]][cc[0]],ar[rr[1]][cc[1]],s['first_product']);b=node('m',ar[rr[0]][cc[1]],ar[rr[1]][cc[0]],s['second_product'])
   cof=node('s',a,b,s['cofactor'])if s['cofactor_sign']==1 else node('s',b,a,s['cofactor']);t=node('m',cof,vr[i],s['term_product']);acc=t if acc is None else node('a',acc,t,s['accumulator'])
  roots.append(acc)
 assert len(N)==54
 users=[0]*54
 for n in N:
  if n['op']!='i':
   for k in n['args']:users[k]+=1
 for k in roots:users[k]+=1
 live=54*20+sum(2*len(n.get('terms',[]))+16 for n in N if n['op']=='i');peak=live+3000;events=[]
 for n in N:
  c=cap(n['bound']);live+=c;peak=max(peak,live+3000)
  if n['op']=='m':peak=max(peak,live+c+3000)
  if n['op']=='i':live-=2*len(n['terms'])+16
  else:
   for k in n['args']:
    users[k]-=1
    if not users[k]:live-=cap(N[k]['bound'])
  events.append(dict(node=n['id'],live_including_numeric_metadata=live+3000))
 assert peak<=73322
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  f.write('CRAMER17D001 54 48\n')
  for n in N:
   b=n['bound'];spec=[v for pair in b['box']for v in pair]+b['degree']+[str(b['D'])];f.write(f"{n['id']} {n['op']} {n['args'][0]} {n['args'][1]}\n"+' '.join(map(str,spec))+'\n')
   if n['op']=='i':
    f.write(str(len(n['terms']))+'\n')
    for m,c in n['terms']:f.write(' '.join(map(str,m))+' '+c+'\n')
  f.write(' '.join(map(str,roots))+'\n')
 dump(out/'literal-source-custody.json',dict(source_normalizations=sourceproof,input_A3_and_v17=dict(A3=A,v17=V),all_input_constant_scales_exact=True,physical_row_order=[2,4,5],physical_columns=[27,31,22],target=17,root_order=roots,unit_identity='detA3=-kappa3; inversekappa3=detA8/(cprod*Saux*T)',kappa17_peer_REQUIRED_before_execution=True))
 dump(out/'native-dense-ownership.json',dict(peak_all_native_numeric_metadata=peak,events=events,source_JSON_preparation='decode only H5 rawmatrix/rho and bottom24 subtree; never whole H5 certificate object',preparation_numeric_bound=40000,readback_numeric_bound=4*6765+15000,owned_cap=73322,integer_temporary_bits=434))
 cpp=ROOT/'code/s61_nonsquare_three_cramer17_dense_20261105_v63.cpp';pins[str(cpp.relative_to(ROOT))]=sha(cpp);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 gate=dict(id=ID+'_WORKER_GATE',all_passed=True,frozen=True,mathematical_job_NOT_run=True,rootcheck_and_kappa17_peer_REQUIRED_before_run=True,
  input_sha256=pins,native_input_sha256=sha(inp),native_nodes=54,source_polynomial_terms=1647,native_all_owned_upper=peak,ALL_phase_owned_upper=max(peak,42060),owned_cap=73322,temporary_bits434=True,TOTAL_RAW_cap=5820000,
  source_root_and_cofactor_signs_fixed=True,no_new_chart_inverse_or_other_coordinate=True,evidence_sha256={p.name:sha(p)for p in out.iterdir()})
 dump(out/'report.json',gate);print(json.dumps({k:v for k,v in gate.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);prepare(p.parse_args().output.resolve())
