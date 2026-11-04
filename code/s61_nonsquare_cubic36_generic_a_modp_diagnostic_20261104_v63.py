#!/usr/bin/env python3
"""One bounded F32003(a)[e,d,z] diagnostic of the safe36 cubic coker."""
import argparse, hashlib, json, math, os, resource, signal, subprocess, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_CUBIC36_GENERIC_A_MODP_DIAGNOSTIC_20261104_v63'
P=32003
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def size_tree(p):return sum(f.stat().st_size for f in p.rglob('*')if f.is_file())
def child_limit():resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
def run_one(wrapper,script,log,deadline,out):
 with log.open('w')as f:
  proc=subprocess.Popen([str(wrapper),'--cpus=0','-q',str(script)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=child_limit)
  reason=None
  while proc.poll()is None:
   if time.monotonic()>=deadline:reason='shared_total_60s_deadline';break
   if size_tree(out)>32*1024*1024:reason='combined_artifact_32MiB_cap';break
   time.sleep(.1)
  if reason:
   os.killpg(proc.pid,signal.SIGKILL);proc.wait()
 return {'returncode':proc.returncode,'stop_reason':reason,'pid':proc.pid}
def main(out):
 started=time.monotonic();deadline=started+60;out.mkdir(parents=True,exist_ok=False)
 pre=ROOT/'reports/S61_NONSQUARE_CUBIC36_FUNCTION_FIELD_PREFLIGHT_20261104_v63/report.json';pr=json.loads(pre.read_text());assert pr['all_read_only_format_provenance_and_statistics_checks_passed']
 for rel,h in pr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 source=ROOT/'work/nonsquare-v62/cubic27-v6/after21';wrapper=ROOT/'work/nonsquare-v62/singular-runtime/run-singular-optimized-v6.sh'
 sanity=out/'field-sanity.sing';sanity.write_text('''ring R=(32003,a),(e,d,z),dp;
poly f=(1/(a+1))*(a+1)-1;
if(f!=0){ERROR("FIELD_RATIONAL_ARITHMETIC_FAILED");quit;}
matrix T[2][3]=a*e,1-a*e,0,0,0,a;
matrix X[3][2]=1,0,1,0,0,1/a;
matrix E[2][2]=1,0,0,1;
if(T*X!=E){ERROR("FIELD_EXPLICIT_MODULE_IDENTITY_FAILED");quit;}
module I=T;module G=std(I);module EE=E;
if(size(reduce(EE,G))!=0){ERROR("FIELD_STD_UNIT_SANITY_FAILED");quit;}
if(vdim(G)!=0){ERROR("FIELD_STD_VDIM_SANITY_FAILED");quit;}
print("EXACT_FIELD_AND_STD_SANITY_PASS");quit;
''')
 result={'id':ID,'ring':'F32003(a)[e,d,z]','module_shape':[36,106],'column_map':'R^106 -> R^36','source_sha256':sha(source),'hard_shared_seconds':60,'child_address_space_limit_bytes':2000000000,'combined_artifact_limit_bytes':32*1024*1024,'no_J_inversion':True,'QQ_promotion':False}
 result['sanity']=run_one(wrapper,sanity,out/'field-sanity.log',deadline,out)
 sanitylog=(out/'field-sanity.log').read_text()
 sanitypass=result['sanity']['returncode']==0 and'EXACT_FIELD_AND_STD_SANITY_PASS'in sanitylog and'?'not in sanitylog
 result['exact_field_std_sanity_passed']=sanitypass
 if not sanitypass:
  result['frontier']='sanity_failed_or_capped_NO_MAIN_GB_LAUNCHED'
 else:
  script=out/'input.sing';inverse={};kept=0;vanished=0;entries=0
  with source.open()as f,script.open('w')as g:
   assert next(f).split()==['36','70','21'];next(f);next(f)
   g.write('ring R=(32003,a),(e,d,z),dp;\nmatrix M[36][106];\n')
   for i in range(36):
    for j in range(70):
     t=next(f).split();l=int(t[0]);assert len(t)==1+4*l;groups={}
     for k in range(l):
      q=t[1+4*k];a,e,d=map(int,t[2+4*k:5+4*k]);num,den=map(int,q.split('/'))if'/'in q else(int(q),1)
      assert den%P
      if den not in inverse:inverse[den]=pow(den,-1,P)
      c=num*inverse[den]%P
      if not c:vanished+=1;continue
      kept+=1;groups.setdefault((e,d),[]).append((a,c))
     terms=[]
     for(e,d),cofs in sorted(groups.items(),reverse=True):
      ct=[]
      for a,c in cofs:
       c=c if c<=P//2 else c-P
       ct.append(str(c)+('*a'if a==1 else '*a^'+str(a)if a else''))
      term='('+ '+'.join(ct)+')'
      if e:term+='*e'+('^'+str(e)if e!=1 else'')
      if d:term+='*d'+('^'+str(d)if d!=1 else'')
      terms.append(term)
     if terms:g.write('M['+str(i+1)+','+str(j+1)+']='+'+'.join(terms)+';\n');entries+=1
   assert not f.read().strip()
   g.write('''poly F=500*a+90*e^2+3*e-24;
poly B=50*a*d+18*d*e^2-9*d*e-2*d+60*e-40;
poly C=100*a+18*e^2-9*e-2;
poly rel=z*d*F*B*C-1;
for(int i=1;i<=36;i++){M[i,70+i]=rel;}
module I=M;
print("EXACT_INPUT36_BY106_READY_NO_J");
module G=std(I);
print("MAIN_STD_COMPLETED");
print("PROFILE_NROWS="+string(nrows(G)));
print("PROFILE_NCOLS="+string(ncols(G)));
print("PROFILE_DIM="+string(dim(G)));
print("PROFILE_VDIM="+string(vdim(G)));
int gt=0;for(int j=1;j<=ncols(G);j++){gt=gt+size(G[j]);}
print("PROFILE_TERMS="+string(gt));
matrix E[36][36];for(int i=1;i<=36;i++){E[i,i]=1;}
module EE=E;
print("PROFILE_QUOTIENT_UNIT="+string(size(reduce(EE,G))==0));
string s=string(G);
print("PROFILE_STD_STRING_BYTES="+string(size(s)));
if(size(s)<8000000){write("'''+str(out/'std-basis.sing')+'''","ring R=(32003,a),(e,d,z),dp; module G="+s+";");}
module L=lead(G);string ls=string(L);
print("PROFILE_LEAD_STRING_BYTES="+string(size(ls)));
if(size(ls)<1000000){write("'''+str(out/'leading-module.sing')+'''","ring R=(32003,a),(e,d,z),dp; module L="+ls+";");}
if(vdim(G)>=0 && vdim(G)<20000){module K=kbase(G);string ks=string(K);print("PROFILE_KBASE_COLUMNS="+string(ncols(K)));if(size(ks)<1000000){write("'''+str(out/'standard-module-basis.sing')+'''","ring R=(32003,a),(e,d,z),dp; module K="+ks+";");}}
print("MAIN_PROFILE_COMPLETE");quit;
''')
  result['literal_reduction']={'source_terms':kept+vanished,'nonzero_reduced_terms':kept,'zero_reduced_terms':vanished,'source_nonzero_entries_written':entries,'all_source36x70_entries_preserved':True,'relation_columns':36,'source_core_order':['a','e','d'],'a_remains_coefficient_field_parameter':True,'source_denominators_checked_32003_units':True,'input_sha256':sha(script),'input_bytes':script.stat().st_size}
  assert kept+vanished==421667 and vanished==14 and entries==1631
  assert size_tree(out)<32*1024*1024
  if time.monotonic()>=deadline:result['frontier']='shared_deadline_during_preparation_MAIN_GB_NOT_LAUNCHED'
  else:
   result['main']=run_one(wrapper,script,out/'main.log',deadline,out)
   text=(out/'main.log').read_text();result['input_ready_marker_observed']='EXACT_INPUT36_BY106_READY_NO_J'in text;result['std_completed_marker_observed']='MAIN_STD_COMPLETED'in text;result['profile_completed_marker_observed']='MAIN_PROFILE_COMPLETE'in text
   result['profile_markers']=[line for line in text.splitlines()if line.startswith('PROFILE_')]
   result['frontier']='PROFILE_COMPLETED_GF_ONLY'if result['profile_completed_marker_observed']else(result['main']['stop_reason']or'MAIN_RETURNED_WITHOUT_COMPLETE_PROFILE')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),pre,source,wrapper]
 result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['artifacts']={str(p.relative_to(out)):{'sha256':sha(p),'bytes':p.stat().st_size}for p in out.rglob('*')if p.is_file()};result['combined_artifact_bytes']=size_tree(out);result['seconds']=time.monotonic()-started;result['no_auto_retry_or_extension']=True;result['frozen']=True
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','artifacts']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);args=p.parse_args();main(args.output.resolve())
