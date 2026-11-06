#!/usr/bin/env python3
"""G: source-bound native support diagnostic at ONE frozen rational beta.

prepare reconstructs all134 QQ equations independently from the original61
encoding, then optionally verifies/adds an exact G6 consequence and adds the CORRECT globalJ=[0,1,2] eta=n0+n34 chart. run invokes
Singular std or slimgb once and, only after detecting1, attempts lift within the SAME
nativeCPUbudget. A modular fixedfiber unit never implies genericchar0 unit.
Diagnostic inputs/logs belong in ignoredwork; any verified mathematical
cofactor certificate is exported as plainJSON to the explicitly given path.
"""
from __future__ import annotations
import argparse,collections,gzip,hashlib,itertools,json,os,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
MODEL='legacy/v02/data/symplectic_slice_reduced_equations.json'
POINT='reports/S61_KNOWN_FIBER_AUDIT_20261031.json'
FRAME='reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz'
P=32003

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):
 raw=p.read_bytes();return json.loads(gzip.decompress(raw)if p.suffix=='.gz'else raw)
def save(p,obj):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x')as f:json.dump(obj,f,ensure_ascii=False,indent=2);f.write('\n')
def put(poly,word,c):
 c=poly.get(word,F(0))+c
 if c:poly[word]=c
 else:poly.pop(word,None)
def encode(poly):return [[list(m),str(c)]for m,c in sorted(poly.items())]
def decode(poly):return {tuple(m):F(c)for m,c in poly}
def expression(poly,field):
 terms=[]
 for mon,c in sorted(poly.items()):
  c=F(c)
  if field=='gf':
   assert c.denominator%P;c=c.numerator*pow(c.denominator,-1,P)%P
   if not c:continue
   text=str(c)
  else:text='('+str(c.numerator)+'/'+str(c.denominator)+')'
  terms.append(text+''.join('*'+('u'if i==45 else'x'+str(i))for i in mon))
 return '+'.join(terms)if terms else'0'

def singular_script(polys,field,algorithm):
 head='ring R='+('32003'if field=='gf'else'0')+',('+','.join(['x'+str(i)for i in range(45)]+['u'])+'),dp;\n'
 head+='ideal inputIdeal=\n'+',\n'.join(expression(p,field)for p in polys)+';\n'
 return head+'''print("__NATIVE_STD_START__");
ideal answerBasis=__ALGORITHM__(inputIdeal);
print("__NATIVE_STD_FINISHED__");
print("__BASIS_COUNT__:"+string(size(answerBasis)));
if(size(answerBasis)==1 && answerBasis[1]==1)
{
 print("__UNIT_DETECTED__");
 print("__LIFT_START__");
 matrix liftmat=lift(inputIdeal,ideal(1));
 matrix liftcheck=matrix(inputIdeal)*liftmat;
 if(liftcheck[1,1]==1)
 {
  print("__NATIVE_LIFT_VERIFIED__");
  int cidx; poly zcoef; intvec zexp;
  for(cidx=1;cidx<=__GENCOUNT__;cidx++)
  {
   print("__COFACTOR__:"+string(cidx-1));
   zcoef=liftmat[cidx,1];
   while(zcoef!=0)
   {
    zexp=leadexp(zcoef);
    print("__TERM__:"+string(leadcoef(zcoef))+":"+string(zexp));
    zcoef=zcoef-lead(zcoef);
   }
  }
  print("__COFACTORS_FINISHED__");
 }
 else {print("__LIFT_IDENTITY_FAILED__");}
}
else
{
 print("__NONUNIT_BASIS_RETURNED__");
 print("__BASIS_DIM__:"+string(dim(answerBasis)));
 print("__BASIS_VDIM__:"+string(vdim(answerBasis)));
 print("__NONUNIT_BASIS_START__");
 print(answerBasis);
 print("__NONUNIT_BASIS_FINISHED__");
}
print("__NATIVE_JOB_FINISHED__");
quit;
'''.replace('__ALGORITHM__',algorithm).replace('__GENCOUNT__',str(len(polys)))

def center_encoded(equation,point,qi,center):
 poly={}
 for word,n,d in equation:
  c=F(n,d);choices=[]
  for i in word:
   if i in qi:choices.append([((),center[qi[i]]),((qi[i],),F(1))])
   else:c*=point[i]
  if not c:continue
  for picked in itertools.product(*choices):
   mon=();v=c
   for w,a in picked:mon+=w;v*=a
   put(poly,tuple(sorted(mon)),v)
 return poly

def G6_consequence(model,point,qids,qi,center,original):
 path='reports/S61_UNIVERSAL_WEIGHT6_20261031.json';data=load(ROOT/path)
 assert data['all_passed']and data['q_to_u_indices_zero_based']==qids
 assert [{sum(model['parameter_weights'][j]for j in mon)for mon,n,d in eq}for eq in model['equations']]==[{4}]*70+[{5}]*64
 raw={};K=[{}for _ in range(70)];recipe=data['selected_non_Bianchi_syzygy'];assert len(recipe)==125
 for slot,c in recipe:
  j,i=divmod(slot,70);c=F(c);assert 0<=j<45 and 0<=i<70
  for mon,n,d in model['equations'][i]:put(raw,tuple(sorted(mon+[qids[j]])),c*F(n,d))
  put(K[i],(j,),c);put(K[i],(),c*center[j])
 target={tuple(mon):F(n,d)for mon,n,d in data['universal_new_equation']}
 assert raw==target and len(raw)==9257
 assert {sum(model['parameter_weights'][j]for j in mon)for mon in raw}=={6}
 g6=center_encoded(data['universal_new_equation'],point,qi,center);back={}
 for coeff,poly in zip(K,original):
  for m,c in coeff.items():
   for n,d in poly.items():put(back,tuple(sorted(m+n)),c*d)
 assert back==g6 and ()not in g6 and all(len(m)<=2 for m in g6)
 return g6,{'source_path':path,'source_sha256':sha(ROOT/path),'universal_all9257_coefficients_QQ_identity_checked':True,'universal_recipe125':recipe,'centered_all_coefficients_QQ_identity_checked':True,'cofactor_K_by_original_F_first70_QQ':list(map(encode,K)),'fixed_G6_QQ':encode(g6),'identity':'G6=sum_i0..69 K_i(deltaq)*F_i(deltaq); K_i=sum_(j,i) c*(deltaq_j+center_j)','degree_atmost2_checked':True}

def prepare(output,field,algorithm,add_g6):
 output.mkdir(parents=True,exist_ok=False);model=load(ROOT/MODEL);pr=load(ROOT/POINT);fr=load(ROOT/FRAME)
 point=list(map(F,next(x for x in pr['points']if x['component_id']=='T01')['S61_coordinates_u1_through_u61']))
 weights=model['parameter_weights'];qids=[i for i,w in enumerate(weights)if w==2];qi={u:i for i,u in enumerate(qids)};center=[point[i]for i in qids]
 assert len(qids)==45 and len(model['equations'])==134
 original=[center_encoded(eq,point,qi,center)for eq in model['equations']]
 assert all(()not in poly and all(len(m)<=2 for m in poly)for poly in original)
 W=[[F(c)for c in row]for row in fr['embedding']]
 assert len(W)==45 and all(len(row)==3 for row in W)and W[:3]==[[F(int(i==j))for j in range(3)]for i in range(3)]
 eta={(3,):F(1),(37,):F(1)}
 for h in range(3):put(eta,(h,),-W[3][h]-W[37][h])
 assert all(sum(eta.get((i,),F(0))*W[i][h]for i in range(45))==0 for h in range(3))
 inverse={():F(1)}
 for mon,c in eta.items():put(inverse,tuple(sorted(mon+(45,))),-c)
 original135=original+[inverse];polys=original135;g6data=None
 if add_g6:
  g6,g6data=G6_consequence(model,point,qids,qi,center,original);polys=[g6]+original135
 if field=='gf':assert all(c.denominator%P for poly in polys for c in poly.values())
 source_script=singular_script(polys,field,algorithm);(output/'native-input.sing').write_text(source_script)
 paths=[MODEL,POINT,FRAME,'code/s61_model.py','code/classification_native_support.py','research/classification/SYMMETRY_FIXED_NORMAL_CHART.md','reports/CLASSIFICATION_SYMMETRY_FIXED_SUM_20261106.json','reports/CLASSIFICATION_SYMMETRY_FIXED_SUM_INDEPENDENT_REVIEW_20261106.json']
 if add_g6:paths+=['reports/S61_UNIVERSAL_WEIGHT6_20261031.json','code/s61_weight6_review.py']
 packet={'task_id':'G','field':'GF(32003)'if field=='gf'else'QQ','field_mode':field,'prime':P if field=='gf'else None,'algorithm':algorithm,'G6_preconditioning':g6data,'native_generator_count':len(polys),'frozen_beta_ONLY':True,'all_original134_retained':True,'variables':['x'+str(i)for i in range(45)]+['u'],'q_to_u_indices_zero_based':qids,'center_q_QQ':list(map(str,center)),'original_full61_frozen_T01':list(map(str,point)),'normal_frame_J':[0,1,2],'eta_normal_indices0':[0,34],'eta_q_indices0':[3,37],'eta_centered_QQ':encode(eta),'eta_rule':'deltaq3+deltaq37-(Wrow3+Wrow37)*deltaqfirst3; GLOBALJ, not sparseK5normals','plane_W45x3_QQ':[[str(c)for c in row]for row in W],'equations_QQ':list(map(encode,polys)),'original135_equations_QQ':list(map(encode,original135)),'equation_term_counts':[len(p)for p in polys],'native_script_sha256':sha(output/'native-input.sing'),'input_sha256':{p:sha(ROOT/p)for p in paths},'generic_or_characteristic_zero_claim':False,'cofactor_certificate_computed':False}
 save(output/'input.json',packet)
 save(output/'prepare-report.json',{'task_id':'G','stage':'prepared_SOURCE_ONLY','all_passed':True,'input_sha256':sha(output/'input.json'),'native_script_sha256':packet['native_script_sha256'],'equations':len(polys),'original134_terms':sum(map(len,original)),'correct_global_J_eta':True,'G6_exact_QQ_identity_checked':bool(g6data),'field':packet['field'],'algorithm':algorithm,'unit_job_UNRUN':True})
 print(json.dumps({'SOURCE_prepared':True,'equations':len(polys),'terms':sum(map(len,polys)),'input_sha256':sha(output/'input.json'),'G6_terms':len(polys[0])if add_g6 else 0}))

def verify_export(log,packet,certificate):
 if '__COFACTORS_FINISHED__'not in log:return {'verified':False,'reason':'no_complete_export'}
 coeff=[{}for _ in packet['equations_QQ']];idx=None
 for line in log.splitlines():
  if line.startswith('__COFACTOR__:'):idx=int(line.split(':',1)[1]);assert 0<=idx<len(coeff)
  elif line.startswith('__TERM__:'):
   _,c,e=line.split(':',2);powers=[int(x)for x in e.replace('(','').replace(')','').split(',')];assert idx is not None and len(powers)==46 and min(powers)>=0
   mon=tuple(i for i,a in enumerate(powers)for _ in range(a));put(coeff[idx],mon,F(c))
 field=packet['field_mode'];ops=0;start=time.monotonic()
 def qmod(c):return c.numerator*pow(c.denominator,-1,P)%P
 def addfield(poly,mon,c):
  if field=='gf':
   z=(poly.get(mon,0)+qmod(c))%P
   if z:poly[mon]=z
   else:poly.pop(mon,None)
  else:put(poly,mon,c)
 def budget(size):
  if ops%10000==0 and(time.monotonic()-start>60 or ops>10000000 or size>1000000):raise TimeoutError('independent_multiplyback_resource_limit')
 def multiplyback(cof,inputs):
  nonlocal ops
  res={}
  for a,b in zip(cof,inputs):
   for m,c in a.items():
    for n,d in b.items():
     addfield(res,tuple(sorted(m+n)),F(c)*d);ops+=1;budget(len(res))
  assert res=={():1 if field=='gf'else F(1)},'cofactor product not1'
 try:multiplyback(coeff,list(map(decode,packet['equations_QQ'])))
 except TimeoutError:return {'verified':False,'reason':'independent_multiplyback_resource_limit','term_products':ops}
 folded=False;used_coeff=coeff;used_inputs=packet['equations_QQ'];fold_reason=None
 if packet['G6_preconditioning']:
  # The already verified136identity plus explicit G6 source identity is itself
  # an exact certificate. Fold to original135 only within the SAME check cap.
  proposed=[dict(x)for x in coeff[1:]]
  try:
   for i,poly in enumerate(packet['G6_preconditioning']['cofactor_K_by_original_F_first70_QQ']):
    for m,c in coeff[0].items():
     for n,d in decode(poly).items():
      addfield(proposed[i],tuple(sorted(m+n)),F(c)*d);ops+=1;budget(len(proposed[i]))
   multiplyback(proposed,list(map(decode,packet['original135_equations_QQ'])))
   used_coeff=proposed;used_inputs=packet['original135_equations_QQ'];folded=True
  except TimeoutError:fold_reason='kept_verified136cofactors_plus_explicitG6identity_due_check_cap'
 body={'task_id':'G','fixed_rational_beta_only':True,'field':packet['field'],'cofactors':list(map(encode,used_coeff)),'generator_input_QQ':used_inputs,'input_sha256':packet['input_sha256'],'independent_sparse_multiplyback_equals1':True,'term_products':ops,'G6_folded_to_original135':folded,'G6_exact_source_identity':packet['G6_preconditioning'],'fold_note':fold_reason,'global_generic_exhaustion_claim':False,'char0_claim':field=='qq'}
 save(certificate,body);return {'verified':True,'certificate_path':str(certificate),'certificate_sha256':sha(certificate),'term_products':ops,'G6_folded_to_original135':folded,'cofactor_count':len(used_coeff)}

def memory_preflight(memory_mb):
 maximum_path=Path('/sys/fs/cgroup/memory.max');current_path=Path('/sys/fs/cgroup/memory.current')
 maximum_raw=maximum_path.read_text().strip();current=int(current_path.read_text())
 maximum=int(maximum_raw)if maximum_raw!='max'else None;requested=memory_mb*1024**2;reserve=2*1024**3
 headroom=None if maximum is None else maximum-current-requested
 return {'cgroup_memory_max_path':str(maximum_path),'cgroup_memory_current_path':str(current_path),'memory_max_bytes':maximum,'memory_current_bytes':current,'requested_native_address_space_bytes':requested,'required_headroom_bytes':reserve,'remaining_headroom_bytes':headroom,'source_budget_feasible':headroom is not None and headroom>=reserve,'note':'Fresh runtime check; stored SOURCE snapshot alone does not authorize future memory use'}

def run(prepared,output,binary,seconds,memory_mb,certificate,expected):
 inputpath=prepared/'input.json';packet=load(inputpath);assert sha(inputpath)==expected
 assert sha(prepared/'native-input.sing')==packet['native_script_sha256']
 for p,h in packet['input_sha256'].items():assert sha(ROOT/p)==h,p
 memory_snapshot=memory_preflight(memory_mb)
 if not memory_snapshot['source_budget_feasible']:raise RuntimeError('Native memory preflight blocked: '+json.dumps(memory_snapshot))
 output.mkdir(parents=True,exist_ok=False);env=os.environ.copy();lib=ROOT/'work/native-cas/root/usr/lib/x86_64-linux-gnu';env['LD_LIBRARY_PATH']=str(lib);env['SINGULAR_PROCS_DIR']=str(lib/'singular/MOD')
 def limits():
  resource.setrlimit(resource.RLIMIT_CPU,(seconds,seconds+1));resource.setrlimit(resource.RLIMIT_AS,(memory_mb*1024**2,memory_mb*1024**2));resource.setrlimit(resource.RLIMIT_FSIZE,(64*1024**2,64*1024**2))
 version=subprocess.run([str(binary),'--version'],env=env,capture_output=True,text=True,timeout=10)
 (output/'native-version.txt').write_text(version.stdout+version.stderr);assert version.returncode==0
 start=time.monotonic();preusage=resource.getrusage(resource.RUSAGE_CHILDREN);wall_stopped=False
 with (output/'native.stdout.txt').open('w')as stdout,(output/'native.stderr.txt').open('w')as stderr:
  p=subprocess.Popen([str(binary),'-q','--no-rc','--no-stdlib','--no-shell',str((prepared/'native-input.sing').resolve())],cwd=ROOT,env=env,stdout=stdout,stderr=stderr,preexec_fn=limits,start_new_session=True)
  try:p.wait(timeout=seconds+15)
  except subprocess.TimeoutExpired:
   wall_stopped=True;os.killpg(p.pid,signal.SIGTERM)
   try:p.wait(timeout=3)
   except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
 usage=resource.getrusage(resource.RUSAGE_CHILDREN)
 log=(output/'native.stdout.txt').read_text();err=(output/'native.stderr.txt').read_text()
 native_errors=[line for line in (log+'\n'+err).splitlines()if line.lstrip().startswith('?')];completed='__NATIVE_JOB_FINISHED__'in log and not native_errors;std='__NATIVE_STD_FINISHED__'in log and not native_errors;unit='__UNIT_DETECTED__'in log and not native_errors
 try:verified=verify_export(log,packet,certificate)if unit else{'verified':False,'reason':'unit_not_detected'}
 except Exception as error:verified={'verified':False,'reason':'invalid_or_unfinished_export','error':repr(error)}
 def marker_int(prefix):
  values=[line.split(':',1)[1].strip()for line in log.splitlines()if line.startswith(prefix+':')]
  return int(values[-1])if values else None
 report={'task_id':'G','backend':'Singular native '+packet['algorithm']+', then optional lift in SAME CPUwindow','field':packet['field'],'fixed_rational_beta_only':True,'normal_frame_J':[0,1,2],'eta':'n0+n34 in GLOBALJframe, qrows3+37','all_original134_retained':True,'equations':packet['native_generator_count'],'G6_preconditioned':bool(packet['G6_preconditioning']),'std_completed':std,'unit_detected_by_backend':unit,'native_lift_identity_checked':'__NATIVE_LIFT_VERIFIED__'in log,'cofactor_independent_check':verified,'job_completed':completed,'nonunit_basis_export_complete':'__NONUNIT_BASIS_FINISHED__'in log,'native_reported_standard_basis_count':marker_int('__BASIS_COUNT__'),'native_reported_dimension':marker_int('__BASIS_DIM__'),'native_reported_vector_dimension':marker_int('__BASIS_VDIM__'),'returncode':p.returncode,'wall_limit_triggered':wall_stopped,'native_cpu_seconds_observed':usage.ru_utime+usage.ru_stime-preusage.ru_utime-preusage.ru_stime,'native_child_peak_rss_kib':usage.ru_maxrss,'native_memory_preflight':memory_snapshot,'native_limits':{'cpu_seconds':seconds,'memory_mb':memory_mb,'wall_seconds':seconds+15,'logfile_cap_bytes':64*1024**2},'seconds':time.monotonic()-start,'input_path':str(inputpath),'input_sha256':sha(inputpath),'singular_binary_sha256':sha(binary),'native_script_sha256':packet['native_script_sha256'],'stdout_sha256':sha(output/'native.stdout.txt'),'stderr_sha256':sha(output/'native.stderr.txt'),'stderr_nonempty':bool(err.strip()),'native_error_lines':native_errors[-20:],'generic_exhaustion_proved':False,'characteristic_zero_inferred_from_modp':False,'status':'completed_fixedfiber_diagnostic'if completed else'resource_limit_or_backend_failure','next':'Use mathematical diagnostics/certificate iffound; no point/prime/chart/backend repetition without new structuralreason'}
 save(output/'report.json',report);print(json.dumps(report,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='mode',required=True)
 a=s.add_parser('prepare');a.add_argument('--output',required=True,type=Path);a.add_argument('--field',choices=['gf','qq'],default='gf');a.add_argument('--algorithm',choices=['std','slimgb'],default='std');a.add_argument('--add-g6',action='store_true')
 b=s.add_parser('run');b.add_argument('--prepared',required=True,type=Path);b.add_argument('--output',required=True,type=Path);b.add_argument('--binary',type=Path,default=ROOT/'work/native-cas/root/usr/bin/Singular');b.add_argument('--seconds',type=int,default=300);b.add_argument('--memory-mb',type=int,default=2048);b.add_argument('--certificate',required=True,type=Path);b.add_argument('--source-sha256',required=True)
 a=p.parse_args()
 if a.mode=='prepare':prepare(a.output,a.field,a.algorithm,a.add_g6)
 else:
  if not 1<=a.seconds<=600 or not 256<=a.memory_mb<=12288:p.error('Scoped limits:1..600CPUseconds,256..12288MiB and at least2GiB fresh cgroup headroom')
  run(a.prepared,a.output,a.binary,a.seconds,a.memory_mb,a.certificate,a.source_sha256)
