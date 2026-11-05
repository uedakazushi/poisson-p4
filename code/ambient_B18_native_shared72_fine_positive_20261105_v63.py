#!/usr/bin/env python3
"""Positive-only native72 fine outputs and sparse605 composition; no responses."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NATIVE_SHARED72_FINE_POSITIVE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;bits=0;prior=None;cursor={};Z=(1,0)
 def read(path):
  nonlocal raw
  p=ROOT/path;pins[path]=sha(p);raw+=2*p.stat().st_size;assert raw<32*1024*1024;return json.loads(p.read_text())
 def ck(n):
  nonlocal bits
  bits=max(bits,abs(n).bit_length());assert bits<=384 and (prior is None or bits<=prior);return n
 def ex(q):q=F(q);return(q.denominator,abs(q.numerator)) if q else Z
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=ck(math.lcm(a[0],b[0]));return d,ck(ck(a[1]*(d//a[0]))+ck(b[1]*(d//b[0])))
 def mul(a,b):return(ck(a[0]*b[0]),ck(a[1]*b[1])) if a[1] and b[1] else Z
 def dump(name,v):(out/name).write_text(json.dumps(v,separators=(',',':'))+'\n')
 def enc(q):return list(map(str,q))
 sp='reports/AMBIENT_B18_NATIVE_SHARED72_INPUT_SOURCE_20261105_v63/'
 sr=read(sp+'report.json');assert sr['all_passed'] and pins[sp+'report.json']=='4f01895013eb59d8d1d5721eff04d8131e9a73f52f2449f9613731e311c704e6'
 se='native-shared75-formal-to72-actual-to605-source-input.json';s=read(sp+se);assert pins[sp+se]==sr['evidence_sha256'][se]
 nrp='reports/AMBIENT_B18_NATIVE_ONE554_NUMERIC_SOURCE_20261105_v63_REVIEW/report.json';nr=read(nrp);assert nr['all_passed'] and pins[nrp]=='f9581e3b3c18c9f98856317d9cf72da89f7330d3149b45ddffb82d4cd9bcb8a9'
 np='reports/AMBIENT_B18_NATIVE_ONE554_POSITIVE_INPUT_20261105_v63/ONE554-positive-four-native-and-nonzero644-outer-input.json';n=read(np)
 assert pins[np]==nr['input_sha256'][np]
 R={tuple(key):{i:(int(d),int(u)) for i,d,u in row} for key,row in n['reused68_positive_covectorL_rows']};assert len(R)==68
 p=s['native_primitive_known_cache'];base=[tuple(m) for m in p['basis18']];bi={m:i for i,m in enumerate(base)};table=[[bi.get(tuple(x+y for x,y in zip(a,b))) for b in base] for a in base]
 G=[[ex(q) for q in row] for row in p['G8_C72']];nu=[[ex(q) for q in row] for row in p['nu_two_J18']]
 program={r['native128']:r for r in s['raw75_original_Koszul_program']};restrictions=s['kernel72_restrictions'];assert len(restrictions)==72
 assert all(table[i][j] is None for a in nu for b in nu for i,x in enumerate(a) if x[1] for j,y in enumerate(b) if y[1])
 # Positive BEFORE-factor book. R envelopes alreadyincludeΓ/D andN/oneν.
 # Any fine target has <=34pieces*48dotrows*3sources*2Koszulentries*
 # 4Gmasks*19optionalνterms <2^20 paths. Onefromeachclasssuffices.
 classes={'R68':[q for row in R.values() for q in row.values()],
  'G8':[q for row in G for q in row], 'nu':[Z]+[q for row in nu for q in row],
  'kernelweights':[ex(q) for rec in restrictions for src,q in rec['original_native_source_combination']]}
 book={};D0=1;U0=2**20
 for name,vals in classes.items():
  d=1
  for den,u in vals:d=ck(math.lcm(d,den))
  u=max([d]+[ck(v*(d//den)) for den,v in vals]);book[name]=[str(d),str(u)];D0=ck(D0*d);U0=ck(U0*u)
 fineprior=max(D0.bit_length(),U0.bit_length());assert 34*48*3*2*4*19<2**20
 # Forconstant605compositionallowonecoordinateLCM/L1rowfactor.
 mapping=s['all605_constant_native_coordinate_recipes'];rowclasses=[]
 for rec in mapping:
  d=1
  for i,q in rec['kernel105_coordinate_pairs']:d=ck(math.lcm(d,F(q).denominator))
  u=sum(abs(F(q).numerator)*(d//F(q).denominator) for i,q in rec['kernel105_coordinate_pairs']);rowclasses.append((d,u))
 composed=max([fineprior]+[max(ck(D0*d).bit_length(),ck(U0*u).bit_length()) for d,u in rowclasses if u]);prior=max(354,composed);assert prior<=384 and 3*prior+1<=1153
 before=count([sr,s,nr,n,R,p,base,bi,table,G,nu,program,restrictions,classes,book,mapping,rowclasses,pins])+20000;assert before<150000
 dump('phase-before-positive-S-and-fine-contractions.json',{'SOURCE_only_no_function_response':True,'positive_classes':book,'positive_common_D':str(D0),'positive_common_U_paths2pow20':str(U0),'fine_source_prior':fineprior,'composition_source_prior':composed,'union_inherited_coarse354_prior':prior,'temporary_prior':3*prior+1,'owned_before_positive':before,'exact_point_S0_override_only_after_fixed_kernel_certificate':True})
 # Exactpointkernelproofisfrozen15b6/faec; onlyPOINTisoverridden.
 kpr='reports/AMBIENT_B18_NATIVE_POINT_KERNEL105_SOURCE_20261105_v63/report.json';kr=read(kpr);assert kr['all_passed'] and pins[kpr]=='15b6ff370efd02939a1fac290ca6c24c1237ab3046f6a82984d22b1071e8d534'
 cubic=[(0,17),(5,16),(6,15),(15,6),(16,5)];pairs=[(0,7),(0,8),(1,7),(1,8),(4,7),(4,8),(5,7),(5,8),(7,7),(7,8)];vf={0:1,1:2,4:3,5:4,7:5,8:6};results=[];points=[];fingerprints=[];fine={}
 def dot(a,b):
  z=Z
  for i,x in a.items():
   if i in b:z=add(z,mul(x,b[i]))
  return z
 for rec in restrictions:
  i=rec['kernel105'];cursor={'stage':'positive_weightedS','kernel105':i};S=[{} for _ in base]
  for src,q in rec['original_native_source_combination']:
   op=program[src];im=op['inputmask'];weight=ex(q)
   for row,g,sign in op['Koszul_entries']:
    for gb in range(18):
     for mask in range(4):
      x=G[g][4*gb+mask]
      if not x[1]:continue
      collision=mask&im;col=4*row+(mask^im)
      if not collision:S[gb][col]=add(S[gb].get(col,Z),mul(weight,x))
      elif collision in [1,2]:
       for vb,y in enumerate(nu[collision-1]):
        dest=table[gb][vb]
        if y[1] and dest is not None:S[dest][col]=add(S[dest].get(col,Z),mul(mul(weight,x),y))
  points.append([i,[[col,*enc(q)] for col,q in sorted(S[0].items())]])
  S[0]={} # actualfullpointidentity, NOTpositivityortruehigherzero.
  encodedS=[[[col,*enc(q)] for col,q in sorted(row.items())] for row in S];fingerprints.append([i,hashlib.sha256(json.dumps(encodedS,separators=(',',':')).encode()).hexdigest()])
  channels=[];piecebook=[]
  for t in range(2):
   a=Z;b=Z;cg=[];ds=[]
   for j,(rb,sb) in enumerate(cubic):q=dot(R[0,t,-1,rb],S[sb]);a=add(a,q);cg.append(q)
   for pair,(x,y) in enumerate(pairs):
    q=[dot(R[1,t,pair,0],S[7+pair])]
    for r,ss in ([(x,y)] if x==y else[(x,y),(y,x)]):q.append(dot(R[1,t,pair,vf[r]],S[vf[ss]]))
    for v in q:b=add(b,v)
    ds.append(q)
   total=add(a,b);channels.append([enc(a),enc(b),enc(total)]);piecebook.append([list(map(enc,cg)),[[enc(q) for q in row] for row in ds]])
   fine[i,t]=total
  results.append({'kernel105':i,'two_channels_Gamma_D_sum_DU':channels,'fine68_piece_DU_sha256':hashlib.sha256(json.dumps(piecebook,separators=(',',':')).encode()).hexdigest()})
 mapped=[]
 for rec in mapping:
  i=rec['case619'];ch=[]
  for t in range(2):
   v=Z
   for k,q in rec['kernel105_coordinate_pairs']:v=add(v,mul(ex(q),fine[k,t]))
   ch.append(enc(v))
  mapped.append([i,ch])
 ep='native72-fine-Gamma-D-and605-constantmap-positive-envelopes.json';data={'positive_only_no_signed_native_or_upper_response':True,'native_kernel72_Gamma_D_sum_fine':results,
  'all72_positive_S18by48_fingerprints_AFTER_point_override':fingerprints,'raw_positive_point_bounds_BEFORE_exact_source_override':points,
  'all605_native_sum_fine_twochannels':mapped,'basis_common_coarse341':True,'basis_coarse_temporary1024':True,
  'mapping_coarse354_temporary1063_union_retained':True,'new_fine_positive_prior':fineprior,'new_coordinate_composition_positive_prior':composed,
  'source_committed_union_prior':prior,'source_temporary_union_prior':3*prior+1,'actual_positive_operation_bits':bits,
  'point_override_scope':'ONLYfullS0. Fixed72constantkernelcertificate; allhighermonic/crosscoefficientskept. DoNOTsetindividualraw4/8/12pointzero.',
  'whole_gate_must_check_fine_outputs_BEFORE_Q':'Useeach72fine sum,compose605fine,thenOLD_FIRST_outerandQ. Nevermultiplycoarse354directlybyQ.',
  'SOURCE_inputs_sameJ18_N_G_nu_GammaD_andcurrentcomparison_frame':True,'OLD_FIRST_outer_honest_or_whole605_values_NOTcomputed':True}
 dump(ep,data);owned=count([sr,s,nr,n,R,p,base,bi,table,G,nu,program,restrictions,classes,book,mapping,rowclasses,kr,S,results,points,fingerprints,fine,mapped,data,pins])+15000;assert owned<150000
 rep={'id':STEM,'all_passed':True,'native_F_FINE_positive_SOURCE_only':True,'no_signed_native_or605upper_values':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep),'phase-before-positive-S-and-fine-contractions.json':sha(out/'phase-before-positive-S-and-fine-contractions.json')},'source_owned_upper':owned,'source_RAW_upper':raw+sum(p.stat().st_size for p in out.iterdir()),'source_union_prior':prior,'source_union_temporary':3*prior+1,'actual_positive_bits':bits,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'fine_prior':fineprior,'mapping_prior':composed,'union':prior,'actual':bits,'owned':owned}))
if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--output',required=True,type=Path);a=pa.parse_args();main(a.output)
