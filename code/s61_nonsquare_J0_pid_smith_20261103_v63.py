#!/usr/bin/env python3
"""Bounded genuine one-variable Smith model for the J0-safe cubic cokernel."""
import argparse,hashlib,json,resource,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',required=True,type=Path);ap.add_argument('--field',choices=['a','e'],required=True)
 ap.add_argument('--output',required=True,type=Path);ap.add_argument('--timeout',type=int,default=180);ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args()
 args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic();names=['a','e'];fi=names.index(args.field);vi=1-fi;variable=names[vi]
 data=args.output/'matrix.sing';run=args.output/'run.sing';proof=args.output/'right-annihilator.sing';diag=args.output/'diagonal.txt';statusfile=args.output/'result.txt'
 with args.source.open()as inp,data.open('w')as out:
  nr,nc,metadata=map(int,inp.readline().split());assert(nr,nc,metadata)==(36,70,0)
  out.write(f'ring R=(0,{args.field}),({variable}),dp;\nmatrix M[{nr}][{nc}];\n')
  for i in range(nr):
   for j in range(nc):
    tokens=inp.readline().split();n=int(tokens[0]);assert len(tokens)==1+4*n;groups={}
    for k in range(n):
     coefficient=tokens[1+4*k];ex=list(map(int,tokens[2+4*k:5+4*k]));assert ex[2]==0
     groups.setdefault(ex[vi],[]).append((ex[fi],coefficient))
    terms=[]
    for power,coeffs in sorted(groups.items()):
     coefficient='+'.join('('+c+')'+(''if k==0 else'*'+args.field+(''if k==1 else'^'+str(k)))for k,c in coeffs)
     terms.append('('+coefficient+')'+(''if power==0 else'*'+variable+(''if power==1 else'^'+str(power))))
    if terms:out.write(f'M[{i+1},{j+1}]='+ '+'.join(terms)+';\n')
  out.write('poly KnownUnit=a*(500*a+90*e^2+3*e-24)*(300*a*e+54*e^3+9*e^2-24*e-4)*(100*a+18*e^2-9*e-2)*(3*e+2);\n')
 script=f'''LIB "jacobson.lib";
execute(read("{data.resolve()}"));
print("EXACT_PID_INPUT_READY_{args.field}");
list S=smith(M,1);
if(S[1]*M*S[3]!=S[2]){{ERROR("Smith U*M*V identity failed");}}
print("EXACT_SMITH_IDENTITY_PASS_{args.field}");
matrix D=S[2];poly ann=1;int full=1;poly di;
for(int i=1;i<={nr};i++){{
 di=D[i,i];write("{diag.resolve()}",string(i)+" "+string(di));
 if(di==0){{full=0;break;}}
 di=di/leadcoef(di);ann=ann*(di/gcd(ann,di));ann=ann/leadcoef(ann);
}}
if(full==1){{
 poly remaining=ann;poly common;
 while(deg(remaining)>0){{common=gcd(remaining,KnownUnit);if(deg(common)<=0){{break;}}remaining=remaining/common;}}
 write("{statusfile.resolve()}","ANN="+string(ann));
 write("{statusfile.resolve()}","REMAINING="+string(remaining));
 matrix E[{nc}][{nr}];for(int i=1;i<={nr};i++){{E[i,i]=ann/D[i,i];}}
 matrix X=S[3]*E*S[1];
 if(M*X!=ann*unitmat({nr})){{ERROR("actual right-annihilator cofactor identity failed");}}
 write("{proof.resolve()}","poly annProof="+string(ann)+";");
 write("{proof.resolve()}","matrix Xproof[{nc}][{nr}];");
 for(int i=1;i<={nc};i++){{for(int j=1;j<={nr};j++){{if(X[i,j]!=0){{write("{proof.resolve()}","Xproof["+string(i)+","+string(j)+"]="+string(X[i,j])+";");}}}}}}
 write("{proof.resolve()}","if(M*Xproof!=annProof*unitmat({nr})){{ERROR(\"right-annihilator replay failed\");}}print(\"EXACT_PID_RIGHT_ANNIHILATOR_REPLAY_PASS\");");
 if(deg(remaining)<=0){{print("EXACT_LOCALIZED_PID_COKER_ZERO_{args.field}");write("{statusfile.resolve()}","LOCALIZED_UNIT=1");}}
 else{{print("EXACT_PID_NONUNIT_FACTOR_{args.field}");write("{statusfile.resolve()}","LOCALIZED_UNIT=0");}}
}}else{{write("{statusfile.resolve()}","GENERIC_RANK_DEFICIENT=1");print("EXACT_PID_GENERIC_RANK_DEFICIENT");}}
quit;
'''
 run.write_text(script)
 receipt={'id':'S61_NONSQUARE_J0_PID_SMITH_20261103_V63','field':f'QQ({args.field})[{variable}]','shape':[nr,nc],
          'known_localization_only':'a,F,Bs,C,3e+2','source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),
          'data_sha256':hashlib.sha256(data.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(run.read_bytes()).hexdigest(),
          'AS_limit_bytes':2200000000,'CPU_limit_seconds':args.timeout}
 if not args.prepare_only:
  def limits():resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));resource.setrlimit(resource.RLIMIT_CPU,(args.timeout,args.timeout+1))
  wrapper=ROOT/'work/nonsquare-v62/singular-runtime/run-singular-optimized-v6.sh'
  with (args.output/'solver.log').open('w')as log:
   try:
    job=subprocess.run([str(wrapper),'--cpus=0','-q',str(run.resolve())],stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout+10,preexec_fn=limits)
    status='completed'if job.returncode==0 else'exit_'+str(job.returncode)
   except subprocess.TimeoutExpired:status='timeout'
  content=(args.output/'solver.log').read_text();receipt.update(status=status,exact_smith_identity_pass='EXACT_SMITH_IDENTITY_PASS_'in content,
      exact_localized_coker_zero=f'EXACT_LOCALIZED_PID_COKER_ZERO_{args.field}'in content,
      result_file_exists=statusfile.exists(),mathematical_unit_result_obtained=f'EXACT_LOCALIZED_PID_COKER_ZERO_{args.field}'in content)
 receipt['seconds']=time.monotonic()-start;(args.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
