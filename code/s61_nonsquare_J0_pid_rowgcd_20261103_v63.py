#!/usr/bin/env python3
"""Localized PID row-gcd/Bezout-column elimination, with exact replay steps."""
import argparse,hashlib,json,re,resource,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
LIBRARY=r'''
proc VerifiedDivRem(poly f,poly g)
{
 if(g==0){ERROR("division by zero");}
 poly q=0;poly r=f;poly term;number c;int previousDegree;
 while((r!=0)&&(deg(r)>=deg(g))){
  previousDegree=deg(r);c=leadcoef(r)/leadcoef(g);
  term=c*var(1)^(deg(r)-deg(g));q=q+term;r=r-term*g;
  if((r!=0)&&(deg(r)>=previousDegree)){ERROR("long division degree failed");}
 }
 if(q*g+r!=f){ERROR("long division identity failed");}
 list result;result[1]=q;result[2]=r;return(result);
}
proc VerifiedExtGCD(poly f,poly g)
{
 poly r0=f;poly r1=g;poly s0=1;poly s1=0;poly t0=0;poly t1=1;
 poly q;poly r;poly s;poly t;number c;list divisionResult;
 while(r1!=0){
  divisionResult=VerifiedDivRem(r0,r1);q=divisionResult[1];r=divisionResult[2];
  if((r!=0)&&(deg(r)>=deg(r1))){ERROR("Euclidean degree failed");}
  s=s0-q*s1;t=t0-q*t1;r0=r1;r1=r;s0=s1;s1=s;t0=t1;t1=t;
 }
 if(r0==0){ERROR("zero gcd input");}
 c=1/leadcoef(r0);r0=c*r0;s0=c*s0;t0=c*t0;
 if(s0*f+t0*g!=r0){ERROR("verified Euclidean Bezout identity failed");}
 list answer;answer[1]=r0;answer[2]=s0;answer[3]=t0;return(answer);
}
proc isKnownUnit(poly g,poly U)
{
 if(g==0){return(0);}
 poly h;list E;list divisionResult;
 while(deg(g)>0){E=VerifiedExtGCD(g,U);h=E[1];if(deg(h)<=0){break;}divisionResult=VerifiedDivRem(g,h);if(divisionResult[2]!=0){ERROR("unit factor division failed");}g=divisionResult[1];}
 if(deg(g)<=0){return(1);}return(0);
}
proc ReplayStep(matrix A,matrix T,poly g,int pivot,poly U)
{
 int n=nrows(A);int m=ncols(A);int i;int j;int ii;
 if(isKnownUnit(g,U)==0){ERROR("pivot not known localized unit");}
 matrix v=A*T;
 if(v[pivot,1]!=g){ERROR("Bezout column identity failed");}
 matrix H[1][m];list divisionResult;
 for(j=1;j<=m;j++){divisionResult=VerifiedDivRem(A[pivot,j],g);if(divisionResult[2]!=0){ERROR("pivot row exact division failed");}H[1,j]=divisionResult[1];}
 list answer;
 if(n==1){matrix Z[1][m];answer[1]=0;answer[2]=Z;return(answer);}
 matrix N[n-1][m];ii=0;
 for(i=1;i<=n;i++){if(i!=pivot){ii++;for(j=1;j<=m;j++){N[ii,j]=A[i,j]-v[i,1]*H[1,j];}}}
 answer[1]=n-1;answer[2]=N;return(answer);
}
'''

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--matrix-input',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--timeout',default=120,type=int);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
 inp=args.matrix_input.resolve();text=inp.read_text();shape=re.search(r'matrix M\[(\d+)\]\[(\d+)\]',text);nr,nc=map(int,shape.groups());field=re.search(r'ring R=\(0,([ae])\)',text)[1]
 library=args.output/'replay-library.sing';library.write_text(LIBRARY);proof=args.output/'proof.sing';run=args.output/'run.sing';statusfile=args.output/'result.txt'
 script=f'''execute(read("{inp}"));
execute(read("{library.resolve()}"));
matrix A=M;int left=nrows(A);int step=0;int m=ncols(A);
write("{proof.resolve()}","execute(read(\\\"{inp}\\\"));");
write("{proof.resolve()}","execute(read(\\\"{library.resolve()}\\\"));matrix Mcur=M;int remainingRows=nrows(Mcur);");
print("LOCALIZED_ROW_GCD_INPUT_READY_{field}");
while(left>0){{
 int pivot=0;int bestCost=2147483647;int cost;int i;int j;int minDegree;int column;int candidate;
 for(i=1;i<=left;i++){{cost=0;for(j=1;j<=m;j++){{if(A[i,j]!=0){{cost=cost+(deg(A[i,j])+1)*(size(A[i,j])+1);}}}}if(cost<bestCost){{pivot=i;bestCost=cost;}}}}
 if(pivot==0){{ERROR("no pivot row selected");}}
 minDegree=2147483647;column=0;
 for(j=1;j<=m;j++){{if(A[pivot,j]!=0){{if(deg(A[pivot,j])<minDegree){{minDegree=deg(A[pivot,j]);column=j;}}}}}}
 if(column==0){{print("EXACT_ACTUAL_COKER_ZERO_ROW");write("{statusfile.resolve()}","ZERO_ROW=1");break;}}
 matrix T[m][1];number inverseLC=1/leadcoef(A[pivot,column]);poly g=inverseLC*A[pivot,column];T[column,1]=inverseLC;
 intvec used;used[m]=0;used[column]=1;poly h;poly gnext;number lc;list E;
 while(deg(g)>0){{
  minDegree=2147483647;candidate=0;
  for(j=1;j<=m;j++){{if((used[j]==0)&&(A[pivot,j]!=0)){{if(deg(A[pivot,j])<minDegree){{minDegree=deg(A[pivot,j]);candidate=j;}}}}}}
  if(candidate==0){{break;}}used[candidate]=1;E=VerifiedExtGCD(g,A[pivot,candidate]);h=E[1];
  if(h!=g){{
   if(E[2]*g+E[3]*A[pivot,candidate]!=E[1]){{ERROR("extended gcd identity failed");}}
   lc=1/leadcoef(E[1]);T=E[2]*T;T[candidate,1]=T[candidate,1]+E[3];
   T=lc*T;g=lc*E[1];
  }}
 }}
 list quotientCheck;for(j=1;j<=m;j++){{quotientCheck=VerifiedDivRem(A[pivot,j],g);if(quotientCheck[2]!=0){{ERROR("full row gcd does not divide entry");}}}}
 if(isKnownUnit(g,KnownUnit)==0){{
  print("EXACT_LOCALIZED_COKER_NONUNIT_ROW_FACTOR_{field}");
  write("{statusfile.resolve()}","NONUNIT_ROW="+string(pivot));write("{statusfile.resolve()}","NONUNIT_FACTOR="+string(g));
  write("{(args.output/'nonunit-current-matrix.txt').resolve()}",A);break;
 }}
 step++;print("EXACT_ROW_GCD_PIVOT_{field} "+string(step)+" LEFT "+string(left)+" ROW "+string(pivot)+" GDEG "+string(deg(g)));
 write("{proof.resolve()}","matrix T"+string(step)+"["+string(m)+"][1];");
 for(j=1;j<=m;j++){{if(T[j,1]!=0){{write("{proof.resolve()}","T"+string(step)+"["+string(j)+",1]="+string(T[j,1])+";");}}}}
 write("{proof.resolve()}","poly g"+string(step)+"="+string(g)+";");
 write("{proof.resolve()}","list step"+string(step)+"=ReplayStep(Mcur,T"+string(step)+",g"+string(step)+","+string(pivot)+",KnownUnit);Mcur=step"+string(step)+"[2];remainingRows=step"+string(step)+"[1];");
 list answer=ReplayStep(A,T,g,pivot,KnownUnit);left=answer[1];A=answer[2];
 kill T;kill used;kill answer;
}}
if(left==0){{
 print("EXACT_LOCALIZED_PID_ROW_COKER_ZERO_{field}");write("{statusfile.resolve()}","LOCALIZED_UNIT=1");
 write("{proof.resolve()}","if(remainingRows!=0){{ERROR(\\\"remaining rows in replay\\\");}}print(\\\"EXACT_LOCALIZED_PID_ROW_REPLAY_PASS_{field}\\\");quit;");
}}
quit;
'''
 # Singular strings use one escaped quote to write a literal quote.
 script=script.replace('\\\\\\"','\\"')
 run.write_text(script)
 def limits():resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));resource.setrlimit(resource.RLIMIT_CPU,(args.timeout,args.timeout+1))
 with (args.output/'solver.log').open('w')as log:
  try:
   result=subprocess.run([str(ROOT/'work/nonsquare-v62/singular-runtime/run-singular-optimized-v6.sh'),'--cpus=0','-q',str(run.resolve())],stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout+10,preexec_fn=limits);status='completed'if result.returncode==0 else'exit_'+str(result.returncode)
  except subprocess.TimeoutExpired:status='timeout'
 content=(args.output/'solver.log').read_text();arithmetic_error=bool(re.search(r'(?m)^\s*\?',content));status='arithmetic_error'if arithmetic_error else status
 report={'id':'S61_NONSQUARE_J0_PID_ROW_GCD_20261103_V63','field':f'QQ({field})[{"e"if field=="a"else"a"}]','shape':[nr,nc],
  'status':status,'arithmetic_error_detected':arithmetic_error,'exact_localized_coker_zero':not arithmetic_error and f'EXACT_LOCALIZED_PID_ROW_COKER_ZERO_{field}'in content,
  'actual_nonunit_row_factor_found':not arithmetic_error and f'EXACT_LOCALIZED_COKER_NONUNIT_ROW_FACTOR_{field}'in content,
  'completed_unit_pivots':content.count(f'EXACT_ROW_GCD_PIVOT_{field} '),'result_exists':statusfile.exists(),
  'AS_limit_bytes':2200000000,'CPU_limit_seconds':args.timeout,'seconds':time.monotonic()-started,
  'input_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),inp,library,run]}}
 (args.output/'receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
