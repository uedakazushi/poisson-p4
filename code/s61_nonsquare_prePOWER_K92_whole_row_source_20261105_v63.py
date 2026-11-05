#!/usr/bin/env python3
"""SOURCE ONLY one whole Schur-row identity, using exact block2+2 algebra."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_WHOLE_ROW_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def support_list(n):
 while n:
  b=n&-n;yield b.bit_length()-1;n-=b
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 last='reports/S61_NONSQUARE_PREPOWER_K92_BORDER4_DIRECT_ENTRY_20261105_v63';r=src(last+'/report.json');assert r['all_passed']and not r['target_NONZERO'];literal=src(last+'/literal-source-and-all-finite-path-identities.json');T={(e['row'],e['column']):dict((int(g),int(v))for g,v in e['integer_polynomial'])for e in literal['all_Tbar']}
 ratios={Q(v,T[1,1][g-68])for g,v in T[2,1].items()};assert ratios=={Q(59375)};kappa=59375;assert all(T[2,j]=={g+68:kappa*v for g,v in T[1,j].items()}for j in range(7))
 files=['reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/ABC72.jsonl','reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/D49.jsonl'];L=1;inputterms=0
 for rel in files:
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for line in f:e=json.loads(line);L=math.lcm(L,Q(e['QQ_content']).denominator)
 P={};source_metadata=[]
 for rel in files:
  with(ROOT/rel).open()as f:
   for ln,line in enumerate(f):
    e=json.loads(line);b=e.get('block','D');i=e['block_row']if b!='D'else e['residual_row'];j=e['block_column']if b!='D'else e['residual_column'];c=Q(e['QQ_content']);fac=c.numerator*(L//c.denominator);S=0;N=0
    for g,v in e['sparse_integer_polynomial']:S|=1<<int(g);N+=abs(int(v)*fac)
    P[b,i,j]=(S,N);inputterms+=len(e['sparse_integer_polynomial']);source_metadata.append(dict(block=b,row=i,column=j,path=rel,line=ln,entry_index=e['entry_index'],content=e['QQ_content'],cleared_scale=str(fac)))
 assert inputterms==23932 and all(P['A',i,j][0]==0 for i in range(2)for j in range(2,4))
 products=[]
 def add(*ps):return(math.prod([1])and __import__('functools').reduce(int.__or__,(p[0]for p in ps),0),sum(p[1]for p in ps))
 def mul(a,b,label):
  S,N=a;T,M=b
  if S.bit_count()>T.bit_count():S,T=T,S
  U=0
  for g in support_list(S):U|=T<<g
  products.append(dict(label=label,actual_support_upper=U.bit_count(),exponent_max=U.bit_length()-1,coefficient_abs_bound=(N*M).bit_length(),naive_integer_pair_upper=S.bit_count()*T.bit_count()))
  assert time.monotonic()-start<45
  return U,N*M
 def scale(p,k):return p[0],abs(k)*p[1]
 dp=add(mul(P['A',0,0],P['A',1,1],'deltaP_plus'),mul(P['A',0,1],P['A',1,0],'deltaP_minus'));dq=add(mul(P['A',2,2],P['A',3,3],'deltaQ_plus'),mul(P['A',2,3],P['A',3,2],'deltaQ_minus'))
 adjP=[[P['A',1,1],P['A',0,1]],[P['A',1,0],P['A',0,0]]];adjQ=[[P['A',3,3],P['A',2,3]],[P['A',3,2],P['A',2,2]]]
 dc=[add(P['C',2,j],scale(P['C',0,j],kappa))for j in range(4)]
 v=[add(*(mul(dc[2+k],adjQ[k][j],f'v{j}term{k}')for k in range(2)))for j in range(2)]
 u=[add(mul(dq,dc[j],f'u{j}_left'),*(mul(v[k],P['A',2+k,j],f'u{j}_cross{k}')for k in range(2)))for j in range(2)]
 w=[add(*(mul(u[k],adjP[k][j],f'w{j}term{k}')for k in range(2)))for j in range(2)]
 vp=[mul(dp,v[j],f'deltaP_v{j}')for j in range(2)];delta=mul(dp,dq,'deltaP_deltaQ');outputs=[];rowmax=0;outslots=0;raw=0
 for j in range(7):
  dd=add(P['D',2,j],scale(P['D',0,j],kappa));terms=[mul(delta,dd,f'output{j}_deltaD')]+[mul(w[k],P['B',k,j],f'output{j}_top{k}')for k in range(2)]+[mul(vp[k],P['B',2+k,j],f'output{j}_bottom{k}')for k in range(2)];q=add(*terms);n=q[0].bit_count();bits=q[1].bit_length();rowmax=max(rowmax,n);outslots+=n;raw+=n*(80+(bits+2)//3);outputs.append(dict(column=j,source='unweighted S_row2-59375*S_row0',support_upper=n,max_exponent=q[0].bit_length()-1,integer_coefficient_bits_upper=bits,integer_abs_upper=str(q[1]),necessary_source_terms=5,
  output_support_hex=hex(q[0])))
 prefix=[dp,dq]+v+u+w+vp+[delta]+dc;retained=sum(q[0].bit_count()for q in prefix);largest_product=max((p['actual_support_upper']for p in products),default=0);owned=2*inputterms+retained+3*rowmax+2*largest_product+15000;allRAW=raw+1500000;naive=sum(p['naive_integer_pair_upper']for p in products)
 dump(out/'whole-row-source-and-operator.json',dict(kappa=59375,source_rows='unweightedR1BQQ rows4+2 minus59375*(4+0)',normalized_prefix_evidence='Trow2=59375*t68*Trow1 through550; alpha0-alpha2=68; R2 rowsnonpivotunit',
  input_source=source_metadata,uniform_integer_denominator=str(L),Ablock='A=[[P,0],[A21,Q]]',deltaP='detP',deltaQ='detQ',adjugates='signed2x2classicaladjugate',
  dc='C2-59375*C0',dd='D2-59375*D0',v='dc_right*adjQ',u='deltaQ*dc_left-v*A21',w='u*adjP',vp='deltaP*v',
  seven_actual_targets='N_j=deltaP*deltaQ*dd_j-w*Btop_j-vp*Bbottom_j; commonintegerdenominatorL^5',
  exact_identity='N=deltaP*deltaQ*(S2-59375*S0); A0=I impliesdeltaPdeltaQ nonzeropolynomial',
  completion_if_ALL7_ZERO='wholeunivariate11detZERO byexact rowdependence; mixedradix injectiveonacceptedwholeNewton=>fixedoriginalminorZERO',
  if_first_nonzero='literalglobalrowidentityfails; no determinantNONZERO/rank conclusion; stopallfurthercolumns',no_unknown_parameter_inverse=True,no_new_weight_precisionwindow_point_prime_source572=True,outputs=outputs))
 dump(out/'resource-gate.json',dict(input_terms=inputterms,two_inputcopies=2*inputterms,prefix_coefficient_support_upper=retained,maximum_output_onecolumn=rowmax,total_all7_output_support_upper=outslots,largest_intermediate_product_support_upper=largest_product,
  whole_owned_upper=owned,new_RAW_upper=allRAW,all_product_operations=products,naive_integer_multiplication_pairs_upper=naive,
  proposed_job='ONE whole-row7identity job, lexcolumns0..6, STOPfirstnonzero; finitepolynomial block2+2adjugate, no fieldGaussian/inverse or series',hard_seconds=60,hard_address_bytes=2000000000,hard_scalar_terms=150000,hard_total_RAW=33554432))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins)
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_new_signed_target_coefficients=True,whole_row_candidate_kappa=59375,whole_owned_upper=owned,new_RAW_upper=allRAW,naive_integer_pairs_upper=naive,
  preflight_fits=owned<150000 and allRAW+external<33554432,NO_whole_row_identity_or_det11_ZERO_NONZERO_claim=True,total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
