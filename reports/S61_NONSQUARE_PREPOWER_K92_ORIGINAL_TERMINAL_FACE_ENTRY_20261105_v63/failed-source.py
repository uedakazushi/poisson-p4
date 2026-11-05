#!/usr/bin/env python3
"""ONE tinyconstant original-terminal-face det/rank; run onlyunderROOTbudget."""
import argparse,hashlib,importlib.util,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_ORIGINAL_TERMINAL_FACE_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_ORIGINAL_TERMINAL_FACE_SOURCE_20261105_v63';r=src(base+'/report.json');assert r['all_passed']and pins[base+'/report.json']=='36f5f68f7da9aa9f6954c7a6b85b7569e0d5d7da7a9fa490d8afecd2f7ec361c';plan=src(base+'/terminal-face-and-normalization.json');nb=src(base+'/terminal-Newton-filter.json');g=src(base+'/resource-gate.json');assert nb['inside_all_Newton_halfspaces']and nb['radix_candidate']==[5,0,0];A=plan['terminal_integer_face_matrix'];assert len(A)==11 and all(len(row)==11 for row in A)
 minorbits=g['integer_det_bit_upper'];productbits=2*minorbits;temporarybits=3*minorbits+1;seen_raw=0;seen_final=0
 def chk(v,b):
  nonlocal seen_raw
  seen_raw=max(seen_raw,abs(v).bit_length());assert abs(v).bit_length()<=b;return v
 def cq(v):
  nonlocal seen_final
  v=Q(v);seen_final=max(seen_final,abs(v.numerator).bit_length(),v.denominator.bit_length());assert abs(v.numerator).bit_length()<=minorbits and v.denominator.bit_length()<=minorbits;return v
 def qm(x,y):
  x=Q(x);y=Q(y);return Q(chk(x.numerator*y.numerator,productbits),chk(x.denominator*y.denominator,productbits))
 def qs(x,y):
  x=Q(x);y=Q(y);a=chk(x.numerator*y.denominator,temporarybits);b=chk(y.numerator*x.denominator,temporarybits);return cq(Q(chk(a-b,temporarybits),chk(x.denominator*y.denominator,temporarybits)))
 def checked_rank(mat):
  m=len(mat);n=len(mat[0]);B=[list(map(Q,row))for row in mat];P=[[Q(i==j)for j in range(m)]for i in range(m)];labels=list(range(m));cols=[];rows=[];ops=[];r=0
  for j in range(n):
   k=next((i for i in range(r,m)if B[i][j]),None)
   if k is None:continue
   if k!=r:B[k],B[r]=B[r],B[k];P[k],P[r]=P[r],P[k];labels[k],labels[r]=labels[r],labels[k];ops.append(['swap',k,r])
   inv=cq(1/B[r][j]);B[r]=[cq(qm(x,inv))for x in B[r]];P[r]=[cq(qm(x,inv))for x in P[r]];ops.append(['scale',r,str(inv)])
   for i in range(m):
    if i!=r and B[i][j]:
     c=B[i][j];B[i]=[qs(x,qm(c,y))for x,y in zip(B[i],B[r])];P[i]=[qs(x,qm(c,y))for x,y in zip(P[i],P[r])];ops.append(['add',i,r,str(-c)])
   rows.append(labels[r]);cols.append(j);r+=1
   if r==m:break
  return B,P,cols,rows,ops
 B=[row[:]for row in A];previous=1;sign=1;steps=[];swaps=[];D=None
 for k in range(10):
  j=next((j for j in range(k,11)if B[j][k]),None)
  if j is None:D=0;steps.append(dict(stage=k,remaining_column_ZERO=True));break
  if j!=k:B[j],B[k]=B[k],B[j];sign=-sign;swaps.append([k,j])
  pivot=B[k][k];cellchecks=[]
  for i in range(k+1,11):
   for j in range(k+1,11):
    x=chk(pivot*B[i][j],productbits);y=chk(B[i][k]*B[k][j],productbits);num=chk(x-y,productbits+1);assert num%previous==0;value=num//previous;chk(value,minorbits);cellchecks.append(dict(row=i,column=j,pivot=str(pivot),oldij=str(B[i][j]),oldik=str(B[i][k]),oldkj=str(B[k][j]),divisor=str(previous),quotient=str(value)));B[i][j]=value
  for i in range(k+1,11):B[i][k]=0
  steps.append(dict(stage=k,pivot=str(pivot),previous=str(previous),all_exact_divisions=cellchecks));previous=pivot
 if D is None:D=sign*B[10][10]
 hp=ROOT/'code/s61_nonsquare_prePOWER_K92_leading_structure_20261105_v63.py';pins[str(hp.relative_to(ROOT))]=sha(hp);sp=importlib.util.spec_from_file_location('lin',hp);lin=importlib.util.module_from_spec(sp);sp.loader.exec_module(lin);rf,R,pc,pr,ops=checked_rank(A);rank=len(pc);free=[j for j in range(11)if j not in pc];kernel=[[Q(0)for j in free]for i in range(11)]
 for k,j in enumerate(free):
  kernel[j][k]=Q(1)
  for i,c in enumerate(pc):kernel[c][k]=-rf[i][j]
 assert bool(D)==(rank==11)
 minor=[[Q(A[i][j])for j in pc]for i in pr];md=lin.determinant(minor)if rank else Q(1);mi=checked_rank(minor)[1]if rank else[];assert md
 factors=[[Q(i==pc.index(j))if j in pc else rf[i][j]for j in range(11)]for i in range(rank)];FD=math.lcm(*(v.denominator for row in factors for v in row));assert FD.bit_length()<=minorbits;FC=[[int(FD*v)for v in row]for row in factors];basis=[[A[i][j]for j in pc]for i in range(11)];assert all(sum(basis[i][k]*FC[k][j]for k in range(rank))==FD*A[i][j]for i in range(11)for j in range(11))
 assert all(sum(A[i][k]*int(FD*kernel[k][j])for k in range(11))==0 for i in range(11)for j in range(len(free)))
 faceQQ=[[Q(v,d)for v in row]for d,row in zip(plan['row_integer_denominator'],A)];assert faceQQ==[[Q(v)for v in row]for row in plan['terminal_QQ_face_matrix']]
 coefficient=D*Q(plan['terminal_integer_det_to_original_coefficient_factor'])
 dump(out/'all-Bareiss-terminal-face-identities.json',dict(original_integer_face=A,row_swaps=swaps,all_steps=steps,exact_integer_det=str(D),all_divisions_verified=True,only_fixed_constant_face=True))
 dump(out/'terminal-face-rank-and-source-coefficient.json',dict(exact_rank=rank,right_kernel=kernel,row_transform=R,row_reduced=rf,row_operations=ops,pivot_rows=pr,pivot_columns=pc,
  unit_minor=minor,unit_minor_inverse=mi,unit_minor_determinant=str(md),all_column_factorization=dict(pivot_column_basis=basis,coefficient_matrix=factors,common_denominator=str(FD),cleared_integer_coefficient_matrix=FC,all121_integer_reproducts_verified=True),original_terminal_QQ_face=faceQQ,original_detM_a5e0d0_coefficient=str(coefficient),original_fixed_det11_coefficient='a7e0d0 withsameQQcoefficient since det11=a²detM',
  NO_otherweight_window_point_prime_or_source572=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new<g['new_RAW_upper']and new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,ONE_constant_terminal_face_job=True,exact_rank=rank,exact_integer_det=str(D),original_detM_a5e0d0_coefficient=str(coefficient),coefficient_NONZERO=bool(coefficient),
  all_fixed_Bareiss_divisions_rank_kernel_unitminor_verified=True,all11_columns_factored=True,normalized_minor_bit_bound=minorbits,uncancelled_product_bit_bound=productbits,uncancelled_rational_temporary_bound=temporarybits,actual_uncancelled_max_bits=seen_raw,actual_normalized_fraction_max_bits=seen_final,no_other_coefficient_or_global57_claim=True,positive_if_any_independent_peer_pending=bool(coefficient),
  total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
