#!/usr/bin/env python3
"""Tiny authorized constant rank; SOURCE ONLY one full-row preconditioning gate."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_LEADING_STRUCTURE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def eye(n):return[[Q(i==j)for j in range(n)]for i in range(n)]
def mm(A,B):return[[sum((a*b for a,b in zip(row,col)),Q(0))for col in zip(*B)]for row in A]
def rank_ops(A):
 B=[list(map(Q,row))for row in A];P=eye(len(A));labels=list(range(len(A)));cols=[];pivrows=[];ops=[];r=0
 for j in range(len(A[0])):
  i=next((i for i in range(r,len(A))if B[i][j]),None)
  if i is None:continue
  if i!=r:B[i],B[r]=B[r],B[i];P[i],P[r]=P[r],P[i];labels[i],labels[r]=labels[r],labels[i];ops.append(['swap',i,r])
  a=B[r][j];B[r]=[v/a for v in B[r]];P[r]=[v/a for v in P[r]];ops.append(['scale',r,str(1/a)])
  for i in range(len(A)):
   if i!=r and B[i][j]:
    a=B[i][j];B[i]=[v-a*w for v,w in zip(B[i],B[r])];P[i]=[v-a*w for v,w in zip(P[i],P[r])];ops.append(['add',i,r,str(-a)])
  pivrows.append(labels[r]);cols.append(j);r+=1
  if r==len(A):break
 assert mm(P,A)==B
 return B,P,cols,pivrows,ops
def inverse(A):
 B,P,c,rows,ops=rank_ops(A);assert B==eye(len(A));assert mm(A,P)==eye(len(A));return P
def determinant(A):
 B=[list(map(Q,r))for r in A];d=Q(1)
 for j in range(len(A)):
  i=next((i for i in range(j,len(A))if B[i][j]),None);assert i is not None
  if i!=j:B[i],B[j]=B[j],B[i];d=-d
  p=B[j][j];d*=p
  for i in range(j+1,len(A)):
   a=B[i][j]/p
   for k in range(j+1,len(A)):B[i][k]-=a*B[j][k]
 return d
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_TROPICAL_PREFIX32_V2_20261105_v63';r=src(base+'/report.json');assert r['all_passed']and r['all33_zero'];selected=src(base+'/selected78-source.json');assert selected['weight']==[1,2244,68];terms=selected['source_terms'];assert len(terms)==78 and all(t['deficit']==0 for t in terms)
 A=[[Q(0)for j in range(11)]for i in range(11)]
 for t in terms:A[t['row']][t['column']]+=Q(t['cleared_integer'])
 B,P,pivotcols,pivotrows,ops=rank_ops(A);rank=len(pivotcols);free=[j for j in range(11)if j not in pivotcols];kernel=[[Q(0)for j in free]for i in range(11)]
 for k,j in enumerate(free):
  kernel[j][k]=Q(1)
  for i,c in enumerate(pivotcols):kernel[c][k]=-B[i][j]
 assert mm(A,kernel)==[[Q(0)for j in free]for i in range(11)]
 minor=[[A[i][j]for j in pivotcols]for i in pivotrows];minv=inverse(minor);mdet=determinant(minor);assert mdet
 order=pivotcols+free;Qmat=[[Q(0)for j in range(11)]for i in range(11)]
 for j,c in enumerate(order):Qmat[c][j]=Q(1)
 U=eye(11)
 for i in range(rank):
  for j,c in enumerate(free):U[i][rank+j]=-B[i][c]
 Qmat=mm(Qmat,U);diag=[[Q(i==j and i<rank)for j in range(11)]for i in range(11)];assert mm(mm(P,A),Qmat)==diag
 assert all(A[1][j]+2*A[2][j]==0 for j in range(11));assert selected['row_LCMs'][1]==selected['row_LCMs'][2]and selected['row_potentials'][1]==selected['row_potentials'][2]
 dump(out/'constant-rank-certificate.json',dict(matrix=A,rank=rank,right_kernel=kernel,pivot_columns=pivotcols,pivot_source_rows=pivotrows,free_columns=free,row_operations=ops,row_transform=P,column_transform=Qmat,R_A_Q=diag,unit_minor=minor,unit_minor_inverse=minv,unit_minor_determinant=str(mdet),
  row1_plus2row2_zero_leading=True,row1and2_LCM_equal=selected['row_LCMs'][1],row1and2_potential_equal=selected['row_potentials'][1],NO_further_detseries=True))
 helper=ROOT/'code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py';pins[str(helper.relative_to(ROOT))]=sha(helper);sp=importlib.util.spec_from_file_location('sourcegate',helper);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
 binding=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02';M,inventory,kept,skipped=h.selected_inputs(state,binding['actual_pure_matrix_roots'])
 w=selected['weight'];rp=selected['row_potentials'];cp=selected['column_potentials'];profiles=[];output_terms=0;raw=0;maxbits=0
 for j in range(11):
  a=M[1][j];b=M[2][j];sa={m for m,v in a['co']};sb={m for m,v in b['co']};union=sa|sb;removed={m for m in union if rp[1]+cp[j]-sum(x*y for x,y in zip(m,w))==0};upper=union-removed;D=math.lcm(a['c'].denominator,b['c'].denominator)
  na=max((abs(v)for m,v in a['co']),default=0)*a['c'].numerator*(D//a['c'].denominator);nb=2*max((abs(v)for m,v in b['co']),default=0)*b['c'].numerator*(D//b['c'].denominator);bits=(na+nb).bit_length();maxbits=max(maxbits,bits);output_terms+=len(upper);raw+=len(upper)*(65+len(str(na+nb))+len(str(D)))
  profiles.append(dict(column=j,source_polynomial_ids=[a['state_id'],b['state_id']],input_terms=[len(a['co']),len(b['co'])],union_terms=len(union),known_leading_zero_monomials=sorted(removed),output_support_upper=len(upper),QQ_denominator=str(D),integer_numerator_bit_upper=bits,
   next_possible_weight_upper=max((sum(x*y for x,y in zip(m,w))for m in upper),default=None),full_coefficient_sum_NOT_computed=True))
 strip=src('reports/S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63/report.json');assert strip['all_passed'];peer=src('reports/S61_NONSQUARE_PREPOWER_H22_FIXED_MINORS_20261104_v63_REVIEW/report.json');assert peer['all_passed'];pure4=src('reports/S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63/native.json');assert pure4['remaining_polynomial']
 dump(out/'one-row-source-gate.json',dict(operation='replace row1 by row1+2row2; all other rows unchanged',actual_rows='normalized G18/a and G39/a, current15 rows1 and3, compact14 rows0 and2',
  potential_alignment='same rowpot and same rowLCM, so NO extra t-monomial/scalar beyond literal +2',determinant_preserved=True,H3_on_bothrows_literalzero=True,columns_original32=[22,11,13,25,15,19,14,18,28,30,29],
  input_coefficient_terms=kept,output_row_support_upper=output_terms,all_input_plus_output_plus15000=kept+output_terms+15000,max_new_numerator_bits=maxbits,artifact_RAW_upper=raw+500000,
  columns=profiles,entire_full_row_zero_NOT_assumed=True,existing_nonzero_P4_precludes_identicallyzero_full_row='P4 originaldet=a^4*C^3*R434; normalizeddet=a^2*C^3*R434; this unit rowop preserves it. Thus at least one of first4 columns is nonzero as polynomial.',
  future_job='ONE11entry exactscalaradd/sourceidentity then samefixedw support-retropicalization; NO newseries/det/otherweight authorized',source572_NOT_replayed=True,NO_polynomial_products=True))
 dump(out/'typed-Schur-recipe.json',dict(constant_rank=rank,residual_size=11-rank,source='B(t)=diag(rowLCM)*diag(t^rowpot)*actualM(t^-w)*diag(t^colpot)',
  exact_identity='R*B(0)*Q=diag(I_rank,0)',block_A='upper rank-square block of R*B(t)*Q; A(0)=I',block_B='upper right',block_C='lower left',block_D='lower right',
  formal_unit_inverse='A^-1=adj(A)/det(A), detA(0)=1',residual='S=D-C*A^-1*B',determinant_identity='detB=(detR*detQ)^-1*detA*detS',
  construction_uses_only_constant_QQ_inverse=True,NO_transformed_polynomial_or_formalseries_coefficients_computed=True,actual_source_constant_unit_minor=str(mdet)))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);total=sum(p.stat().st_size for p in out.iterdir());assert kept+output_terms+15000<150000 and raw+500000<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,leading_constant_exact_rank=rank,right_kernel_dimension=11-rank,unit_minor_nonzero=True,all_constant_R_A_Q_and_kernel_identities_verified=True,NO_further_detseries_or_fullpolyrowadd=True,
  one_row_output_support_upper=output_terms,one_row_live_upper=kept+output_terms+15000,one_row_RAW_upper=raw+500000,full_row_nonzero_from_existingP4_source=True,no_det11_nonzero_or_global57_claim=True,
  total_RAW_before_report=total,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
