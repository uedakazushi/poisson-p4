#!/usr/bin/env python3
"""ONE fixed constant leading face: exact rank/minor/kernel/all column factors."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import json,hashlib,time,signal,resource,gc
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_TIED_EXACT_20261105_v64';ADD='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_TIED_EXACT_STORAGE_GATE_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def bits(n):return abs(n).bit_length()
def countq(A):return 2*sum(len(row)for row in A)
class Arithmetic:
 def __init__(self,H):self.H=H;self.uncancelled=0;self.stored=0
 def bounded(self,n,d,limit):
  b=max(bits(n),bits(d));self.uncancelled=max(self.uncancelled,b);assert d and b<=limit;return Q(n,d)
 def div(self,a,b):assert b;return self.bounded(a.numerator*b.denominator,a.denominator*b.numerator,2*self.H+2)
 def mul(self,a,b):return self.bounded(a.numerator*b.numerator,a.denominator*b.denominator,2*self.H+2)
 def sub(self,a,b):return self.bounded(a.numerator*b.denominator-b.numerator*a.denominator,a.denominator*b.denominator,3*self.H+4)
 def store(self,q):b=max(bits(q.numerator),bits(q.denominator));self.stored=max(self.stored,b);assert b<=self.H;return q

def det_bareiss(M,H):
 B=[row[:]for row in M];n=len(B);prev=1;sign=1;piv=[];peak=0
 if n==0:return 1,piv,peak
 for k in range(n-1):
  p=next((i for i in range(k,n)if B[i][k]),None);assert p is not None,'selected rank minor unexpectedly singular'
  if p!=k:B[k],B[p]=B[p],B[k];sign=-sign
  q=B[k][k];piv.append(str(q))
  for i in range(k+1,n):
   for j in range(k+1,n):
    a=B[i][j]*q;b=B[i][k]*B[k][j];z=a-b;peak=max(peak,bits(a),bits(b),bits(z));assert peak<=2*H+2 and z%prev==0;B[i][j]=z//prev;assert bits(B[i][j])<=H
   B[i][k]=0
  prev=q
 return sign*B[-1][-1],piv,peak

def verify_integer_dot(Z,coeff,target,H):
 d=1
 for q in coeff:d=d*q.denominator//gcd(d,q.denominator)
 assert bits(d)<=H;iv=[q.numerator*(d//q.denominator)for q in coeff];assert max(map(bits,iv),default=0)<=2*H+1
 for row,want in zip(Z,target):
  prods=[a*b for a,b in zip(row,iv)];s=sum(prods);assert max([bits(x)for x in prods]+[bits(s)])<=2*H+80;assert s==want*d
 return str(d)

def main():
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 ap=ROOT/f'reports/{ADD}/report.json';a=read(ap);assert a['all_passed']and a['ONE_constant_job_safe'];bp=ap.parent/'implementation-storage-addendum.json';assert sha(bp)==a['evidence_sha256'][bp.name];b=read(bp)
 for p,h in b['actual_READ_input_sha256'].items():assert sha(ROOT/p)==h
 assert sha(Path(__file__))==b['worker_sha256'];prepared=ap.stat().st_size+bp.stat().st_size;assert prepared<=65536
 ev=read(ROOT/b['source_face']);Z=[[int(x)for x in row]for row in ev['integer_cleared_leading_matrix24']];D=list(ev['row_LCM24']);H=ev['Hadamard_all_minor_integer_bit_upper'];assert H==1618 and max(bits(x)for row in Z for x in row)==64;del ev;gc.collect()
 A=[[Q(x)for x in row]for row in Z];permutation=list(range(24));pcols=[];r=0;ar=Arithmetic(H);peak=25000+576+countq(A)+200;assert peak<=29608
 events=out/'exact-fixed-row-operations.jsonl'
 with events.open('wb')as f:
  for col in range(24):
   p=next((i for i in range(r,24)if A[i][col]),None)
   if p is None:continue
   if p!=r:A[p],A[r]=A[r],A[p];permutation[p],permutation[r]=permutation[r],permutation[p];f.write(enc(['swap',p,r]))
   pivot=A[r][col];f.write(enc(['divide_row',r,str(pivot)]));A[r]=[ar.store(ar.div(q,pivot))for q in A[r]]
   for i in range(24):
    if i==r or not A[i][col]:continue
    factor=A[i][col];f.write(enc(['subtract_row_multiple',i,r,str(factor)]))
    for j in range(24):A[i][j]=ar.store(ar.sub(A[i][j],ar.mul(factor,A[r][j])))
   pcols.append(col);r+=1
   if r==24:break
 assert all(not q for row in A[r:]for q in row);rows=permutation[:r];C=[row[:]for row in A[:r]];peak=max(peak,25000+576+countq(A)+countq(C)+200);assert peak<=29608
 # Every original column is exactly in the selected original column span.
 basis=[[Z[i][j]for j in pcols]for i in range(24)];column_LCD=[]
 for j in range(24):column_LCD.append(verify_integer_dot(basis,[C[i][j]for i in range(r)],[Z[i][j]for i in range(24)],H))
 minor,minor_pivots,ffbits=det_bareiss([[Z[i][j]for j in pcols]for i in rows],H);assert minor;peak=max(peak,25000+576+countq(A)+countq(C)+r*r+200);assert peak<=29608
 free=[j for j in range(24)if j not in pcols];kernels=[];kernel_LCD=[]
 for j in free:
  v=[Q(0)]*24;v[j]=Q(1)
  for i,k in enumerate(pcols):v[k]=-C[i][j]
  kernel_LCD.append(verify_integer_dot(Z,v,[0]*24,H));kernels.append({'free_column':j,'coordinates24':[str(q)for q in v]})
 # Kernel columns have identity on all free rows, so they are independent.
 assert all(Q(z['coordinates24'][j])==int(j==z['free_column'])for z in kernels for j in free)
 factor={'all_passed':True,'rank':r,'selected_basis_column_positions':pcols,'selected_nonzero_minor_row_positions':rows,'selected_nonzero_minor_column_positions':pcols,'selected_minor_determinant':str(minor),'selected_minor_Bareiss_pivots':minor_pivots,'all24_column_factors_r_by24':[[str(q)for q in row]for row in C],'all24_integer_dot_identity_LCDs':column_LCD,'kernel_basis':kernels,'kernel_integer_dot_identity_LCDs':kernel_LCD,'all_column_factor_identities_checked':True,'all_kernel_dots_checked':True,'kernel_independent_free_identity':True,'rank_upper_from_all_column_factors_and_lower_from_minor':True,'only_this_fixed_leading_face':True}
 positive=r==24;actual_det=None
 if positive:
  # Replay the frozen literal row program on I only after full rank is known.
  del A,C;gc.collect();U=[[Q(int(i==j))for j in range(24)]for i in range(24)]
  with events.open()as f:
   for line in f:
    op=json.loads(line)
    if op[0]=='swap':_,i,j=op;U[i],U[j]=U[j],U[i]
    elif op[0]=='divide_row':_,i,q=op;q=Q(q);U[i]=[ar.store(ar.div(x,q))for x in U[i]]
    else:
     _,i,j,q=op;q=Q(q)
     for k in range(24):U[i][k]=ar.store(ar.sub(U[i][k],ar.mul(q,U[j][k])))
  for i in range(24):verify_integer_dot([list(row)for row in zip(*Z)],U[i],[int(i==j)for j in range(24)],H)
  inversions=sum(rows[i]>rows[j]for i in range(24)for j in range(i+1,24));actual_det=minor*(-1)**inversions
  factor.update({'full_integer_determinant':str(actual_det),'constant_inverse_of_integer_cleared_face':[[str(q)for q in row]for row in U],'actual_fixed24_minor_nonzero':True})
  # Typed fixed finite program, not expanded scalar inverse or Gram recipe.
  assignment=read(ROOT/b['assignment_data']);leads=assignment['all_actual_leading_entries'];u=assignment['row_dual24'];v=assignment['column_dual24'];refmap=[[None]*24 for _ in range(24)]
  for z in leads:refmap[z['i']][z['j']]=z['source']
  program={'all_passed':True,'actual_entry_source_matrix24':refmap,'fixed_row_order':assignment['fixed_rows24'],'fixed_column_order':assignment['fixed_columns24'],'row_t0_exponents24':[-x for x in u],'column_t0_exponents24':[-x for x in v],'constant_preconditioner':'Z_inverse * diag(row_LCM24), using exact savedU and sourceD','actual_initial_RF_matrix':'G0=constant_preconditioner * diag(t0^(-u)) * ActualA * diag(t0^(-v)); every actual source entry is retained.','all24_actual_pivot_RF_definitions':[{'k':k,'pivot':'G[k][k,k]','G_next':'for i,j>k: G[k+1][i,j]=G[k][i,j]-G[k][i,k]*G[k][k,j]/G[k][k,k]','nonzero':'Formal weighted leading matrix after exact constant preconditioner isI; all prefix minors have leading1, hence every actual K5 pivot is nonzero.'}for k in range(24)],'RHS_forward':'h0=constant_preconditioner*diag(t0^(-u))*selected_rhs; h[k+1]_i=h[k]_i-G[k]_ik*h[k]_k/pivot_k. Backsolve the final upper triangular actualRF system; actual_unknown=diag(t0^(-v))*solution.','normal42_transport':'Frozen49 unit RHS/backsubstitution, ten constraint particular RHS solution, selected24 solve, undo one P290 column scale by old_x_free=P290*y_free, ten pivot restoration h_part+N_alias*y_free, then49 restoration. Return original83 normal columns41..82.','all87_other_top_RHS_and66_zero_constraint_RHS_kept':True,'cokernel':'full177 coker=top111 coker direct_sumK5^66. No zero summands removed.','chart':'Old base monomial units*P290*detA*detB and the finite set of proven nonzero actualRF pivot numerators/denominators. Complements remain open.','no_scalar_inverse_expansion_or_all111_claim':True}
  (out/'actual-fixed24-finite-RF-inverse-program.json').write_bytes(enc(program))
 else:factor.update({'full_integer_determinant':'0','actual_fixed24_determinant_still_unresolved':True,'no_adjugates_computed':True,'next_SOURCE_only':'Compare complete cached actualrows11/44 for literalRF dependence, not just the leading face. No new row/weight search.'})
 (out/'fixed24-exact-rank-minor-kernel-column-factors.json').write_bytes(enc(factor));assert ar.uncancelled<=3*H+4 and ffbits<=2*H+2
 readbase=b['actual_READ_input_RAW_bytes']+prepared;raw=readbase+sum(p.stat().st_size for p in out.rglob('*')if p.is_file());assert raw<=b['fresh_RAW_upper']
 report={'id':ID,'all_passed':True,'rank_of_this_fixed_leading_face':r,'determinant_nonzero':positive,'all24_column_factor_identities':True,'kernel_dimension':24-r,'nonzero_rank_minor_verified':True,'actual_fixed24_minor_nonzero_proved':positive,'maximum_stored_rational_bits':ar.stored,'maximum_uncancelled_rational_bits':ar.uncancelled,'allowed_uncancelled_rational_bits':3*H+4,'maximum_fraction_free_temporary_bits':ffbits,'whole_owned_numeric_upper':29608,'checked_owned_numeric_peak':peak,'actual_fresh_RAW_before_small_report':raw,'no_new_source_values_rowset_weight_or_points':True,'input_sha256':{**b['actual_READ_input_sha256'],str(ap.relative_to(ROOT)):sha(ap),str(bp.relative_to(ROOT)):sha(bp)},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));assert readbase+sum(p.stat().st_size for p in out.rglob('*')if p.is_file())<=b['fresh_RAW_upper'];signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
