#!/usr/bin/env python3
"""ONE32factor/65product leading rank profile; then typed residual SOURCE only."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63';W=(1,2244,68)
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 def mod(rel,name):
  p=ROOT/rel;pins[rel]=sha(p);s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
 gate='reports/S61_NONSQUARE_PREPOWER_K92_CORRECTION_FACTOR_SOURCE_V2_20261105_v63';rpt=src(gate+'/report.json');assert rpt['all_passed']and pins[gate+'/report.json']=='aa6a12c4910285b94caf96446ab801c80e5a753d5af362b290a05caebb4b866b';plan=src(gate+'/next32-factor-profile-source.json');proof=src(gate+'/rank-support-zero-proof.json');assert len(plan['physical_atoms'])==37 and len(plan['requests32'])==32
 base='reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63';rt=src(base+'/retropicalized-source.json');row=src(base+'/actual-row11.json');R=[[Q(v)for v in rr]for rr in rt['row_transform']];C=[[Q(v)for v in rr]for rr in rt['column_transform']];DR=math.lcm(*(v.denominator for rr in R for v in rr));DQ=math.lcm(*(v.denominator for rr in C for v in rr));DD=DR*DQ;assert DD==int(plan['common_linear_denominator']);Ds=list(map(int,rt['constant_clearing_rowLCMs']))
 helper=mod('code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py','sourcegate');lin=mod('code/s61_nonsquare_prePOWER_K92_leading_structure_20261105_v63.py','constantQQ');bind=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02';M,inventory,kept,skip=helper.selected_inputs(state,bind['actual_pure_matrix_roots'])
 for p in row['columns']:M[1][p['column']]=dict(state_id=-1001-p['column'],c=Q(p['QQ_common_factor']),co=[(tuple(m),int(v))for m,v in p['primitive_integer_polynomial']])
 atoms={}
 for a in plan['physical_atoms']:
  p=M[a['row']][a['column']];assert p['state_id']==a['state_polynomial_id'];mon,z=p['co'][a['source_index']];assert mon==tuple(a['monomial'])and z==int(a['source_integer'])and p['c']==Q(a['source_QQfactor']);gap=rt['row_potentials'][a['row']]+rt['column_potentials'][a['column']]-sum(u*v for u,v in zip(mon,W));assert gap==a['gap'];v=p['c']*z*Ds[a['row']];assert v.denominator==1 and v.numerator==int(a['cleared_integer']);atoms[a['id']]=v.numerator
 FZ=[[0]*4 for i in range(7)];BZ=[[0]*7 for i in range(4)];entries=[];products=0
 for q in plan['requests32']:
  u=q['global_row'];v=q['global_column'];total=0
  for c in q['contributors']:
   a=plan['physical_atoms'][c['physical_id']];assert a['id']==c['physical_id']and a['gap']==q['gap'];rr=int(c['R_integer']);cc=int(c['Q_integer']);assert R[u][a['row']]*DR==rr and C[a['column']][v]*DQ==cc;total+=rr*atoms[a['id']]*cc;products+=2;assert abs(total)<=int(q['integer_numerator_abs_bound'])
  (FZ if q['kind']=='F'else BZ)[q['row']][q['column']]=total;entries.append(dict(q,exact_integer_numerator=str(total),QQ_value=str(Q(total,DD))))
 assert products==188;DZ=plan['known_D_leading_integer_numerators'];PZ=[[0]*7 for i in range(7)]
 for i,m,j in plan['F0B0_products65']:PZ[i][j]+=FZ[i][m]*BZ[m][j]
 HZ=[[DD*DZ[i][j]-PZ[i][j]for j in range(7)]for i in range(7)];assert all(abs(v).bit_length()<=plan['H_integer_bit_upper']for rr in HZ for v in rr)
 DQmat=[[Q(v,DD)for v in rr]for rr in DZ];FQ=[[Q(v,DD)for v in rr]for rr in FZ];BQ=[[Q(v,DD)for v in rr]for rr in BZ];HQ=[[Q(v,DD**2)for v in rr]for rr in HZ];FP=lin.mm(FQ,BQ);assert HQ==[[DQmat[i][j]-FP[i][j]for j in range(7)]for i in range(7)]
 rf,P,pcols,prows,ops=lin.rank_ops(HZ);rank=len(pcols);assert rank<=6;free=[j for j in range(7)if j not in pcols];N=[[Q(0)for j in free]for i in range(7)]
 for k,j in enumerate(free):
  N[j][k]=Q(1)
  for i,c in enumerate(pcols):N[c][k]=-rf[i][j]
 assert lin.mm(HZ,N)==[[Q(0)for j in free]for i in range(7)]
 order=pcols+free;Q2=[[Q(0)for j in range(7)]for i in range(7)]
 for j,c in enumerate(order):Q2[c][j]=Q(1)
 U=lin.eye(7)
 for i in range(rank):
  for j,c in enumerate(free):U[i][rank+j]=-rf[i][c]
 Q2=lin.mm(Q2,U);R2=[[v*DD**2 for v in rr]for rr in P];diag=[[Q(i==j and i<rank)for j in range(7)]for i in range(7)];assert lin.mm(lin.mm(R2,HQ),Q2)==diag
 minor=[[Q(HZ[i][j])for j in pcols]for i in prows];mdet=lin.determinant(minor)if rank else Q(1);minv=lin.inverse(minor)if rank else[]
 dump(out/'all37-and32-factor-identities.json',dict(all37_source_atoms_bound=True,physical_atoms=plan['physical_atoms'],all32_integer_linear_identities=entries,common_linear_denominator=str(DD),F_integer=FZ,B_integer=BZ,F_QQ=FQ,B_QQ=BQ,linear_products=188,nonconstant_Ainverse_coefficient_queries=0))
 dump(out/'actual-H-leading.json',dict(D_integer=DZ,D_denominator=str(DD),F0B0_products65=plan['F0B0_products65'],product_integer=PZ,H_integer=HZ,H_common_denominator=str(DD**2),H_QQ=HQ,whole_integer_and_QQ_product_identity_verified=True,
  H_determinant_known_zero_from_source_support_NOT_new_det_attempt=True))
 dump(out/'rank-kernel-unitminor.json',dict(exact_rank=rank,right_kernel=N,unit_minor_rows=prows,unit_minor_columns=pcols,integer_unit_minor=minor,integer_unit_minor_inverse=minv,integer_unit_minor_determinant=str(mdet),
  rank_operations_on_integerH=ops,integer_row_reduced=rf,row_transform_for_actualQQ_H=R2,column_transform_for_actualQQ_H=Q2,R2_H_Q2=diag,all_rank_kernel_and_unitidentities_verified=True))
 n=7-rank;alpha=proof['outer_row_potentials'];beta=proof['outer_column_potentials'];residualrefs=[]
 for a in range(n):
  for b in range(n):
   coeffs=[dict(row=i,column=j,coefficient=str(R2[rank+a][i]*Q2[j][rank+b]),original_S_degree=alpha[i]+beta[j]+1)for i in range(7)for j in range(7)if R2[rank+a][i]and Q2[j][rank+b]];residualrefs.append(dict(row=a,column=b,first_possible_coefficient='E1=(R2*H1*Q2)_residual because newoffdiagonalconstant0',H1_source_terms=coeffs))
 dump(out/'typed-next-residual-Schur.json',dict(fixed_weight=W,original_rank4_block_A0='I4',original_S='D-C*A^-1*B retained',outer_scaling='H=t^-alpha*S*t^-beta',alpha=alpha,beta=beta,
  new_leading_rank=rank,new_residual_size=n,row_transform_R2=R2,column_transform_Q2=Q2,typed_identity='R2*H*Q2=[[U,V],[W,J]], U0=Irank, V0=W0=J0=0; E=J-W*U^-1*V; E0=0',
  new_inverse_formal_unit_only=True,whole_determinant_order_known_atleast=26579+n,NO_further_coefficients_evaluated=True,first_possible_entry_recipe=residualrefs,
  next_SOURCE_goal='For each referenced H1 coefficient, retain original C*A^-1*B; positiveA gaps>=1 and requiredroute slack<=1, so enumerate source near-shortest A walks (at mostsimplepath+one gap1selfloop) before a finite coefficient budget.',
  no_other_weight_K_point_prime_or_source572=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);actual=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert actual<plan['artifact_RAW_upper']and actual+external<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,ONE_finite_factor_rank_profile=True,all37_source_atoms_and32_linear_identities_verified=True,linear_products=188,quadratic_products=65,H_integer_maxbits=max(abs(v).bit_length()for rr in HZ for v in rr),
  exact_H_leading_rank=rank,remaining_typed_residual_size=n,kernel_and_nonzero_unitminor_verified=True,new_formal_residual_constant_zero=True,next_source_only_order_lower_bound=26579+n,no_further_series_coefficient_det_weight_K_or_source572=True,
  no_det11_zero_nonzero_or_global57_claim=True,actual_new_RAW_before_report=actual,actual_external_read_RAW=external,actual_new_plus_read_RAW=actual+external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
