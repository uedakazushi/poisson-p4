#!/usr/bin/env python3
"""ONE fixed finite leading7 calculation; no series/weight/layer fallback."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_FINITE_LEADING7_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def bareiss(A):
 B=[list(map(int,row))for row in A];n=len(B);previous=1;sign=1;events=[]
 for k in range(n-1):
  i=next((i for i in range(k,n)if B[i][k]),None)
  if i is None:return 0,events,dict(zero_column=k,after_nonzero_pivot_prefix=True)
  if i!=k:B[i],B[k]=B[k],B[i];sign=-sign;events.append(dict(op='swap',rows=[i,k]))
  p=B[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    z=B[i][j]*p-B[i][k]*B[k][j];assert z%previous==0;q=z//previous;events.append(dict(op='exactdivision',k=k,row=i,column=j,old_ij=str(B[i][j]),old_ik=str(B[i][k]),old_kj=str(B[k][j]),pivot=str(p),divisor=str(previous),quotient=str(q)));B[i][j]=q
   B[i][k]=0
  previous=p
 return sign*B[-1][-1],events,dict(last_entry=str(B[-1][-1]),row_swap_sign=sign)
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 def mod(rel,name):
  p=ROOT/rel;pins[rel]=sha(p);s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
 gate='reports/S61_NONSQUARE_PREPOWER_K92_SCHUR_GAP_SOURCE_PREFLIGHT_20261105_v63';gr=src(gate+'/report.json');assert gr['all_passed']and pins[gate+'/report.json']=='53846a4145284e794555fb843b4a4d3a6d8f81bf4201aa6fde0dc946c339d119';plan=src(gate+'/finite-leading7-source.json');support=src(gate+'/minplus-proof.json');assert len(plan['physical_atoms'])==66 and len(plan['entry_recipes'])==41 and plan['first_possible_determinant_gap']==13255
 base='reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63';rt=src(base+'/retropicalized-source.json');row=src(base+'/actual-row11.json');assert rt['new_leading_constant_rank']==4 and rt['new_tropical_degree']==232227
 helper=mod('code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py','sourcegate');lin=mod('code/s61_nonsquare_prePOWER_K92_leading_structure_20261105_v63.py','constantQQ')
 bind=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02';M,inventory,kept,skip=helper.selected_inputs(state,bind['actual_pure_matrix_roots'])
 for p in row['columns']:M[1][p['column']]=dict(state_id=-1001-p['column'],c=Q(p['QQ_common_factor']),co=[(tuple(m),int(v))for m,v in p['primitive_integer_polynomial']])
 R=[[Q(v)for v in rr]for rr in rt['row_transform']];C=[[Q(v)for v in rr]for rr in rt['column_transform']];A0=[[Q(v)for v in rr]for rr in rt['new_leading_constant_matrix']];assert lin.mm(lin.mm(R,A0),C)==[[Q(i==j and i<4)for j in range(11)]for i in range(11)]
 DR=int(plan['R_common_denominator']);DQ=int(plan['Q_common_denominator']);DD=int(plan['common_denominator']);assert DR*DQ==DD;Ds=list(map(int,rt['constant_clearing_rowLCMs']));physical={};atom_receipts=[]
 for a in plan['physical_atoms']:
  i=a['row'];j=a['column'];p=M[i][j];assert p['state_id']==a['source_state_id'];mon,v=p['co'][a['source_index']];assert mon==tuple(a['monomial'])and v==int(a['source_integer'])and p['c']==Q(a['source_QQfactor']);gap=rt['row_potentials'][i]+rt['column_potentials'][j]-sum(x*y for x,y in zip(mon,(1,2244,68)));assert gap==a['gap'];z=p['c']*v*Ds[i];assert z.denominator==1 and z.numerator==int(a['cleared_integer']);physical[a['id']]=z.numerator;atom_receipts.append(dict(a,whole_source_atom_checked=True))
 HZ=[[0]*7 for i in range(7)];entry_cert=[];products=0
 tight={tuple(v)for v in support['all_tight_entries']};assert tight=={(q['residual_row'],q['residual_column'])for q in plan['entry_recipes']}
 for q in plan['entry_recipes']:
  i=q['residual_row'];j=q['residual_column'];u=q['global_row'];v=q['global_column'];assert u==i+4 and v==j+4 and q['gap']==support['Schur_min'][i][j]and support['correction_C_Ainv_B_min'][i][j]>q['gap'];total=0
  for c in q['contributors']:
   a=plan['physical_atoms'][c['physical_id']];assert a['id']==c['physical_id']and a['gap']==q['gap'];rc=int(c['R_cleared_coefficient']);cc=int(c['Q_cleared_coefficient']);assert R[u][a['row']]*DR==rc and C[a['column']][v]*DQ==cc;total+=rc*physical[a['id']]*cc;products+=2;assert abs(total).bit_length()<=plan['max_integer_entry_numerator_bits']
  HZ[i][j]=total;entry_cert.append(dict(residual_row=i,residual_column=j,gap=q['gap'],integer_numerator=str(total),common_denominator=str(DD),contributors=q['contributors']))
 assert products==260 and all(HZ[i][j]==0 for i in range(7)for j in range(7)if(i,j)not in tight)
 detZ,divisions,end=bareiss(HZ);assert abs(detZ).bit_length()<=plan['max_integer_det_bits'];rf,P,pcols,prows,ops=lin.rank_ops(HZ);rank=len(pcols);free=[j for j in range(7)if j not in pcols];kernel=[[Q(0)for j in free]for i in range(7)]
 for k,j in enumerate(free):
  kernel[j][k]=Q(1)
  for i,c in enumerate(pcols):kernel[c][k]=-rf[i][j]
 assert lin.mm(HZ,kernel)==[[Q(0)for j in free]for i in range(7)]
 if rank==7:assert detZ==lin.determinant(HZ)
 else:assert detZ==0
 minor=[[Q(HZ[i][j])for j in pcols]for i in prows];mi=lin.inverse(minor)if rank else[];md=lin.determinant(minor)if rank else Q(1);detR=lin.determinant(R);detQ=lin.determinant(C);assert detR and detQ
 coefficient=Q(detZ,DD**7)/(detR*detQ*math.prod(Ds));norm=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/normalization.json');Icoef=coefficient/Q(norm['content_lower_divisor']);assert Icoef.denominator==1
 dump(out/'all66-source-atoms.json',dict(atoms=atom_receipts,source_state_sha256=pins[str(state.relative_to(ROOT))],all66_checked=True,no_source572_replay=True))
 dump(out/'all41-entry-identities.json',dict(integer_matrix7=HZ,QQ_matrix7=[[str(Q(v,DD))for v in rr]for rr in HZ],common_denominator=str(DD),entries=entry_cert,all_other8_entries_zero_by_dual=True,integer_products=products))
 dump(out/'det-rank-certificate.json',dict(integer_determinant=str(detZ),exact_Bareiss_divisions=divisions,final_Bareiss_state=end,exact_rank=rank,right_kernel=kernel,row_transform=P,row_reduced_matrix=rf,row_operations=ops,pivot_columns=pcols,pivot_rows=prows,unit_minor=minor,unit_minor_inverse=mi,unit_minor_determinant=str(md),all_matrix_and_kernel_identities_checked=True))
 dump(out/'leading-coefficient-normalization.json',dict(weight=[1,2244,68],source_matrix_tropical_degree=232227,Schur_minimum_order=13255,original_weight=218972,normalized_source_matrix_monomial=[12,97,19],before_known_a_squared_source_row_normalization_monomial=[14,97,19],
  determinant_R=str(detR),determinant_Q=str(detQ),all11_row_LCMs=list(map(str,Ds)),row_LCM_product=str(math.prod(Ds)),fourblock_A0=lin.eye(4),fourblock_constant_determinant='1',residual_common_denominator=str(DD),
  determinant_identity='coefficient(detM,a12e97d19)=det(integerLeading7)/(DD^7 * detR * detQ * product(rowLCM)); all formalSchur corrections strictlylater on tight edges',
  exact_QQ_coefficient=str(coefficient),normalized_integer_scalar_I_monomial=[7,97,19],normalized_integer_scalar_I_coefficient=str(Icoef.numerator),normalization_cglobal=norm['content_lower_divisor'],
  row_preconditioner_determinant='1',all_original_row_column_order_signs_preserved=True,nonzero_primary=bool(detZ),zero_means_ONLY_this_leading_layer=True))
 # Source-only next structure; no new coefficient request is evaluated.
 dump(out/'next-source-structure.json',dict(leading7_rank=rank,leading7_kernel_dimension=7-rank,leading7_unit_minor_rows=prows,leading7_unit_minor_columns=pcols,
  possible_next='Use this constant rank/profile after t^minplus row/column scaling to a typed unit Schur of size7-rank; derive next physical/source-gap request BEFORE any further coefficient job',
  important='ResidualS still contains original C A^-1 B; firstlayer correctionwasstrictlylater, notidenticallyzero. It must be retained at the next source gap.',
  NO_new_series_weight_K_or_fallback=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);totalbytes=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert totalbytes<677280 and totalbytes+external<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,ONE_finite_leading7_job=True,all66_source_atoms_bound=True,all41_integer_entry_identities=True,integer_products=products,exact_integer_determinant=str(detZ),exact_leading7_rank=rank,
  exact_target_QQ_coefficient=str(coefficient),normalized_source_monomial=[12,97,19],actual_det11_nonzero_primary=bool(detZ),independent_peer_pending=True,zero_ONLY_thislayer=not bool(detZ),no_new_series_weight_K_or_source572=True,no_global57_claim=True,
  actual_new_RAW_before_report=totalbytes,actual_external_read_RAW=external,actual_new_plus_read_RAW=totalbytes+external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
