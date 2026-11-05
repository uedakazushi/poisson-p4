#!/usr/bin/env python3
"""ONE residual6 finite integer DAG at a10e91d20; no further window or series."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_RESIDUAL6_FINITE_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 def mod(rel,name):
  p=ROOT/rel;pins[rel]=sha(p);s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
 gate='reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL6_FINITE_SOURCE_20261105_v63';gr=src(gate+'/report.json');assert gr['all_passed']and pins[gate+'/report.json']=='3c49242955f122800f8d48eeb5c842bc6e564e652d6e927d30a5514038b4e13f';dag=src(gate+'/finite-residual6-DAG.json');budget=src(gate+'/resource-gate.json');assert len(dag['literals'])==95 and len(dag['H_recipes'])==31 and len(dag['residual6_tight_entries'])==25
 lin=mod('code/s61_nonsquare_prePOWER_K92_leading_structure_20261105_v63.py','constantQQ');dethelper=mod('code/s61_nonsquare_prePOWER_K92_finite_leading7_20261105_v63.py','bareiss')
 first='reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63';rt=src(first+'/retropicalized-source.json');R1=[[Q(v)for v in rr]for rr in rt['row_transform']];Q1=[[Q(v)for v in rr]for rr in rt['column_transform']];Ds=list(map(int,rt['constant_clearing_rowLCMs']))
 second='reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63';typed=src(second+'/typed-next-residual-Schur.json');R2=[[Q(v)for v in rr]for rr in typed['row_transform_R2']];Q2=[[Q(v)for v in rr]for rr in typed['column_transform_Q2']]
 values={};requests={};source_receipts=[]
 for l in dag['literals']:
  if l['kind']=='saved_coefficient':requests.setdefault(l['source']['path'],{}).setdefault(l['source']['line'],[]).append(l)
  else:
   s=l['source'];assert s['kind']=='R2Q2';i,j=s['R2_index'];k,n=s['Q2_index'];value=R2[i][j]*Q2[k][n];assert value==Q(l['value'])and R2[i][j]==Q(s['R2_value'])and Q2[k][n]==Q(s['Q2_value']);values[l['id']]=value;source_receipts.append(dict(literal_id=l['id'],kind='R2Q2',exact_value=str(value),verified=True))
 aloop=None
 for rel in set(requests)|{'reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/ABC72.jsonl'}:
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for ln,line in enumerate(f):
    wanted=requests.get(rel,{}).get(ln,[]);e=json.loads(line)if wanted or 'ABC72' in rel else None
    if e is not None and e.get('block')=='A'and e['block_row']==1 and e['block_column']==1:
     co={int(k):Q(e['QQ_content'])*int(v)for k,v in e['sparse_integer_polynomial']};assert co.get(0)==1 and co.get(1)==Q(77,6625);aloop=co[1]
    if not wanted:continue
    pco={int(k):int(v)for k,v in e['sparse_integer_polynomial']};c=Q(e['QQ_content'])
    for l in wanted:
     assert e['entry_index']==l['source']['entry_index']and c==Q(l['source']['QQ_content']);value=c*pco.get(l['gap'],0);assert value==Q(l['value'])and value;values[l['id']]=value;source_receipts.append(dict(literal_id=l['id'],path=rel,line=ln,entry_index=e['entry_index'],gap=l['gap'],source_content=str(c),source_integer=str(pco[l['gap']]),exact_value=str(value),verified=True))
 assert len(values)==95 and aloop==Q(77,6625)
 DL=int(dag['common_literal_denominator']);DEN=int(dag['coefficient_integer_basis_denominator']);DG=int(dag['scalar_inverse_G_formula']['common_denominator']);DH=DL**4*DG;assert DEN==DL*DH and dag['scalar_inverse_G_formula']['indices_used']==[0,1]
 lifted={i:int(v*DL)for i,v in values.items()};assert all(Q(z,DL)==values[i]for i,z in lifted.items());GQ={0:Q(1),1:-aloop};GZ={i:int(v*DG)for i,v in GQ.items()};assert all(Q(z,DG)==GQ[i]for i,z in GZ.items())and GQ[1]+aloop*GQ[0]==0
 dump(out/'all95-source-and-lift.json',dict(literals=dag['literals'],source_receipts=source_receipts,literal_common_denominator=str(DL),literal_integer_lifts={str(i):str(v)for i,v in lifted.items()},G_ONLY0and1_QQ={str(i):str(v)for i,v in GQ.items()},G_integer_lifts={str(i):str(v)for i,v in GZ.items()},G_common_denominator=str(DG),scalar_inverse_identity_checked_ONLY_through1=True,NO_G69_or_inverse_series=True))
 Hnum={};Hcert=[];termcount=0
 for h in dag['H_recipes']:
  total=0;trace=[];bound=math.ceil(DH*Q(h['absolute_bound']))
  for ti,t in enumerate(h['terms']):
   assert len(t['literals'])<=3 and t['G']in[0,1];z=t['sign']*GZ[t['G']]*DL**(4-len(t['literals']))
   for i in t['literals']:z*=lifted[i]
   total+=z;assert abs(z)<=bound and abs(total)<=bound;trace.append(dict(term_index=ti,sign=t['sign'],literals=t['literals'],G=t['G'],integer_term=str(z)));termcount+=1
  Hnum[h['id']]=total;Hcert.append(dict(id=h['id'],Hrow=h['Hrow'],Hcolumn=h['Hcolumn'],relative_gap=h['relative_gap'],original_S_degree=h['original_S_degree'],integer_numerator=str(total),common_denominator=str(DH),exact_QQ_coefficient=str(Q(total,DH)),all_term_integer_identities=trace))
 assert len(Hnum)==31 and termcount==62
 EZ=[[0]*6 for i in range(6)];Ecert=[]
 for e in dag['residual6_tight_entries']:
  total=0;trace=[];bound=math.ceil(DEN*Q(e['absolute_bound']))
  for t in e['terms']:
   z=lifted[t['outer_literal']]*Hnum[t['H_recipe']];total+=z;assert abs(z)<=bound and abs(total)<=bound;trace.append(dict(outer_literal=t['outer_literal'],H_recipe=t['H_recipe'],integer_term=str(z)))
  assert abs(total).bit_length()<=budget['coefficient_integer_bit_upper'];EZ[e['row']][e['column']]=total;Ecert.append(dict(row=e['row'],column=e['column'],relative_gap=e['relative_gap'],integer_numerator=str(total),exact_QQ_coefficient=str(Q(total,DEN)),common_denominator=str(DEN),all_outer_integer_identities=trace))
 detZ,divisions,end=dethelper.bareiss(EZ);assert abs(detZ).bit_length()<=budget['determinant_integer_bit_upper'];rf,P,pcols,prows,ops=lin.rank_ops(EZ);rank=len(pcols);free=[j for j in range(6)if j not in pcols];kernel=[[Q(0)for j in free]for i in range(6)]
 for k,j in enumerate(free):
  kernel[j][k]=Q(1)
  for i,c in enumerate(pcols):kernel[c][k]=-rf[i][j]
 assert lin.mm(EZ,kernel)==[[Q(0)for j in free]for i in range(6)]
 if rank==6:assert detZ==lin.determinant(EZ)
 else:assert detZ==0
 minor=[[Q(EZ[i][j])for j in pcols]for i in prows];mdet=lin.determinant(minor)if rank else Q(1);minv=lin.inverse(minor)if rank else[]
 detR1=lin.determinant(R1);detQ1=lin.determinant(Q1);detR2=lin.determinant(R2);detQ2=lin.determinant(Q2);assert detR1*detQ1*detR2*detQ2
 coef=Q(detZ,DEN**6)/(detR1*detQ1*detR2*detQ2*math.prod(Ds));norm=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/normalization.json');icoef=coef/Q(norm['content_lower_divisor']);assert icoef.denominator==1
 dump(out/'all31-H-coefficients.json',dict(H_denominator=str(DH),all31=Hcert,all62_integer_term_products_verified=True,NO_polynomial_convolution=True))
 dump(out/'all25-E6-coefficients.json',dict(integer_leading6=EZ,common_denominator=str(DEN),all25=Ecert,all_other11_entries_zero_by_source_dual=True))
 dump(out/'det-rank-certificate.json',dict(integer_determinant=str(detZ),exact_Bareiss_divisions=divisions,Bareiss_final=end,exact_rank=rank,right_kernel=kernel,row_transform_integerE=P,row_reduced_matrix=rf,row_operations=ops,pivot_rows=prows,pivot_columns=pcols,
  unit_minor=minor,unit_minor_inverse=minv,unit_minor_determinant=str(mdet),all_rank_and_kernel_identities_checked=True))
 dump(out/'fixed-leading-normalization.json',dict(fixed_weight=[1,2244,68],original_tropical_degree=232227,firstSchur_outer_order=26579,nextSchur_order=74,total_order=26653,original_weight=205574,normalized11_source_monomial=[10,91,20],before_known_a_squared_source_rowclear_monomial=[12,91,20],
  coefficient_integer_basis_denominator=str(DEN),determinant_R1=str(detR1),determinant_Q1=str(detQ1),determinant_R2=str(detR2),determinant_Q2=str(detQ2),all11_rowLCMs=list(map(str,Ds)),rowLCM_product=str(math.prod(Ds)),
  first_fourblock_constant_determinant='1',second_oneblock_constant_determinant='1',row_preconditioner_determinant='1',all_original_row_column_order_signs_preserved=True,
  exact_QQ_coefficient=str(coef),normalized_scalar_I_monomial=[5,91,20],normalized_scalar_I_integer_coefficient=str(icoef.numerator),coefficient_nonzero_PRIMARY=bool(detZ),independent_peer_pending=True,
  identity='coef(detM,a10e91d20)=det(Eint)/(DEN^6*detR1*detQ1*detR2*detQ2*product(rowLCM))'))
 n=6-rank
 dump(out/'next-SOURCE-only.json',dict(leading6_rank=rank,remaining_size=n,ifzero_next='Use constant rank/kernel/unitminor to normalize the scaled residual6 leading block before any new sparse-path source budget; retain BOTH previous Schur corrections.',ifnonzero_next='Independent source/finitecoefficient peer then typed actualpivot completion on existing chart; no global57 inference',NO_further_coefficients_weight_K_point_prime_or_source572=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);actual=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert actual<budget['new_RAW_upper']and actual+external<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,ONE_exact_finite_DAG_job=True,all95_source_literals_bound=True,ONLY_G0_G1_computed=True,all31_H_and25_E_coefficients_saved=True,H_integer_product_terms=62,
  exact_integer_determinant=str(detZ),exact_leading6_rank=rank,exact_target_QQ_coefficient=str(coef),normalized11_source_monomial=[10,91,20],actual_det11_nonzero_PRIMARY=bool(detZ),independent_peer_pending=True,zero_ONLY_this_coefficient=not bool(detZ),
  no_inverse_series_otherweight_precisionwindow_or_source572=True,no_global57_claim=True,actual_new_RAW_before_report=actual,actual_external_read_RAW=external,actual_new_plus_read_RAW=actual+external,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
