#!/usr/bin/env python3
"""ONE full literal row addition, followed by same-weight SOURCE gap/Schur gate."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63';W=(1,2244,68)
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 def mod(rel,name):
  p=ROOT/rel;pins[rel]=sha(p);s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
 previous='reports/S61_NONSQUARE_PREPOWER_K92_LEADING_STRUCTURE_20261105_v63';gate=src(previous+'/report.json');assert gate['all_passed']and pins[previous+'/report.json']=='11e88eb715e0c3ae36b2ff9ed1062e337ab4ce8631d581f6879ff99c56e79c32';bound=src(previous+'/one-row-source-gate.json');assert bound['output_row_support_upper']==627
 old=src('reports/S61_NONSQUARE_PREPOWER_K92_TROPICAL_PREFIX32_V2_20261105_v63/selected78-source.json');assert old['row_LCMs'][1]==old['row_LCMs'][2]and old['row_potentials'][1]==old['row_potentials'][2]
 helper=mod('code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py','profile');lin=mod('code/s61_nonsquare_prePOWER_K92_leading_structure_20261105_v63.py','constantQQ')
 binding=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02';M,inventory,kept,skipped=helper.selected_inputs(state,binding['actual_pure_matrix_roots']);assert kept==21690
 outputs=[];all_support=0;peak=kept;maxbits=0
 dump(out/'phase.json',dict(phase='ONE literalrow add starts; no polynomial products'))
 for j in range(11):
  a=M[1][j];b=M[2][j];D=math.lcm(a['c'].denominator,b['c'].denominator);s=a['c'].numerator*(D//a['c'].denominator);t=2*b['c'].numerator*(D//b['c'].denominator);poly={}
  for p,scale in[(a,s),(b,t)]:
   for mon,v in p['co']:
    value=poly.get(mon,0)+scale*v;assert abs(value).bit_length()<=bound['columns'][j]['integer_numerator_bit_upper'];before=kept+all_support+len(poly)+15000;assert before<150000
    if value:poly[mon]=value
    elif mon in poly:del poly[mon]
    peak=max(peak,kept+all_support+len(poly))
  assert len(poly)<=bound['columns'][j]['output_support_upper'];G=math.gcd(*(abs(v)for v in poly.values()))if poly else 1;content=Q(G,D);normal=[(m,v//G)for m,v in sorted(poly.items())];assert content>0;pa=dict(a['co']);pb=dict(b['co']);pn=dict(normal)
  # Different normalization-based whole-coefficient replay, no convolution.
  for mon in set(pa)|set(pb)|set(pn):assert a['c']*pa.get(mon,0)+2*b['c']*pb.get(mon,0)==content*pn.get(mon,0)
  for mon,v in normal:maxbits=max(maxbits,abs(v).bit_length())
  result=dict(column=j,source_state_ids=[a['state_id'],b['state_id']],source_row_labels=['G18/a','G39/a'],QQ_common_factor=str(content),primitive_integer_polynomial=[[list(m),str(v)]for m,v in normal],source_common_denominator=str(D),source_integer_scales=[str(s),str(t)],all_coefficient_identity_checked=True)
  outputs.append(result);all_support+=len(normal);M[1][j]=dict(state_id=-1001-j,c=content,co=normal)
 assert all_support>0 and all(not M[1][j]['co']for j in[8,9,10]);assert all_support<=627
 dump(out/'actual-row11.json',dict(operation='row1 <- row1+2row2',fixed_weight=W,no_exponent_shift=True,old_rowpot_equal=old['row_potentials'][1],old_rowLCM_equal=old['row_LCMs'][1],columns=outputs,actual_total_terms=all_support,full_row_nonzero=True,no_new_parameter_inverse=True,all_source572_NOT_replayed=True))
 vals=[]
 for row in M:
  rr=[]
  for p in row:
   values=[sum(a*b for a,b in zip(mon,W))for mon,v in p['co']];assert len(values)==len(set(values));rr.append(max(values,default=helper.NEG))
  vals.append(rr)
 T,count=helper.assignment(vals);rp,cp,match=helper.dual_potentials(vals);assert sum(rp)+sum(cp)==T
 Ds=[math.lcm(*(p['c'].denominator for p in row))for row in M];A=[[Q(0)for j in range(11)]for i in range(11)];hist={};records=[];delta=None
 for i,row in enumerate(M):
  for j,p in enumerate(row):
   for idx,(mon,v)in enumerate(p['co']):
    gap=rp[i]+cp[j]-sum(a*b for a,b in zip(mon,W));assert gap>=0;hist[gap]=hist.get(gap,0)+1
    if gap==0:A[i][j]=p['c']*v*Ds[i];assert A[i][j].denominator==1
    elif delta is None or gap<delta:delta=gap
 assert delta is not None
 B,P,pcols,prows,ops=lin.rank_ops(A);r=len(pcols);free=[j for j in range(11)if j not in pcols];order=pcols+free;Qmat=[[Q(0)for j in range(11)]for i in range(11)]
 for j,c in enumerate(order):Qmat[c][j]=Q(1)
 U=lin.eye(11)
 for i in range(r):
  for j,c in enumerate(free):U[i][r+j]=-B[i][c]
 Qmat=lin.mm(Qmat,U);diag=[[Q(i==j and i<r)for j in range(11)]for i in range(11)];assert lin.mm(lin.mm(P,A),Qmat)==diag
 minor=[[A[i][j]for j in pcols]for i in prows];mdet=lin.determinant(minor);minv=lin.inverse(minor);assert mdet
 DR=math.lcm(*(v.denominator for row in P for v in row));DQ=math.lcm(*(v.denominator for row in Qmat for v in row));RZ=[[int(v*DR)for v in row]for row in P];QZ=[[int(v*DQ)for v in row]for row in Qmat];delta_max=0
 for i,row in enumerate(M):
  for j,p in enumerate(row):
   for idx,(mon,v)in enumerate(p['co']):
    gap=rp[i]+cp[j]-sum(a*b for a,b in zip(mon,W))
    if gap==delta:
     value=p['c']*v*Ds[i];assert value.denominator==1;delta_max=max(delta_max,abs(value.numerator));records.append(dict(row=i,column=j,source_state_id=p['state_id'],source_index=idx,monomial=mon,QQfactor=str(p['c']),integer=str(v),row_LCM=str(Ds[i]),cleared_integer=str(value.numerator)))
 n=11-r;boundnum=0
 for i in range(r,11):
  for j in range(r,11):boundnum=max(boundnum,delta_max*sum(abs(v)for v in RZ[i])*sum(abs(QZ[k][j])for k in range(11)))
 target=T-n*delta;z=target-5;mon=(5+z%68,(z//68)//33,(z//68)%33);box_inside=5<=mon[0]<=72 and 0<=mon[1]<=107 and 0<=mon[2]<=32
 oldwb=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/weighted-output.json');newton_inside=box_inside and all(lo<=sum(a*b for a,b in zip(mon,w))<=hi for w,(lo,hi)in zip(oldwb['weights'],oldwb['output_metadata']['linear']))
 dump(out/'retropicalized-source.json',dict(fixed_weight=W,new_exact_valuation_matrix=vals,new_tropical_degree=T,old_tropical_degree=232363,tropical_drop=232363-T,optimal_matching_count=count,row_potentials=rp,column_potentials=cp,matching=match,
  ALL_source_gaps_histogram=sorted(hist.items()),smallest_positive_gap=delta,constant_clearing_rowLCMs=list(map(str,Ds)),new_leading_constant_matrix=A,new_leading_constant_rank=r,right_nullity=n,
  row_transform=P,column_transform=Qmat,R_A0_Q=diag,unit_minor_rows=prows,unit_minor_columns=pcols,unit_minor_determinant=str(mdet),unit_minor_inverse=minv,constant_row_operations=ops,NO_further_determinantseries=True))
 dump(out/'next-delta-block-source-gate.json',dict(target=f'ONE{n}x{n}matrix S_delta',fixed_weight=W,delta=delta,rank=r,residual_size=n,typed_identity='R*B(t)*Q has constant diag(Ir,0); every positive entrydegree>=delta; Schur leading S_delta=(R*B_delta*Q)_residual, correction starts2delta',
  conditional_nonzero='det(S_delta)!=0 =>ord detB=(11-rank)*delta; fixed det11 NONZERO',known_order_lower_bound=n*delta,target_original_weight=target,target_unique_original_monomial=mon,monomial_in_known_Newton_bounds=newton_inside,
  physical_delta_source_records=records,number_physical_delta_entries=len(records),R_common_denominator=str(DR),Q_common_denominator=str(DQ),R_cleared_integer=RZ,Q_cleared_integer=QZ,
  maximum_residual_integer_numerator_bits=boundnum.bit_length(),common_residual_denominator=str(DR*DQ),maximum_scalar_products=n*n*121,source_fixed_terms=kept+all_support,
  future_residual_support_slots=n*n,upper_initialized_integer_slots=kept+all_support+121+n*n+15000,determinant_bit_upper=n*boundnum.bit_length()+math.factorial(max(1,n)).bit_length(),
  NO_delta_residual_matrix_OR_determinant_evaluated=True,NO_new_formalcoeffquery_or_otherweight_or_K=True))
 dump(out/'phase.json',dict(phase='FINAL ONEfullrow exact; sameweight source gaps/leadingconstant only'))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);actual=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert actual+external<32*1024*1024 and peak+15000<150000
 report=dict(id=ID,all_passed=True,frozen=True,ONE_literal_fullrow_add=True,all11_wholecoefficient_identities_checked=True,actual_newrow_terms=all_support,actual_newrow_nonzero=True,primitive_integer_coefficient_maxbits=maxbits,
  same_weight=W,new_tropical_degree=T,tropical_degree_drop=232363-T,new_constant_rank=r,residual_size=n,smallest_positive_gap=delta,next_delta_physical_entries=len(records),next_delta_residual_NOT_evaluated=True,
  no_determinantseries_otherweight_or_K=True,no_det11_nonzero_or_global57_claim=True,peak_source_plus15000scratch=peak+15000,actual_new_RAW_before_report=actual,actual_external_read_RAW=external,actual_new_plus_read_RAW=actual+external,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
