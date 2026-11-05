#!/usr/bin/env python3
"""SOURCE ONLY: one fixed saved-row comparison and typed residual5 next target."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_SOURCE_STRUCTURE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL6_FINITE_ENTRY_20261105_v63';r=src(base+'/report.json');assert r['all_passed']and r['exact_leading6_rank']==1 and r['zero_ONLY_this_coefficient'];cert=src(base+'/det-rank-certificate.json');E=src(base+'/all25-E6-coefficients.json');P=[[Q(v)for v in rr]for rr in cert['row_transform_integerE']];EZ=[[int(v)for v in rr]for rr in E['integer_leading6']];DEN=int(E['common_denominator']);pcols=cert['pivot_columns'];free=[j for j in range(6)if j not in pcols];assert len(pcols)==1
 hp=ROOT/'code/s61_nonsquare_prePOWER_K92_leading_structure_20261105_v63.py';pins[str(hp.relative_to(ROOT))]=sha(hp);sp=importlib.util.spec_from_file_location('lin',hp);lin=importlib.util.module_from_spec(sp);sp.loader.exec_module(lin)
 order=pcols+free;Q3=[[Q(0)for j in range(6)]for i in range(6)]
 for j,c in enumerate(order):Q3[c][j]=Q(1)
 U=lin.eye(6);rf=[[Q(v)for v in rr]for rr in cert['row_reduced_matrix']]
 for j,c in enumerate(free):U[0][1+j]=-rf[0][c]
 Q3=lin.mm(Q3,U);R3=[[DEN*v for v in rr]for rr in P];EQ=[[Q(v,DEN)for v in rr]for rr in EZ];assert lin.mm(lin.mm(R3,EQ),Q3)==[[Q(i==j and i==0)for j in range(6)]for i in range(6)]
 # This one fixed pair is selected from the already-observed identical
 # valuation patterns, not a row-pair search or a new polynomial evaluation.
 pair={};read_terms=0
 for rel in['reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/ABC72.jsonl','reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/D49.jsonl']:
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for ln,line in enumerate(f):
    e=json.loads(line);b=e.get('block','D');i=e['block_row']if b!='D'else e['residual_row'];j=e['block_column']if b!='D'else e['residual_column']
    if b not in['C','D']or i not in[2,4]:continue
    content=Q(e['QQ_content']);co={int(k):int(v)for k,v in e['sparse_integer_polynomial']};pair[b,i,j]=dict(path=rel,line=ln,content=content,co=co);read_terms+=len(co)
 assert len(pair)==22;ratio=None;counter=None;checked=0;checks=[]
 for b,width in[('C',4),('D',7)]:
  for j in range(width):
   a=pair[b,2,j];z=pair[b,4,j];degrees=sorted(set(a['co'])|set(z['co']))
   for g in degrees:
    va=a['content']*a['co'].get(g,0);vb=z['content']*z['co'].get(g,0)
    if ratio is None:
     if not va and not vb:continue
     if not va or not vb:counter=dict(block=b,column=j,exponent=g,left=str(va),right=str(vb),reason='onezero');break
     ratio=va/vb
    checked+=1
    if va!=ratio*vb:counter=dict(block=b,column=j,exponent=g,left=str(va),right=str(vb),candidate_constant_ratio=str(ratio),scaled_right=str(ratio*vb),left_source=dict(path=a['path'],line=a['line']),right_source=dict(path=z['path'],line=z['line']));break
   checks.append(dict(block=b,column=j,comparison_complete=counter is None))
   if counter:break
  if counter:break
 full_equal=counter is None and ratio is not None
 dump(out/'fixed-full-row-comparison.json',dict(pair_global_transformed_rows=[6,8],pair_blocks=['Crows2,4','Drows2,4'],one_pair_ONLY=True,saved_terms_read=read_terms,candidate_constant_ratio=None if ratio is None else str(ratio),coefficient_comparisons=checked,all11_entries_proportional=full_equal,
  first_countercoefficient=counter,checks=checks,NO_new_polynomial_coefficients_or_products=True,if_all_equal_scope='formal11detzero, then fixedactualdetzero only via proven mixedradix injection; peerrequired',NO_full_detzero_claim_if_counter=True))
 nxt=src('reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/next-residual6-source.json');aa=nxt['row_dual_potentials'];bb=nxt['column_dual_potentials'];Emin=nxt['residual6_min'];Tmin=nxt['transformed_positive_part_min'];refs=[];queries=set();dropped=0
 for i in range(5):
  for j in range(5):
   terms=[]
   for p in range(6):
    if not R3[1+i][p]:continue
    for q in range(6):
     if not Q3[q][1+j]:continue
     g=aa[p]+bb[q]+1
     if Emin[p][q]>g:dropped+=1;continue
     assert g>=1 and Emin[p][q]>=aa[p]+bb[q];queries.add((p,q,g));terms.append(dict(source_E_row=p,source_E_column=q,source_E_degree=g,QQ_coefficient=str(R3[1+i][p]*Q3[q][1+j]),new1block_inverse_slack_upper=g-(Tmin[1+p][0]+Tmin[0][1+q])))
   refs.append(dict(residual5_row=i,residual5_column=j,first_positive_entry='Z1=(R3*scaledE6_1*Q3)_residual; newoffdiagcorrectionstarts2',terms=terms))
 assert all(t['new1block_inverse_slack_upper']<=1 for e in refs for t in e['terms']);maxg=max((g for p,q,g in queries),default=0)
 dump(out/'typed-residual5-and-next-source.json',dict(R3=R3,Q3=Q3,actual_leading6_denominator=str(DEN),R3_E0_Q3='diag(I1,0_5)',scaled_source='L6=t^-aa*previousResidual6*t^-bb',aa=aa,bb=bb,
  typed_residual='R3*L6*Q3=[[V,W],[X,Y]],V0=1,W0=X0=Y0=0; Z5=Y-X*V^-1*W; Z0=0; Z1=(R3*L6_1*Q3)_residual',remaining_size=5,whole_order_known_atleast=26658,
  first_positive_matrix_source_refs=refs,unique_previousE_coefficient_queries=sorted(queries),source_impossible_terms_removed=dropped,maximum_derived_previousE_degree=maxg,
  retained_previous_correction='previousE=J-Wold*Uold^-1*Vold; dualbounds imply ONLYUold_inverse0and1 needed at these queries, Uinv1=-U1',original4block_correction='H=D-C*A^-1*B remains; derive finite nearshortest window frommaxg BEFORE coefficients',
  no_new_coefficient_or_determinant_or_precisionwindow=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);total=sum(p.stat().st_size for p in out.iterdir());assert read_terms+15000<150000 and total<1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,one_fixed_saved_row_comparison=True,full_row_proportional=full_equal,countercoefficient_found=counter is not None,
  new_typed_residual_size=5,constant_unit_normalization_verified=True,known_order_atleast26658=True,next_E_source_queries=len(queries),maximum_derived_E_degree=maxg,previous_scalar_inverse_slack_atmost1=True,
  no_new_polynomial_or_series_coefficients_or_determinant=True,no_full_det11_zero_nonzero_or_global57_claim=True,total_RAW_before_report=total,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
