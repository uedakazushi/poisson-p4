#!/usr/bin/env python3
"""SOURCE fulla-axis numeratorface +fixedrank6 kernel propagation cost, no rankjob."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_WHOLE_AXIS_SOURCE_V2_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_ORIGINAL_TERMINAL_FACE_SOURCE_20261105_v63';p=src(base+'/terminal-face-and-normalization.json');src(base+'/report.json');pr=src('reports/S61_NONSQUARE_PREPOWER_K92_ORIGINAL_TERMINAL_FACE_ENTRY_V2_20261105_v63/report.json');assert pr['all_passed']and pr['exact_rank']==6;k=src('reports/S61_NONSQUARE_PREPOWER_K92_ORIGINAL_TERMINAL_FACE_ENTRY_V2_20261105_v63/terminal-face-rank-and-source-coefficient.json')
 hp=ROOT/'code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py';pins[str(hp.relative_to(ROOT))]=sha(hp);sp=importlib.util.spec_from_file_location('source',hp);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h);binding=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);M,inventory,kept,skipped=h.selected_inputs(state,binding['actual_pure_matrix_roots']);row=src('reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63/actual-row11.json')
 for e in row['columns']:M[1][e['column']]=dict(state_id=-1001-e['column'],c=Q(e['QQ_common_factor']),co=[(tuple(m),int(v))for m,v in e['primitive_integer_polynomial']])
 axisD=[math.lcm(p['row_integer_denominator'][i],*(M[i][j]['c'].denominator for j in range(11)if any(mon[1]==mon[2]==0 for mon,z in M[i][j]['co'])))for i in range(11)];rowS=[d//old for d,old in zip(axisD,p['row_integer_denominator'])];assert all(d==s*old for d,s,old in zip(axisD,rowS,p['row_integer_denominator']))
 axis={};records=[];total=0;raw_axis=0;min_exp=[10**9]*3;degrees=[]
 for i in range(11):
  for j in range(11):
   atom=M[i][j];pp={};shift=p['row_dual'][i]+p['column_dual'][j];scale=axisD[i]*atom['c'];assert scale.denominator==1
   for mon,z in atom['co']:
    for d,x in enumerate(mon):min_exp[d]=min(min_exp[d],x)
    if mon[1]or mon[2]:continue
    g=mon[0]+shift;assert g>=0;pp[g]=int(scale)*z
   axis[i,j]=pp;total+=len(pp);degrees.append(max(pp,default=0));raw_axis+=len(pp)*(80+max((len(str(v))for v in pp.values()),default=1));records.append(dict(row=i,column=j,state_polynomial_id=atom['state_id'],original_QQ_content=str(atom['c']),axis_row_integer_clear=axisD[i],known_extra_QQ_row_scale=rowS[i],a_rowcol_shift=shift,axis_integer_polynomial=sorted(pp.items())))
 const=[[axis[i,j].get(0,0)for j in range(11)]for i in range(11)];constant_equal=const==[[rowS[i]*v for v in rr]for i,rr in enumerate(p['terminal_integer_face_matrix'])];assert constant_equal
 pc=k['pivot_columns'];prows=k['pivot_rows'];FC=[[Q(x)for x in rr]for rr in k['all_column_factorization']['coefficient_matrix']];R=[[Q(x)/rowS[j]for j,x in enumerate(rr)]for rr in k['row_transform']];free=[j for j in range(11)if j not in pc];RD=math.lcm(*(x.denominator for rr in FC for x in rr));LD=math.lcm(*(x.denominator for rr in R[6:]for x in rr));assert RD.bit_length()<=521 and LD.bit_length()<=521
 jobs=[];scalar_adds=0;slots=0;bits=0
 for side in['right_kernel','left_kernel']:
  for i in range(11 if side=='right_kernel'else 5):
   for j in range(5 if side=='right_kernel'else 11):
    terms=[]
    if side=='right_kernel':
     col=free[j];terms=[(i,col,RD)]+[(i,c,-int(RD*FC[z][col]))for z,c in enumerate(pc)if FC[z][col]]
    else:terms=[(z,j,int(LD*R[6+i][z]))for z in range(11)if R[6+i][z]]
    support=set();bound=0;adds=0
    for u,v,s in terms:support.update(axis[u,v]);bound+=abs(s)*sum(abs(x)for x in axis[u,v].values());adds+=len(axis[u,v])
    jobs.append(dict(side=side,row=i,column=j,source_terms=[dict(row=u,column=v,integer_scale=str(s))for u,v,s in terms],support_upper=len(support),coefficient_bits_upper=bound.bit_length(),integer_abs_upper=str(bound),scalar_additions_upper=adds));scalar_adds+=adds;slots+=len(support);bits=max(bits,bound.bit_length())
 live=2*total+2*max((e['support_upper']for e in jobs),default=0)+len(jobs)*12+15000;raw=raw_axis+slots*(100+(bits+2)//3)+1200000
 dump(out/'whole-axis-literal-source.json',dict(axis_is_coefficient_face_e0d0=True,NOT_numeric_beta_point=True,same_original_fixedweight=[1,2244,68],all_source_variables_minimum_exponents=min_exp,
  no_negative_e_or_d_exponents=min_exp[1]>=0 and min_exp[2]>=0,axis_coefficients=records,axis_total_integer_coefficients=total,max_entry_degree=max(degrees),
  normalization='Caxis(a)=diag(axis_rowLCD*a^terminal_rowdual)*Mprime(a,0,0)*diag(a^terminal_coldual); detCaxis=(prodaxisLCD)*a^-5*detMaxis',
  Caxis0_EQUALS_knownQQrowscale_of_verified_terminal_rank6_matrix=True,extra_QQ_row_scales=rowS,axis_row_denominators=axisD,fixed6minor_rows=prows,fixed6minor_columns=pc,fixed6minor_constant_determinant=str(Q(k['unit_minor_determinant'])*math.prod(rowS[i]for i in prows)),fixed6minor_NONZERO_over_QQ_a_without_ranksearch=True,
  FULLU_boundary_scope='d=0 makesoriginalFULLUlocalizationundefined; thisis ONLY the polynomialclearednumeratorcoefficientface, NOT actualPoissonrankataboundarypoint',
  other_named_unit_axis_images=dict(a='a',d='0',F='500a-24',B='-40',C='100a-2',J='-4'),
  no_inference_of_fullU_failure_or_fullrank_from_axis_rank=True))
 dump(out/'fixed-kernel-propagation-source.json',dict(right_factorization_FC=FC,right_clear=str(RD),left_bottom5_R=R[6:],left_clear=str(LD),all110_finite_polynomial_linear_checks=jobs,
  mathematical_completion='ifall55rightORall55left identitiesZERO, rankwholeaxis<=6; fixednonzero6minor givesexactrank6. IfbothsetsZERO, constantterminalnormalizers splitwholeaxisinto6unitblockpluszero5block; transverse(e,d) initial5x5then finiteSOURCEtarget.',
  if_nonzero='recordconcretepolynomialobstruction; no assumedwholeaxisrank6',no_polynomial_fieldGaussian=True,no_new_coefficient_values_evaluated=True))
 dump(out/'resource-gate.json',dict(job='ONE complete110 polynomialLINEARkernel-propagation checks, no det/ranksearch; dependsROOTbudget',input_axis_coefficients=total,max_entry_degree=max(degrees),all_output_coefficient_slots_upper=slots,scalar_additions_upper=scalar_adds,
  maximum_new_integer_bits=bits,whole_owned_slots_upper=live,new_RAW_upper=raw,hard_seconds=60,hard_address_bytes=2000000000,hard_scalar_terms=150000,hard_total_RAW=33554432))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,whole_axis_literal_source_complete=True,axis_integer_coefficient_terms=total,max_axis_entry_degree=max(degrees),axis_constant_EQUALS_knownQQrowscaled_terminal_rank6=True,fixed6minor_NONZERO_sourcebound=True,
  d0_boundary_not_FULLU_point=True,no_kernel_propagation_or_axis_rank_newarithmetic=True,proposed110_linear_checks_owned=live,proposed_checks_bits=bits,proposed_checks_additions=scalar_adds,proposed_new_RAW=raw,
  no_full_detM_or_global57_claim=True,total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
