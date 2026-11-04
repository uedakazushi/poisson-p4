#!/usr/bin/env python3
"""Read-only statistics and provenance for the safe36 cubic coker input."""
import argparse, hashlib, json, math, resource, signal, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_CUBIC36_FUNCTION_FIELD_PREFLIGHT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(500000000,500000000))
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('30s read-only preflight cap')));signal.alarm(30)
 out.mkdir(parents=True,exist_ok=False)
 src=ROOT/'work/nonsquare-v62/cubic27-v6/after21';h=sha(src)
 assert h=='c65b3d80b2207e75ba794821bb41ebc3ffc864959dab27307af77cdeab9d8fd1'
 man=ROOT/'research/handoff/20261103_v62/v6_cubic27_coker/manifest.json';manifest=json.loads(man.read_text())
 pinned=next(x for x in manifest['omitted_large_generated_files']if x['path']==str(src.relative_to(ROOT)))
 assert pinned['sha256']==h and pinned['bytes']==src.stat().st_size
 phase=ROOT/'work/nonsquare-v62/cubic27-v6/phase21-check.json';ph=json.loads(phase.read_text());assert ph['all_passed']and ph['known_unit_pivots']==21
 code=ROOT/'work/nonsquare-v62/cubic27-v6/forced21.cpp';codeh=next(x['sha256']for x in manifest['members']if x['path']==str(code.relative_to(ROOT)));assert sha(code)==codeh
 basis=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json';br=json.loads(basis.read_text());assert br['all_passed']and br['alternative_ad_triangular_minor']=='a^3*d/10'
 assert all(t in code.read_text()for t in ['"a"','"500*a+90*e^2+3*e-24"','"50*a*d+18*d*e^2-9*d*e-2*d+60*e-40"','"100*a+18*e^2-9*e-2"'])
 prime=32003;rows=[];colterms=[0]*70;colnonzero=[0]*70;degrees=[0,0,0];totaldegree=0;eddegree=0;terms=0;nonzero=0;maxbits=[0,0];denoms=set();constantden=True;bad_den=0;div_num=0;grouped=0;largest_entry={'terms':0};hist={}
 with src.open()as f:
  n,m,k=map(int,next(f).split());assert(n,m,k)==(36,70,21)
  used=list(map(int,next(f).split()));ids=list(map(int,next(f).split()));assert used==ph['used_original_columns']and ids==ph['remaining_original_rows'];assert len(used)==21 and len(set(used))==21 and len(ids)==36
  originalcols=[j for j in range(91)if j not in used]
  for i in range(n):
   rs={'row_index':i,'original_row_index':ids[i],'nonzero_entries':0,'coefficient_terms':0,'max_degree_a_e_d':[0,0,0],'max_total_degree':0,'max_degree_e_d':0,'constant_rational_denominator_lcm':1,'coefficient_numerator_gcd':0,'all_numerators_divisible_by_32003':True,'grouped_function_field_terms_e_d':0}
   for j in range(m):
    token=next(f).split();l=int(token[0]);assert len(token)==1+4*l and l>=0
    hist[l]=hist.get(l,0)+1
    if l:nonzero+=1;rs['nonzero_entries']+=1;colnonzero[j]+=1
    terms+=l;rs['coefficient_terms']+=l;colterms[j]+=l;edmons=set();last=None
    if l>largest_entry['terms']:largest_entry={'row_index':i,'column_index':j,'original_row_index':ids[i],'original_column_index':originalcols[j],'terms':l}
    for t in range(l):
     q=token[1+4*t];ex=tuple(map(int,token[2+4*t:5+4*t]));assert min(ex)>=0 and(last is None or last>ex);last=ex
     if '/'in q:num,den=map(int,q.split('/'))
     else:num,den=int(q),1
     assert num and den>0 and math.gcd(abs(num),den)==1
     denoms.add(den);bad_den+=den%prime==0;div_num+=num%prime==0
     rs['all_numerators_divisible_by_32003']&=num%prime==0
     rs['constant_rational_denominator_lcm']=math.lcm(rs['constant_rational_denominator_lcm'],den)
     rs['coefficient_numerator_gcd']=math.gcd(rs['coefficient_numerator_gcd'],abs(num))
     maxbits=[max(maxbits[0],abs(num).bit_length()),max(maxbits[1],den.bit_length())]
     for v in range(3):degrees[v]=max(degrees[v],ex[v]);rs['max_degree_a_e_d'][v]=max(rs['max_degree_a_e_d'][v],ex[v])
     totaldegree=max(totaldegree,sum(ex));eddegree=max(eddegree,ex[1]+ex[2]);rs['max_total_degree']=max(rs['max_total_degree'],sum(ex));rs['max_degree_e_d']=max(rs['max_degree_e_d'],ex[1]+ex[2]);edmons.add(ex[1:])
    grouped+=len(edmons);rs['grouped_function_field_terms_e_d']+=len(edmons)
   for field in ['constant_rational_denominator_lcm','coefficient_numerator_gcd']:rs[field]=str(rs[field])
   assert rs['nonzero_entries']>0;rows.append(rs)
  assert not f.read().strip()
 assert bad_den==0 and all(not r['all_numerators_divisible_by_32003']for r in rows)
 stats={'source_sha256':h,'source_bytes':src.stat().st_size,'source_format':{'header':[n,m,k],'header_third_integer_is_completed_pivots_not_rank':True,'next_two_lines':['21 used original column indices','36 retained original row indices'],'entry_order':'row-major; exactly2520 records','polynomial_record':'termcount, then rational_coefficient a_exponent e_exponent d_exponent for each term','variable_order':['a','e','d'],'all_exponents_nonnegative':True,'terms_strict_descending_lex_no_duplicates':True,'rational_coefficients_nonzero_canonical':True},
 'shape':[n,m],'orientation':'COLUMN map R^70 -> R^36; coker generically zero over QQ(a,e,d)',
 'used_original_columns':used,'retained_original_rows':ids,'remaining_original_column_order':originalcols,
 'nonzero_entries':nonzero,'zero_entries':n*m-nonzero,'total_coefficient_terms':terms,'max_degree_a_e_d':degrees,'max_total_degree':totaldegree,'max_degree_e_d_over_field_QQ_a':eddegree,
 'max_coefficient_numerator_denominator_bits':maxbits,'distinct_positive_constant_denominators':len(denoms),'maximum_constant_denominator':str(max(denoms)),
 'all_input_denominators_32003_units':bad_den==0,'coefficient_numerators_divisible_by_32003':div_num,'no_row_identically_zero_by_scalar_32003_content':True,
 'function_field_grouped_terms_upper_bound_e_d':grouped,'largest_entry':largest_entry,'row_statistics':rows,
 'column_statistics':[{'column_index':j,'original_column_index':originalcols[j],'nonzero_entries':colnonzero[j],'terms':colterms[j]}for j in range(m)],
 'entry_termcount_histogram':[[k,v]for k,v in sorted(hist.items())],
 'provenance':{'alternative_ad_57_row_basis_required_units':['a','d','F','B','C'],'basis_triangular_minor':'a^3*d/10','first21_schur_pivot_and_clean_factors':['a','F','B','C'],'no_J_inversion':True,'no_parameter_denominators_in_saved_entries':True,'constant_QQ_scalars_in_prefix_not_reaudited_mod32003':True,'raw42x70_coker_has_generic_free_rank6_and_is_not_the_target':True},
 'proposed_field':'F32003(a)','proposed_polynomial_variables':['e','d','z'],'proposed_relation':'z*d*F*B*C-1 (a is already an invertible coefficient-field element)',
 'proposed_module_columns':'[M36x70 | (z*d*F*B*C-1)*I36], shape36x106; the relation must be appended in EVERY target component',
 'alternative_requested_relation_caveat':'z*a*F*B*C-1 defines saved coker without removing d=0; not an exhaustive actual cubic-defect model there because alternative ad basis needs d unit',
 'no_GB_field_point_or_source_regeneration_launched':True,
 'scope':'read-only exact source statistics and plan only; no generic-Fp support, finite-basis, annihilator or QQ-promotion theorem'}
 (out/'statistics.json').write_text(json.dumps(stats,indent=2)+'\n')
 wrapper=ROOT/'work/nonsquare-v62/singular-runtime/run-singular-optimized-v6.sh'
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),src,man,phase,code,basis,wrapper]
 report={'id':ID,'all_read_only_format_provenance_and_statistics_checks_passed':True,'shape':[n,m],'total_terms':terms,'max_degree_a_e_d':degrees,'max_total_degree':totaldegree,'max_e_d_degree':eddegree,'all_denominators_32003_units':True,'d_localization_required_for_actual_rank_comparison':True,'no_GB_or_field_job_launched':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'statistics_sha256':sha(out/'statistics.json'),'statistics_bytes':(out/'statistics.json').stat().st_size,'seconds':time.monotonic()-start,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args()
 try:main(a.output.resolve())
 except(TimeoutError,MemoryError)as error:
  a.output.mkdir(parents=True,exist_ok=True);(a.output/'BLOCKED.json').write_text(json.dumps({'completed':False,'constraint':type(error).__name__,'message':str(error),'no_auto_retry':True},indent=2)+'\n');raise
