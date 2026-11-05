#!/usr/bin/env python3
"""SOURCE ONLY: precise P01 continuation, then bottom-tail sum before P4 product."""
import argparse,hashlib,importlib.util,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_K92_70_TAIL_ASSOCIATION_PREFLIGHT_20261105_v63'
P01='S61_NONSQUARE_PREPOWER_K92_70_LAPLACE_GUARDED_ATTEMPT_20261105_v63'
STATE_SHA='b979ac727295fe6689bd8cd6805189899271ce7c124103e4479728268dce9aa7'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def module(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def source(rel):
  p=ROOT/rel;pins[rel]=sha(p);return read(p)
 def loadmodule(rel,name):
  p=ROOT/rel;pins[rel]=sha(p);return module(p,name)
 prof=loadmodule('code/s61_nonsquare_prePOWER_K92_70_laplace_source_preflight_20261104_v63.py','bounds70')
 parser=loadmodule('code/s61_nonsquare_prePOWER_K92_restore70_sparse_state_20261105_v63.py','p01reader')
 pr=source(f'reports/{P01}/report.json');assert pr['safety_and_restore_checks_passed']and not pr['full_scalar_complete']
 statepath=ROOT/'reports'/P01/'state.bin';assert sha(statepath)==STATE_SHA;pins[str(statepath.relative_to(ROOT))]=STATE_SHA
 dump(out/'phase.json',dict(phase='read P01 exact state once; no coefficient products'))
 state=parser.read_state(statepath);h=state['header'];pending=state['pending'];polys=state['polynomials'];assert h['stage']=='global_accumulate'and h['group']==0 and h['term']==3 and h['groups_completed']==0 and h['terms_completed']==3
 assert pending is not None and not pending['value_known']and (pending['dst'],pending['lhs'],pending['rhs'])==(h['acc_id'],h['top_id'],h['weighted_id'])
 acc=polys[h['acc_id']];top=polys[h['top_id']];wb=polys[h['weighted_id']];i=pending['left'];j=pending['right'];a=top['coefficients'][i][0];b=wb['coefficients'][j][0];shift=h['known_output_a_shift'];assert shift==5
 assert pending['key']==(a[0]+b[0]-shift,a[1]+b[1],a[2]+b[2]);scale=top['factor']*wb['factor']/acc['factor'];assert scale.denominator==1 and pending['scale']==-scale.numerator
 remaining_pairs=(len(top['coefficients'])-i-1)*len(wb['coefficients'])+(len(wb['coefficients'])-j)
 matrixids=h['matrix_ids'];weightids=h['weight_ids'];fixedids={x for r in matrixids for x in r}|set(weightids);fixed=sum(len(polys[x]['coefficients'])for x in fixedids);assert fixed==25266
 M=[[prof.atom_meta(polys[x])for x in r]for r in matrixids];weights=[prof.atom_meta(polys[x])for x in weightids];T0=prof.atom_meta(top);A0=prof.atom_meta(acc);WB0=prof.atom_meta(wb)
 # Future exact source support of the two contiguous pieces of the pending product.
 def slice_meta(p,begin,end):
  return prof.atom_meta(dict(coefficients=p['coefficients'][begin:end],factor=p['factor']))
 def ashift(meta,k):
  if meta is None:return None
  return dict(meta,box=((meta['box'][0][0]+k,meta['box'][0][1]+k),)+tuple(meta['box'][1:]),degree=tuple(x+k for x in meta['degree']),linear=tuple((x+k*w[0],y+k*w[0])for(x,y),w in zip(meta['linear'],prof.W)))
 oneleft=slice_meta(top,i,i+1);righttail=slice_meta(wb,j,len(wb['coefficients']));lefttail=slice_meta(top,i+1,len(top['coefficients']));remaining=prof.add(prof.mul(oneleft,righttail),prof.mul(lefttail,WB0));resume_bound=prof.add(A0,ashift(remaining,-shift))
 resume=dict(origin_state_sha256=STATE_SHA,state_version='K9270P01',stage='global_accumulate',group=0,term=3,completed_H3_terms=3,completed_groups=0,dst=pending['dst'],lhs=pending['lhs'],rhs=pending['rhs'],left_index=i,right_index=j,pending_key=pending['key'],integer_scale=pending['scale'],
  pending_value_known=False,serialized_scratch_integer_MUST_NOT_be_committed=True,first_recomputed_update_exactly_once=True,original_top_terms=len(top['coefficients']),original_weighted_bottom_terms=len(wb['coefficients']),original_accumulator_terms=len(acc['coefficients']),remaining_coefficient_pairs=remaining_pairs,
  restored_owned_terms=state['total_terms'],restored_active_polys=len(polys),resume_to_term_boundary_accumulator_support_upper=resume_bound['support'],resume_to_boundary_L1_bits=resume_bound['bits'],NO_undo_or_division_of_partial_accumulator=True,completed_first_three_terms_NOT_repeated=True)
 dump(out/'P01-exact-resume.json',resume)
 # All further operations are metadata bounds. The exact old coefficient payload is released.
 del state,polys,acc,top,wb,oneleft,righttail,lefttail,A0,WB0
 binding=source('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');schedule=binding['normative_sign_schedule'];assert len(schedule)==10
 globalc=Q(source('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/normalization.json')['content_lower_divisor']);assert globalc==resume_bound['c']
 output=None;peak_work=0;updates=remaining_pairs;profiles=[];maxbits=resume_bound['bits'];maxfactor=0;gap_before_product=[]
 for gi,cc in enumerate(itertools.combinations(range(8),4)):
  other=[j for j in range(8)if j not in cc]
  if gi==0:T=T0;tp=dict(peak_cells_upper=0,pair_updates_upper=0,reused_actual_P4=True);first=4
  else:T,tp=prof.det_profile([[M[i][j]for j in cc]for i in range(4)]);first=0
  assert T is not None;tail=None;terms=[];work=tp['peak_cells_upper'];ops=tp['pair_updates_upper'];previous_support=0
  for ti in range(first,10):
   t=schedule[ti];rr=[4]+[5+r for r in t['pure_remaining_lower6']];B,bp=prof.det_profile([[M[r][c]for c in other]for r in rr]);assert B is not None
   weighted=prof.mul(weights[ti],B);tail=prof.add(tail,weighted);stream_scale=weighted['c']/tail['c'];assert stream_scale.denominator==1
   # The final common tail factor is selected BEFORE any coefficient work below;
   # the transient metadata tail factors here only supply source bounds.
   stage=max(T['support']+previous_support+bp['peak_cells_upper'],T['support']+tail['support']+B['support'])
   work=max(work,stage);ops+=bp['pair_updates_upper']+weights[ti]['support']*B['support'];previous_support=tail['support']
   terms.append(dict(term=ti,H3_selected=t['selected_H3_lower6'],bottom_rows=rr,H3_sign=t['combined_multiplier'],bottom4=B,bottom4_DP=bp,weighted_bottom=weighted,tail_after_term=tail,working_without_fixed_and_global_upper=stage));maxbits=max(maxbits,tail['bits'],B['bits'],weighted['bits'])
  assert tail is not None
  # Exact finite arithmetic common-factor recipe: gcd_Q of all products of the
  # predetermined DP-content factors, independent of polynomial coefficients.
  tailc=tail['c']
  for t in terms:
   t['integer_stream_scale_to_FIXED_tail_content']=str((t['weighted_bottom']['c']/tailc).numerator);assert(t['weighted_bottom']['c']/tailc).denominator==1
  product=prof.mul(T,tail);outer_scale=product['c']/globalc;assert outer_scale.denominator==1
  ops+=T['support']*tail['support'];work=max(work,T['support']+tail['support']);peak_work=max(peak_work,work);updates+=ops;output=prof.add(output,ashift(product,-shift));maxbits=max(maxbits,product['bits']);maxfactor=max(maxfactor,max(product['c'].numerator.bit_length(),product['c'].denominator.bit_length()))
  profiles.append(dict(group=gi,top4_columns=list(cc),bottom4_columns=other,first_tail_H3_term=first,terms=terms,tail_common_QQfactor=str(tailc),tail_support_upper=tail['support'],top4=T,top4_profile=tp,column_Laplace_sign=(-1)**sum(cc),
   integer_outer_scale=str(outer_scale.numerator),working_without_fixed_global_upper=work,coefficient_pairs_upper=ops,tail_NOT_multiplied_by_P4_until_all_terms_complete=True))
  gap_before_product.append(dict(group=gi,tail_support_upper=tail['support'],available_global_cells_at_loose_peak=135000-fixed-work))
  dump(out/'phase.json',dict(phase='SOURCE ONLY tail-association profiling',groups_profiled=gi+1))
 assert len(profiles)==70
 finalbound=prof.add(resume_bound,output);assert finalbound['c']==globalc
 # Guarded correctness does not assert the loose peak or full output fits.
 CELL=135000;WORK=15000;POLYS=512;safe_bytes=CELL*144+POLYS*1216+4096+6*1024*1024
 assert safe_bytes<32*1024*1024 and maxbits<1024 and maxfactor<4096
 gate=dict(source_only=True,coefficient_products_NOT_performed=True,old572sourceoperations_NOT_repeated=True,old_first3_complete_terms_NOT_repeated=True,old_partial_fourth_term_NO_undo=True,
  P01_exact_pending_remaining_pairs=remaining_pairs,P01_starting_owned_cells=resume['restored_owned_terms'],P01_starting_plus_work=resume['restored_owned_terms']+WORK,
  reused_fixed_input_terms=fixed,new_top_level_group_products=70,old_top_level_group_products=697,group0_tail_terms='4..9 only AFTER originalterm3 reaches its boundary',other_group_tail_terms='0..9',
  maximum_unconstrained_temporary_work_upper=peak_work,work_fixed_reserve_WITHOUT_global_accumulator=fixed+peak_work+WORK,whole_completion_support_upper=finalbound['support'],ordinary_completion_guaranteed=False,
  remaining_coefficient_pairs_upper=updates,integer_L1_bit_upper=maxbits,QQfactor_bit_upper=maxfactor,owned_numeric_hardcap=CELL,metadata_and_scratch_reserve=WORK,total_hardcap=CELL+WORK,
  common_QQfactor_bits_stop=4096,integer_coefficient_bits_stop=1024,pending_integer_bits_stop=16384,clock_seconds=45,active_poly_cap=POLYS,RAW_prefix_upper=safe_bytes,total_RAW_hardcap=33554432,
  all_active_coefficients_and_cursor_saved_BEFORE_overflow=True,guarded_safety_gate_passes=True,CPU_and_completion_NOT_promised=True,loose_peak_exceeds_pool=fixed+peak_work+WORK>150000)
 dump(out/'groups70-tail-source.json',dict(groups=profiles,all_original_H3_signs_retained=True,no_extra_chart=True,g_remains_factored=True));dump(out/'resource-gate.json',gate)
 recipe=dict(version='K9270P02',origin_version='K9270P01',origin_state_sha256=STATE_SHA,
  phase0='Restore all101 exact P01 activepolys with originalIDs/commonfactors. DO NOT commit scratch with value_known=false. Recompute pending top434*weighted16099 once, then only the lexicographic remaining product pairs; preserve every prior acc coefficient.',
  phase1='At exact term3 boundary increment completedterms to4/releaseweighted182. Reuse P4id114. Sum signed H3weight[t]*bottom4[t] for group0 t4..9 into ONE fixed-content tail accumulator. No P4 multiplication per term. Then stream P4*tail into existing global accumulator once, with a^-5 and exact commonQQscale.',
  phase2='For groups1..69 compute one P4; sum all10 signed H3weight*bottom4 into tail FIRST, then multiply P4 once with column-Laplace sign. Retain all ten terms and all original source/tails.',
  tail_content='Before starting any group choose gcd_Q of all positive H3factor*bottomDPfactor products; per-term integer ratios recorded in groups70-tail-source.json. No coefficient gcd/unknown division.',
  required_state_fields='P02 stage=resume_P01_outer/build_tail/bottom_DP/new_global_accumulate; inherited P01origin hash,cursor and originalIDs; currentgroup,tail term,tailid,topid,bottomid,accid; DP row/mask/order; operation dst/lhs/rhs/left/right/key/scale/value_known; all active integerpolys/commonQQfactors. No auto continuation.',
  checkpoint='One K9270P02 binary, every write prechecked to totalRAW32MiB; no duplicated JSON coefficients/atomic doublebuffer. New parser must restore all exact state but not resume.',
  coefficient_guards='Every term insertion,coefficient or factor growth BEFORE commit:135000numeric+15000allmetadata/scratch,1024bitcoeff,4096bitfactor,16384bitpending,int32exponents,45secondclock. Customcheckedops only.',
  exactness='P01prefix + exact remaining term3 + P4*sum(t4..9) + sum_{groups1..69} P4*sum(t0..9) is algebraically identical to original70×10 expression. No old term/source replay,undo,division of acc,stale update or omitted term.',
  proof_scope='Prefix is not det/nonzero. ALL70groups and all original700 H3terms must be accounted for unless a separate source noncancellation certificate is proved.',scalar_run_authorized_by_this_packet=False)
 dump(out/'P02-operator-recipe.json',recipe)
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
 dump(out/'phase.json',dict(phase='FINAL source-only P02 association gate; no scalar continuation'))
 total=sum(p.stat().st_size for p in out.iterdir());assert total<4*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_new_coefficient_arithmetic=True,exact_P01_boundary_resume_pinned=True,all70groups_and_all10signs_preserved=True,guarded_safety_gate_passes=True,ordinary_completion_guaranteed=False,resource_gate=gate,
  total_RAW_bytes_before_report=total,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-began)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
