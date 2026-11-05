#!/usr/bin/env python3
"""ONE exact D49 linear job; entry-streamed; then SOURCE ONLY Schur valuations."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63';W=(1,2244,68);INF=10**12
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
class Guard(RuntimeError):pass
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={};cursor={};current={};verify={};done=0;adds=0;peak=0
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 def mod(rel,name):
  p=ROOT/rel;pins[rel]=sha(p);s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
 gate='reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D_SOURCE_PREFLIGHT_20261105_v63';g=src(gate+'/report.json');assert g['all_passed']and pins[gate+'/report.json']=='52e68cdc4811791a4b4b60d3fbaced40a851197bdc370395c3fa8e812d91e8aa';operator=src(gate+'/all49-linear-operator.json');budget=src(gate+'/resource-gate.json');entries=operator['entries'];assert len(entries)==49
 base='reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63';rt=src(base+'/retropicalized-source.json');row=src(base+'/actual-row11.json');R=[[Q(v)for v in rr]for rr in rt['row_transform']];C=[[Q(v)for v in rr]for rr in rt['column_transform']];DD=int(operator['QQ_common_denominator']);Ds=list(map(int,operator['row_LCMs']));rp=operator['row_potentials'];cp=operator['column_potentials']
 helper=mod('code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py','sourcegate');mp=mod('code/s61_nonsquare_prePOWER_K92_schur_gap_source_preflight_20261105_v63.py','minplus');support=src('reports/S61_NONSQUARE_PREPOWER_K92_SCHUR_GAP_SOURCE_PREFLIGHT_20261105_v63/minplus-proof.json')
 bind=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02';M,inventory,kept,skip=helper.selected_inputs(state,bind['actual_pure_matrix_roots'])
 for p in row['columns']:M[1][p['column']]=dict(state_id=-1001-p['column'],c=Q(p['QQ_common_factor']),co=[(tuple(m),int(v))for m,v in p['primitive_integer_polynomial']])
 source_fixed=kept+row['actual_total_terms'];metas=[];total_actual=0;maxbits=0;linefile=out/'D49.jsonl';scalechecks=[]
 def guard():
  nonlocal peak
  owned=source_fixed+len(current)+len(verify)+15000;peak=max(peak,owned)
  if owned>150000:raise Guard('BEFORE_COMMIT_owned')
  if time.monotonic()-start>45:raise Guard('45s_before_commit')
 try:
  with linefile.open('w')as stream:
   for e in entries:
    ei=e['entry_index'];u=e['global_row'];v=e['global_column'];current={};verify={};current_contributors=[]
    dump(out/'phase.json',dict(phase='ONEexactD entry-stream calculation',next_entry=ei,completed_entries=done))
    for ci,c in enumerate(e['source_contributors']):
     i=c['source_row'];j=c['source_column'];p=M[i][j];assert p['state_id']==c['state_polynomial_id'];scale=int(c['integer_scale']);qscale=R[u][i]*C[j][v]*Ds[i]*p['c'];assert qscale*DD==scale;current_contributors.append((c,p,qscale));scalechecks.append(dict(entry=ei,contributor=ci,state_polynomial_id=p['state_id'],integer_scale=str(scale),exact_QQscale=str(qscale)))
     for si,(mon,z)in enumerate(p['co']):
      gap=rp[i]+cp[j]-sum(a*b for a,b in zip(mon,W));assert gap>=0;cursor=dict(entry=ei,stage='primary',contributor=ci,source_index=si,gap=gap,next_update_uncommitted=True);guard();value=current.get(gap,0)+scale*z;assert abs(value).bit_length()<=e['integer_coefficient_L1_bits']
      if value:current[gap]=value
      elif gap in current:del current[gap]
      adds+=1
    # Reverse-order QQ replay is independent of the integer accumulation association.
    for ci,(c,p,qscale)in reversed(list(enumerate(current_contributors))):
     i=c['source_row'];j=c['source_column']
     for si in range(len(p['co'])-1,-1,-1):
      mon,z=p['co'][si];gap=rp[i]+cp[j]-sum(a*b for a,b in zip(mon,W));cursor=dict(entry=ei,stage='reverseQQverification',contributor=ci,source_index=si,gap=gap,next_update_uncommitted=True);guard();value=verify.get(gap,Q(0))+qscale*z
      if value:verify[gap]=value
      elif gap in verify:del verify[gap]
    assert set(current)==set(verify)and all(verify[k]==Q(z,DD)for k,z in current.items());assert len(current)<=e['union_exponent_support_upper'];assert all(current.get(k,0)==0 for k in e['certified_zero_coefficients_to_recheck']);assert current.get(0,0)==0
    G=math.gcd(*(abs(z)for z in current.values()))if current else 1;content=Q(G,DD)if current else Q(1);co=[[k,str(z//G)]for k,z in sorted(current.items())];valuation=min(current,default=INF);bits=max((abs(z).bit_length()for z in current.values()),default=0);maxbits=max(maxbits,bits);total_actual+=len(current)
    obj=dict(entry_index=ei,residual_row=e['residual_row'],residual_column=e['residual_column'],QQ_content=str(content),sparse_integer_polynomial=co,formal_variable='t',exact_valuation=None if valuation==INF else valuation,source_contributors=e['source_contributors'],all_coefficient_integer_and_reverseQQ_identity_checked=True,
     D0_zero_checked=True,previous_known_zero_coefficients_checked=e['certified_zero_coefficients_to_recheck']);line=json.dumps(obj,separators=(',',':'))+'\n';assert len(line.encode())<=e['raw_serialization_upper'];assert stream.tell()+len(line.encode())<budget['all49_coefficients_RAW_upper'];offset=stream.tell();stream.write(line);stream.flush();metas.append(dict(entry_index=ei,residual_row=e['residual_row'],residual_column=e['residual_column'],terms=len(current),valuation=valuation,leading_QQ_coefficient=None if valuation==INF else str(Q(current[valuation],DD)),line_byte_offset=offset,line_bytes=len(line.encode()),line_sha256=hashlib.sha256(line.encode()).hexdigest()));done+=1;current={};verify={};cursor=dict(stage='entry_boundary',completed_entries=done,next_update_uncommitted=False)
  assert done==49 and adds==36717 and maxbits<=168
 except Guard as ex:
  dump(out/'partial-entry-prefix.json',dict(version='K92D49P1',complete=False,reason=str(ex),cursor=cursor,completed_entries=done,current_integer_polynomial=[[k,str(z)]for k,z in sorted(current.items())],partial_reverseQQ=[[k,str(z)]for k,z in sorted(verify.items())],source_input_sha256=pins,do_NOT_use_partial_as_full_D=True))
  dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,guard_reason=str(ex),completed_entries=done,partial_entry_and_prior_stream_preserved=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()}));signal.alarm(0);return
 dump(out/'all-scale-identities.json',dict(QQ_common_denominator=str(DD),all_contributor_integer_vs_QQscales_verified=True,identities=scalechecks));dump(out/'D49-index.json',dict(entries=metas,complete49=True,total_terms=total_actual,coefficient_maxbits=maxbits,formalD_scope_only=True))
 # SOURCE ONLY: exact D valuations replace the loose Boolean values; do not
 # evaluate new Schur coefficients, determinant, or inverse series.
 D=[[INF]*7 for i in range(7)]
 for m in metas:D[m['residual_row']][m['residual_column']]=m['valuation']
 route=support['correction_C_Ainv_B_min'];S=[[min(D[i][j],route[i][j])for j in range(7)]for i in range(7)];order,count=mp.min_assignment(S);alpha,beta,match=mp.dual_min(S);tight=[]
 for i in range(7):
  for j in range(7):
   if S[i][j]<INF and alpha[i]+beta[j]==S[i][j]:tight.append(dict(row=i,column=j,gap=S[i][j],kind='knownDonly'if D[i][j]<route[i][j]else'correctionOnly'if route[i][j]<D[i][j]else'knownD_and_correction'))
 target=rt['new_tropical_degree']-order;z=target-5;mon=(5+z%68,(z//68)//33,(z//68)%33);inside=5<=mon[0]<=72 and 0<=mon[1]<=107 and 0<=mon[2]<=32
 dump(out/'next-source-Schur-valuation.json',dict(EXACT_D_valuations=D,C_Ainverse_B_support_lower_bounds_unchanged=route,new_Schur_support_lower_bounds=S,new_determinant_order_lower_bound=order,old_lower_bound=13255,new_minimum_assignment_count=count,row_dual_potentials=alpha,column_dual_potentials=beta,matching=match,tight_entries=tight,
  next_original_weight=target,next_unique_original_monomial=mon,monomial_in_box=inside,Schur_correction_RETAINED=True,NO_new_Schur_coefficient_OR_determinant_evaluated=True,next_goal='Build finite shortest-Ainverse-path coefficient DAG only for tight edges; use exactsavedD where present. Fresh product/height/source gate BEFORE any next arithmetic.'))
 dump(out/'phase.json',dict(phase='FINAL ALL49exactD complete; updatedSchurSOURCEonly'))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);actual=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert actual<budget['new_artifacts_plus500k_upper']and actual+external<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,ONE_all49_exactD_linear_job=True,all49_integer_and_reverseQQ_whole_coefficient_identities=True,all49_D0_and41_prior_zero_checks=True,entry_stream_release=True,
  actual_D_total_terms=total_actual,nonzero_D_entries=sum(m['terms']>0 for m in metas),max_integer_numerator_bits=maxbits,scalar_additions_primary=adds,reverseQQ_additions=adds,peak_source_output_and_scratch=peak,
  source_only_next_order_lower_bound=order,source_only_next_unique_monomial=mon,tight_next_edges=len(tight),correction_involved_next_edges=sum(t['kind']!='knownDonly'for t in tight),no_new_Schurcoef_inverse_series_det_weight_K_or_source572=True,no_det11_zero_nonzero_or_global57_claim=True,
  actual_new_RAW_before_report=actual,actual_external_read_RAW=external,actual_new_plus_read_RAW=actual+external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
