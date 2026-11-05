#!/usr/bin/env python3
"""ONE exactABC72 linear job, then SOURCE ONLY two-distance/residual6 gap proof."""
import argparse,hashlib,heapq,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63';W=(1,2244,68);INF=10**12
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def two_distances(edges):
 dist=[[]for _ in range(6)];queue=[(0,4)]
 while queue:
  d,u=heapq.heappop(queue)
  if d in dist[u]or len(dist[u])==2:continue
  dist[u].append(d)
  if u==5:continue
  for v,w in edges[u]:
   if len(dist[v])<2 and d+w not in dist[v]:heapq.heappush(queue,(d+w,v))
 return dist[5],dist
class Guard(RuntimeError):pass
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={};current={};verify={};cursor={};done=adds=peak=0
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 def mod(rel,name):
  p=ROOT/rel;pins[rel]=sha(p);s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
 gate='reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC_SOURCE_PREFLIGHT_20261105_v63';gr=src(gate+'/report.json');assert gr['all_passed']and pins[gate+'/report.json']=='0815fa83d5bc481ee068c7da5e8e87d9292ddd27d6cbc64024cd72665d0d31d5';op=src(gate+'/all72-linear-operator.json');budget=src(gate+'/resource-gate.json');entries=op['entries'];assert len(entries)==72
 base='reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63';rt=src(base+'/retropicalized-source.json');row=src(base+'/actual-row11.json');R=[[Q(v)for v in rr]for rr in rt['row_transform']];C=[[Q(v)for v in rr]for rr in rt['column_transform']];DD=int(op['QQ_common_denominator']);Ds=list(map(int,op['row_LCMs']));rp=op['row_potentials'];cp=op['column_potentials']
 helper=mod('code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py','sourcegate');mp=mod('code/s61_nonsquare_prePOWER_K92_schur_gap_source_preflight_20261105_v63.py','minplus');bind=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02';M,inventory,kept,skip=helper.selected_inputs(state,bind['actual_pure_matrix_roots'])
 for p in row['columns']:M[1][p['column']]=dict(state_id=-1001-p['column'],c=Q(p['QQ_common_factor']),co=[(tuple(m),int(v))for m,v in p['primitive_integer_polynomial']])
 source_fixed=kept+row['actual_total_terms'];metas=[];actual_terms=0;maxbits=0;scalechecks=[];linefile=out/'ABC72.jsonl'
 def guard():
  nonlocal peak
  peak=max(peak,source_fixed+len(current)+len(verify)+15000)
  if peak>150000:raise Guard('BEFORE_COMMIT_owned')
  if time.monotonic()-start>45:raise Guard('45s_before_commit')
 try:
  with linefile.open('w')as stream:
   for e in entries:
    ei=e['entry_index'];u=e['global_row'];v=e['global_column'];current={};verify={};contributors=[];dump(out/'phase.json',dict(phase='ONEABC72 entry-stream',completed=done,next_entry=ei))
    for ci,c in enumerate(e['source_contributors']):
     i=c['source_row'];j=c['source_column'];p=M[i][j];assert p['state_id']==c['state_polynomial_id'];scale=int(c['integer_scale']);qscale=R[u][i]*C[j][v]*Ds[i]*p['c'];assert qscale*DD==scale;contributors.append((c,p,qscale));scalechecks.append(dict(entry=ei,contributor=ci,source_polynomial=p['state_id'],integer_scale=str(scale),exact_QQscale=str(qscale)))
     for si,(mon,z)in enumerate(p['co']):
      gap=rp[i]+cp[j]-sum(a*b for a,b in zip(mon,W));assert gap>=0;cursor=dict(entry=ei,stage='primary',contributor=ci,source_index=si,gap=gap,next_update_uncommitted=True);guard();value=current.get(gap,0)+scale*z;assert abs(value).bit_length()<=e['integer_coefficient_L1_bits']
      if value:current[gap]=value
      elif gap in current:del current[gap]
      adds+=1
    for ci,(c,p,qscale)in reversed(list(enumerate(contributors))):
     i=c['source_row'];j=c['source_column']
     for si in range(len(p['co'])-1,-1,-1):
      mon,z=p['co'][si];gap=rp[i]+cp[j]-sum(a*b for a,b in zip(mon,W));cursor=dict(entry=ei,stage='reverseQQ',contributor=ci,source_index=si,gap=gap,next_update_uncommitted=True);guard();value=verify.get(gap,Q(0))+qscale*z
      if value:verify[gap]=value
      elif gap in verify:del verify[gap]
    assert set(current)==set(verify)and all(verify[k]==Q(z,DD)for k,z in current.items());assert len(current)<=e['union_exponent_support_upper'];assert all(current.get(c['gap'],0)==int(c['integer_numerator'])for c in e['known_coefficients_to_recheck'])
    G=math.gcd(*(abs(z)for z in current.values()))if current else 1;content=Q(G,DD)if current else Q(1);co=[[k,str(z//G)]for k,z in sorted(current.items())];positive=sorted(k for k in current if k>0);val=positive[0]if positive else INF;maxbits=max(maxbits,max((abs(z).bit_length()for z in current.values()),default=0));actual_terms+=len(current)
    obj=dict(entry_index=ei,block=e['block'],block_row=e['block_row'],block_column=e['block_column'],global_row=u,global_column=v,QQ_content=str(content),sparse_integer_polynomial=co,formal_variable='t',exact_min_positive_gap=None if val==INF else val,known_coefficient_checks=e['known_coefficients_to_recheck'],source_contributors=e['source_contributors'],all_integer_and_reverseQQ_identities_verified=True);line=json.dumps(obj,separators=(',',':'))+'\n';assert len(line.encode())<=e['raw_serialization_upper'];assert stream.tell()+len(line.encode())<budget['all72_coefficient_RAW_upper'];offset=stream.tell();stream.write(line);stream.flush();metas.append(dict(entry_index=ei,block=e['block'],block_row=e['block_row'],block_column=e['block_column'],terms=len(current),first_two_positive_gaps=positive[:2],exact_min_positive_gap=val,constant_QQ_value=str(Q(current.get(0,0),DD)),line_byte_offset=offset,line_bytes=len(line.encode()),line_sha256=hashlib.sha256(line.encode()).hexdigest()));done+=1;current={};verify={};cursor=dict(stage='entry_boundary',completed=done,next_update_uncommitted=False)
  assert done==72 and adds==26556 and maxbits<=164
 except Guard as ex:
  dump(out/'partial-entry-prefix.json',dict(version='K92ABC72P1',complete=False,reason=str(ex),cursor=cursor,completed_entries=done,current_integer=[[k,str(z)]for k,z in sorted(current.items())],partial_reverseQQ=[[k,str(z)]for k,z in sorted(verify.items())],input_sha256=pins));dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,guard_reason=str(ex),completed_entries=done,partial_state_preserved=True,input_sha256=pins));signal.alarm(0);return
 dump(out/'all-scale-identities.json',dict(QQ_common_denominator=str(DD),all_integer_vs_QQscales_checked=True,identities=scalechecks));dump(out/'ABC72-index.json',dict(complete72=True,entries=metas,total_terms=actual_terms,max_integer_bits=maxbits))
 # SOURCE ONLY after the exact linear job: at most two positive edge weights
 # suffice for the two shortest DISTINCT walk costs (a third parallel edge has
 # two cheaper alternatives with the same prefix/suffix).
 gaps={'A':[[[]for j in range(4)]for i in range(4)],'B':[[[]for j in range(7)]for i in range(4)],'C':[[[]for j in range(4)]for i in range(7)]}
 for e in metas:gaps[e['block']][e['block_row']][e['block_column']]=e['first_two_positive_gaps']
 profbase='reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63';profile=src(profbase+'/report.json');assert profile['exact_H_leading_rank']==1;kcert=src(profbase+'/rank-kernel-unitminor.json');Hknown=src(profbase+'/actual-H-leading.json');R2=[[Q(v)for v in rr]for rr in kcert['row_transform_for_actualQQ_H']];Q2=[[Q(v)for v in rr]for rr in kcert['column_transform_for_actualQQ_H']];alpha,beta=src('reports/S61_NONSQUARE_PREPOWER_K92_CORRECTION_FACTOR_SOURCE_V2_20261105_v63/rank-support-zero-proof.json')['outer_row_potentials'],src('reports/S61_NONSQUARE_PREPOWER_K92_CORRECTION_FACTOR_SOURCE_V2_20261105_v63/rank-support-zero-proof.json')['outer_column_potentials']
 dbase='reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63';src(dbase+'/report.json');dfile=ROOT/dbase/'D49.jsonl';pins[str(dfile.relative_to(ROOT))]=sha(dfile);Drem=[[INF]*7 for i in range(7)]
 with dfile.open()as f:
  for line in f:
   e=json.loads(line);i=e['residual_row'];j=e['residual_column'];oldg=alpha[i]+beta[j];Drem[i][j]=min((int(k)for k,z in e['sparse_integer_polynomial']if int(k)!=oldg),default=INF)
 Pnext=[[INF]*7 for i in range(7)];pathproof=[]
 for i in range(7):
  for j in range(7):
   edges=[[]for k in range(6)]
   for k in range(4):
    edges[4]+=[(k,g)for g in gaps['C'][i][k]];edges[k]+=[(5,g)for g in gaps['B'][k][j]]
    for l in range(4):edges[k]+=[(l,g)for g in gaps['A'][k][l]]
   costs,labels=two_distances(edges);oldg=alpha[i]+beta[j];assert not costs or costs[0]>=oldg
   Pnext[i][j]=next((c for c in costs if c!=oldg),INF);pathproof.append(dict(row=i,column=j,old_degree=oldg,two_shortest_distinct_support_costs=costs,all_vertex_labels=labels,first_after_removed_P0=Pnext[i][j]))
 Hplus=[[min(Drem[i][j],Pnext[i][j])-(alpha[i]+beta[j])if min(Drem[i][j],Pnext[i][j])<INF else INF for j in range(7)]for i in range(7)];assert all(g>=1 for rr in Hplus for g in rr)
 Tplus=[[min((Hplus[i][j]for i in range(7)for j in range(7)if R2[u][i]and Q2[j][v]),default=INF)for v in range(7)]for u in range(7)];E=[[min(Tplus[1+i][1+j],Tplus[1+i][0]+Tplus[0][1+j])for j in range(6)]for i in range(6)];order,count=mp.min_assignment(E);aa,bb,match=mp.dual_min(E);total_order=26579+order;target=rt['new_tropical_degree']-total_order;z=target-5;mon=(5+z%68,(z//68)//33,(z//68)%33)
 dump(out/'two-distance-support-proof.json',dict(Aplus_first_two_positive_gaps=gaps['A'],B_first_two_positive_gaps=gaps['B'],C_first_two_positive_gaps=gaps['C'],parallel_edge_truncation_proof='anythirdparallelweight has two smaller distinct costs with sameprefix/suffix, so cannot contribute top2walkcosts',original_P0_removed_EXACTLY_from_previousfactorprofile=True,walks=pathproof,D_after_removed_D0_min=Drem,P_after_removed_P0_support_min=Pnext,Hplus_support_min=Hplus,NO_inverse_coefficients_computed=True))
 dump(out/'next-residual6-source.json',dict(typed_source='H=t^-alpha*(D-C*A^-1*B)*t^-beta, R2*H*Q2 hasconstantdiag(I1,0)',R2=R2,Q2=Q2,transformed_positive_part_min=Tplus,
  scalar_Uinverse_constant_unit='U0=1; leading inversepath cost0, positivecycles never lower cost',residual6_min=E,minimum_residual6_determinant_order=order,outer_order=26579,total_original_order_lower_bound=total_order,minimum_matching_count=count,row_dual_potentials=aa,column_dual_potentials=bb,matching=match,
  next_original_weight=target,next_unique_monomial_in_older_box=mon,NO_next_residual_coefficient_or_determinant_evaluated=True,next_goal='derive finite shortestpath coefficient requests at these exactA/B/C support costs, preserving both original4block correction and new1block Schur correction; no blindK/weight'))
 dump(out/'phase.json',dict(phase='FINAL ALL72exactABC complete; all121formalinput available; nextresidual6SOURCEonly'))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);actual=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert actual<budget['new_artifacts_plus500k_upper']and actual+external<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,ONE_ABC72_exactlinear_job=True,A16_B28_C28_complete=True,all72_integer_reverseQQ_identities_verified=True,A0I4_B0C0_and32previous_coefficients_verified=True,entry_stream_release=True,
  actual_ABC_terms=actual_terms,max_integer_numerator_bits=maxbits,scalar_additions_primary=adds,reverseQQ_additions=adds,peak_source_output_scratch=peak,all121_formal_input_complete=True,
  source_only_next_residual6_order=order,source_only_total_order_lower_bound=total_order,source_only_next_monomial=mon,NO_inverse_series_det_newweight_K_or_source572=True,no_det11_zero_nonzero_or_global57_claim=True,
  actual_new_RAW_before_report=actual,actual_external_read_RAW=external,actual_new_plus_read_RAW=actual+external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
