#!/usr/bin/env python3
"""SOURCE ONLY: 70 fixed P4 Laplace groups, reuse K92SP001, no scalar products."""
import argparse, hashlib, itertools, json, math, resource, signal, struct, time
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_K92_70_LAPLACE_SOURCE_PREFLIGHT_20261104_v63'
PREV='S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63'
EXPECTED_STATE='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02'
W=((1,1,1),(2,1,0),(2,1,1),(1,1,2),(2,1,2),(1,2,1),(1,0,1),(0,1,1),(1,1,0))
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def read(p):return json.loads(p.read_text())
def qgcd(a,b):return Q(math.gcd(a.numerator,b.numerator),math.lcm(a.denominator,b.denominator))

def selective_state(path,wanted):
 """No DP coefficient is allocated; skip its variable-length integer byte record."""
 selected={};inventory=[];kept=skipped=0;peak_record=0
 with path.open('rb') as f:
  def get(fmt):
   n=struct.calcsize('<'+fmt);b=f.read(n);assert len(b)==n;return struct.unpack('<'+fmt,b)[0]
  def z(limit,keep=True):
   n=get('H');sgn=get('B');assert sgn in (0,1) and n<=(limit+7)//8
   if not keep:f.seek(n,1);return None
   b=f.read(n);assert len(b)==n and (not n or b[0]);v=int.from_bytes(b,'big');assert v.bit_length()<=limit
   return -v if sgn else v
  assert f.read(8)==b'K92SP001';n=get('I');assert n<65536;header=json.loads(f.read(n));np=get('I');assert np<=1024
  assert header['source_completed']==572 and not header['complete']
  for _ in range(np):
   pid=get('I');role=get('B');src=get('i');term=get('i');row=get('i');mask=get('i');complete=get('B')
   retain=src in wanted or (role==2 and term in range(4))
   num=z(4096,retain);den=z(4096,retain);nt=get('I');peak_record=max(peak_record,int(retain))
   co=[] if retain else None
   for j in range(nt):
    m=(get('i'),get('i'),get('i'));v=z(1024,retain)
    if retain:assert v and (not co or co[-1][0]<m);co.append((m,v))
   inventory.append(dict(id=pid,role=role,source=src,term=term,row=row,mask=mask,complete=bool(complete),terms=nt,retained=retain))
   if retain:
    assert complete;factor=Q(num,den);assert num>0 and den>0 and factor.numerator==num and factor.denominator==den
    selected[pid]=dict(source=src,role=role,term=term,coefficients=co,factor=factor);kept+=nt
   else:skipped+=nt
  pending=get('B');assert pending in (0,1)
  if pending:
   for _ in range(3):get('i')
   for _ in range(2):get('Q')
   for _ in range(3):get('q')
   z(16384,False);get('B');z(16384,False)
  assert not f.read(1)
 return header,selected,inventory,kept,skipped,peak_record

@lru_cache(None)
def count_lattice(box,degree):
 lo,hi=degree;shift=sum(a for a,b in box);width=[b-a+1 for a,b in box]
 def upto(t):
  t-=shift;s=0
  for mask in range(8):
   n=t-sum(width[j] for j in range(3) if mask>>j&1)
   if n>=0:s+=(-1 if mask.bit_count()%2 else 1)*math.comb(n+3,3)
  return s
 return upto(hi)-upto(lo-1)

def atom_meta(p):
 co=p['coefficients']
 if not co:return None
 return dict(support=len(co),box=tuple((min(m[j] for m,v in co),max(m[j] for m,v in co))for j in range(3)),degree=(min(sum(m)for m,v in co),max(sum(m)for m,v in co)),
  linear=tuple((min(sum(w[j]*m[j]for j in range(3))for m,v in co),max(sum(w[j]*m[j]for j in range(3))for m,v in co))for w in W),c=p['factor'],L=sum(abs(v)for m,v in co),bits=max(abs(v).bit_length()for m,v in co))

ONE=dict(support=1,box=((0,0),)*3,degree=(0,0),linear=((0,0),)*9,c=Q(1),L=1,bits=1)
def mul(a,b):
 if a is None or b is None:return None
 box=tuple((x[0]+y[0],x[1]+y[1])for x,y in zip(a['box'],b['box']));degree=tuple(x+y for x,y in zip(a['degree'],b['degree']));L=a['L']*b['L']
 return dict(support=min(a['support']*b['support'],L,count_lattice(box,degree)),box=box,degree=degree,linear=tuple((x[0]+y[0],x[1]+y[1])for x,y in zip(a['linear'],b['linear'])),c=a['c']*b['c'],L=L,bits=L.bit_length())
def add(a,b):
 if a is None:return b
 if b is None:return a
 c=qgcd(a['c'],b['c']);r=a['c']/c;s=b['c']/c;assert r.denominator==s.denominator==1
 box=tuple((min(x[0],y[0]),max(x[1],y[1]))for x,y in zip(a['box'],b['box']));degree=(min(a['degree'][0],b['degree'][0]),max(a['degree'][1],b['degree'][1]));L=a['L']*r.numerator+b['L']*s.numerator
 return dict(support=min(a['support']+b['support'],L,count_lattice(box,degree)),box=box,degree=degree,linear=tuple((min(x[0],y[0]),max(x[1],y[1]))for x,y in zip(a['linear'],b['linear'])),c=c,L=L,bits=L.bit_length())
def det_profile(M):
 """Only support/content/height arithmetic, never coefficient multiplication."""
 k=len(M);prev={0:ONE};peak=1;pair_ops=0;stages=[]
 for row in range(k):
  nex={};old_live=sum(x['support']for x in prev.values())
  for mask,a in sorted(prev.items()):
   for j,b in enumerate(M[row]):
    if mask>>j&1 or b is None:continue
    pair_ops+=a['support']*b['support'];nm=mask|1<<j;nex[nm]=add(nex.get(nm),mul(a,b));peak=max(peak,old_live+sum(x['support']for x in nex.values()))
   old_live-=a['support']
  prev=nex;stages.append(dict(row=row+1,layer_support=sum(x['support']for x in prev.values()),states=len(prev)))
 result=prev.get((1<<k)-1)
 return result,dict(peak_cells_upper=peak,pair_updates_upper=pair_ops,stages=stages)

def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):
  p=ROOT/rel;pins[rel]=sha(p);return read(p)
 prior=src(f'reports/{PREV}/report.json');assert prior['safety_and_restore_checks_passed'] and not prior['full_scalar_complete']
 binding=src(f'reports/{PREV}/source-binding.json');roots=binding['actual_pure_matrix_roots'];schedule=binding['normative_sign_schedule'];assert len(schedule)==10
 norm=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/normalization.json');globalc=Q(norm['content_lower_divisor']);ashift=norm['known_output_a_power']
 frozen=ROOT/'reports'/PREV/'state.bin';assert sha(frozen)==EXPECTED_STATE;pins[str(frozen.relative_to(ROOT))]=EXPECTED_STATE
 wanted={i for r in roots for i in r}|set(range(593,599));dump(out/'phase.json',dict(phase='selective state scan; no DP allocation/no source arithmetic'))
 header,kept,inventory,keptcells,skipcells,record=selective_state(frozen,wanted);by_source={x['source']:pid for pid,x in kept.items()if x['source']>=0};literal={x['term']:pid for pid,x in kept.items()if x['role']==2}
 assert len(literal)==4 and all(i in by_source for i in wanted)
 meta={pid:atom_meta(p)for pid,p in kept.items()};matrixids=[[by_source[i]for i in r]for r in roots];M=[[meta[i]for i in r]for r in matrixids]
 weightids=[literal[0],literal[1],by_source[593],literal[2],literal[3],by_source[594],by_source[597],by_source[596],by_source[595],by_source[598]];weights=[meta[i]for i in weightids]
 matrix_unique={i for r in matrixids for i in r};assert matrix_unique|set(weightids)==set(kept)
 fixed=keptcells;peakbits=max(x['bits']for x in meta.values()if x);input_factorbits=max(max(x['c'].numerator.bit_length(),x['c'].denominator.bit_length())for x in meta.values()if x)
 dump(out/'state-reuse.json',dict(state_sha256=EXPECTED_STATE,original_header=header,retained_polynomials=len(kept),retained_terms=fixed,skipped_terms_without_allocation=skipcells,one_record_scratch=record,all_DP_data_skipped=True,matrix_state_ids=matrixids,weight_state_ids=weightids,source_nodes_NEVER_reevaluated=True,inventory=inventory,retained_metadata=meta))
 # Do not retain exact coefficient payloads during metadata propagation.
 del kept
 groups=[];workmax=0;updates=0;outputs=None;integer_scalars=[];maxbits=peakbits;maxfactorbits=input_factorbits
 for gi,cc in enumerate(itertools.combinations(range(8),4)):
  other=[j for j in range(8)if j not in cc];top,pt=det_profile([[M[i][j]for j in cc]for i in range(4)]);assert top is not None
  terms=[];gwork=pt['peak_cells_upper'];gupdates=pt['pair_updates_upper'];combined=None
  for ti,t in enumerate(schedule):
   rr=[4]+[5+i for i in t['pure_remaining_lower6']];bottom,pb=det_profile([[M[i][j]for j in other]for i in rr]);assert bottom is not None
   weighted=mul(weights[ti],bottom);allterm=mul(top,weighted);scale=allterm['c']/globalc;assert scale.denominator==1
   # Buffer strategy: P4 stored, bottom4 DP streamed; then weighted-bottom stored,
   # then streamed P4*weighted-bottom directly to one checked global accumulator.
   stage=max(pt['peak_cells_upper'],top['support']+pb['peak_cells_upper'],top['support']+bottom['support']+weighted['support'])
   gwork=max(gwork,stage);gupdates+=pb['pair_updates_upper']+weights[ti]['support']*bottom['support']+top['support']*weighted['support'];outputs=add(outputs,allterm);combined=add(combined,weighted)
   maxbits=max(maxbits,allterm['bits'],weighted['bits'],bottom['bits'],top['bits']);maxfactorbits=max(maxfactorbits,max(allterm['c'].numerator.bit_length(),allterm['c'].denominator.bit_length()))
   terms.append(dict(term=ti,H3_selected=t['selected_H3_lower6'],pure_bottom_state_rows=rr,sign=t['combined_multiplier']*((-1)**sum(cc)),bottom4=bottom,bottom4_DP=pb,weighted_bottom=weighted,temporary_cells_upper=stage,global_integer_scale=str(scale.numerator)))
  workmax=max(workmax,gwork);updates+=gupdates;groups.append(dict(group=gi,top4_columns=list(cc),bottom4_columns=other,Laplace_sign=(-1)**sum(cc),top4=top,top4_DP=pt,ten_H3_terms=terms,working_without_fixed_or_accumulator_upper=gwork,pair_updates_upper=gupdates))
  dump(out/'phase.json',dict(phase='SOURCE ONLY 70 Laplace-group support',completed_groups=gi+1))
 assert len(groups)==70 and ashift==5
 # Full output remains a loose bound; dynamic before-commit guards promise only safety.
 fullupper=outputs['support'];CELL=135000;WORK=15000;POLYS=512;coefbits=1024;factorbits=4096
 safe_bytes=CELL*144+POLYS*1216+4096+3*1024*1024+2*1024*1024+1024*1024
 assert safe_bytes<32*1024*1024 and fixed<CELL and maxbits<coefbits and maxfactorbits<factorbits
 dump(out/'groups70.json',dict(all70_groups=groups,all10_H3_terms_retained=True,source_pure8_order=[22,11,13,25,15,19,14,18],source_lower6_order=['A','B','D','C','E','F'],no_coefficient_products_evaluated=True))
 recipe=dict(expression='sum_C (-1)^sum(C) det(P4[:,C])*sum_t sign[t]*H3weight[t]*det([G20;remaining_lower3[t]][:,C_complement])',
  execution_order='C lexicographic70; compute one P4 cofactor; t=0..9, compute one bottom4 minor with consumed previous-layer release; multiply weight*bottom4 into one temporary; stream P4*temporary into global accumulator; release bottom4/weighted temporary before next t; release P4 before next C',
  state_reuse='Select state.bin only pinned matrix88 roots and ten weights; do not allocate DP records or rerun572sourceoperations. Source-factor normalization from original state is preserved.',
  signs='normative source-binding ten signs times (-1)^sum(top4columns); all rows/columns retain original order',
  normalization=dict(common_QQfactor=str(globalc),known_output_a_shift=ashift,g_kept_factored=True,known_source_a_squared_already_removed=True),
  checkpoints='Every addmul coefficient BEFORE-COMMIT: ownedpool135000/work15000,coef1024bits,commonQQ4096bits,exponentint32,pending16384bits,clock45s; store all active source/cofactor/bottomDP/temp/accumulator polynomials and group/term/row/mask/operandcursor/pendingupdate ONCE as binary, not duplicateJSON',
  exact_prefix_restore='All active factors/integers/cursors. A partial sum or nonzero prefix is NOT det11/nonzero proof. Need ALL70groups and ALL10terms.',
  no_scalar_run_in_this_packet=True,no_completion_promise=True,no_new_chart_or_candidate=True)
 dump(out/'operator-recipe.json',recipe)
 resource_gate=dict(fixed_input_actual_terms=fixed,fixed_matrix_unique_polynomials=len(matrix_unique),fixed_H3_weight_polynomials=len(weightids),old_saved_DP_terms_NOT_allocated=skipcells,
  maximum_unconstrained_working_cells_upper=workmax,maximum_owned_WITHOUT_global_accumulator_upper=fixed+workmax+WORK,remaining_global_accumulator_at_worst_work=CELL-fixed-workmax,
  ordinary_output_support_upper=fullupper,ordinary_completion_memory_guaranteed=fixed+workmax+fullupper+WORK<=150000,all_pair_updates_combinatorial_upper=updates,
  all_full_group_polynomial_bounds_metadata_ONLY=True,coefficient_L1_bit_upper=maxbits,QQfactor_bit_upper=maxfactorbits,
  proposed_hard_owned_numeric_cell_cap=CELL,work_reserve=WORK,total_live_hardcap=CELL+WORK,proposed_poly_cap=POLYS,binary_total_RAW_upper=safe_bytes,
  input_binary_referenced_unchanged_NOT_copied=True,no_large_new_source_input=True,guarded_safety_gate_passes=True,completion_and_runtime_NOT_guaranteed=True)
 dump(out/'resource-gate.json',resource_gate)
 note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 dump(out/'phase.json',dict(phase='FINAL SOURCE ONLY 70-group gate; NO scalar run'))
 total=sum(p.stat().st_size for p in out.iterdir());assert total<3*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,all572source_and_four_H3cofactors_reused=True,old_DP_skipped_without_coefficient_allocation=True,new_coefficient_products_NOT_performed=True,groups=70,H3terms_per_group=10,
  resource_gate=resource_gate,total_RAW_bytes_before_report=total,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-began)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2,default=str))

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
