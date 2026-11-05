#!/usr/bin/env python3
"""ONE whole seven-component row identity, lexicographic stopfirstnonzero."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_WHOLE_ROW_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
class GuardStop(Exception):pass
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={};pool={};literal_names=set();peak=0;paircount=0;ops=[];current=None;combine_context=None;results=[]
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_ROW_SOURCE_20261105_v63';r=src(base+'/report.json');assert r['all_passed']and pins[base+'/report.json']=='c2f2216d2597d5bfeedac89a939fce55706e74800093effcf6e8d0cca2224064';gate=src(base+'/resource-gate.json');plan=src(base+'/whole-row-source-and-operator.json');receipt=src('reports/S61_NONSQUARE_PREPOWER_K92_WHOLE_ROW_WORKER_GATE_20261105_v63/report.json');assert receipt['all_passed']and receipt['whole_owned_scalar_upper']==121965
 L=int(plan['uniform_integer_denominator']);kappa=59375;source_records=[];input_terms=0
 for rel in['reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/ABC72.jsonl','reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/D49.jsonl']:
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for ln,line in enumerate(f):
    e=json.loads(line);b=e.get('block','D');i=e['block_row']if b!='D'else e['residual_row'];j=e['block_column']if b!='D'else e['residual_column'];c=Q(e['QQ_content']);assert L%c.denominator==0;fac=c.numerator*(L//c.denominator);key=f'{b}{i}_{j}';pp={int(g):int(v)*fac for g,v in e['sparse_integer_polynomial']};pool[key]=pp;literal_names.add(key);input_terms+=len(pp);source_records.append(dict(name=key,path=rel,line=ln,entry_index=e['entry_index'],QQ_content=str(c),cleared_integer_factor=str(fac),terms=len(pp)))
 assert input_terms==23932
 typ=src('reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63/typed-next-residual-Schur.json');assert typ['alpha'][0]-typ['alpha'][2]==68
 norm=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER4_DIRECT_SOURCE_20261105_v63/direct-border4-source.json')['R2_normalized'];assert list(map(Q,norm[1]))==[Q(i==0)for i in range(7)]and list(map(Q,norm[2]))==[Q(i==2)for i in range(7)]
 actual=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER4_DIRECT_ENTRY_20261105_v63/literal-source-and-all-finite-path-identities.json');T={(e['row'],e['column']):{int(g):int(v)for g,v in e['integer_polynomial']}for e in actual['all_Tbar']};assert all(T[2,j]=={g+68:kappa*v for g,v in T[1,j].items()}for j in range(7));assert all(not pool[f'A{i}_{j}']for i in range(2)for j in range(2,4))
 assert all(pool[f'A{i}_{j}'].get(0,0)==(L if i==j else 0)for i in range(4)for j in range(4))
 bound_products={e['label']:e for e in gate['all_product_operations']};operation_id=0
 def guard(extra=0):
  nonlocal peak
  n=sum(len(p)for p in pool.values())+extra+15000;peak=max(peak,n)
  if n>121965:raise GuardStop('owned_scalar_guard')
  if time.monotonic()-start>45:raise GuardStop('45s_guard')
 def mul(a,b,label):
  nonlocal paircount,current,operation_id
  p=pool[a];q=pool[b];ap=sorted(p);bp=sorted(q);dst={};bound=bound_products[label];operation_id+=1
  current=dict(operation=operation_id,label=label,left=a,right=b,lhs_index=0,rhs_index=0,partial_product=dst,next_pair_not_committed=True,parent_combine=combine_context)
  for ii,g in enumerate(ap):
   x=p[g]
   for jj,h in enumerate(bp):
    if paircount%4096==0:
     current['lhs_index']=ii;current['rhs_index']=jj;guard(len(dst)+2)
    z=g+h;v=dst.get(z,0)+x*q[h]
    assert abs(v).bit_length()<=bound['coefficient_abs_bound']
    if v:dst[z]=v
    else:dst.pop(z,None)
    paircount+=1
   assert len(dst)<=bound['actual_support_upper']
  guard(len(dst)+2);current=None;return dst
 def add_into(dst,pp,scale):
  for g,x in pp.items():
   v=dst.get(g,0)+scale*x
   if v:dst[g]=v
   else:dst.pop(g,None)
 def combine(name,terms):
  nonlocal current,combine_context
  pool[name]={};records=[];combine_context=dict(output=name,full_terms=terms,completed_terms=0)
  for a,b,sgn,label in terms:
   tmp=mul(a,b,label);add_into(pool[name],tmp,sgn);records.append(dict(left=a,right=b,sign=sgn,source_label=label));combine_context['completed_terms']=len(records);current=dict(kind='completed_combine_term',parent_combine=combine_context,no_pending_product=True);guard(len(tmp)+2);del tmp
  ops.append(dict(output=name,op='sum_products',terms=records));current=None;combine_context=None;return name
 def diff_literal(name,b,i,j,k,scale):
  pool[name]=dict(pool[f'{b}{i}_{j}']);add_into(pool[name],pool[f'{b}{k}_{j}'],-scale);ops.append(dict(output=name,op='literal_linear_combination',left=f'{b}{i}_{j}',right=f'{b}{k}_{j}',right_scalar=-scale));guard()
 try:
  combine('deltaP',[('A0_0','A1_1',1,'deltaP_plus'),('A0_1','A1_0',-1,'deltaP_minus')]);combine('deltaQ',[('A2_2','A3_3',1,'deltaQ_plus'),('A2_3','A3_2',-1,'deltaQ_minus')]);assert pool['deltaP'].get(0)==pool['deltaQ'].get(0)==L**2
  adjP=[[('A1_1',1),('A0_1',-1)],[('A1_0',-1),('A0_0',1)]];adjQ=[[('A3_3',1),('A2_3',-1)],[('A3_2',-1),('A2_2',1)]]
  for j in range(4):diff_literal('dc'+str(j),'C',2,j,0,kappa)
  for j in range(2):combine('v'+str(j),[('dc'+str(2+k),adjQ[k][j][0],adjQ[k][j][1],f'v{j}term{k}')for k in range(2)])
  for j in range(2):combine('u'+str(j),[('deltaQ','dc'+str(j),1,f'u{j}_left')]+[('v'+str(k),f'A{2+k}_{j}',-1,f'u{j}_cross{k}')for k in range(2)])
  for j in range(2):combine('w'+str(j),[('u'+str(k),adjP[k][j][0],adjP[k][j][1],f'w{j}term{k}')for k in range(2)])
  for j in range(2):combine('vp'+str(j),[('deltaP','v'+str(j),1,f'deltaP_v{j}')])
  combine('delta',[('deltaP','deltaQ',1,'deltaP_deltaQ')]);assert pool['delta'].get(0)==L**4
  dump(out/'all-source-and-prefix-identities.json',dict(uniform_integer_denominator=str(L),source=source_records,ALL_source23932_bound=True,actual_prefix_ratio59375_all7polys=True,alpha0_minusalpha2=68,
   all_prefix_operations=ops,all_prefix_integer_polynomials={k:sorted(v.items())for k,v in pool.items()if k not in literal_names},exact_A_topright2x2_ZERO=True,deltaP_anddeltaQ_constant=L**2,scalar_denominator_for_outputs=str(L**5)))
  for j in range(7):
   name='dd'+str(j);diff_literal(name,'D',2,j,0,kappa);outname='N'+str(j)
   combine(outname,[('delta',name,1,f'output{j}_deltaD')]+[('w'+str(k),f'B{k}_{j}',-1,f'output{j}_top{k}')for k in range(2)]+[('vp'+str(k),f'B{2+k}_{j}',-1,f'output{j}_bottom{k}')for k in range(2)])
   pp=pool[outname];ub=plan['outputs'][j];assert len(pp)<=ub['support_upper']and all(abs(v)<=int(ub['integer_abs_upper'])for v in pp.values());iszero=not pp
   dump(out/f'whole-row-column{j}.json',dict(column=j,integer_polynomial=sorted(pp.items()),common_denominator=str(L**5),is_ZERO=iszero,exact_operation=ops[-1],literal_difference=ops[-2],sourcebound=ub))
   results.append(dict(column=j,is_ZERO=iszero,integer_coefficient_terms=len(pp),sha256=sha(out/f'whole-row-column{j}.json')));del pool[name];del pool[outname]
   if not iszero:break
  allzero=len(results)==7 and all(e['is_ZERO']for e in results);scope=dict(candidate_whole_row_identicallyZERO_primary=allzero,all7_columns_verified=len(results)==7,first_nonzero_source_obstruction=None if allzero else results[-1]['column'],
   exact_mathematical_identity='deltaPdeltaQ*(S2-59375*S0)=N/L^5; A(0)=I so deltaPdeltaQ hasconstant1 and is nonzero',
   ifallzero='formalrowrelation=>fixed11detM(t^-w)ZERO; acceptedmixedradixNewtoninjectivity=>originalfixedminorZERO overQQ(a,e,d)',
   only_this_candidate=True,independent_peer_required_before_accepted_global_zero=True,NO_fullcubic57_or_otherminor_rank_inference=True)
  dump(out/'exact-source-row-conclusion.json',scope)
 except GuardStop as err:
  signal.alarm(0);state=dict(version='K92WR001',reason=str(err),current_multiplication=current,all_active_polynomials={k:sorted(p.items())for k,p in pool.items()},all_completed_operations=ops,completed_columns=results,source_pins=pins,uniform_integer_denominator=str(L),source572_replayed=False,restore='JSON1 exactintegers; inspectonly untilnewROOTbudget; current pair atcursorNOTcommitted')
  dump(out/'state-K92WR001.json',state);dump(out/'report.json',dict(id=ID,all_passed=False,hard_guard_preserved=True,reason=str(err),no_row_or_determinant_claim=True,state_sha256=sha(out/'state-K92WR001.json'),paircount=paircount,owned_peak=peak,seconds=time.monotonic()-start));return
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new<gate['new_RAW_upper']and new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,ONE_whole_row_identity_job=True,source23932_and_prefix_ratio_bound=True,whole_row_kappa=59375,checked_columns=results,all7_ZERO=allzero,candidate_det11_identicallyZERO_PRIMARY=allzero,
  first_nonzero_obstruction_column=None if allzero else results[-1]['column'],no_det_NONZERO_or_global57_claim=True,independent_peer_pending=True,whole_owned_scalar_peak=peak,actual_integer_pairs=paircount,
  no_new_weight_window_point_prime_or_source572=True,total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
