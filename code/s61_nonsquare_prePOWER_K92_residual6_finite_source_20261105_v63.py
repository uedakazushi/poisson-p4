#!/usr/bin/env python3
"""SOURCE ONLY: finite residual6 leading DAG, closed scalar powers, no coefficients."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_RESIDUAL6_FINITE_SOURCE_20261105_v63';W=(1,2244,68);INF=10**12
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 abase='reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63';ar=src(abase+'/report.json');assert ar['all_passed']and pins[abase+'/report.json']=='abf17b55f023c7b87b535383e95e3173499b0ae28b9165f3fd5e9cde8ca806f7';nexts=src(abase+'/next-residual6-source.json');assert nexts['minimum_residual6_determinant_order']==74
 dbase='reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63';src(dbase+'/report.json');profile='reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63';pr=src(profile+'/report.json');assert pr['exact_H_leading_rank']==1;typed=src(profile+'/typed-next-residual-Schur.json');R2=[[Q(v)for v in rr]for rr in typed['row_transform_R2']];Q2=[[Q(v)for v in rr]for rr in typed['column_transform_Q2']];alpha=typed['alpha'];beta=typed['beta']
 polys={};prov={};terms_count=0
 for rel in[abase+'/ABC72.jsonl',dbase+'/D49.jsonl']:
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for ln,line in enumerate(f):
    e=json.loads(line);block=e.get('block','D');i=e['block_row']if block!='D'else e['residual_row'];j=e['block_column']if block!='D'else e['residual_column'];c=Q(e['QQ_content']);vals={int(g):c*int(v)for g,v in e['sparse_integer_polynomial']};assert all(vals.values());polys[block,i,j]=vals;prov[block,i,j]=dict(path=rel,line=ln,entry_index=e['entry_index'],QQ_content=e['QQ_content']);terms_count+=len(vals)
 assert len(polys)==121 and terms_count==23932
 A=[[polys['A',i,j]for j in range(4)]for i in range(4)];assert all(not A[i][j]for i in range(2)for j in range(2,4));assert all(A[i][j].get(0,Q(0))==Q(i==j)for i in range(4)for j in range(4))
 amin=[[min((g for g in A[i][j]if g>0),default=INF)for j in range(4)]for i in range(4)];star=[[0 if i==j else amin[i][j]for j in range(4)]for i in range(4)]
 for k in range(4):
  for i in range(4):
   for j in range(4):star[i][j]=min(star[i][j],star[i][k]+star[k][j])
 cycleproof=[];pathproof=[]
 for k in range(4):
  for length in range(1,5):
   for mids in itertools.permutations([i for i in range(4)if i!=k],length-1):
    walk=(k,)+mids+(k,);cost=sum(amin[a][b]for a,b in zip(walk,walk[1:]))
    if cost>=INF:continue
    if length==1 and k==1:continue
    assert cost>69;cycleproof.append(dict(walk=walk,min_cost=cost))
 for k in range(4):
  for l in range(4):
   if k==l:continue
   alternatives=[]
   for count in[1,2]:
    for mids in itertools.permutations([i for i in range(4)if i not in[k,l]],count):
     walk=(k,)+mids+(l,);cost=sum(amin[a][b]for a,b in zip(walk,walk[1:]));
     if cost<INF:alternatives.append(dict(walk=walk,cost=cost))
   if amin[k][l]>=INF:assert star[k][l]>=INF
   else:assert star[k][l]==amin[k][l]and all(p['cost']>amin[k][l]+69 for p in alternatives)
   pathproof.append(dict(start=k,end=l,direct_min=amin[k][l],shortest=star[k][l],alternative_simplepaths=alternatives))
 looplow={g:v for g,v in A[1][1].items()if 0<g<=69};assert looplow=={1:Q(77,6625),69:Q(-97,175)};aloop=Q(77,6625);bloop=Q(-97,175);DG=aloop.denominator**69*bloop.denominator;assert abs(aloop)<1 and abs(bloop)<1
 dump(out/'closed-inverse-source-proof.json',dict(fixed_weight=W,extra_window_derived_from_residual6_targets_max69=True,exact_A_topright2x2_ZERO=True,A_min_positive=amin,Ainverse_minplus=star,
  all_nonloop1_cycles_cost_gt69=cycleproof,all_alternative_simplepaths_gt_directplus69=pathproof,A11_through69={'constant':'1','t1':str(aloop),'t69':str(bloop)},
  G_recipe='G(k)=(-77/6625)^k for0<=k<69; G(69)=(-77/6625)^69+97/175',G_common_denominator=str(DG),G_absolute_bound='2',
  inverse_near_shortest='diag!=1:delta_a0; diag1:G(a); offdiagonal touching1:-sum Akl[d]*G(a-d); offdiagnotouch1:-Akl[a]; top-rightzero',
  no_Gpowers_or_inversecoefficients_evaluated=True))
 literals=[];litkeys={};Gused=set();inv_cache={};Hrecipes=[];Hkeys={};Habs=[]
 def literal(block,i,j,g):
  key=(block,i,j,g)
  if key in litkeys:return litkeys[key]
  value=polys[block,i,j].get(g,Q(0))
  if not value:return None
  q=len(literals);literals.append(dict(id=q,kind='saved_coefficient',block=block,row=i,column=j,gap=g,value=str(value),source=prov[block,i,j]));litkeys[key]=q;return q
 def scalar(value,reference):
  key=('scalar',str(value),json.dumps(reference,sort_keys=True))
  if key in litkeys:return litkeys[key]
  assert value;q=len(literals);literals.append(dict(id=q,kind='constant_QQ',value=str(value),source=reference));litkeys[key]=q;return q
 def invterms(k,l,a):
  key=(k,l,a)
  if key in inv_cache:return inv_cache[key]
  result=[]
  if a<0 or star[k][l]>=INF or a<star[k][l]:pass
  else:
   assert a<=star[k][l]+69
   if k==l:
    if k==1:assert a<=69;result=[dict(sign=1,literals=[],G=a)]
    elif a==0:result=[dict(sign=1,literals=[],G=0)]
   elif k==1 or l==1:
    for d in sorted(A[k][l]):
     if 0<d<=a:
      assert a-d<=69;result.append(dict(sign=-1,literals=[literal('A',k,l,d)],G=a-d))
   else:
    q=literal('A',k,l,a)
    if q is not None:result=[dict(sign=-1,literals=[q],G=0)]
  for p in result:Gused.add(p['G'])
  inv_cache[key]=result;return result
 def Hrecipe(i,j,g):
  key=(i,j,g)
  if key in Hkeys:return Hkeys[key]
  assert 1<=g<=69;theta=alpha[i]+beta[j]+g;terms=[];d=literal('D',i,j,theta)
  if d is not None:terms.append(dict(sign=1,literals=[d],G=0))
  for k in range(4):
   cp=polys['C',i,k];cm=min(cp,default=INF)
   if cm>=INF:continue
   for l in range(4):
    bp=polys['B',l,j];bm=min(bp,default=INF)
    if bm>=INF or star[k][l]>=INF:continue
    baseline=cm+star[k][l]+bm;assert baseline>=alpha[i]+beta[j];slack=theta-baseline
    if slack<0:continue
    assert slack<=g<=69
    citems=[c for c in cp if c<=cm+slack];bitems=[b for b in bp if b<=bm+slack]
    for c in sorted(citems):
     for b in sorted(bitems):
      aa=theta-c-b
      for inv in invterms(k,l,aa):terms.append(dict(sign=-inv['sign'],literals=[literal('C',i,k,c)]+inv['literals']+[literal('B',l,j,b)],G=inv['G'],path_source=dict(C=[i,k,c],Ainverse=[k,l,aa],B=[l,j,b])))
  bound=Q(0)
  for t in terms:
   product=Q(2)
   for q in t['literals']:product*=abs(Q(literals[q]['value']))
   bound+=product
  q=len(Hrecipes);Hkeys[key]=q;Hrecipes.append(dict(id=q,Hrow=i,Hcolumn=j,relative_gap=g,original_S_degree=theta,terms=terms,absolute_bound=str(bound)));Habs.append(bound);return q
 E=nexts['residual6_min'];ar=nexts['row_dual_potentials'];bc=nexts['column_dual_potentials'];Tplus=nexts['transformed_positive_part_min'];entries=[];Eabs=[];expanded=0
 for i in range(6):
  for j in range(6):
   g=ar[i]+bc[j]
   if E[i][j]!=g:continue
   assert 1<=g<=69 and Tplus[1+i][1+j]==g<Tplus[1+i][0]+Tplus[0][1+j];refs=[];bound=Q(0)
   for p in range(7):
    if not R2[1+i][p]:continue
    for q in range(7):
     if not Q2[q][1+j]:continue
     w=R2[1+i][p]*Q2[q][1+j];wid=scalar(w,dict(kind='R2Q2',R2_index=[1+i,p],Q2_index=[q,1+j],R2_value=str(R2[1+i][p]),Q2_value=str(Q2[q][1+j])));hid=Hrecipe(p,q,g);refs.append(dict(outer_literal=wid,H_recipe=hid));bound+=abs(w)*Habs[hid];expanded+=len(Hrecipes[hid]['terms'])
   entries.append(dict(row=i,column=j,relative_gap=g,known_new1block_correction_STRICTLY_later=True,terms=refs,absolute_bound=str(bound)));Eabs.append(bound)
 DL=math.lcm(*(Q(q['value']).denominator for q in literals));DEN=DL**5*DG;Nbound=max((math.ceil(DEN*b)for b in Eabs),default=0);bits=Nbound.bit_length();detbits=6*bits+math.factorial(6).bit_length();nodes=len(literals)+len(Hrecipes)+sum(len(h['terms'])for h in Hrecipes)+len(entries)+sum(len(e['terms'])for e in entries)+len(inv_cache)+len(Gused);live=terms_count+nodes+6*6*4+15000;raw_upper=(len(literals)+sum(len(h['terms'])for h in Hrecipes)+len(entries)+70)*(detbits//3+600)+500000
 assert all(len(t['literals'])+1<=5 for h in Hrecipes for t in h['terms']);assert live<150000
 dump(out/'finite-residual6-DAG.json',dict(fixed_weight=W,derived_residual6_order=74,total_original_order=26653,unique_normalized11_monomial=[10,91,20],
  all_tight_new1block_corrections_STRICTLY_later=True,scalar_inverse_G_formula=dict(a=str(aloop),b=str(bloop),common_denominator=str(DG),indices_used=sorted(Gused)),
  literals=literals,H_recipes=Hrecipes,residual6_tight_entries=entries,all_other36_entries_ZERO_by_integer_dual=True,ordinary6det_sign='unchanged residualrow andcolumn order',
  common_literal_denominator=str(DL),coefficient_integer_basis_denominator=str(DEN),integer_lift='literal n_i=DL*value_i; Gint=DG*G(k). Each Eterm withm<=5 literals contributes sign*product(n_i)*DL^(5-m)*Gint.',
  source_algebra='H=D-C*A^-1*B; inverseoffdiagonal andG signs retained in each Hterm; E firstlayer isR2*H_g*Q2 because newUcorrectionstrictlylater',
  NO_inverse_Gpowers_H_or_E_coefficients_evaluated=True))
 dump(out/'resource-gate.json',dict(input_exact_formal_terms=terms_count,literal_constants=len(literals),H_coefficient_requests=len(Hrecipes),H_product_terms=sum(len(h['terms'])for h in Hrecipes),outer_H_uses=sum(len(e['terms'])for e in entries),expanded_product_terms=expanded,
  inverse_request_keys=len(inv_cache),G_closedpower_indices=len(Gused),largest_G_power=max(Gused,default=0),tight_residual6_entries=len(entries),DAG_nodes_upper=nodes,
  whole_live_upper=live,coefficient_integer_bit_upper=bits,determinant_integer_bit_upper=detbits,integer_basis_denominator_bits=DEN.bit_length(),new_RAW_upper=raw_upper,
  proposed_job='ONE finite closedGpowers+integer-DAG coefficients+6det/rank; scalarunitidentitymodt70 maybechecked; NO generalinverse series/source572/otherweight orblindK',
  hard_seconds=60,hard_address_bytes=2000000000,hard_scalar_terms=150000,hard_total_RAW=33554432))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);total=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert raw_upper+external<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,same_fixed_weight=W,derived_maximum_local_window=69,no_blind_K_or_weight_change=True,
  all_inverse_near_shortest_cases_proved_by_paths_andcycles=True,all_new1block_corrections_strictly_later_on_tight_edges=True,literal_constants=len(literals),H_requests=len(Hrecipes),H_product_terms=sum(len(h['terms'])for h in Hrecipes),tight_E6_entries=len(entries),G_indices=len(Gused),max_Gpower=max(Gused,default=0),
  whole_live_upper=live,coefficient_integer_bit_upper=bits,determinant_integer_bit_upper=detbits,new_RAW_upper=raw_upper,no_target_coefficients_or_Gpowers_or_det_computed=True,no_det11_zero_nonzero_or_global57_claim=True,
  total_RAW_before_report=total,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
