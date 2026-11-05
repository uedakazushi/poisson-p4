#!/usr/bin/env python3
"""SOURCE ONLY residual5 first-positive finite DAG, caches and both corrections."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_FINITE_SOURCE_20261105_v63';INF=10**12
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base5='reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_SOURCE_STRUCTURE_20261105_v63';r5=src(base5+'/report.json');assert r5['all_passed'];typ5=src(base5+'/typed-residual5-and-next-source.json');assert typ5['maximum_derived_previousE_degree']==70;R3=[[Q(v)for v in rr]for rr in typ5['R3']];Q3=[[Q(v)for v in rr]for rr in typ5['Q3']]
 ar='reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63';src(ar+'/report.json');nxt=src(ar+'/next-residual6-source.json');Emin=nxt['residual6_min'];Tp=nxt['transformed_positive_part_min'];aa=nxt['row_dual_potentials'];bb=nxt['column_dual_potentials']
 fp='reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63';typ=src(fp+'/typed-next-residual-Schur.json');R2=[[Q(v)for v in rr]for rr in typ['row_transform_R2']];Q2=[[Q(v)for v in rr]for rr in typ['column_transform_Q2']];alpha=typ['alpha'];beta=typ['beta']
 closed=src('reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL6_FINITE_SOURCE_20261105_v63/closed-inverse-source-proof.json');assert all(p['min_cost']>70 for p in closed['all_nonloop1_cycles_cost_gt69']);assert all(p['cost']>q['direct_min']+70 for q in closed['all_alternative_simplepaths_gt_directplus69']for p in q['alternative_simplepaths'])
 polys={};prov={};input_terms=0
 for rel in[ar+'/ABC72.jsonl','reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/D49.jsonl']:
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for ln,line in enumerate(f):
    e=json.loads(line);b=e.get('block','D');i=e['block_row']if b!='D'else e['residual_row'];j=e['block_column']if b!='D'else e['residual_column'];c=Q(e['QQ_content']);polys[b,i,j]={int(g):c*int(v)for g,v in e['sparse_integer_polynomial']};prov[b,i,j]=dict(path=rel,line=ln,entry_index=e['entry_index'],QQ_content=e['QQ_content']);input_terms+=len(polys[b,i,j])
 assert input_terms==23932 and {g:v for g,v in polys['A',1,1].items()if 0<g<=70}=={1:Q(77,6625),69:Q(-97,175)}
 Ast=closed['Ainverse_minplus'];Hcache={};Ecache={};prior='reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL6_FINITE_ENTRY_20261105_v63'
 hc=src(prior+'/all31-H-coefficients.json');ec=src(prior+'/all25-E6-coefficients.json')
 for h in hc['all31']:Hcache[h['Hrow'],h['Hcolumn'],h['relative_gap']]=Q(h['exact_QQ_coefficient'])
 for e in ec['all25']:Ecache[e['row'],e['column'],e['relative_gap']]=Q(e['exact_QQ_coefficient'])
 nodes=[];memo={};cache_receipts=[];source_queries=[];Gindices=set();Hmemo={};Tmemo={};Ememo={};invMemo={};maxHgap=0
 def node(op,args,D,N,**more):
  assert D>0 and N>=0;key=(op,tuple(args),json.dumps(more,sort_keys=True))
  if key in memo:return memo[key]
  k=len(nodes);nodes.append(dict(id=k,op=op,args=args,denominator_upper=str(D),integer_abs_upper=str(N),**more));memo[key]=k;return k
 def lit(value,source):
  v=Q(value);return node('literal',[],v.denominator,abs(v.numerator),value=str(v),source=source)
 zero=lit(0,dict(constant='0'));one=lit(1,dict(constant='1'))
 def add(terms):
  merged={}
  for sign,k in terms:
   if k!=zero and sign:merged[k]=merged.get(k,0)+sign
  items=sorted((k,s)for k,s in merged.items()if s)
  if not items:return zero
  if items==[(one,1)]:return one
  if len(items)==1 and items[0][1]==1:return items[0][0]
  D=math.lcm(*(int(nodes[k]['denominator_upper'])for k,s in items));N=sum(abs(s)*int(nodes[k]['integer_abs_upper'])*(D//int(nodes[k]['denominator_upper']))for k,s in items)
  return node('sum',[k for k,s in items],D,N,weights=[s for k,s in items])
 def mul(args):
  if zero in args:return zero
  args=sorted(k for k in args if k!=one)
  if not args:return one
  if len(args)==1:return args[0]
  return node('product',args,math.prod(int(nodes[k]['denominator_upper'])for k in args),math.prod(int(nodes[k]['integer_abs_upper'])for k in args))
 def neg(k):return add([(-1,k)])
 def power(k,n):
  assert 0<=n<=70
  if n==0:return one
  if n==1:return k
  return node('power',[k],int(nodes[k]['denominator_upper'])**n,int(nodes[k]['integer_abs_upper'])**n,exponent=n)
 Gbase=lit(Q(-77,6625),dict(saved_scalar_inverse_G1=prior+'/all95-source-and-lift.json'));Blo=lit(Q(-97,175),dict(source_A11_gap69=prov['A',1,1]))
 def G(k):
  assert 0<=k<=70;Gindices.add(k)
  if k<69:return power(Gbase,k)
  if k==69:return add([(1,power(Gbase,69)),(-1,Blo)])
  return add([(1,power(Gbase,70)),(1,mul([lit(2,dict(constant='2')),neg(Gbase),Blo]))])
 def coefficient(b,i,j,g):
  value=polys[b,i,j].get(g,Q(0));source_queries.append(dict(block=b,row=i,column=j,gap=g,present=bool(value),source=prov[b,i,j]))
  if not value:return zero
  return lit(value,dict(kind='saved_formal_coefficient',block=b,row=i,column=j,gap=g,**prov[b,i,j]))
 def invterms(k,l,a):
  key=(k,l,a)
  if key in invMemo:return invMemo[key]
  terms=[]
  if a<0 or Ast[k][l]>=INF or a<Ast[k][l]:pass
  else:
   assert a<=Ast[k][l]+70
   if k==l:
    if k==1:terms=[(1,G(a))]
    elif a==0:terms=[(1,one)]
   elif k==1 or l==1:
    for d in sorted(polys['A',k,l]):
     if 0<d<=a:terms.append((-1,mul([coefficient('A',k,l,d),G(a-d)])))
   else:terms=[(-1,coefficient('A',k,l,a))]
  result=add(terms);invMemo[key]=result;return result
 def H(i,j,g):
  nonlocal maxHgap
  key=(i,j,g)
  if key in Hmemo:return Hmemo[key]
  assert 1<=g<=70;maxHgap=max(maxHgap,g)
  if key in Hcache:
   v=Hcache[key];cache_receipts.append(dict(kind='completed_H',key=key,value=str(v),source=prior+'/all31-H-coefficients.json'));result=zero if not v else lit(v,dict(kind='completed_H',key=key,path=prior+'/all31-H-coefficients.json'));Hmemo[key]=result;return result
  theta=alpha[i]+beta[j]+g;terms=[(1,coefficient('D',i,j,theta))]
  for k in range(4):
   C=polys['C',i,k];cm=min(C,default=INF)
   if cm>=INF:continue
   for l in range(4):
    B=polys['B',l,j];bm=min(B,default=INF)
    if bm>=INF or Ast[k][l]>=INF:continue
    baseline=cm+Ast[k][l]+bm;assert baseline>=alpha[i]+beta[j];slack=theta-baseline
    if slack<0:continue
    assert slack<=g<=70
    for c in sorted(c for c in C if c<=cm+slack):
     for b in sorted(b for b in B if b<=bm+slack):terms.append((-1,mul([coefficient('C',i,k,c),invterms(k,l,theta-c-b),coefficient('B',l,j,b)])))
  result=add(terms);Hmemo[key]=result;return result
 def T(u,v,g):
  key=(u,v,g)
  if key in Tmemo:return Tmemo[key]
  if g==0:result=one if u==v==0 else zero
  else:
   assert 1<=g<=70;terms=[]
   for i in range(7):
    if not R2[u][i]:continue
    for j in range(7):
     if Q2[j][v]:terms.append((1,mul([lit(R2[u][i]*Q2[j][v],dict(kind='R2Q2',indices=[u,i,j,v])),H(i,j,g)])))
   result=add(terms)
  Tmemo[key]=result;return result
 def E(i,j,g):
  key=(i,j,g)
  if key in Ememo:return Ememo[key]
  if g<Emin[i][j]:result=zero
  elif key in Ecache:
   v=Ecache[key];cache_receipts.append(dict(kind='completed_E',key=key,value=str(v),source=prior+'/all25-E6-coefficients.json'));result=zero if not v else lit(v,dict(kind='completed_E',key=key,path=prior+'/all25-E6-coefficients.json'))
  else:
   terms=[(1,T(1+i,1+j,g))];wm=Tp[1+i][0];vm=Tp[0][1+j];slack=g-wm-vm
   assert slack<=1
   if slack>=0:
    for w in range(wm,wm+slack+1):
     for v in range(vm,vm+slack+1):
      u=g-w-v
      if u not in[0,1]:continue
      inv=one if u==0 else neg(T(0,0,1));terms.append((-1,mul([T(1+i,0,w),inv,T(0,1+j,v)])))
   result=add(terms)
  Ememo[key]=result;return result
 roots=[];maxE=0
 for i in range(5):
  for j in range(5):
   terms=[]
   for p in range(6):
    if not R3[1+i][p]:continue
    for q in range(6):
     if not Q3[q][1+j]:continue
     g=aa[p]+bb[q]+1;maxE=max(maxE,g);terms.append((1,mul([lit(R3[1+i][p]*Q3[q][1+j],dict(kind='R3Q3',indices=[1+i,p,q,1+j])),E(p,q,g)])))
   roots.append(dict(row=i,column=j,node=add(terms)))
 assert maxE==70 and maxHgap<=70
 RD=math.lcm(*(int(nodes[r['node']]['denominator_upper'])for r in roots));bounds=[int(nodes[r['node']]['integer_abs_upper'])*(RD//int(nodes[r['node']]['denominator_upper']))for r in roots];bits=max(bounds,default=0).bit_length();detbits=5*bits+math.factorial(5).bit_length();allbits=max(max(int(n['integer_abs_upper']).bit_length(),int(n['denominator_upper']).bit_length())for n in nodes)
 # Structural maximum matching of nonzero symbolic roots, not coefficient rank.
 dp={0}
 for i in range(5):dp={mask|1<<j for mask in dp for j in range(5)if not(mask>>j&1)and int(nodes[roots[5*i+j]['node']]['integer_abs_upper'])}
 structural_full=31 in dp;live=input_terms+len(nodes)+len(source_queries)+len(cache_receipts)+25*4+15000;raw_upper=(len(nodes)+len(roots)+100)*((max(allbits,detbits)+2)//3+600)+500000
 dump(out/'finite-residual5-DAG.json',dict(fixed_weight=[1,2244,68],source_target='firstpositive residual5 matrix Z1, NOT blind widerprecision determinant',known_total_order_beforeZ=26653,coefficient_matrix_degree=1,residual5_roots=roots,nodes=nodes,
  cached_completed_coefficients_used=cache_receipts,source_literal_coefficient_queries=source_queries,source_relations='H=D-C*A^-1*B; T=R2*H*Q2; E=T22-W*U^-1*V with Uinv0=1,Uinv1=-U1; L=t^-aa*E*t^-bb; Z1=(R3*L1*Q3)_residual',
  both_previous_Schur_corrections_retained=True,old_scalar_Uinverse_indices_ONLY=[0,1],original_G_indices=sorted(Gindices),new_A11_window70_has_NOextra_term=True,G70_formula='(-a)^70+2*a*b',
  maximum_derived_E_degree=maxE,maximum_derived_H_degree=maxHgap,root_integer_common_denominator=str(RD),no_node_target_value_evaluated=True))
 dump(out/'resource-gate.json',dict(input_exact_formal_terms=input_terms,DAG_nodes=len(nodes),source_queries=len(source_queries),completed_H_or_E_cache_uses=len(cache_receipts),H_request_keys=len(Hmemo),E_request_keys=len(Ememo),T_request_keys=len(Tmemo),inverse_request_keys=len(invMemo),
  original_G_indices=sorted(Gindices),maximum_original_G_power=max(Gindices,default=0),whole_live_upper=live,max_node_numerator_or_denominator_bits=allbits,root_integer_bit_upper=bits,determinant5_integer_bit_upper=detbits,
  new_RAW_upper=raw_upper,structural_nonzero_roots=sum(bool(int(nodes[r['node']]['integer_abs_upper']))for r in roots),structural_perfect_matching=structural_full,
  proposed_JOB='ONE finite DAG evaluation and5x5rank/unitminor profile; determinantnonzero possibleONLYif structuralmatching; no nextwindow/weight/source572',hard_seconds=60,hard_address_bytes=2000000000,hard_scalar_terms=150000,hard_total_RAW=33554432))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);total=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert live<150000 and raw_upper+external<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,both_previous_corrections_and_all_signs_retained=True,old_scalar_inverse_ONLY0and1=True,original_G_indices=sorted(Gindices),DAG_nodes=len(nodes),H_request_keys=len(Hmemo),E_request_keys=len(Ememo),completed_cache_uses=len(cache_receipts),
  symbolic_nonzero_roots=sum(bool(int(nodes[r['node']]['integer_abs_upper']))for r in roots),structural_perfect_matching=structural_full,whole_live_upper=live,max_node_bits=allbits,root_integer_bit_upper=bits,determinant_bit_upper=detbits,new_RAW_upper=raw_upper,
  no_new_H_E_Z_Gpower_coefficients_or_determinant_evaluated=True,no_det11_zero_nonzero_or_global57_claim=True,total_RAW_before_report=total,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
