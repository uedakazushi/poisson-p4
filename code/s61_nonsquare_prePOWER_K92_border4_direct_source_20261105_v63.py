#!/usr/bin/env python3
"""SOURCE ONLY one allowed coefficient via polynomial 3x3 bordered Schur."""
import argparse,hashlib,importlib.util,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_BORDER4_DIRECT_SOURCE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 profile='reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_RANK_ENTRY_20261105_v63';pr=src(profile+'/report.json');assert pr['all_passed']and pr['exact_firstpositive5_rank']==0;skip=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_CLOSED_GAP_V2_20261105_v63/report.json');assert skip['new_whole_order_lower_bound']==26726;intersection=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_SUPPORT_INTERSECTION_20261105_v63/report.json');assert intersection['first_possible_leading_Newton_match']==dict(border_det_degree=141,whole_order=26794,original_monomial=[5,91,18])
 src(profile+'/rank-kernel-unitminor.json');known=src('reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_SUPPORT_ZERO_20261105_v63/report.json');assert known['Newton_sharpened_order_lower_bound']==26705
 t2=src('reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63/typed-next-residual-Schur.json');oldR=[[Q(x)for x in r]for r in t2['row_transform_R2']];QQ=[[Q(x)for x in r]for r in t2['column_transform_Q2']];alpha=t2['alpha'];beta=t2['beta'];K2=oldR[1][0];assert K2
 RR=[oldR[0]]+[[x/K2 for x in r]for r in oldR[1:]];assert all(sum(bool(x)for x in r)==1 for r in RR)
 r3=src('reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_SOURCE_STRUCTURE_20261105_v63/typed-residual5-and-next-source.json');aa=r3['aa'];bb=r3['bb'];freeP=[1,2,3,4,0];freeQ=list(range(5));shifts=[[aa[p]+bb[q]+1 for q in freeQ]for p in freeP]
 target=345;known_min=1;entry_max=342;maxgap=550
 closed=src('reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL6_FINITE_SOURCE_20261105_v63/closed-inverse-source-proof.json');alt=min(p['cost']-r['direct_min']for r in closed['all_alternative_simplepaths_gt_directplus69']for p in r['alternative_simplepaths']);cycle=min(r['min_cost']for r in closed['all_nonloop1_cycles_cost_gt69']);assert alt==135 and maxgap<cycle
 polys={};prov={};input_terms=0
 for rel in['reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/ABC72.jsonl','reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/D49.jsonl']:
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for ln,line in enumerate(f):
    e=json.loads(line);b=e.get('block','D');i=e['block_row']if b!='D'else e['residual_row'];j=e['block_column']if b!='D'else e['residual_column'];c=Q(e['QQ_content']);polys[b,i,j]={int(g):c*int(v)for g,v in e['sparse_integer_polynomial']};prov[b,i,j]=dict(path=rel,line=ln,entry_index=e['entry_index'],QQ_content=e['QQ_content']);input_terms+=len(polys[b,i,j])
 assert input_terms==23932
 L=math.lcm(*(x.denominator for p in polys.values()for x in p.values()));I={k:{g:int(L*x)for g,x in p.items()}for k,p in polys.items()};assert all(Q(v,L)==polys[k][g]for k,p in I.items()for g,v in p.items())
 ap=polys['A',1,1];assert {g:x for g,x in ap.items()if 0<g<=maxgap}=={1:Q(77,6625),69:Q(-97,175)}
 assert all(min((g for g in polys['A',i,i]if g>0),default=10**12)>maxgap for i in[0,2,3])
 DA={0:Q(1),1:Q(77,6625),69:Q(-97,175)};Dm=math.lcm(*(x.denominator for x in DA.values()));DAint={g:int(Dm*x)for g,x in DA.items()};HD=L**5*Dm;Rd=math.lcm(*(x.denominator for r in RR for x in r));Qd=math.lcm(*(x.denominator for r in QQ for x in r));Ri=[[int(Rd*x)for x in r]for r in RR];Qi=[[int(Qd*x)for x in r]for r in QQ];TD=HD*Rd*Qd
 paths={};pathproof=[]
 for k in range(4):
  for l in range(4):
   found=[]
   if k==l:found=[dict(vertices=[k],edges=[],minimum=0,sign=1,through1=k==1)]
   else:
    for n in range(3):
     for mid in itertools.permutations([v for v in range(4)if v not in[k,l]],n):
      vs=[k]+list(mid)+[l];edges=list(zip(vs,vs[1:]));cost=sum(min(I['A',a,b],default=10**12)for a,b in edges)
      if cost<=closed['Ainverse_minplus'][k][l]+maxgap:found.append(dict(vertices=vs,edges=edges,minimum=cost,sign=(-1)**len(edges),through1=1 in vs))
   paths[k,l]=found;pathproof.append(dict(start=k,end=l,paths=found,nearshortest_extrawindow=maxgap))
 counts={'Hbar_contributions':0,'cofactor_bound_products':0,'border_bound_products':0,'det5_bound_products':0};Ast=closed['Ainverse_minplus'];HB={};HBrecipes=[]
 # Absolute upper bounds only. No signed Hbar/Tbar/border coefficient is summed.
 for i in range(7):
  for j in range(7):
   shift=alpha[i]+beta[j];bound={};recipes=[]
   def put(g,v,r):
    if 0<=g<=maxgap and v:bound[g]=bound.get(g,0)+v;recipes.append(r);counts['Hbar_contributions']+=1
   for d,dv in I['D',i,j].items():
    for z,zv in DAint.items():put(d+z-shift,abs(dv*zv)*L**4,dict(kind='DA_D',Dgap=d,DAgap=z))
   for k in range(4):
    for l in range(4):
     if Ast[k][l]>=10**12:continue
     C=I['C',i,k];B=I['B',l,j]
     if not C or not B or min(C)+Ast[k][l]+min(B)-shift>maxgap:continue
     Abar={};source_kind='finite_simple_paths'
     for path in paths[k,l]:
      pp={0:1}
      for u,v in path['edges']:
       qq={}
       for g,x in pp.items():
        for h,y in I['A',u,v].items():
         if g+h<=Ast[k][l]+maxgap:qq[g+h]=qq.get(g+h,0)+abs(x*y)
       pp=qq
      n=len(path['edges']);scale=L**(3-n)
      for g,x in pp.items():
       for h,y in ({0:Dm}if path['through1']else DAint).items():
        if g+h<=Ast[k][l]+maxgap:Abar[g+h]=Abar.get(g+h,0)+abs(x*y)*scale
     # Abar is a bound-source polynomial in known literal coefficients, not a new target.
     for c,cv in C.items():
      if c+min(B)+min(Abar,default=10**12)-shift>maxgap:continue
      for b,bv in B.items():
       if c+b+min(Abar,default=10**12)-shift>maxgap:continue
       for a,av in Abar.items():put(c+a+b-shift,abs(cv*av*bv),dict(kind=source_kind,k=k,l=l,Cgap=c,Agap=a,Bgap=b))
   HB[i,j]=bound;HBrecipes.append(dict(row=i,column=j,recipes=recipes,coefficient_upper=[[g,str(v)]for g,v in sorted(bound.items())]))
 TB={}
 for u in range(7):
  for v in range(7):
   p={}
   for i in range(7):
    for j in range(7):
     if Ri[u][i]and Qi[j][v]:
      for g,h in HB[i,j].items():p[g]=p.get(g,0)+abs(Ri[u][i]*Qi[j][v])*h
   p.pop(0,None)
   if u==v==0:p[0]=TD
   TB[u,v]=p
 actualbase='reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_FOUR_MATRIX_ENTRY_20261105_v63'
 actual=src(actualbase+'/literal-source-and-all-finite-path-identities.json');actualrank=src('reports/S61_NONSQUARE_PREPOWER_K92_BORDER5_ACTUAL_SUPPORT_RANK_20261105_v63/actual-prefix-support-and-rankone-proof.json');assert actualrank['rank68']==1
 p=Q(actualrank['exact_nonzero_entry']['QQ']);assert p==Q(50416319635200000000000000,49)
 assert TD==int(actual['Tbar_common_denominator'])
 for e in actual['all_Tbar']:
  k=e['row'],e['column'];TB[k]={g:v for g,v in TB[k].items()if g>138}|{int(g):abs(int(v))for g,v in e['integer_polynomial']}
 def add(*pp):
  r={}
  for p in pp:
   for g,v in p.items():r[g]=r.get(g,0)+v
  return r
 def conv(p,q,kind,cap=maxgap):
  r={}
  for a,x in p.items():
   for b,y in q.items():
    if a+b<=cap:r[a+b]=r.get(a+b,0)+x*y;counts[kind]+=1
  return r
 P=[0,6,1];free=list(range(2,6));cof={};detp={}
 for i in range(3):
  for j in range(3):
   rr=[P[k]for k in range(3)if k!=i];cc=[P[k]for k in range(3)if k!=j]
   cof[i,j]=add(conv(TB[rr[0],cc[0]],TB[rr[1],cc[1]],'cofactor_bound_products'),conv(TB[rr[0],cc[1]],TB[rr[1],cc[0]],'cofactor_bound_products'))
 for j in range(3):detp=add(detp,conv(TB[P[0],P[j]],cof[0,j],'cofactor_bound_products'))
 NB={};records=[]
 for u,i in enumerate(free):
  for v,j in enumerate(free):
   terms=[conv(detp,TB[i,j],'border_bound_products')]
   for k in range(3):
    for l in range(3):terms.append(conv(conv(TB[i,P[k]],cof[l,k],'border_bound_products'),TB[P[l],j],'border_bound_products'))
   full=add(*terms);shift=139+aa[i-1]+bb[j-1];tr={g-shift:n for g,n in full.items()if 1<=g-shift<=entry_max};NB[u,v]=tr
   records.append(dict(row=u,column=v,Tbar_row=i,Tbar_column=j,removed_t_power=shift,known_constant_ZERO_from_rankone68=True,needed_coefficient_upper=[[g,str(n)]for g,n in sorted(tr.items())]))
 dp={0:{0:1}};stage_slots=[]
 for i in range(4):
  nxt={}
  for mask,pb in dp.items():
   for j in range(4):
    if mask>>j&1:continue
    k=mask|1<<j;v=conv(pb,NB[i,j],'det5_bound_products',target-(3-i));nxt[k]=add(nxt.get(k,{}),v)
  dp=nxt;stage_slots.append(sum(map(len,dp.values())))
 possible=target in dp.get(15,{});targetbound=dp.get(15,{}).get(target,0);allbits=max([TD.bit_length(),targetbound.bit_length()]+[n.bit_length()for pb in NB.values()for n in pb.values()]);Tslots=sum(map(len,TB.values()));Nslots=sum(map(len,NB.values()));cofslots=len(detp)+sum(map(len,cof.values()));DPslots=max(stage_slots,default=0)
 live=2*input_terms+Tslots+Nslots+cofslots+2*DPslots+max(map(len,HB.values()),default=0)+15000;raw_upper=(Tslots+Nslots+2*DPslots+1000)*((allbits+2)//3+100)+3000000
 hp=ROOT/'code/s61_nonsquare_prePOWER_K92_leading_structure_20261105_v63.py';pins[str(hp.relative_to(ROOT))]=sha(hp);sp=importlib.util.spec_from_file_location('lin',hp);lin=importlib.util.module_from_spec(sp);sp.loader.exec_module(lin)
 r1=src('reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63/retropicalized-source.json');R1=[[Q(x)for x in r]for r in r1['row_transform']];Q1=[[Q(x)for x in r]for r in r1['column_transform']]
 unitfactor=1/(lin.determinant(R1)*lin.determinant(Q1)*lin.determinant(RR)*lin.determinant(QQ)*TD**16*p**3*math.prod(map(int,r1['constant_clearing_rowLCMs'])))
 dump(out/'direct-border4-source.json',dict(fixed_weight=[1,2244,68],whole_coefficient_order=27338,original_monomial=[5,91,10],required_border4_coefficient_degree=target,entry_maximum_needed=entry_max,maximum_Hbar_gap=maxgap,
  DA={g:str(v)for g,v in DA.items()},R2_normalized=RR,Q2=QQ,source_uniform_denominator=str(L),DA_denominator=str(Dm),Hbar_denominator=str(HD),Tbar_denominator=str(TD),R2_denominator=str(Rd),Q2_denominator=str(Qd),
  finite_simple_path_proof=pathproof,Hbar_bound_recipes=HBrecipes,Tbar_actual_prefix138_reused=True,Tbar_absolute_bounds=[dict(row=i,column=j,coefficient_upper=[[g,str(v)]for g,v in sorted(pb.items())])for(i,j),pb in TB.items()],
  pivot_rows=P,pivot_columns=P,pivot_determinant_order=71,pivot_determinant_leading=str(p),free_rows=free,free_columns=free,entries=records,
  exact_identity='det4border=det3pivot*Schur4entry; normalizedborder4=t^-(139+aa_p+bb_q)*detTbar[[0,6,1,i],[0,6,1,j]]; differsfromtypedresidual4 byformalunit withconstant p/c',
  full_determinant_identity='det(all4x4border)=det(pivot3)^3*detTbar; sumrowcolumnshifts627; wholeouter26579; target345+627-3*71+26579=27338',
  original_detM_target_from_integer_det4_target_factor=str(unitfactor),integer_entry_denominator='TD^4',integer_det4_denominator='TD^16',
  determinant_target_has_source_support=possible,all_previous_corrections_and_all16borderterms_retained=True,no_scalar_inverse_or_series=True,no_signed_new_target_coefficients=True))
 dump(out/'resource-gate.json',dict(input_terms=input_terms,Tbar_slots_upper=Tslots,border16_slots_upper=Nslots,pivot_and9cofactor_slots_upper=cofslots,DP_stage_slots_upper=stage_slots,whole_owned_scalar_upper=live,
  source_copies_counted=2,operation_pairs_upper=counts,max_integer_or_determinant_coefficient_bits=allbits,target_integer_bound=str(targetbound),new_RAW_upper=raw_upper,
  one_proposed_job='ONE integer DIRECT4border coefficient[t345]det of16finitepolynomialentries; requiredsource supports only, no inverse/windowsearch',hard_seconds=60,hard_address_bytes=2000000000,hard_scalar_terms=150000,hard_total_RAW=33554432))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins)
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_new_signed_target_coefficients=True,direct4border_inverse_free=True,actual_required_target_degree=target,required_Hgap=maxgap,
  source_support_possible=possible,whole_owned_scalar_upper=live,new_RAW_upper=raw_upper,preflight_fits=live<150000 and raw_upper+external<33554432,
  no_det11_zero_nonzero_or_global57_claim=True,total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
