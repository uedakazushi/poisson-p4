#!/usr/bin/env python3
"""SOURCE support theorem at the first alternative-path boundary, no new values."""
import argparse,hashlib,itertools,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_BORDER5_CLOSED_GAP_20261105_v63';W=(1,2244,68)
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 base='reports/S61_NONSQUARE_PREPOWER_K92_POLYNOMIAL_BORDER5_SOURCE_20261105_v63';r=src(base+'/report.json');assert r['all_passed']and not r['source_support_possible'];old=src(base+'/polynomial-border5-source.json');pr=src('reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL5_RANK_ENTRY_20261105_v63/report.json');assert pr['exact_firstpositive5_rank']==0
 closed=src('reports/S61_NONSQUARE_PREPOWER_K92_RESIDUAL6_FINITE_SOURCE_20261105_v63/closed-inverse-source-proof.json');alt=min(p['cost']-r['direct_min']for r in closed['all_alternative_simplepaths_gt_directplus69']for p in r['alternative_simplepaths']);horizon=alt-1;assert alt==135 and horizon==134
 polys={};terms=0
 for rel in['reports/S61_NONSQUARE_PREPOWER_K92_EXACT_ABC72_20261105_v63/ABC72.jsonl','reports/S61_NONSQUARE_PREPOWER_K92_EXACT_D49_20261105_v63/D49.jsonl']:
  p=ROOT/rel;pins[rel]=sha(p)
  with p.open()as f:
   for line in f:
    e=json.loads(line);b=e.get('block','D');i=e['block_row']if b!='D'else e['residual_row'];j=e['block_column']if b!='D'else e['residual_column'];polys[b,i,j]={int(g)for g,v in e['sparse_integer_polynomial']};terms+=len(polys[b,i,j])
 assert terms==23932
 typ=src('reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63/typed-next-residual-Schur.json');alpha=typ['alpha'];beta=typ['beta'];Ast=closed['Ainverse_minplus'];DA={0,1,69};HB={};checks=[]
 for i in range(7):
  for j in range(7):
   sh=alpha[i]+beta[j];S={d+z-sh for d in polys['D',i,j]for z in DA if 0<=d+z-sh<=horizon}
   for k in range(4):
    for l in range(4):
     if Ast[k][l]>=10**12:continue
     C=polys['C',i,k];B=polys['B',l,j]
     if not C or not B or min(C)+Ast[k][l]+min(B)-sh>horizon:continue
     A={0}if k==l==1 else DA if k==l else polys['A',k,l]if(k==1 or l==1)else{a+z for a in polys['A',k,l]for z in DA}
     for c in C:
      if c+min(B)+min(A,default=10**12)-sh>horizon:continue
      for b in B:
       if c+b+min(A,default=10**12)-sh>horizon:continue
       S.update(c+a+b-sh for a in A if 0<=c+a+b-sh<=horizon)
   HB[i,j]=S
 RR=old['R2_normalized'];QQ=old['Q2'];TB={}
 for u in range(7):
  for v in range(7):
   S=set().union(*(HB[i,j]for i in range(7)for j in range(7)if RR[u][i]!='0'and QQ[j][v]!='0'));S.discard(0)
   if u==v==0:S.add(0)
   prior=next(e for e in old['Tbar_positive_absolute_bounds']if(e['row'],e['column'])==(u,v));oldS={g for g,n in prior['coefficient_upper']};assert {g for g in S if g<=114}==oldS
   TB[u,v]=S;checks.append(dict(row=u,column=v,all_support_through134=sorted(S),new_support_after114=sorted(g for g in S if g>114)))
 def conv(p,q):return{a+b for a in p for b in q if a+b<=horizon}
 cof0=conv(TB[0,0],TB[6,6])|conv(TB[0,6],TB[6,0]);rows=[];maxneeded=0;supports={}
 for u,p in enumerate(old['freeP']):
  for v,q in enumerate(old['freeQ']):
   i=1+p;j=1+q;shift=next(e['removed_t_power']for e in old['entries']if(e['row'],e['column'])==(u,v));N=conv(cof0,TB[i,j])|conv(conv(TB[0,6],TB[6,j])|conv(TB[6,6],TB[0,j]),TB[i,0])|conv(conv(TB[6,0],TB[0,j])|conv(TB[0,0],TB[6,j]),TB[i,6]);upto=horizon-shift;S={g-shift for g in N if 2<=g-shift<=upto};supports[u,v]=S;maxneeded=max(maxneeded,max(S,default=0));rows.append(dict(row=u,column=v,shift=shift,provable_support_through_degree=upto,known_minimum=2,all_possible_low_support=sorted(S),uncontrolled_remainder_order=upto+1))
 assert maxneeded==6 and all(not e['new_support_after114']for e in checks)
 lowdet=set();assignments=[]
 for pp in itertools.permutations(range(5)):
  S={0}
  for i,j in enumerate(pp):S={a+b for a in S for b in supports[i,j]}
  lowdet|=S
  if S:assignments.append(dict(columns=pp,possible_degrees=sorted(S)))
 maxlow=max(lowdet);assert maxlow<52
 # One uncontrolled entry plus four entries of order>=2 is the first unknown contribution.
 remainders=[e['uncontrolled_remainder_order']for e in rows];unknown_order=min(remainders)+8;assert unknown_order==73
 oldwhole=26705;outer=26653;newformal=outer+unknown_order;assert newformal==26726
 wb=src('reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/weighted-output.json');box=wb['output_metadata']['box'];weights=wb['weights'];linear=wb['output_metadata']['linear'];limit=232227-newformal;best=None;bestmon=None
 for a in range(box[0][0],box[0][1]+1):
  for e in range(box[1][0],box[1][1]+1):
   lo,hi=box[2]
   for w,(lower,upper)in zip(weights,linear):
    b=w[0]*a+w[1]*e
    if w[2]:lo=max(lo,(lower-b+w[2]-1)//w[2]);hi=min(hi,(upper-b)//w[2])
    elif not lower<=b<=upper:hi=lo-1;break
   hi=min(hi,(limit-a-2244*e)//68)
   if lo<=hi:
    v=a+2244*e+68*hi
    if best is None or v>best:best=v;bestmon=(a,e,hi)
 newton=232227-best;assert newton>=newformal
 dump(out/'all-source-supports-through-first-path-break.json',dict(first_alternative_path_slack=alt,strict_closed_horizon=horizon,Tbar=checks,border5=rows,
  low_determinant_support=sorted(lowdet),all120_permutation_supports=assignments,no_new_signed_target_coefficients=True))
 dump(out/'zero-and-next-gap-proof.json',dict(known_whole_order_before=oldwhole,known_low_border_det_order=52,maximum_low_polynomial_det_degree=maxlow,
  exact_low_polynomial_determinant_IDENTICALLY_ZERO='detlow degree<52, but detB5 order>=52 and all remainders start73; thus everydetlow coefficientzero',
  remainder_of_each_border_entry_atleast=min(remainders),other_four_entries_order_atleast8=True,whole_order_after_support=newformal,
  newton_sharpened_order=newton,next_possible_original_monomial=bestmon,next_possible_original_weight=best,
  coefficient_a26e91d19_ZERO_by_source_support=True,no_full_detM_zero_or_nonzero_claim=True,no_new_coefficient_window_weight_point_prime_or_source572=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins)
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_new_signed_target_values=True,closed_horizon=134,all_border_low_support_max6=True,low_det_polynomial_max_degree=maxlow,
  allowed_a26e91d19_coefficient_ZERO=True,new_whole_order_lower_bound=newton,next_possible_monomial=bestmon,no_full_det11_zero_nonzero_or_global57_claim=True,
  total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 assert new+external<33554432;dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
