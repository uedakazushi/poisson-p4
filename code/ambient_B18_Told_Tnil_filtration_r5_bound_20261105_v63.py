#!/usr/bin/env python3
"""Fixed16 homogeneous products and exact natural-filtration theorem.
No higher natural coefficient, jet, Groebner search or cubic value.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63';OUT=ROOT/'reports'/ID
ZERO=(0,)*8
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def unit(var,n=1):return tuple(n if i==var else 0 for i in range(8))
def plusmon(a,b):return tuple(x+y for x,y in zip(a,b))
def lead(p):return max(p,key=lambda k:(-k[0],k[1]))
def termdict(v):return {(t,tuple(e)):F(c)for t,e,c in v if F(c)}
def monmul(p,e,c=F(1)):return {(t,plusmon(m,e)):a*c for(t,m),a in p.items()if a*c}
def terms(p):return [[t,list(e),c]for(t,e),c in sorted(p.items())]
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 exact=load('AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63');assert exact['actual_T_length']==11 and exact['actual_Tnil_length']==16
 old=load('AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63','exact-Loewy11-and-source-basis-corollary.json');G=load('AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63','known-initial40-and-length-upper17.json')
 assert load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW')['exact_natural_injection']
 natural=load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63','actual-natural-Told-to-Tnil-source-map-and-origin.json');N=natural['actual_origin_data'];assert natural['old_minimal_generator_free56_rows']==[22,23,29,30,55]
 z=load('AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63','actual-three-first-columns-and-fixed-quotient-proof.json');assert z['actual_cokernel_HF']==[3,1,1]
 assert load('AMBIENT_B18_TOLD_TNIL_ZERO3_SOURCE_20261105_v63_REVIEW')['all_passed']
 assert load('AMBIENT_B18_HONEST_NIL_T_QUADRATIC_REMAINING_20261105_v63_REVIEW')['actual_HF_Tnil_through2']==[5,5,3]
 red=load('AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63_REVIEW');assert red['actual_m3_equals_rA_m2']
 cub=load('AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63','literal-known-cubic-input-and-exact24-row-reduction.json')
 save('phase.json',dict(phase='actual source and filtered bases pinned; fixed tiny coefficient gate only'))
 origins=[Q.sparse(v)for v in N['actual_origin_map_columns5']];assert origins[4][1]==1 and all(not origins[i]for i in[1,2,3])
 F5=[Q.pairs(v.items())for v in z['quotient5_by40_after_actual_rank35_first_relations']]
 first3=[Q.pairs(v.items())for v in z['actual_three_minimal_first_columns40']]
 # Actual oldgr1 basis: e3*r, e4*z5, e4*r.
 old1=[(3,7),(4,5),(4,7)]
 map1=[Q.mvec(F5,{8*t+v:a for t,a in origins[j].items()})for j,v in old1]
 assert not map1[0]and not map1[1]and map1[2]
 zero3=[Q.mvec(F5,v)for v in first3];columns=zero3+[map1[2]];rows=[0,2,3,4]
 minor=[{j:col[r]for j,col in enumerate(columns)if col.get(r,F(0))}for r in rows];det4=Q.det_fixed(minor);assert det4
 assert all(not c.get(1,F(0))for c in columns) # remainingQ firstdual is r*e42
 # Whole F1 has oldm1 plus the three zero-origin top lifts. Their first
 # images are independent modulo im(oldm1), so F2 is contained in oldm1.
 # The above rank1 map on the exact oldgr1 basis then identifies F2.
 gens=[termdict(v)for v in G['homogeneous_module_generators40']];lms=[lead(p)for p in gens]
 def divide(p):
  work=p.copy();nf={};co={};steps=0
  while work:
   t,e=lead(work);a=work[t,e];hit=None
   for j,(row,m)in enumerate(lms):
    if t==row and all(x>=y for x,y in zip(e,m)):hit=j;break
   if hit is None:Q.put(nf,(t,e),a);work.pop((t,e));continue
   j=hit;delta=tuple(x-y for x,y in zip(e,lms[j][1]));c=a/gens[j][lms[j]];Q.put(co,(j,delta),c);work=Q.plus(work,Q.scaled(monmul(gens[j],delta),-c));steps+=1;assert steps<1000
  rebuilt=Q.plus(nf,*(monmul(gens[j],m,c)for(j,m),c in co.items()));assert rebuilt==p
  return nf,co
 u={(3,unit(7)):F(1)};v={(4,unit(5)):F(1)};q1=(4,plusmon(unit(5),unit(7)));q2=(4,unit(7,2));products=[];nfcols=[]
 for label,p in[('u=e3*r',u),('v=e4*z5',v)]:
  for var in range(8):
   raw=monmul(p,unit(var));nf,co=divide(raw);assert set(nf)<={q1,q2};nfcols.append({i:nf[k]for i,k in enumerate([q1,q2])if nf.get(k,F(0))})
   products.append(dict(source=label,variable=var,input=terms(raw),normal_form=terms(nf),exact_initial_ideal_cofactors=[[j,list(m),a]for(j,m),a in co.items()]))
 assert nfcols[15]=={0:F(1)}
 extra=next((j for j,col in enumerate(nfcols)if col.get(1,F(0))),None)
 productrank=2 if extra is not None else 1
 if extra is None:
  assert all(col==Q.scaled(nfcols[15],col.get(0,F(0)))for col in nfcols)
 else:assert Q.det_fixed([{0:nfcols[15].get(i,F(0)),1:nfcols[extra].get(i,F(0))}for i in range(2)])
 r2cert=[]
 for label,p in[('u=e3*r',u),('v=e4*z5',v)]:
  raw=monmul(p,unit(7,2));nf,co=divide(raw);assert not nf
  r2cert.append(dict(source=label,input=terms(raw),zero_initial_NF=True,cofactors=[[j,list(m),a]for(j,m),a in co.items()]))
 # Degree3 initial-zero means the ACTUAL class lies in oldm4=0.
 assert G['homogeneous_module_generators40'][-2]==[[3,list(unit(7,2)),'1'],[4,list(plusmon(unit(5),unit(7))),'11767/3652']]
 assert G['homogeneous_module_generators40'][-1]==[[4,list(plusmon(unit(5),unit(7,2))),'1']]
 # Additional parent-authorized2-column test uses only ZERO3 head3 firsts.
 pi3=[Q.pairs(w)for w in cub['quadratic_projection3_by180_sparse']];m2=list(itertools.combinations_with_replacement(range(8),2));mi={m:i for i,m in enumerate(m2)}
 rfirst={36*t+mi[tuple(sorted((i,7)))]:a for j,a in first3[2].items()for t,i in[divmod(j,8)]}
 head3r=Q.mvec(pi3,rfirst);head4r2=Q.mvec(pi3,{36*t+mi[(7,7)]:a for t,a in origins[4].items()})
 pairdet=Q.det_fixed([{0:head3r.get(i,F(0)),1:head4r2.get(i,F(0))}for i in[1,2]])
 pairrank=2 if pairdet else 1
 if not pairdet:
  k=next(i for i,a in head4r2.items()if a);ratio=head3r.get(k,F(0))/head4r2[k];assert head3r==Q.scaled(head4r2,ratio)
 else:ratio=None
 full21=productrank==2 or pairrank==1
 proof=dict(old_top_source_rows=[22,23,29,30,55],actual_old_gr1_basis=['e3*rA','e4*zA5','e4*rA'],actual_old_gr1_to_new_gr1_columns5=map1,
  actual_zero3_top_first_columns5=zero3,fixed_four_column_minor_rows=rows,fixed_four_column_minor_determinant=det4,
  exact_F2='i^-1(m²Tnil)=m²Told + k*(e3*rA) + k*(e4*zA5)',
  F2_proof='Originrank2 identifies F1 modulo oldm1 with the three ZERO3 top classes. Their first images are independent modulo image(oldm1), by the4-column minor. Thus F2 is contained in oldm1. Its quotient by oldm2 is exactly the kernel of the oldgr1 map, spanned by the two literal displayed classes.',
  all16_products_with_exact_initial_ideal_cofactors=products,degree2_product_rank=productrank,mF2_length=1+productrank,
  mF2_contains_old_m3_and_is_contained_old_m2=True,
  all_r2_u_v_degree3_zero_NF_cofactors=r2cert,
  actual_r2_F2_zero_proof='r²oldm²=0 by oldm4=0. For u=e3*r multiply actual quadratic initial relation e3*r²+c*e4*z5*r by r; its omitted terms have order>=4 and vanish. The actual cubic e4*z5*r² also has remainder in oldm4=0. Thus r²u=r²v=0 in the ACTUAL old module.',
  exact_F3='i^-1(m³Tnil), identified with m³Tnil because m³Q=0; length3 and contained in F2',
  actual_new_m5_zero='m⁵Tnil=r²*m³Tnil=i(r²F3)=0 by r-reduction and F3 subset F2',
  remaining_actual_HF_options=[[5,5,3,3],[5,5,3,2,1]],complete917_rank_options=[0,1],HF553111_EXCLUDED=True,
  m2_k_rA_Jordan_options=[[2,2,2],[3,2,1]],rA_NOT_DVR_ell=True,
  extra_actual_gr2_pair_columns=[head3r,head4r2],extra_pair_minor_rows=[1,2],extra_pair_minor_determinant=pairdet,extra_pair_rank=pairrank,extra_pair_ratio_if_dependent=ratio,
  extra_pair_criterion='If the two columns are dependent, u-d*(e4*r²) maps to newm³, while its r-image has nonzero oldq1 coefficient -11767/3652. Injectivity forces newm4 nonzero; together with newm5=0 this forcesHF55321. If independent the criterion is negative.',
  FULL_HF55321_forced=full21,actual_tail_not_chosen_unless_forced=not full21,
  no_new_jet_source_higher_natural_coefficients_or_cubic_values=True)
 save('actual-F2-source-products-and-m5-zero-proof.json',proof)
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 result=dict(id=ID,all_passed=True,frozen=True,actual_m5_Tnil_zero=True,actual_HF553111_excluded=True,known16_product_rank=productrank,extra_two_column_rank=pairrank,FULL_HF55321_forced=full21,actual_complete917_rank_options=[0,1],no_new_source_jet_or_higher_map_coefficients=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20
 save('report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
