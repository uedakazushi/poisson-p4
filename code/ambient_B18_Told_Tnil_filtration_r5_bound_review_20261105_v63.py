#!/usr/bin/env python3
"""Independent fixed source/filtration audit, no producer import or new jet."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63_REVIEW'
AUTHOR='AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def put(p,k,c):
 if c:
  p[k]=p.get(k,Q(0))+c
  if not p[k]:del p[k]
def add(p,q,c=Q(1)):
 p=dict(p)
 for k,a in q.items():put(p,k,c*a)
 return p
def dot(row,v):return sum((c*v.get(i,Q(0))for i,c in row.items()),Q(0))
def mv(A,v):return {i:c for i,row in enumerate(A)if(c:=dot(row,v))}
def determinant(A):
 A=[r[:]for r in A];n=len(A);d=Q(1)
 for j in range(n):
  k=next((i for i in range(j,n)if A[i][j]),None)
  if k is None:return Q(0)
  if k!=j:A[k],A[j]=A[j],A[k];d=-d
  p=A[j][j];d*=p
  for i in range(j+1,n):
   c=A[i][j]/p
   for k in range(j+1,n):A[i][k]-=c*A[j][k]
 return d
def poly(rows):return {(t,tuple(e)):Q(c)for t,e,c in rows if Q(c)}
def shift(p,e,c=Q(1)):return {(t,tuple(a+b for a,b in zip(m,e))):v*c for(t,m),v in p.items()if v*c}
def unit(i,n=1):return tuple(n if k==i else 0 for k in range(8))
def pterms(p):return [[t,list(e),c]for(t,e),c in sorted(p.items())]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True);assert not(out/'report.json').exists()
 ap=ROOT/'reports'/AUTHOR/'report.json';ar=json.loads(ap.read_text());assert ar['all_passed']and sha(ap)=='6f74542bb99048ce026ffb340a457ace53819a866d57be5b3384c41331560faa'
 pins={str(ap.relative_to(ROOT)):sha(ap)}
 for p,h in ar['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 E=read(f'reports/{AUTHOR}/actual-F2-source-products-and-m5-zero-proof.json')
 old=read('reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/exact-Loewy11-and-source-basis-corollary.json');G=read('reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/known-initial40-and-length-upper17.json')
 nat=read('reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json');Z=read('reports/AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63/actual-three-first-columns-and-fixed-quotient-proof.json')
 F=read('reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json');quad=read('reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63/literal-known-cubic-input-and-exact24-row-reduction.json')
 red=read('reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63_REVIEW/report.json');assert red['all_passed']and red['actual_m3_equals_rA_m2']
 assert Z['actual_cokernel_HF']==[3,1,1]and Z['actual_cokernel_length']==5
 assert old['exact_length']==11 and old['dimensions_m_power_T'][-1]==0
 # Source degree1 basis is actual filtered basis, not merely a guessed normal form.
 basis=[(p['target'],p['monomial_original8'])for p in old['actual_filtered_monomial_basis11']if p['degree']==1]
 assert basis==[(3,list(unit(7))),(4,list(unit(5))),(4,list(unit(7)))]
 origin=[{i:Q(c)for i,c in enumerate(v)if Q(c)}for v in nat['actual_origin_data']['actual_origin_map_columns5']]
 assert not any(origin[i]for i in [1,2,3])
 P5=[{int(i):Q(c)for i,c in row.items()}for row in Z['quotient5_by40_after_actual_rank35_first_relations']]
 for col in F['all849_actual_minimal_first_columns40_sparse']:assert not mv(P5,dict((i,Q(c))for i,c in col))
 images=[]
 for head,var in [(3,7),(4,5),(4,7)]:images.append(mv(P5,{8*i+var:c for i,c in origin[head].items()}))
 assert not images[0]and not images[1]and images[2]
 new3=[{int(i):Q(c)for i,c in row.items()}for row in Z['actual_three_minimal_first_columns40']];zeroimages=[mv(P5,v)for v in new3];columns=zeroimages+[images[2]];rr=[0,2,3,4]
 det4=determinant([[c.get(i,Q(0))for c in columns]for i in rr]);assert det4==Q(86403125,7746048)
 # Own reversed-generator fixed reduction, and independent verification of saved cofs.
 generators=[poly(v)for v in G['homogeneous_module_generators40']];lm=[max(g,key=lambda x:(-x[0],x[1]))for g in generators]
 def divide(p):
  v=dict(p);nf={};co={};steps=0
  while v:
   term=max(v,key=lambda x:(-x[0],x[1]));target,e=term;c=v[term];chosen=None
   for j in range(39,-1,-1):
    tt,ee=lm[j]
    if target==tt and all(a>=b for a,b in zip(e,ee)):chosen=j;break
   if chosen is None:put(nf,term,c);del v[term];continue
   j=chosen;ee=tuple(a-b for a,b in zip(e,lm[j][1]));a=c/generators[j][lm[j]];put(co,(j,ee),a);v=add(v,shift(generators[j],ee),-a);steps+=1;assert steps<2000
  rebuilt=dict(nf)
  for(j,e),c in co.items():rebuilt=add(rebuilt,shift(generators[j],e),c)
  assert rebuilt==p;return nf,co
 checks=[];nf2=[]
 for index,case in enumerate(E['all16_products_with_exact_initial_ideal_cofactors']):
  p=poly(case['input']);assert case['variable']==index%8
  seed={(3,unit(7)):Q(1)}if index<8 else{(4,unit(5)):Q(1)}
  assert p==shift(seed,unit(index%8))
  nf,co=divide(p);assert nf==poly(case['normal_form'])
  saved=poly(case['normal_form'])
  for j,e,c in case['exact_initial_ideal_cofactors']:saved=add(saved,shift(generators[j],tuple(e)),Q(c))
  assert saved==p
  nf2.append(nf);checks.append({'input':pterms(p),'own_reversed_normal_form':pterms(nf),'own_reversed_cofactors':[[j,list(e),c]for(j,e),c in co.items()]})
 q1=(4,tuple(a+b for a,b in zip(unit(5),unit(7))));q2=(4,unit(7,2))
 assert any(v.get(q1,Q(0))for v in nf2)and all(not v.get(q2,Q(0))for v in nf2)
 for case in E['all_r2_u_v_degree3_zero_NF_cofactors']:
  p=poly(case['input']);nf,co=divide(p);assert not nf
  rebuilt={}
  for j,e,c in case['cofactors']:rebuilt=add(rebuilt,shift(generators[j],tuple(e)),Q(c))
  assert rebuilt==p
 g38={(3,unit(7,2)):Q(1),(4,q1[1]):Q(11767,3652)};qC={(4,tuple(a+b for a,b in zip(unit(5),unit(7,2)))):Q(1)}
 assert generators[38]==g38 and generators[39]==qC
 # Extra2 diagnostic in the ACTUAL degree2 quotient; source paths are independent.
 pi3=[dict((i,Q(c))for i,c in row)for row in quad['quadratic_projection3_by180_sparse']];m2=list(itertools.combinations_with_replacement(range(8),2));mi={m:i for i,m in enumerate(m2)};u2={}
 for row,c in new3[2].items():target,var=divmod(row,8);put(u2,36*target+mi[tuple(sorted((var,7)))],c)
 v2={36*target+mi[(7,7)]:c for target,c in origin[4].items()};a,b=mv(pi3,u2),mv(pi3,v2);det2=determinant([[a.get(i,Q(0)),b.get(i,Q(0))]for i in [1,2]])
 assert det2==Q(-60595625,1871074032)
 proof={'same_actual_S_and_rA_not_ell':True,'all_input_hashes_verified':True,'own_F2_source_minor':det4,'own_oldgr1_images':images,'own_ZERO3_first_images':zeroimages,
  'all16_saved_and_own_reverse_cofactors_verified':True,'own16_checks':checks,'own16_degree2_rank':1,'own_extra_pair_determinant':det2,'both_stronger_forcing_tests_negative':True,
  'F2_actual_equality':'F2=m²Told + k e3r + k e4z5, length5; no strict filtration assumed',
  'mF2_actual_length':2,'mF2_justification':'mF2 contains m³Told automatically since m²Told subsetF2; its quotient modulo m³ has fixed rank1',
  'r2_F2_actual_zero':'r²m²Told subsetm4Told=0. G38 lifts with remainder>=3; r*remainder liesm4=0. qC lifts with remainder>=4=0. Thus r²e3r=r²e4z5=0.',
  'full_torsion_module_m5_zero':'m³Q=0 identifies m³Tnil with F3 subsetF2 insideTold; r-reduction givesm5Tnil=r²m³Tnil=0',
  'actual_HF_options':[[5,5,3,3],[5,5,3,2,1]],'actual_complete917_rank_options':[0,1],
  'core_scope':'finite module Tnil only; no nil-ideal power or cube conclusion',
  'frozen_note_typo_clarification':'Opening HF331 is a typo; code/input and proof use actual Q HF311 oflength5. Author bytes unchanged.',
  'no_author_import_or_new_source_jet_ranksearch':True}
 p=out/'independent-fixed-filtration-and-r5-proof.json';blob=json.dumps(enc(proof),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 r={'id':STEM,'all_passed':True,'actual_m5_Tnil_zero_accepted':True,'actual_HF553111_excluded':True,'actual_HF_options':[[5,5,3,3],[5,5,3,2,1]],'actual_complete917_rank_options':[0,1],'new_source_jet_or_higher_map_coefficient':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()
