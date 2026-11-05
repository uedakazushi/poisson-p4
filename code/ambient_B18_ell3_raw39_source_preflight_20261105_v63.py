#!/usr/bin/env python3
"""Pure exact cubic-truncation theorem and 75-column input for the ell³ base."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL3_RAW39_SOURCE_PREFLIGHT_20261105_v63'
PATHS={
 'source':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'sourcepeer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json',
 'raw2':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'raw2peer':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json',
 'table2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json',
 'trunc2peer':'reports/AMBIENT_B18_ELL2_RAW45_TRUNCATION_PREFLIGHT_20261105_v63_REVIEW/report.json',
 'trunc2note':'research/notes/AMBIENT_B18_ELL2_RAW45_TRUNCATION_PREFLIGHT_20261105_v63.md',
 'DVRnote':'research/notes/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63.md'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(p,q,c=Q(1)):
 r=dict(p)
 for m,v in q.items():
  x=r.get(m,Q(0))+c*v
  if x:r[m]=x
  else:r.pop(m,None)
 return r
def mul(p,q):
 r={}
 for (i,j),a in p.items():
  for (k,l),b in q.items():
   m=i+k,j+l
   if m[1]<3 and sum(m)<=4:r=add(r,{m:a*b})
 return r
def compose(v,r,s):
 out={}
 for mon,c in v:
  p={(0,0):Q(c)}
  for _ in range(mon[0]):p=mul(p,r)
  for _ in range(mon[1]):p=mul(p,s)
  out=add(out,p)
 return out
def enc(p):return [[list(m)if isinstance(m,tuple)else m,str(v)]for m,v in sorted(p.items())]
def main():
 t=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={p:sha(ROOT/p)for p in PATHS.values()}
 D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()if not k.endswith('note')}
 assert all(v['all_passed']for k,v in D.items()if k.endswith('peer'))
 assert len(D['raw2']['actual33_basis'])==33 and D['table2']['m5_zero']
 S=D['source'];assert len(S['quotient_basis_monomials'])==20
 r={(1,0):Q(1)};ell={(0,1):Q(1)};s=add(r,ell,Q(1,5));beta=[]
 for a,b in zip(S['beta2_before_ell'],S['beta3_before_ell']):
  aa=[[m[7:],c]for m,c in a];bb=[[m[7:],c]for m,c in b if not any(m[:7])]
  assert all(not any(m[:7])for m,c in a)
  beta.append(add(compose(aa,r,s),compose(bb,r,s)))
 gA=compose(S['gA_plane3'],r,s);bot=[row[7:]for row in S['YB_in_YA_linear'][7:]]
 rb=add({m:Q(bot[0][0])*c for m,c in r.items()},{m:Q(bot[0][1])*c for m,c in s.items()});sb=add({m:Q(bot[1][0])*c for m,c in r.items()},{m:Q(bot[1][1])*c for m,c in s.items()})
 gB=compose(S['gB_own_plane3'],rb,sb);assert gA[(3,0)]==Q(125,24)and gB[(3,0)]==Q(-15125,144)
 ambient=[(u,a,l)for u in range(8)for l in range(3)for a in range(5-l-(u!=0))];ambient.sort(key=lambda m:(m[1]+m[2]+(m[0]!=0),-m[0],m[2],m[1]));assert len(ambient)==75;ix={m:i for i,m in enumerate(ambient)}
 def term(u,p):return {ix[u,a,l]:v for(a,l),v in p.items()if(u,a,l)in ix}
 rel=[]
 for i in range(7):
  for j in range(i+1,7):rel.append((['pair',i,j],add(term(i+1,beta[j]),term(j+1,beta[i]),Q(-1))))
 for i in range(7):rel.append((['gA',i],add(term(i+1,gA),term(0,mul(gA,beta[i])))))
 for i in range(7):rel.append((['gB',i],term(i+1,gB)))
 rows=[]
 for label,v in rel:
  multipliers=[(0,0),(1,0),(0,1)]if label[0]=='pair'else[(0,0)]
  for a,l in multipliers:
   row={}
   for k,c in v.items():
    u,b,h=ambient[k]
    if(u,a+b,l+h)in ix:row[ix[u,a+b,l+h]]=c
   if row:rows.append([label,[a,l],enc(row)])
 assert len(rows)<=77 and sum(len(v)for lab,mul,v in rows)<=938
 data={
  'rings':'B=R0/(IA*IB), U=R0/(IA intersectIB), T=(IA intersectIB)/(IA*IB)',
  'same_actual_parameter':'ell=5*(sA-rA)',
  'actual_branch_DVR_ranks':[3,3],'actual_union_DVR_free_rank':6,
  'torsion_proof':'U embeds in the directsum of the two Dfree3 conductor rings, hence ellregular/free6. T is annihilated by J=IA+IB, for f∈IA andx∈intersection, fx∈IAIB (and similarlyIB). ActualJ containsell². Thus exact0→T→B→U→0 remainsleft-exact moduloell^n and T injects underB3→B2.',
  'corrected_actual_module_counts':{'B_mod_ell_length':20,'B_mod_ell2_length':33,'B_DVR_free_rank':6,'T_mu':14,'T_length':21,'T_ell_exponent1_count':7,'T_ell_exponent2_count':7,'B_mod_ell3_length':39,'raw_coker_mod_ell3_length':57},
  'union_basis6':['(1,1)','(r,r)','(r²,0)','(0,r²)','(ell²,0)','(ell*r,0)'],
  'union_weights6':[0,1,2,2,2,2],
  'union_order_proof':'Weierstrass r³=A2*r²+A1*r+A0 withordell(A2,A1,A0)≥(1,2,3). On the sixbasis, r raisesweighteddegree≥1 and eachzi (graphbetaorder≥2) raisesit≥2. Inparticularm^5 U⊂ell³ U; thisusesnoquarticcoefficient values.',
  'new_actual_m5_B3_zero_proof':'m5 B3 maps0 toU3, so liesT. ReductiontoacceptedB2 annihilatesm5 becauseactualraw33 hasm5zero. T injects underreductionB3→B2 asell²T=0. Thusm5B3 liesT intersectkerreduction=0.',
  'same_proof_for_cubic_truncation':'Truncatedpairs have thesameJ quadraticspan andDfree3 branches; acceptedell²truncation identifiestheirB2 withactualraw33. Hence theirB3 alsom5zero. Generatoromissions haveorder≥4 andproducterrors≥5, zero inbothB3 rings. Thus(IAIB,ell³)=(IA_le3*IB_le3,ell³) exactly.',
  'beta3_z_unit_scope':'AcceptedI+M2 generatornormalization removeszdependentbeta3; introducedplaneorder4terms do notchangeproductideal becausem5B3zero.',
  'ambient75_basis':ambient,'beta7_cubic_plane_after_s_equals_r_plus_ell_over5':list(map(enc,beta)),'gA_cubic':enc(gA),'gB_cubic':enc(gB),
  'source_plane_module_relations35':[[lab,enc(v)]for lab,v in rel],'all_atmost77_plane_multiple_rows75':rows,
  'input_nonzero_count':sum(len(v)for lab,mul,v in rows),'quotient_rank_expected_from_source_theorem':36,
  'ring_presentation_completeness':'Eliminate49monic zi(zj+beta_j). The21pair and14g rows areclosed underallzi actions; gAgB degree6zero. Degree≥5 iszero bytheorem. Thus these75ambient/≤77rows exactlypresentB3, with noGaussian/computedquotienttableyet.',
  'individual_branch_warning':'AlthoughproductB3 ignoresquarticconductorterms, actual R0/IA andR0/IB moduloell³ maydependonthem. Do nottakebranchquotientsusingonlycubicg. A laterexactsource105evaluation onB39 canrecoveractualrawbA/bB ideals andtheir9dimquotients withoutguessingquartics.',
  'future_base_builder_profile':{'ambient':75,'rows_cap':77,'input_terms_cap':938,'QQ_and_numeric_metadata_envelope':30000,'seconds':60,'memory_bytes':2<<30,'input_plus_output_bytes':2<<20,'no_physical_source_solve':True},
  'new_base_quotient_rank_table_or_physical_source_values_computed':False}
 ev=out/'exact-m5-ell3-cubic-truncation-and-75-column-source-input.json';ev.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'PROVED_exact_ell3_product39_cubic_source_reduction','actual_product_ell3_length':39,'actual_raw_coker_ell3_length':57,'actual_product_ell3_m5_zero':True,'actual_product_ell3_exact_cubic_truncation':True,'new_quotient_rank_or_physical_source_values_computed':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
