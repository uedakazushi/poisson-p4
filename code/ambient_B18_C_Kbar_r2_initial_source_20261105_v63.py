#!/usr/bin/env python3
"""Tiny cached SOURCE: r²Kbar projected to old gr², higher tails formal."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63'
PATHS={
 'old40':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/known-initial40-and-length-upper17.json',
 'old40_peer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63_REVIEW/report.json',
 'linear':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/linear-module-rewrite-certificate.json',
 'exact':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/exact-Loewy11-and-source-basis-corollary.json',
 'exact_peer':'reports/AMBIENT_B18_T_COMBINED619_QUARTIC_WITNESS_20261105_v63_REVIEW/report.json',
 'seven':'reports/AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63/actual-seven-barbasis-typed56-cycle-operator-and-full11-unit-complement.json',
 'values_peer':'reports/AMBIENT_B18_C_CONNECTING_HEADS2930_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'R3_peer':'reports/AMBIENT_B18_C_CONNECTING_R3_FACTORIZATION_SOURCE_20261105_v63_REVIEW/report.json',
}
MAXPRE=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def h(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def pre(n):
 global MAXPRE
 MAXPRE=max(MAXPRE,n);assert n<=1153
def addq(a,b):pre(2*max(h(a),h(b))+1);q=a+b;assert h(q)<=384;return q
def mulq(a,b):pre(h(a)+h(b));q=a*b;assert h(q)<=384;return q
def divq(a,b):assert b;pre(h(a)+h(b));q=a/b;assert h(q)<=384;return q
def parse(s):
 a=str(s).split('/');n=int(a[0]);d=int(a[1])if len(a)>1 else 1
 assert d;pre(max(abs(n).bit_length(),abs(d).bit_length()));q=Q(n,d);assert h(q)<=384;return q
def exponent(*ii):return tuple(ii.count(i)for i in range(8))
def enc(p):return[[t,list(m),str(q)]for(t,m),q in sorted(p.items())]
def polyadd(p,q,c=Q(1),shift=(0,)*8):
 out=p.copy()
 for(t,m),v in q.items():
  key=(t,tuple(a+b for a,b in zip(m,shift)));out[key]=addq(out.get(key,Q(0)),mulq(c,v))
  if not out[key]:del out[key]
 return out
def order(k):return(-k[0],)+k[1]
def leading(p):return max(p,key=order)
def reduce(p,gens):
 source=p.copy();rem={};tr=[]
 while p:
  key=leading(p)
  for i,g in enumerate(gens):
   l=leading(g)
   if l[0]==key[0]and all(a<=b for a,b in zip(l[1],key[1])):
    q=divq(p[key],g[l]);sh=tuple(b-a for a,b in zip(l[1],key[1]))
    p=polyadd(p,g,-q,sh);tr.append([i,list(sh),str(q)]);assert len(tr)<1000;break
  else:rem[key]=p.pop(key)
 # Exact membership cofactor multiplyback, independently from the order.
 rebuilt=rem.copy()
 for i,sh,q in tr:rebuilt=polyadd(rebuilt,gens[i],parse(q),tuple(sh))
 assert rebuilt==source
 return rem,tr
def main():
 start=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={'old40':['homogeneous_module_generators40','standard_degree1','standard_degree2'],
  'linear':['linear_basis37'],'exact':['exact_length','exact_Loewy_HF','actual_filtered_monomial_basis11'],
  'seven':['actual_Kbar11_by7','actual_correction_coefficients','actual_four_head_parameters']}
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 prior=S.count(d)+20000;assert prior<150000
 raw=sum((ROOT/p).stat().st_size for p in PATHS.values());assert raw+200000<33554432
 assert d['exact']['exact_length']==11 and d['exact']['exact_Loewy_HF']==[5,3,2,1]
 assert d['old40']['standard_degree1']==['e3*r','e4*z5','e4*r']
 assert d['old40']['standard_degree2']==['e4*z5*r','e4*r²']
 gens=[{(t,tuple(m)):parse(q)for t,m,q in g}for g in d['old40']['homogeneous_module_generators40']]
 l37=[{(ix//8,exponent(ix%8)):parse(q)for ix,q in row}for row in d['linear']['linear_basis37']]
 assert gens[:37]==l37
 leads={leading(g)for g in l37};expected={(t,exponent(i))for t in range(5)for i in range(8)}-{(3,exponent(7)),(4,exponent(5)),(4,exponent(7))}
 assert leads==expected and len(leads)==37
 q1=(4,exponent(5,7));q2=(4,exponent(7,7));v=(4,exponent(5));u=(3,exponent(7))
 linear_results=[];alphas=[]
 for t in[1,2]:
  inp={(t,exponent(7)):Q(1)};nf,tr=reduce(inp,l37)
  if t==1:assert set(nf)<={v}
  else:assert nf.get(u)==Q(1,4)and set(nf)<={u,v}
  alphas.append(nf.get(v,Q(0)));linear_results.append({'head':t,'input':enc(inp),'normal_form':enc(nf),'exact_linear_source_cofactors':tr})
 K=[[parse(x)for x in row]for row in d['seven']['actual_Kbar11_by7']]
 assert len(K)==11 and len(K[0])==7
 projected=[];receipts7=[]
 for j in range(7):
  # ONLY constant-head rows can contribute modulo oldm³. The six higher
  # old11 rows begin in mTold, and r² sends them into m³Told.
  inp={(t,exponent(7,7)):K[t][j]for t in range(5)if K[t][j]}
  nf,tr=reduce(inp,gens)
  assert set(nf)<={q1,q2}
  projected.append([nf.get(q1,Q(0)),nf.get(q2,Q(0))]);receipts7.append({'Kbar_basis_column':j,'degree2_input':enc(inp),'normal_form':enc(nf),'exact_old40_source_cofactors':tr})
 assert all(row==[0,0]for row in projected[:5])and all(row[1]==0 for row in projected)
 a3=parse(d['seven']['actual_correction_coefficients']['a3=b3/b1']);as_=parse(d['seven']['actual_correction_coefficients']['aS=(b2-b3/4)/b1'])
 c=Q(11767,3652)
 expectedU=addq(-c,-mulq(a3,alphas[0]));expectedS=addq(alphas[1],-mulq(as_,alphas[0]))
 assert projected[5][0]==expectedU and projected[6][0]==expectedS
 rank=int(any(row[0]for row in projected));witness=next((j for j,row in enumerate(projected)if row[0]),None)
 data={'scope':'Tiny cached exact initial reductions only: r²Kbar→m²Told/m³Told. NOT a full action on Kbar, not ell-DVR torsion data.',
  'same_original8_variables':['z0','z1','z2','z3','z4','z5','z6','r'],
  'before_arithmetic_SOURCE_preflight':{'owned_upper':prior,'RAW_with_output_reserve':raw+200000,'committed_bits':384,'temporary_bits':1153,'seconds':60,'memory_bytes':2000000000},
  'all37_linear_generators_literal_match_cached_source':True,'leading37_complement_exactly_known_three_standards':True,
  'two_exact_linear_reductions_and_membership_cofactors':linear_results,
  'actual_alpha1_alpha2':list(map(str,alphas)),
  'actual_gr2_basis':['q1=z5*r*h4','q2=r²*h4'],
  'all7_exact_degree2_reductions_and_membership_cofactors':receipts7,
  'actual_projected_r2_matrix2_by7':[[str(row[i])for row in projected]for i in range(2)],
  'projected_rank':rank,'fixed_nonzero_witness_column':witness,
  'formula_U_S_q1':['-11767/3652-a3*alpha1','alpha2-aS*alpha1'],
  'actual_initial_membership_and_completeness_proof':'The37 literal source linear initials and G38 are actual initial relations by accepted source/peer. Exact sourcecofactor multiplybacks certify everyreduction. Their degree2 quotient has two displayed standards; accepted exactTold HF2=2 proves these ARE the actual gr²basis. No unknown higher initial relation is discarded in degree2.',
  'higher_tails_retained':'For actual r²Kbar, each displayed vector may additionally have a component along oldm³Told=QQ*r³h4. Those components remain FORMAL/uncomputed. No full r² matrix/r-Jordan choice is asserted.',
  'one_valid_consequence':'If projectedrank1 then actualr²Kbar!=0; together with independently provedr³Kbar0, at leastone r-block has length3. This does not determine fullr-Jordan partition or anyell-DVR/H0invariant.',
  'no_new_root_Physical_IFT_gauge_jet_or_full_action_query':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),
  'maximum_exact_QQ_preop_bits':MAXPRE,'selective_receipts':receipts,'all_passed':True}
 owned=S.count(d)+S.count(data)+S.count(gens)+S.count(l37)+S.count(K)+10000;assert owned<150000;data['whole_SOURCE_owned_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'two-cached-alpha-reductions-and-projected-r2-seven-source-certificates.json';p.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_CACHED_INITIAL_ONLY_PROJECTED_R2_KBAR_NO_HIGHER_ACTIONS',
  'actual_alpha1_alpha2':list(map(str,alphas)),'projected_r2_rank':rank,'input_sha256':pins,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{p.name:sha(p)},'whole_SOURCE_owned_upper':owned,'maximum_exact_QQ_preop_bits':MAXPRE,'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)
if __name__=='__main__':main()
