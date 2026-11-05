#!/usr/bin/env python3
"""70 point-C4 kernel checks and complete raw33->A18 ring-map source gate."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count,add,image
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_NATIVE70_A18_SOURCE_20261105_v63'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 R=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json')
 assert read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json')['all_passed']
 m=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['native_S0_columns'])
 C=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations']
 records=[]
 for i,c in enumerate(C):
  cn={j-640:Q(x) for j,x in c['original768_combination'] if j>=640}
  if cn:records.append({'case619':i,'original768_constant_combination':c['original768_combination'],'native128_constant_projection':[[j,str(x)] for j,x in sorted(cn.items())],'actual_a1_27_times9':c['actual_a1_27_times9']})
 assert [r['case619'] for r in records]==list(range(549,619)) and len(records)==70
 S0=list(map(sp,m['native_S0_columns']));assert len(S0)==128 and all(all(0<=i<48 for i in v) for v in S0)
 A=list(map(lambda M:list(map(sp,M)),R['coordinate_actions9_sparse_columns']));basis=R['actual33_basis'];assert len(basis)==33 and R['coordinate_order9']==['zA_0','zA_1','zA_2','zA_3','zA_4','zA_5','zA_6','rA','ell']
 pairs=[(7,7),(7,8),(0,7),(0,8),(4,7),(4,8),(5,7),(5,8)]
 mon=[()]+[(i,) for i in range(9)]+[tuple(sorted(p)) for p in pairs];mi={v:i for i,v in enumerate(mon)};assert len(mon)==18
 def word(triple):
  sector,r,e=triple;return ([sector-1] if sector else [])+[7]*r+[8]*e
 Pi=[{mi[tuple(sorted(word(v)))]:Q(1)} if tuple(sorted(word(v))) in mi else {} for v in basis]
 assert Pi[0]=={0:Q(1)} and {i for v in Pi for i in v}==set(range(18))
 def mul18(v,w):
  q={}
  for i,a in v.items():
   for j,b in w.items():
    k=mi.get(tuple(sorted(mon[i]+mon[j])))
    if k is not None:q=add(q,{k:a*b})
  return q
 positive_bits=0
 def mark(v):
  nonlocal positive_bits
  positive_bits=max(positive_bits,v[0].bit_length(),v[1].bit_length());return v
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def applypos(cols,v):
  q={}
  for j,b in v.items():
   for i,a in cols[j].items():q[i]=plus(q.get(i,(1,0)),mul(ex(a),b))
  return q
 # First derive priors for every point kernel and all raw multiplication
 # partial sums, before doing the exact source products.
 for rec in records:applypos(S0,{i:ex(x) for i,x in rec['native128_constant_projection']})
 for b in basis:
  for j in range(33):
   v={j:(1,1)}
   for k in word(b):v=applypos(A[k],v)
   applypos(Pi,v)
 for cols in A:
  for v in cols:applypos(Pi,{i:ex(x) for i,x in v.items()})
 assert positive_bits<=384
 # Now the bounded, independent full48-coordinate point identities.
 for rec in records:
  cn=sp(rec['native128_constant_projection']);assert not image(S0,cn)
  rec['all48_intrinsic_native_S0_components_exact_zero']=True
 for i,b in enumerate(basis):
  for j in range(33):
   v={j:Q(1)}
   for k in word(b):v=image(A[k],v)
   assert image(Pi,v)==mul18(Pi[i],Pi[j])
 for k,cols in enumerate(A):
  for j,v in enumerate(cols):assert image(Pi,v)==mul18({mi[(k,)]:Q(1)},Pi[j])
 owned=count([R,m,C,records,A,Pi,mon])+10000;assert owned<60000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'all70-C4-kernels-and-full-raw33-A18-ring-map-source.json'
 ev={'A18_basis_monomials':mon,'A18_coordinate_order9':R['coordinate_order9'],'selected_second_pairs8':pairs,'raw33_to_A18_columns':[[[i,str(x)] for i,x in sorted(v.items())] for v in Pi],
 'all1089_raw33_ring_products_exact':True,'all297_coordinate_intertwiners_exact':True,'unital_surjection_rank18':True,'all70_post27_native_point_kernel_records':records,
 'all70_native_S0_kernel_is_fullC4_not_scalar_origin':True,'source_positive_prior_bits':positive_bits,'source_owned_bound':owned,
 'conditional_second_extension_theorem':'LetC72=A18[tA,tB]/(tA^2-nuA,tB^2-nuB),free4overA18withm_A18^3=0. All48componentsofS0*cNare0,soactualS*cN liesm_A18*C72^12. IfNformalhasTHE_SAMEactualcoefficient0and1asN,DeltaNinm_A18^2*C72;DeltaN*(S*cN)=0. ThusweightednativeX2isindependentofN2. ThisdoesNOTapplytoindividualpivots;Xp1retainsLeS0.',
 'remaining_source_gates_before_actual_native70':['actualnuA/Braw33pushandC72moniclaw','actualkappa9_C132pushmatchesSAMEcanonicalbranch/C36andC16firsts','deriveREGULAR_CI_G8(notradical9)andKoszul28plusfree4usingliteralproperformula','actualN0/N1sourcebindingsall9directionsbeforeanysecondextension','Gamma/sourcepivotcarrysameframe','freshwholeheight_owned_RAW_andindependentpeerbefore70values'],
 'no_native70_second_responses_or_actual_G8_values':True,'no_new_physical_gauge_or_sourcejet':True}
 ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'all70_fullC4_point_kernel_verified':True,'all1089_ring_products_verified':True,'native_second_extension_theorem_conditional_on_Nfirst_bindings':True,'no_native70_response_values':True,'source_prior_bits':positive_bits,'source_owned_bound':owned,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'kernel70':True,'ringproducts':1089,'prior':positive_bits,'owned':owned,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
