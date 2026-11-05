#!/usr/bin/env python3
"""Complete first-linear-syzygy block, with cached-plane-only cubic adjoints."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_CPRIME_RIGHT85_SOURCE_PREFLIGHT_20261105_v63'
PATHS={
 'first':'reports/AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63/phase1-actual-fullnine-first-V-columns.json',
 'quadratic':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase4-all-twenty-actual-quadratic-values.json',
 'quadratic_peer':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63_REVIEW/report.json',
 'dual':'reports/AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_SOURCE_PREFLIGHT_20261105_v63/actual-curve-source-and-one-residual-cubic-functional.json',
 'dual_peer':'reports/AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_SOURCE_20261105_v63_REVIEW/report.json',
 'closure_peer':'reports/AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def sp(v):return [[i,c]for i,c in sorted(v.items())if c]
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def add(v,w,c=Q(1)):
 z=dict(v)
 for i,a in w.items():
  z[i]=z.get(i,Q(0))+c*a
  if not z[i]:del z[i]
 return z
def count(v):
 if isinstance(v,(Q,int)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(count(x)for x in v.values())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 return 0
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(D[k]['all_passed']for k in ['quadratic_peer','dual_peer','closure_peer'])
 first=D['first'];L=[rec(v)for v in first['all44_actual_linear_columns27_sparse']];cp=first['fixed_pivot_columns'];assert len(cp)==24
 mons2=list(itertools.combinations_with_replacement(range(9),2));mi2={m:i for i,m in enumerate(mons2)};mons3=list(itertools.combinations_with_replacement(range(9),3));mi3={m:i for i,m in enumerate(mons3)}
 columns=[]
 for j in cp:
  for var in range(9):
   v={}
   for i,c in L[j].items():target,x=divmod(i,9);v=add(v,{45*target+mi2[tuple(sorted((x,var)))]:c})
   columns.append(v)
 piv={};pc=[];pr=[];piv_values=[];kernel=[];free=[];input_count=count(D);peak=0
 def check(live):
  nonlocal peak
  peak=max(peak,input_count+count(live));assert peak<=60000,peak
  for group in [piv_values,kernel]:
   for v in group:
    for c in v.values():assert max(abs(c.numerator).bit_length(),c.denominator.bit_length())<=2048
 for j,source in enumerate(columns):
  v=dict(source);a={j:Q(1)}
  for row,(w,b)in sorted(piv.items()):
   c=v.get(row,Q(0))
   if c:v=add(v,w,-c);a=add(a,b,-c)
  if v:
   row=min(v);c=v[row];w={i:z/c for i,z in v.items()};b={i:z/c for i,z in a.items()};piv[row]=(w,b);pc.append(j);pr.append(row);piv_values.append(w)
  else:kernel.append(a);free.append(j)
  check([L,columns,piv,pc,pr,piv_values,kernel,free,v,a])
 assert len(pc)==131 and len(kernel)==85
 digest=hashlib.sha256()
 for k,v in enumerate(kernel):
  result={}
  for j,c in v.items():result=add(result,columns[j],c)
  assert not result and all(v.get(j,Q(0))==Q(k==i)for i,j in enumerate(free));digest.update(json.dumps(enc(sp(v)),separators=(',',':')).encode()+b'\n')
 # Witness the rank without any new131-squareinverse: each normalizedpivot
 # is an explicit combination of the216 original columns and has distinct
 # leading row; its selected-column transformation is triangular with units.
 witnesses=[]
 for row,(w,a)in sorted(piv.items()):
  result={}
  for j,c in a.items():result=add(result,columns[j],c)
  assert result==w and min(w)==row and w[row]==1;witnesses.append([row,sp(w),sp(a)])
 dual=rec(D['dual']['residual_actual_cubic_dual495_sparse']);assert all(mons3[i%165][0]>=7 for i in dual);folded=[];zero_z=0
 for syz in kernel:
  a={}
  # Only pureplane entries can hit thisliteraldual. Everyz-variablefactor
  # killsallterms, regardless of unavailable nonplane Q2 coefficients.
  for source,c in syz.items():
   j,var=divmod(source,9);old=cp[j]
   if var<7:zero_z+=1;continue
   for target in range(3):
    for k,mon in enumerate([(7,7),(7,8),(8,8)]):
     coefficient=c*dual.get(165*target+mi3[tuple(sorted(mon+(var,)))],Q(0))
     if coefficient:a=add(a,{9*old+3*target+k:coefficient})
  folded.append(a)
 assert all((key//9)in cp for a in folded for key in a)
 # Do not pair theseadjoints with actualQ2 values inthisSOURCEbuilder.
 data={'coordinate_order9':first['coordinate_order9'],'actual_V1_fixed24_source_columns':cp,'all216_first_product_columns135_sparse':list(map(sp,columns)),
  'rank131_pivot_original_columns':pc,'rank131_pivot_rows_creation_order':pr,'rank131_sparse_normalized_pivot_and_source_witnesses':witnesses,
  'all85_syzygy_columns216_sparse':list(map(sp,kernel)),'free_source_rows85':free,'free_tail_I85_verified':True,'all85_source_multiplybacks_zero':True,'kernel_sha256':digest.hexdigest(),
  'all85_cached_plane_Q2_adjoints396_sparse':list(map(sp,folded)),'adjoint_leaf_order':'9*originalVcolumn +3*GammaTarget +planePairIndex; pairs(rr,rell,ellell)',
  'all_required_actual_leaves':sorted({i for a in folded for i in a}),'z_factor_callback_terms_proved_zero':zero_z,
  'all_higher_source_or_gauge_callbacks_absent':True,'no_V3_needed':'Eachu1 is anactualhomogeneouslinearcoefficientvector onliteral24originalVcolumns. V1*u1=0 identically, soactualrelationV*u1 hasinitialcubic V2*u1. No source-column higherchange orconstant-kernel correction occurs.',
  'complete_firstlinear_syzygy_block':85,'not_full_constant19_actual_cubic_block':True,'future_value_shape':[1,85],
  'future_input_plus_saved_live_numeric_cap':60000,'future_derived_bit_guard':2048,'future_seconds':60,'future_memory_bytes':2<<30,'future_input_plus_output_bytes':32<<20,
  'source_builder_peak_numeric_slots':peak,'actual85_values_not_evaluated':True,'positive_criterion':'Anynonzero addsactualinitial toknownrank26, provingupper9/e2; allzeroonlypartialnegative.'}
 ep=out/'complete-right85-and-cached-plane-source-adjoints.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert sum((ROOT/p).stat().st_size for p in PATHS.values())+len(blob.encode())<32<<20;ep.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_complete_right85_cached_plane_adjoint','input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'actual85_values_evaluated':False,'peak_numeric_slots':peak,'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps(report))
if __name__=='__main__':main()
