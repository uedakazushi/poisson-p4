#!/usr/bin/env python3
"""Source-only cache and callback gate for the three zero-origin natural columns."""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_SOURCE_PREFLIGHT_20261104_v63';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 assert load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW')['all_passed']
 nat=load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63','actual-natural-Told-to-Tnil-source-map-and-origin.json');orig=nat['actual_origin_data'];O=load('AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63','actual-fullV56-m0-and-typed-comparison.json')
 X=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase2-actual-radical-X-firsts.json');M=load('AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63','actual-minimal-Tnil5-849-and-d3-linear-operator.json');m0=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase1-actual-m0-and-section.json');B=load('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63','actual-first640-X-mod-m2.json');At=load('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63','phase1-actual-atom-firsts.json');chain=load('AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63','actual-UP-QB-640-typed-operator.json');R=load('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json')
 oldcols=[23,29,30]
 for j,k in enumerate(oldcols,1):assert all(not F(row[k])for row in O['all24_by56_m0'])and Q.sparse(orig['all_five_old_cycle_origins56'][j])=={k:F(1)}
 src=load('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63','actual-source-coefficients.json');J=load('AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63','full9-firstjet-certificate.json');nodes=src['DAG_scalar_nodes'];rootfamilies=[(20153,'CA'),(40931,'CB'),(20144,'CA'),(40922,'CB')];profiles=[]
 for root,family in rootfamilies:
  seen=set();todo=[root]
  while todo:
   i=todo.pop()
   if i in seen:continue
   seen.add(i);a=nodes[i]
   if a[0]in['a','m','n','i']:todo.extend(a[1:])
   else:assert a[0]in['v','c']
  leaves={i for i in seen if nodes[i][0]=='v'};assert leaves<=set(J['both57_source2'][family]['unknown57'])
  uses=collections.Counter([root])
  for i in seen:
   if nodes[i][0]in['a','m','n','i']:uses.update(nodes[i][1:])
  live=set();peak=0
  for i in sorted(seen):
   live.add(i);peak=max(peak,len(live))
   if nodes[i][0]in['a','m','n','i']:
    for k in nodes[i][1:]:
     uses[k]-=1
     if not uses[k]:live.remove(k)
  profiles.append(dict(root=root,family=family,node_count=len(seen),leaf_count=len(leaves),all_leaves_in_cached57=True,closure=sorted(seen),maximum_live_degree2_vectors=peak))
 assert [(x['node_count'],x['leaf_count'])for x in profiles]==[(2189,48),(2194,48),(2197,48),(2202,48)]
 assert chain['actual_F12_C4'][5][0]==246 and chain['actual_F12_C4'][7][0]==250
 assert chain['operator_scalar_nodes'][246]==['sum',[[194,1,1],[245,1,1]]]and chain['operator_scalar_nodes'][250]==['sum',[[198,1,1],[249,1,1]]]
 assert [chain['operator_scalar_nodes'][i]for i in[194,195,198,199]]==[['actual_C0_scalar_function',i]for i in[20153,40931,20144,40922]]
 cA=[r for r in At['all35_ordered_CA_prefix_firsts']if r['actual_source_function']==7790];assert len(cA)==7
 cb=next(x for x in chain['accepted_actual_regular_function_identities']if x['kind']=='accepted_actual_CB_f63_tail');atom={i:Q.pairs(v)for i,v in At['all79_atom_values']};assert all(i in atom and not atom[i].get(0,F(0))for i in cb['coordinate_quotient_nodes'])
 # Cached exact pure-r restrictions provide the nine last-coordinate firsts.
 curves=[]
 for oldcol in oldcols:
  cc=[]
  for sec,power in enumerate([3,3,2]):
   v={int(k):F(c)for k,c in O['all56_sector_curve_restrictions'][oldcol][sec]};assert not v.get(power,F(0));cc.append(v.get(power+1,F(0)))
  curves.append(cc)
 # First18 origin cleanup: expose, do NOT erase, remaining A2 and native2 callbacks.
 free=[i for i in range(44)if i not in m0['selected_source_columns24']];I=[free[i]for i in M['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[Q.sparse(v)for v in M['fixed_minor_data']['inverse']];xp=X['X0_pivot_columns'];unknown_native=[{}for _ in range(780)]
 for j,col in enumerate(R['native_tensor_Koszul6_HB2_resolution'][1]['columns']):
  for row,pol in col:
   if row==5:
    for e,c in pol:
     if sum(e)==1 and any(e[:5]):Q.put(unknown_native[640+4*j],e.index(1),F(c))
 residual=[]
 for j,oldcol in enumerate(oldcols,1):
  c=Q.sparse(orig['all_five_new_cycle_origins44'][j]);a=Q.mvec(Ji,{k:c[i]for k,i in enumerate(I)if c.get(i,F(0))});ev={};sens={}
  for k,x in a.items():
   col=xp[k]
   if col<640:
    for row,sp in B['source_columns640_rows56_sparse_dual9'][col]:
     if row<48 and row%4==0:Q.put(ev,row,-x*Q.pairs(sp).get(0,F(0)))
   sens=Q.plus(sens,Q.scaled(unknown_native[col],-x))
  assert set(ev)<= {32};residual.append(dict(old_free56=oldcol,new18_pivot_coefficients=[[xp[k],str(v)]for k,v in sorted(a.items())],unknown_A_i5_even2_source_weights=[[i//4,str(v)]for i,v in sorted(ev.items())],unknown_native_u5_even2_weights=[[i,str(v)]for i,v in sorted(sens.items())],NOT_final81_cleanup=True))
 assert [F(x['unknown_A_i5_even2_source_weights'][0][1])for x in residual]==[F(11520,59821),F(40320,59821),F(-161280,59821)]
 maxlive=max(p['maximum_live_degree2_vectors']for p in profiles)
 cost=105*45+4*45+maxlive*45+324+3*44*8+3*18+3*81+15000
 assert cost<150000
 ev=dict(old_zero_origin_generators=oldcols,old_m0_columns_zero=True,old_n1_cancels='cold1=−nold0*mold1_ej because mold0_ej=0',old_second_derivatives_not_needed_for_cnew20='Zrow20_0=0',old_m1_possible_slots=324,old_last_coordinate_z_derivatives_zero_by_definition=True,nine_r_first_values_cached=curves,
  sector_function_matrix={'23':['−bB5*nuB','bA5*nuA','delta5'],'29':['delta7','0','−bB7'],'30':['0','delta7','bA7']},
  delta5='source20153−40931',delta7='source20144−40922',four_cached_point_function_profiles=profiles,CA_rawf63_first_prefixes=cA,CB_rawf63_prefix_nodes=cb['coordinate_quotient_nodes'],CB_prefix_origins_allzero=True,
  old_m1_Hadamard_rule='Hi firsts are degree2 f_(zi zq) with q<=i orq=r; firstf=0 for these columns so triangular z+c r² shift does not alter these coefficients. Lastcoordinatefirst is saved curve coefficient r^(power+1).',
  source_new_Z20_second_direct_columns='old trace masks3,1,2 respectively. even(mu(A)*tau_mask) carries nu; so these seconds use A_C36 firsts pluscachednu2, no A2.',
  remaining_after_first18_cleanup=residual,
  missing_gate='Do full81 zero-origin cleanup only after recovering oldm1. Then test exact A2/native-u5-even2 sensitivity in the fixed4 cotangent quotient. The one A8,5 callback cannot be guessed0; its weights may cancel but are not yet tested.',
  cache_only_oldm1_value_cost_upper=cost,all_4_point_AST_last_use_covered_by_saved57=True,prospective_value_job_not_launched=True,MS661_job_distinct=True,no_actual_zero3_first_map_or_HF_cokernel_decided=True)
 p=OUT/'three-oldm1-cache-and-cleanup-callback-gate.json';p.write_text(json.dumps(enc(ev),indent=2)+'\n')
 for p2 in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p2.relative_to(ROOT))]=sha(p2)
 rep=dict(id=ID,all_passed=True,frozen=True,source_only=True,three_old_n1_callbacks_cancel=True,old_m1_cache_recovery_source_gate_closed=True,one_A8_5_even2_callback_still_present_after18=True,final81_and_fixed4_sensitivity_gate_OPEN=True,no_new_source_or_jet_value=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st)
 (OUT/'report.json').write_text(json.dumps(rep,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
