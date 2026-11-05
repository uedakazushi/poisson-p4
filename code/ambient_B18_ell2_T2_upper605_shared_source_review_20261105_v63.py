#!/usr/bin/env python3
"""Independent complete605 inventory and atom-first OLD2 SOURCE contraction."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_ELL2_T2_UPPER605_SHARED_SOURCE_20261105_v63';ID=BASE+'_REVIEW'
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 pins={};raw=0;pre=0;bits=0;cursor={}
 def read(p,h=None,parse=True):
  nonlocal raw
  p=ROOT/p;k=str(p.relative_to(ROOT));z=sha(p)
  if h:assert z==h,(k,h,z)
  if k not in pins:raw+=p.stat().st_size
  pins[k]=z;return json.loads(p.read_text())if parse else None
 rp=f'reports/{BASE}/report.json';r=read(rp,'af0728e250dda8979857b60e1d414e19399f311f46bd756cb1382deb4df3cc89');assert r['all_passed']
 en=next(iter(r['evidence_sha256']));E=read(f'reports/{BASE}/{en}',r['evidence_sha256'][en])
 inputs={p:read(p,h)for p,h in r['input_sha256'].items()}
 read('code/ambient_B18_ell2_T2_upper605_shared_source_20261105_v63.py',r['code_sha256'],False)
 def suffix(s):return inputs[next(p for p in inputs if p.endswith(s))]
 j=suffix('617-direct-old-known-first-formal-Jacobian-positive-input.json')
 s=suffix('combined-OLD2-D-and-four-pivot-source-factor-inventory.json')
 jf=suffix('617-formal-direct-old-third-and-second-symbol-sensitivities.json')
 of=suffix('actual617-formal-combined-OLD2-second-symbol-coefficients.json')
 phys=suffix('actual-source105-ell2-raw33-values.json')
 native=suffix('all617-native-constant-sources-factor-through-fixed105-kernel.json')
 for p,v in inputs.items():
  if p.endswith('/report.json'):assert v['all_passed']
 # Bind the formal certificates to the same original graph, not a name-only
 # citation. Accepted independent formal proofs remain lemmas, without any
 # assigning actual second/third values to their formal symbol variables.
 for st,data in [('AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63',jf),('AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63',of)]:
  rr=read(f'reports/{st}/report.json');assert rr['all_passed']
  for p,h in rr['input_sha256'].items():read(p,h,False)
  nn=next(iter(rr['evidence_sha256']));assert sha(ROOT/'reports'/st/nn)==rr['evidence_sha256'][nn]
  primary=next(p for p in pins if p.endswith(nn));assert pins[primary]==rr['evidence_sha256'][nn]
 globalpeer=read('reports/AMBIENT_B18_ELL2_CONDUCTORS_GLOBAL_M3_SOURCE_20261105_v63_REVIEW/report.json','fd0d61dabb2751bfe54eb858aa3ba3774e35da46f1b3fb36e0b461a4dc07b239');assert globalpeer['all_passed']
 gt=dict(s['typed_program']);jt=dict(j['typed_source_program']);assert len(gt)==757 and len(jt)==629 and all(gt[k]==v for k,v in jt.items())
 assert len(j['literal_C0_program'])==532 and len(gt)==len(s['typed_program'])
 boundaries={i:(kind,v)for i,kind,v in j['literal_boundary_firsts_point_r_ell_and_membership']};assert len(boundaries)==96 and all(v[0]=='0'for kind,v in boundaries.values())
 leaf67={n[1]for n in gt.values()if n[0]=='actual_C0_scalar_function'};assert len(leaf67)==67 and leaf67<=set(boundaries)
 seen=set();external=set()
 for i,n in j['literal_C0_program']:
  assert i not in seen and n[0]in['c','n','a','m'];deps=[]if n[0]=='c'else n[1:2]if n[0]=='n'else n[1:3]
  external.update(k for k in deps if k not in seen);seen.add(i)
 assert external<=set(boundaries)
 # All typed atoms vanish at the origin. Consequently the derivative of a
 # multi-product in independent second coefficients has no linear term.
 # Keep its full original graph; this derivative is not a value evaluator.
 support={}
 for i,n in s['typed_program']:
  if n[0]=='actual_C0_scalar_function':v={n[1]}
  elif n[0]=='sum':v=set().union(*(support[k]for k,a,b in n[1]if a))if n[1]else set()
  else:
   assert n[0]=='product'and n[3]and all(k in support for k in n[3]);v=set(support[n[3][0]])if len(n[3])==1 and n[1]else set()
  support[i]=v
 Dsupport=set().union(*(support[k]for c,roots in s['D_old_root_columns']for row,k in roots))
 opaque={330164,330184,330289,330309,330524,330544,330594,330614,566400,566401,566402,566403};conductors={505938,299650}
 assert not Dsupport&opaque
 outer=s['outer_four_pivot_source'];assert not outer['inverse_nodes']
 assert len(outer['typed_program'])==8 and all(n[0]=='actual_C0_scalar_function'and n[1]in opaque for i,n in outer['typed_program'])
 assert all(n[0]in['g','n','a']for i,n in outer['literal_C0_program'])
 assert {i for i,n in outer['literal_C0_program']if n[0]=='g'}==set(range(566400,566404))
 assert outer['original4_old_pivots']==[604,613,622,631]
 recipes=s['kernel617_source_recipes'];free={x['free_case619']for x in recipes};assert len(recipes)==617 and free==set(range(619))-{563,609}
 assert j['original619_source768_combinations']==s['original619_c0_source_combinations']
 f1={x['free_case619']:x['two_formal_J0_Jr_Jell_symbol_rows']for x in jf['all617_formal_directOLD_J0_Jr_Jell']};f2=dict(of['all617_combined_OLD2_independent_symbol_coefficients']);assert set(f1)==set(f2)==free
 parsed=count([r,E,inputs,j,s,jf,of,phys,native,globalpeer,gt,jt,boundaries,leaf67,external,support,Dsupport,recipes,f1,f2,pins])
 # Sparse derivative maps/cached column contractions +617 result copies and
 # all dictionary labels/controllers are charged before exact arithmetic.
 prior=parsed+45000;assert prior<150000 and raw+4_000_000<32*1024*1024
 save(out/'positive-before-atom-first-formal-contractions.json',{'all_decoded_numeric_labels_and_alias_maps':parsed,'whole_owned_numeric_upper':prior,'RAW_input_plus_output_upper':raw+4000000,'rational_committed_hard':384,'temporary_hard':1153,'known_positive_formal2_prior120':True,'no_actual_second_third_605_response':True})
 def ck(x):
  nonlocal bits
  b=max(abs(x.numerator).bit_length(),x.denominator.bit_length());bits=max(bits,b);assert b<=384,(cursor,b);return x
 def mul(x,y):
  nonlocal pre
  b=max(abs(x.numerator).bit_length()+abs(y.numerator).bit_length(),x.denominator.bit_length()+y.denominator.bit_length());pre=max(pre,b);assert b<=1153;return ck(x*y)
 def add(x,y):
  nonlocal pre
  b=max(abs(x.numerator).bit_length()+y.denominator.bit_length(),abs(y.numerator).bit_length()+x.denominator.bit_length())+1;pre=max(pre,b,x.denominator.bit_length()+y.denominator.bit_length());assert pre<=1153;return ck(x+y)
 def madd(v,k,x):
  z=add(v.get(k,Q(0)),x)
  if z:v[k]=z
  else:v.pop(k,None)
 def axpy(v,w,a):
  for k,x in w.items():madd(v,k,mul(a,x))
 # Independent atom-first derivative, in reverse child summation order.
 # This differs from the producer's per-case identical-root fold first.
 grad={}
 for i,n in s['typed_program']:
  z={}
  if n[0]=='actual_C0_scalar_function':z={n[1]:Q(1)}
  elif n[0]=='sum':
   for k,a,b in reversed(n[1]):axpy(z,grad[k],Q(a,b))
  elif len(n[3])==1:axpy(z,grad[n[3][0]],Q(n[1],n[2]))
  grad[i]=z;assert set(z)<=support[i]
 D=[[[Q(x)for x in row]for row in ch]for ch in s['D_two_by10_by56']]
 Dcols={}
 for old,roots in reversed(s['D_old_root_columns']):
  channels=[]
  for ch in range(2):
   z={}
   for row,k in reversed(roots):
    for fn,a in grad[k].items():
     for pair in range(10):
      if D[ch][pair][row]:madd(z,(pair,fn),mul(a,D[ch][pair][row]))
   channels.append(z)
  Dcols[old]=channels
 assert not {fn for cols in Dcols.values()for z in cols for pair,fn in z}&opaque
 og={fn:{fn:Q(1)}for fn in opaque}
 for i,n in outer['literal_C0_program']:
  if n[0]=='g':assert i in og
  elif n[0]=='n':og[i]={k:-x for k,x in og[n[1]].items()}
  else:og[i]=dict(og[n[2]]);axpy(og[i],og[n[1]],Q(1))
 G=[{i:Q(x)for i,x in row}for row in outer['Gamma_two']];ot={i:n[1]for i,n in outer['typed_program']};Ocols={}
 for old,roots in outer['four_source_root_columns']:
  Ocols[old]=[]
  for ch in range(2):
   z={}
   for row,(kind,k)in reversed(roots):
    if row in G[ch]:axpy(z,og[ot[k]]if kind=='typed'else og[k],G[ch][row])
   Ocols[old].append(z)
 assert set().union(*(set(z)for cols in Ocols.values()for z in cols))==opaque
 combos=[{k:Q(x)for k,x in row if k<640}for row in s['original619_c0_source_combinations']]
 a619={i:{k:Q(x)for k,x in row}for i,row in s['raw619_old4_plane_a1_nonempty28coefficients']}
 assert [[i,str(x)]for i,x in sorted(a619[609].items())]==[[16,'-500/29'],[17,'640/29']]
 rebuilt=[];outer605zero=[];affected=[];outer_unknown={};total_slots=0
 for rec in reversed(recipes):
  case=rec['free_case619'];cursor={'case':case};c=dict(combos[case]);aa=dict(a619.get(case,{}))
  for b,a in zip([563,609],rec['subtract_actual_basis563609']):
   if Q(a):axpy(c,combos[b],-Q(a));axpy(aa,a619.get(b,{}),-Q(a))
  channels=[];outerchannels=[];hasopaque=False
  for ch in range(2):
   z={}
   for old,x in sorted(c.items(),reverse=True):axpy(z,Dcols.get(old,[{},{}])[ch],x)
   oz={}
   for slot,x in sorted(aa.items(),reverse=True):
    pi,var=divmod(slot,9);assert pi<4 and var in[7,8];pair=8 if var==8 else 9
    for fn,a in Ocols[outer['original4_old_pivots'][pi]][ch].items():madd(oz,(pair,fn),mul(-x,a))
   axpy(z,oz,Q(1));expected={(p,fn):Q(x)for p,fn,x in f2[case][ch]};assert z==expected,('formal2',case,ch)
   total_slots+=len(z);channels.append(z);outerchannels.append(oz)
   assert {(p,fn):x for(p,fn),x in z.items()if fn in opaque}==oz
   if oz:hasopaque=True
   for p,fn in oz:outer_unknown.setdefault(fn,set()).add(p)
  if hasopaque:affected.append(case)
  else:outer605zero.append(case);assert not any(outerchannels)
  rebuilt.append((case,channels))
 assert total_slots==1909 and sorted(affected)==[579,580,593,594,597,598,605,606,611,613,615,616]==E['opaque_cases12']
 good=sorted(free-set(affected));assert len(good)==605 and good==E['cases605']and sorted(outer605zero)==good
 assert set(outer_unknown)==opaque and all(v=={8,9}for v in outer_unknown.values())
 # Apply the already independently proved full formal DIRECTOLD certificate
 # to the SAME catalog/recipes; inspect every surviving variable, including
 # genuine conductor thirds that a full-m3 theorem does not kill.
 physical={i for i,v in phys['all_requested164_root_values']};assert len(physical)==164 and not physical&(opaque|conductors)
 thirds=set();seconds=set();global3=set();global2=set();g3cases=set();unprovided=set();third_case_counts=0
 for case in sorted(free):
  thirdflag=False
  for ch in range(2):
   for fn,vs in f1[case][ch]:
    assert fn in boundaries and len(vs)==3
    for component,x in enumerate(vs):
     if not Q(x):continue
     assert fn in physical|conductors,('honest uncancelled J1',case,fn)
     (global3 if component==0 else global2).add(fn)
     if component==0:
      thirdflag=True
      if fn in conductors:g3cases.add(case)
      if case in good:assert fn in physical;thirds.add(fn)
   for pair,fn,x in f2[case][ch]:
    assert Q(x)and fn in physical|conductors|opaque
    if fn not in physical|conductors:unprovided.add((fn,pair))
    if case in good:
     assert fn not in opaque
     if fn in physical:seconds.add(fn)
  third_case_counts+=thirdflag
 assert len(thirds)==30 and sorted(thirds)==E['corrected605_directOLD_third_physical164_ids']
 assert len(seconds)==26 and sorted(seconds)==E['corrected605_combinedOLD2_physical164_ids']
 assert g3cases=={579,613,616}and sorted(g3cases)==E['conductor_third_cases_only_in_opaque12']and not g3cases&set(good)
 assert sorted(global3)==E['all617_directOLD_third_remaining_functions']and sorted(global2)==E['all617_directOLD_second_remaining_functions']
 assert unprovided=={(fn,p)for fn in opaque for p in[8,9]}and len(unprovided)==24
 assert E['all24_opaque_second_symbols_retained']==[list(x)for x in sorted(unprovided)]
 assert native['all617_literal_full48S0cN_zero_verified']and native['all617_full128_reconstructions_verified']
 nm={x['free_case619']:x for x in native['all617_source_maps']};assert set(nm)==free
 indices=sorted({k for case in good for k,x in nm[case]['kernel105_coordinates']})
 assert len(indices)==72 and sum(len(nm[case]['kernel105_coordinates'])for case in good)==120
 assert indices==E['native105_basis_indices_needed_for605']and E['native_coordinate_nonzeros_for605']==120
 # Both complete formal proof reports are pinned. Actual SECOND conductors
 # vanish only by the independent global physical-R0 m3 theorem; no third
 # value and no genuine honest higher jet was assigned zero here.
 assert E['gA_gB_global_m3_SECONDzero_source']in inputs
 assert E['full757_multiproducts_FIRSTcross_kept_not_pruned104']and E['corrected_source_basis563609_upper_NOT_set0']
 actual=count([grad,D,Dcols,og,G,ot,Ocols,combos,a619,rebuilt,outer605zero,affected,outer_unknown,nm,indices])+parsed+12000
 assert actual<=prior and all(sha(ROOT/p)==h for p,h in pins.items())
 data={'all_passed':True,'exact_typed629_subset_of757_and_literal532_shared96boundary_programs':True,'all67_typed_leaves_in_same96catalog':True,
  'all_typed_and_actual_boundary_origins0_for_formal_SECOND_linearization_only':True,'whole757_multiproducts_and_FIRSTcross_untouched':True,
  'independent_method':'Atom-first derivative maps and D column contractions, then literal source619 subtraction and four exact linear outer functions; reverse iteration order. No producer imports and no actual second/third coefficients.',
  'all617_combined_OLD2_formal_rows_rebuilt_exactly':True,'formal2_nonzero_slots':total_slots,
  'outer_exactly_linear_in12_actual_functions_not_just_first_order':True,'D_linear_opaque_support_empty':True,
  'all605_complete_corrected_outer_SECOND_contractions0':True,'raw609_plane_a1_kept':[[16,'-500/29'],[17,'640/29']],
  'cases605':good,'opaque_cases12':sorted(affected),'all24_opaque_symbols':E['all24_opaque_second_symbols_retained'],
  'DIRECTOLD_full_formal_certificate_applied_to_same_original_sources':True,'DIRECTOLD_third_physical_IDs30':sorted(thirds),'combinedOLD2_physical_IDs26':sorted(seconds),
  'genuine_gA_gB_third_cases_NOT_zero_filled':sorted(g3cases),'gA_gB_SECOND0_only_by_independent_global_full_m3':True,
  'all_genuine_honest_higher_formal_sensitivities0_with_conductor_exceptions_retained':True,
  'native72_basis120coordinates_inventory_only_with_full_kernel_source_proof_pinned':True,
  'future_scope':'Every FIRST cross, corrected basis563/609, fullPsi/comparison frame, individual native L2S0 and accepted native92 outer remain. Separate native72 proof/whole positive gate/worker/ROOT GO required. No actual605 responses, zero/rank/HF/length/H0 claim.',
  'whole_owned_numeric_upper':prior,'observed_held_numeric_upper':actual,'RAW_input_plus_output_upper':raw+4000000,'actual_known_QQ_committed_bits':bits,'actual_preoperation_bits':pre}
 ep=out/'independent-complete605-inventory-and-atom-first617-linear-outer-proof.json';save(ep,data)
 rep={'id':ID,'all_passed':True,'frozen':True,'SOURCE_ONLY_no_actual605_responses':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t}
 save(out/'report.json',rep);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':prior,'actual_owned':actual,'bits':bits,'pre':pre,'seconds':rep['seconds']}))
if __name__=='__main__':main()
