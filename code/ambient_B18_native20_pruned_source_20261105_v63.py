#!/usr/bin/env python3
"""Exact native20 SOURCE prune; no native functional or compiler responses."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NATIVE20_PRUNED_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;pointbits=0
 def read(path):
  nonlocal raw
  p=ROOT/path;pins[path]=sha(p);raw+=2*p.stat().st_size;assert raw<32*1024*1024;return json.loads(p.read_text())
 def dump(name,v):(out/name).write_text(json.dumps(v,separators=(',',':'))+'\n')
 def ck(q):
  nonlocal pointbits
  pointbits=max(pointbits,abs(q.numerator).bit_length(),q.denominator.bit_length());assert pointbits<=32;return q
 sp='reports/AMBIENT_B18_NATIVE_SHARED72_INPUT_SOURCE_20261105_v63/'
 sr=read(sp+'report.json');assert sr['all_passed'] and pins[sp+'report.json']=='4f01895013eb59d8d1d5721eff04d8131e9a73f52f2449f9613731e311c704e6'
 se='native-shared75-formal-to72-actual-to605-source-input.json';s=read(sp+se);assert pins[sp+se]==sr['evidence_sha256'][se]
 fp='reports/AMBIENT_B18_NATIVE_SHARED72_FINE_POSITIVE_20261105_v63/'
 fr=read(fp+'report.json');assert fr['all_passed'] and pins[fp+'report.json']=='072c20b07c37929a97f49da4b0f8235f846d30a4d05f9ebada0d5c858edffa76'
 fe='native72-fine-Gamma-D-and605-constantmap-positive-envelopes.json';f=read(fp+fe);assert pins[fp+fe]==fr['evidence_sha256'][fe]
 prpath='reports/AMBIENT_B18_NATIVE_SHARED72_SOURCE_FINE_20261105_v63_REVIEW/report.json';peer=read(prpath);assert peer['all_passed'] and pins[prpath]=='0384d2623fab486eb0dbd1a4df239bd04e6a7ff778f7464c454ec324c541683e'
 op='reports/AMBIENT_B18_ELL2_T2_UPPER20_PRUNED_SOURCE_20261105_v63/'
 orr=read(op+'report.json');assert orr['all_passed'] and pins[op+'report.json']=='34e180e5161344f1aee520c86c53d2691e05a19827b3c2533f8ecc71e2d70702'
 oe='pruned20-exceptional-OLD-FIRST-outer-source-and605-coverage.json';o=read(op+oe);assert pins[op+oe]==orr['evidence_sha256'][oe]
 cases=o['exceptional_cases20'];assert cases==[525,529,530,550,552,558,560,562,566,567,569,573,576,577,589,595,599,603,612,617]
 recipes=dict(o['corrected20_source_recipes']);assert recipes[617]==['-20','0'] and all(v==['0','0'] for i,v in recipes.items() if i!=617)
 maps={r['case619']:r for r in s['all605_constant_native_coordinate_recipes']};fine={r['kernel105']:r for r in f['native_kernel72_Gamma_D_sum_fine']};restr={r['kernel105']:r for r in s['kernel72_restrictions']}
 used=sorted({i for case in cases for i,q in maps[case]['kernel105_coordinate_pairs']});assert len(used)==25
 zero=[i for i in used if all(int(ch[2][1])==0 for ch in fine[i]['two_channels_Gamma_D_sum_DU'])];active=sorted(set(used)-set(zero));assert len(zero)==10 and len(active)==15
 # Exact zero SOURCE certificates from nonnegative envelopes, not uncomputed
 # response defaults. Sumzero implies BOTH Gamma andD envelopezero aswell.
 certificates=[]
 for i in zero:
  ch=fine[i]['two_channels_Gamma_D_sum_DU'];assert all(q==['1','0'] for row in ch for q in row)
  certificates.append({'kernel105':i,'two_channels_Gamma_D_sum_DU':ch,'actual_two_channels_sum_for_fixed_kernel':['0','0'],
   'proof':'all six nonnegative D/U bounds have numerator0; sourcepeer0384 verifies exact fullS0 restriction before pointoverride.',
   'original_native_source_combination':restr[i]['original_native_source_combination']})
 selectedrestr=[restr[i] for i in active];rawids=sorted({j for rec in selectedrestr for j,q in rec['original_native_source_combination']});assert rawids==[3,5,8,19,23,29,32,43,47,50,62,65,66,80,103,110]
 program={r['native128']:r for r in s['raw75_original_Koszul_program']};ops=[program[j] for j in rawids];assert len(ops)==16
 p=s['native_primitive_known_cache'];G0=[[ck(F(q)) for q in row[:4]] for row in p['G8_C72']];nu0=[ck(F(row[0])) for row in p['nu_two_J18']]
 assert nu0==[F(0),F(0)] and all(row[:3]==[F(0)]*3 for row in G0)
 before=count([sr,s,fr,f,peer,orr,o,maps,fine,restr,used,zero,active,certificates,selectedrestr,rawids,program,ops,p,G0,nu0,pins])+15000;assert before<150000
 dump('phase-before-constant-point-restriction-checks.json',{'source_point_only_no_response':True,'known_point_prior32_temp97':True,'inherited_higher_intermediate_prior354_temp1063_UNCHANGED':True,'source_owned_before_point':before})
 points={}
 for op in ops:
  v={};im=op['inputmask']
  for row,gen,sign in op['Koszul_entries']:
   for mask,q in enumerate(G0[gen]):
    if q and not(mask&im):
     col=4*row+(mask^im);z=ck(v.get(col,F(0))+ck(F(sign)*q))
     if z:v[col]=z
     else:v.pop(col,None)
  points[op['native128']]=v
 assert {j for j,v in points.items() if v}=={8,32,80}
 for rec in selectedrestr:
  image={}
  for j,q in rec['original_native_source_combination']:
   for row,a in points[j].items():
    z=ck(image.get(row,F(0))+ck(F(q)*a))
    if z:image[row]=z
    else:image.pop(row,None)
  assert not image
 prunedmaps=[]
 for case in cases:
  full=maps[case]['kernel105_coordinate_pairs'];live=[[i,q] for i,q in full if i in active];z=[[i,q] for i,q in full if i in zero]
  assert sorted(live+z)==full
  prunedmaps.append({'case619':case,'source_subtract563609':recipes[case],
   'full_kernel105_coordinate_pairs':full,'kernel105_coordinate_pairs':live,'explicit_proven_zero_kernel_terms':z})
 assert next(r for r in prunedmaps if r['case619']==617)['full_kernel105_coordinate_pairs']==[[2,'1200'],[3,'4800'],[37,'20'],[82,'1']]
 cfr='reports/AMBIENT_B18_NATIVE_KERNEL617_TO105_SOURCE_20261105_v63/all617-native-constant-sources-factor-through-fixed105-kernel.json';cm=read(cfr);cmap={r['free_case619']:r for r in cm['all617_source_maps']}
 for case in cases:assert list(map(str,cmap[case]['subtract563609']))==recipes[case]
 envmap=dict(f['all605_native_sum_fine_twochannels'])
 data={'SOURCE_metadata_and_constant_point_only_no_functional_responses':True,'exceptional_cases20':cases,
  'native_primitive_known_cache':p,'raw_auxiliary_program':ops,'raw_auxiliary_source_indices':rawids,
  'actual_kernel_restrictions':selectedrestr,'actual_kernel_original105_indices':active,
  'actual_kernel_fine_Gamma_D_sum_DU':[fine[i] for i in active],
  'explicit_zero_kernel_certificates':certificates,'referenced_kernel_indices_before_zero_prune':used,
  'all20_constant_native_coordinate_recipes':prunedmaps,
  'all20_native_sum_fine_twochannels':[[j,envmap[j]] for j in cases],
  'raw_point_S0_columns48':[[j,[[i,str(q)] for i,q in sorted(v.items())]] for j,v in sorted(points.items())],
  'actual_point_G8_eight_by4':[[str(q) for q in row] for row in G0],
  'all15_retained_kernel_full48_point_restrictions_zero':True,
  'pointnonzero_raw_indices_retained':[8,32,80],
  'raw4_and12_not_needed_after_prune':'Omitted by exact source dependency/zero-kernel certificates; NOT declared actualpointzero.',
  'formal_auxiliary_rule':'Allraw16 partialN Gamma/DcoefficientsremainFORMAL untilactual15fullpointkernelrestrictions; raw8/32/80haveNONZEROfullS0.',
  'native_maps_ALREADY_corrected_no_second_lower_subtraction':True,
  'case617_rule':'Corrected native coefficients ALREADYinclude+20*source563. OLD/FIRST raw21 lane applies its own source617+20source563. DoNOTadd563again to native map; doNOTsetraw563upper0.',
  'zero_kernel_rule':'Eachomittedkernel has an explicit SOURCE nonnegative bound0 in bothchannels Gamma/D/sum. No silent .get(...,0) for anyotheractualkernel.',
  'inherit_basis_coarse341_temp1024_native_map354_temp1063':True,
  'fine_bounds_checked_BEFORE_native20_mapping_OLD_FIRST_outer_Q':True,
  'no_new_physical_gauge_honest_or_native_response':True,
  'whole20_immutableworker_custody_freshpeer_and_ROOTGO_REQUIRED':True}
 ep='pruned-native20-formal16-toactual15-and-explicit10zero-kernel-source.json';dump(ep,data)
 owned=count([sr,s,fr,f,peer,orr,o,maps,fine,restr,used,zero,active,certificates,selectedrestr,rawids,program,ops,p,G0,nu0,points,prunedmaps,cm,cmap,envmap,data,pins])+10000;assert owned<150000
 rep={'id':STEM,'all_passed':True,'SOURCE_pruning_and_known_point_only':True,'all_native_or20_response_values_UNRUN':True,
  'referenced_kernels25':len(used),'proven_zero_kernels10':len(zero),'retained_actual_kernels15':len(active),'formal_raw_auxiliaries16':len(rawids),
  'native20_recipe_coordinate_records':sum(len(r['kernel105_coordinate_pairs']) for r in prunedmaps),'input_numeric_count':count(data),
  'source_owned_upper':owned,'source_RAW_upper':raw+sum(p.stat().st_size for p in out.iterdir()),'constant_point_prior32':True,'actual_constant_point_bits':pointbits,
  'higher_native_prior354_temp1063_unchanged':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep),'phase-before-constant-point-restriction-checks.json':sha(out/'phase-before-constant-point-restriction-checks.json')},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'input_numeric':count(data),'source_owned':owned,'native15_formal16':True}))
if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--output',required=True,type=Path);a=pa.parse_args();main(a.output)
