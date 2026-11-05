#!/usr/bin/env python3
"""Boolean support/lifetime and literal integer-height classes only."""
import argparse,json,math,resource,signal,time
from pathlib import Path
from fractions import Fraction as F
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER605_BOOLEAN_PROFILE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={};raw=0
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024;return json.loads((ROOT/p).read_text())
 prefix='reports/AMBIENT_B18_ELL2_T2_UPPER605_OLD_FIRST_CACHE_SOURCE_20261105_v63/'
 g=read(prefix+'report.json');assert g['all_passed'];s=read(prefix+'pruned605-OLD757-C0532-B21-FIRST162-and-native92-source-input.json')
 n=read('reports/AMBIENT_B18_NATIVE_SHARED72_INPUT_SOURCE_20261105_v63/native-shared75-formal-to72-actual-to605-source-input.json')
 fine=read('reports/AMBIENT_B18_NATIVE_SHARED72_FINE_POSITIVE_20261105_v63/native72-fine-Gamma-D-and605-constantmap-positive-envelopes.json')
 basis=[tuple(x) for x in s['ordinary_B21_basis']];bi={x:i for i,x in enumerate(basis)};assert len(basis)==21
 table=[[bi.get(tuple(x+y for x,y in zip(a,b))) for b in basis] for a in basis]
 def product(a,b):return {table[i][j] for i in a for j in b if table[i][j] is not None}
 leaves={i:{k for k,q in enumerate(v) if F(q)} for i,mode,v in s['actual_or_formal_B21_boundary96']};assert all(0 not in v for v in leaves.values())
 C=dict(leaves)
 for i,d in s['shared_literal_C0_program532']:
  if d[0]=='c':v={0} if d[1] else set()
  elif d[0]=='n':v=set(C[d[1]])
  elif d[0]=='a':v=C[d[1]]|C[d[2]]
  else:assert d[0]=='m';v=product(C[d[1]],C[d[2]])
  C[i]=v;assert 0 not in v
 C={i:C[i] for i,d in s['shared_literal_C0_program532']};T={}
 for i,d in s['shared_full_typed_program757']:
  if d[0]=='actual_C0_scalar_function':v=set(leaves[d[1]])
  elif d[0]=='sum':v=set().union(*(T[k] for k,a,b in d[1] if a)) if d[1] else set()
  else:
   assert d[0]=='product';v={0} if d[1] else set()
   for k in d[3]:v=product(v,T[k])
  T[i]=v;assert 0 not in v
 L0={j:{i for i,q in row if F(q)} for j,row in s['native_L0_selected_columns']};Le={j:{i for i,q in row if F(q)} for j,row in s['native_Le_selected_columns']};first={}
 for rec in s['individual_FIRST_column_recipes162']:
  d=rec['first_column_recipe']
  if d['kind']=='literal56':v={i for i,q in d['column'] if F(q)}
  else:
   assert d['kind']=='native_ell_SUM';v=set()
   for name,mat in [('Se',L0),('S0',Le)]:
    for k,q in d[name]:
     if F(q):v|=mat[k]
  first[rec['pivot27'],rec['variable9']]=v
 pairs=s['fixed_pair_order10'];D=[[{i for i,q in enumerate(row) if F(q)} for row in ch] for ch in s['D_two_by10_by56']]
 first_terms=[]
 for t in range(2):
  for p in range(27):
   for var in [0,1,4,5,7,8]:
    z=0
    for pair,(a,b) in enumerate(pairs):
     if var==a:z+=len(D[t][pair]&first[p,b])
     if var==b and a!=b:z+=len(D[t][pair]&first[p,a])
    first_terms.append([t,p,var,z])
 # Integer-only LCM/cleared absolute numerator classes of LITERAL operands.
 # No signed polynomial or response contraction is performed here.
 def nums(v):
  if isinstance(v,(list,tuple)):
   for x in v:yield from nums(x)
  else:yield F(v)
 operands={
  'boundary96_21':[v for fn,mode,v in s['actual_or_formal_B21_boundary96']],
  'GammaD':[q for row in s['Gamma_two'] for i,q in row]+[q for ch in s['D_two_by10_by56'] for row in ch for q in row],
  'typed_scalar_coefficients':[F(d[1],d[2]) for i,d in s['shared_full_typed_program757'] if d[0]=='product']+[F(a,b) for i,d in s['shared_full_typed_program757'] if d[0]=='sum' for k,a,b in d[1]],
  'raw619_constant_sources':[q for row in s['raw619_original768_constant_sources'] for i,q in row],
  'raw619_a1':[q for row in s['raw619_actual_a1_27times9'] for i,q in row],
  'correct605_weights':[q for i,row in s['corrected605_recipes_basis563609'] for q in row],
  'FIRST_literal_columns':[q for r in s['individual_FIRST_column_recipes162'] if r['first_column_recipe']['kind']=='literal56' for i,q in r['first_column_recipe']['column']],
  'FIRST_native_ell_SeS0':[q for r in s['individual_FIRST_column_recipes162'] if r['first_column_recipe']['kind']=='native_ell_SUM' for name in ['Se','S0'] for i,q in r['first_column_recipe'][name]],
  'FIRST_L0_Le':[q for name in ['native_L0_selected_columns','native_Le_selected_columns'] for j,row in s[name] for i,q in row],
  'accepted_native92':s['accepted_native92_Gamma_plane_targets23_twochannels_rr_rell']}
 classes={};inputbits=0
 for name,v in operands.items():
  den=1;xs=list(nums(v))
  for q in xs:den=math.lcm(den,q.denominator);assert den.bit_length()<=384
  U=max([den]+[abs(q.numerator)*(den//q.denominator) for q in xs]);assert U.bit_length()<=384
  classes[name]={'D':str(den),'U_including_one':str(U),'count':len(xs),'D_bits':den.bit_length(),'U_bits':U.bit_length()};inputbits=max(inputbits,den.bit_length(),U.bit_length())
 cs=sum(map(len,C.values()));ts=sum(map(len,T.values()));fs=sum(map(len,first.values()));ls=sum(map(len,leaves.values()))
 stage_source=count([g,s,n,fine,C,T,first,leaves,L0,Le,D,table,operands,classes,first_terms,pins])+15000;assert stage_source<150000
 data={'no_QQ_jet_convolution_or_signed_response_or_native_functional_values':True,'shared_input_numeric_counts':{'OLD_FIRST_cache':count(s),'native75_72_cache':count(n),'native_FINE':count(fine)},'actual_support_counts':{'boundary96_scalar_cells':ls,'C0_532_scalar_cells':cs,'typed757_scalar_cells':ts,'FIRST162_scalar_cells':fs,'FIRST324_scalar_dot_terms':sum(r[3] for r in first_terms)},'all324_FIRST_table_operational_term_counts':first_terms,'Boolean_program_fingerprints':{'C0':H.sha_text(json.dumps([[i,sorted(v)] for i,v in sorted(C.items())])) if hasattr(H,'sha_text') else __import__('hashlib').sha256(json.dumps([[i,sorted(v)] for i,v in sorted(C.items())]).encode()).hexdigest(),'typed':__import__('hashlib').sha256(json.dumps([[i,sorted(v)] for i,v in sorted(T.items())]).encode()).hexdigest()},'literal_scalar_operand_classes':classes,'only_input_LCM_integer_class_height':inputbits,'future_lifetime_recipe':'Hash/bind input union, keep immutable source once. ComputeC0jets then compileC0Gamma-table and deleteC0jets BEFOREtyped757. Deleteall OLD polynomial jets before native75stream. Only72kernel restrictions and605component/class arrayspersist. FIRST162 may be discarded after324table. Guardcodec serializes onlycurrentphaseNEWstates plus completedtables/maps, immutable input byhash.','future_positive_OLD_workspace_DU_numeric_upper':3*max(cs,ts)+3*ls+3*fs+10000,'future_positive_OLD_temp_integer_proof':'EverypositiveDU mul/add/lcm/resultmustcommitguard384beforeuse; individualbinarytemporaries<=769<1153. STOPonlyifheightfails; noneofthisprofilepromisesOLDwholeheight.','native354_is_intermediate_only_USEfine254_BEFORE_Q':True,'whole_input_pruning_and_dedup_and_positive_assembly_and_codec_gate_pending':True}
 ep=out/'Boolean-shared-program-supports-FIRST324-and-literal-LCM-classes.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'SOURCE_Boolean_and_literal_integer_only':True,'owned_source_numeric_upper':stage_source,'literal_operand_class_max_bits':inputbits,'RAW_upper':raw,'no_605_upper_or_component_or_native_response_values':True,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','owned_source_numeric_upper','literal_operand_class_max_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
