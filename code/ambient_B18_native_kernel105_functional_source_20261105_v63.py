#!/usr/bin/env python3
"""Native105 factored linear SOURCE and ownership only; no response values."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NATIVE_KERNEL105_FUNCTIONAL_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def values(v):
 if isinstance(v,(list,tuple)):
  for x in v:yield from values(x)
 else:yield F(v)
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;bits=0
 def read(p):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size;assert raw<32*1024*1024;return json.loads(q.read_text())
 def ck(n):
  nonlocal bits
  bits=max(bits,abs(n).bit_length());assert bits<=384 and 3*bits+1<=1153;return n
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 kpre='reports/AMBIENT_B18_NATIVE_POINT_KERNEL105_SOURCE_20261105_v63/'
 kr=read(kpre+'report.json');assert kr['all_passed'] and pins[kpre+'report.json']=='15b6ff370efd02939a1fac290ca6c24c1237ab3046f6a82984d22b1071e8d534'
 ke='constant-kernel105-basis-unit23-minor-and-native-functional-determinacy.json';k=read(kpre+ke);assert pins[kpre+ke]==kr['evidence_sha256'][ke]
 ip='reports/AMBIENT_B18_NATIVE_ONE554_POSITIVE_INPUT_20261105_v63/ONE554-positive-four-native-and-nonzero644-outer-input.json';d=read(ip);p=d['full_primitive_input_J18_N_G_nu_GammaD']
 framep='reports/AMBIENT_B18_NATIVE_ONE554_FRAME_SOURCE_20261105_v63/report.json';assert read(framep)['all_passed']
 basis=[tuple(m) for m in p['basis18']];bi={m:i for i,m in enumerate(basis)};table=[[bi.get(tuple(x+y for x,y in zip(a,b))) for b in basis] for a in basis]
 assert len(basis)==18 and p['N_known_positions']==[0,1,2,3,4,5,6,15,16]
 nus=[[F(q) for q in row] for row in p['nu_two_J18']]
 assert all(not row[i] for row in nus for i in [0,5,6])
 assert all(table[i][j] is None for a in nus for b in nus for i,x in enumerate(a) if x for j,y in enumerate(b) if y)
 cols=k['explicit105_kernel_basis_columns'];assert len(cols)==105
 free=[x['free_native_source'] for x in cols];assert len(set(free))==105
 pairs=[(i,j) for i in range(8) for j in range(i+1,8)]
 intrinsic=[]
 for src in range(112):
  j,mask=divmod(src,4);a,b=pairs[j];intrinsic.append({'native128':src,'inputmask':mask,'Koszul_entries':[[a,b,1],[b,a,-1]]})
 used=sorted({src for col in cols for src,q in col['column']});assert set(used)<=set(range(112))
 assert all(max(len(x['column']),1)<=3 for x in cols)
 # Independent symbolic linearity map; AUX values outside kerS0 are FORMAL
 # partial-N extensions. Only each105 restricted target is an actual jet.
 restrictions=[]
 for i,col in enumerate(cols):
  restrictions.append({'kernel105':i,'original_native_source_combination':col['column'],
   'two_actual_outputs':[['sum_auxiliary_coefficients',t,col['column']] for t in range(2)]})
 A=[q for row in p['Gamma_two56'] for i,q in row]+[q for ch in p['D_two_by10_by56'] for row in ch for q in row]
 classes={'known_partial_N9':p['N_known12_by12_by9_by4'],'canonical_G8':p['G8_C72'],'monics':p['nu_two_J18'],
  'Gamma_and_D_scalar_coefficients':A,'native_kernel105_weights':[q for col in cols for src,q in col['column']]}
 h={}
 for name,v in classes.items():
  xs=list(values(v));den=1
  for q in xs:den=ck(math.lcm(den,q.denominator))
  u=max([den]+[ck(abs(q.numerator)*(den//q.denominator)) for q in xs]);h[name]={'D':str(den),'U_with_one':str(u),'operand_count':len(xs)}
 # Any one requested target has <=3 source terms and two Koszul entries.
 # Twelve physical outputs,16 G/N mask pairs,9 N bases,18 G bases,
 # <=19 optional-monic choices and11 Gamma/D selector choices cover ALL
 # intermediate subsets. Double nu histories are removed by support first.
 paths=3*2*12*16*9*18*19*11;assert paths<2**28
 den=1;u=2**28
 for v in h.values():den=ck(den*int(v['D']));u=ck(u*int(v['U_with_one']))
 prior=max(den.bit_length(),u.bit_length());assert prior<=384 and 3*prior+1<=1153
 pruned={'basis18':p['basis18'],'N_known_positions':p['N_known_positions'],'N_known12_by12_by9_by4':p['N_known12_by12_by9_by4'],
  'G8_C72':p['G8_C72'],'nu_two_J18':p['nu_two_J18'],'Gamma_two56':p['Gamma_two56'],'D_two_by10_by56':p['D_two_by10_by56']}
 # Whole future ownership is explicit, but is not an execution/runtime gate.
 n=count([pruned,cols,intrinsic,restrictions,table,h]);workspace={
  'parsed_N_first8_by12_by9_by4':8*12*9*4,
  'parsed_G576_nu36_Gamma18_D1120':576+36+18+1120,
  'streamed_weightedS_atmost3rows18base4mask_index_value':3*18*4*2,
  'streamed_partialNSc_atmost12rows18base4mask_index_value':12*18*4*2,
  'optional_all112_auxiliary_coefficients_with_indices':112*2*2,
  'all210_restricted_targets_index_channel_value':210*3,
  'kernel_basis_and_recipe_parse':count([cols,intrinsic,restrictions]),
  'table_and_monomial_route_metadata':count(table)+1000,
  'guard_scalar_wrappers_report_serialization':8000}
 live=2*(n+1000)+sum(workspace.values());assert live<60000 and live<150000
 rawfuture=2_000_000+live*(2*math.ceil(prior*math.log10(2))+42);assert rawfuture<32*1024*1024
 ep='factored-two-by105-native-kernel-functional-source-and-prior.json'
 ev={'source_only_no_functional_values':True,'same_J18_with_r2ell_primitive_input':pruned,'intrinsic112_original_Koszul_columns':intrinsic,
  'kernel105_restrictions':restrictions,'fixed23_complement':k['fixed23_source_complement'],'fixed23_minor':'-1/128',
  'linearity_theorem':'Finite extraction, Res4 and N^t are R0-linear in constant source coefficients. Restrict each partial-N auxiliary row to the fixed105 kernel columns; unknown L3*S0 and L2_nonplane contributions vanish only after this restriction/pureplane extraction.',
  'functional_formula':'Gamma[L0S3_r2ell+LrS2_rell+LeS2_rr+LrrSe+LreSr] +sum_nonzz Dpair[L0S2pair+unique Lfirst*Sfirst]; ordinary rr once.',
  'all_point_kernels_bound_to_full_Res4S0_not_M0cancellation':True,
  'auxiliary112_values_FORMAL_outside_kerS0_NOT_actual_unrequested_nativeX':True,
  'unknown_N2mixed_N3_NOT_set_to_true_zero':True,
  'all_varying_nu_tracemask_carries_retained':True,
  'no_honest79_value_or_OLD_FIRST_outer_source_carry_supplied':True,
  'no_605_or617_HF_length_annihilation_inference':True,
  'positive_operand_classes':h,'common_D':str(den),'common_U_with_paths2pow28':str(u),'direct_paths_per_target_upper':paths,
  'committed_prior_bits':prior,'temporary_prior_bits':3*prior+1,'zero_double_nu_support_checked_before_products':True,
  'future_live_formula':workspace,'source_numeric_input_count':n,'future_numeric_forecast':live,'future_RAW_forecast':rawfuture,
  'future_immutable_worker_custody_runtime_guard_and_fresh_parent_GO_REQUIRED':True}
 dump(ep,ev);owned=count([kr,k,d,p,basis,bi,table,nus,cols,intrinsic,used,restrictions,classes,h,pruned,workspace,ev,pins])+10000
 assert owned<150000 and raw+(out/ep).stat().st_size<32*1024*1024
 rep={'id':STEM,'all_passed':True,'SOURCE_linearity_ownership_positive_only':True,'no_response_or_value_worker':True,
  'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'source_numeric_owned_upper':owned,'source_RAW_upper':raw+(out/ep).stat().st_size,
  'committed_prior':prior,'temporary_prior':3*prior+1,'future_forecast_only_numeric':live,'future_forecast_only_RAW':rawfuture,
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'prior':prior,'temp':3*prior+1,'owned':owned,'forecast':live}))
if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--output',required=True,type=Path);a=pa.parse_args();main(a.output)
