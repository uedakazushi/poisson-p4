#!/usr/bin/env python3
"""Independent corrected-span source/budget check; no signed quadratics."""
import itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_G3_quadratic_source_review_20261105_v63 import sha,sp,add,scale,mv
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_G3_PLANE_QUADRATIC_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
def main():
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 pp=ROOT/'reports'/A/'report.json';pr=json.loads(pp.read_text());assert pr['all_passed'];inputs={'primary':{'path':str(pp.relative_to(ROOT)),'sha256':sha(pp)}};D={}
 for p,h in pr['input_sha256'].items():
  path=ROOT/p;assert sha(path)==h;inputs[p]={'path':p,'sha256':h}
  if path.name!='actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json':D[path.name]=json.loads(path.read_text())
 ep=pp.parent/'actual-enlarged-G3quadratic-factored-source-and-known-NF-input.json';assert sha(ep)==pr['evidence_sha256'][ep.name]
 E=json.loads(ep.read_text());inputs['evidence']={'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)}
 assert sha(ROOT/'code/ambient_B18_ell2_T2_G3_plane_quadratic_source_20261105_v63.py')==pr['code_sha256']
 assert sha(ROOT/'research/notes'/f'{A}.md')==pr['note_sha256']
 for stem in ['AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_20261105_v63_REVIEW','AMBIENT_B18_ELL2_T2_G3_PLANE_KERNEL_INPUT_GATE_20261105_v63_REVIEW','AMBIENT_B18_ELL2_T2_G3_QUADRATIC_INPUT_GATE_20261105_v63_REVIEW']:
  path=ROOT/'reports'/stem/'report.json';assert json.loads(path.read_text())['all_passed'];inputs[stem]={'path':str(path.relative_to(ROOT)),'sha256':sha(path)}
 # Previously independently rebuilt P and NF are unchanged source maps.
 oldp=ROOT/'reports/AMBIENT_B18_ELL2_T2_G3_QUADRATIC_INPUT_GATE_20261105_v63/actual176-G3quadratic-factored-source-and-known-NF-input.json'
 oldreport=json.loads((oldp.parent/'report.json').read_text());assert sha(oldp)==oldreport['evidence_sha256'][oldp.name]
 old=json.loads(oldp.read_text());inputs['previous_verified_P_NF']={'path':str(oldp.relative_to(ROOT)),'sha256':sha(oldp)}
 for key in ['actual_selected126_rawD31_Pi42_linear_columns','PiY_L48_degree2_NF42_times20_columns','active_variable_global_indices','degree2_active_pair_order20']:
  assert E[key]==old[key]
 prefix=next(p for p in pr['input_sha256'] if p.endswith('actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'))
 F=fields(ROOT/prefix,['actual_active_I3_origin_columns','PiY3_from102_columns'])
 I3=list(map(sp,F['actual_active_I3_origin_columns']));PiY=list(map(sp,F['PiY3_from102_columns']))
 t=D['actual-active-unit-bases-and26-frame-transport.json'];pref=D['active-five-d3-source-split-and-budget.json']
 Pi42=list(map(sp,t['active_Pi42_columns']));fixed=list(map(sp,t['active_F2_fixed182_columns']));D30=list(map(sp,pref['active_d3_columns']))
 norms={};reps={}
 for receipt in t['active_F2_fixed182_receipts']:
  j=receipt['source_column'];p=receipt['pivot'];w={j:Q(1)}
  for k,c in receipt['subtract']:w=add(w,reps[k],-Q(c))
  w=scale(w,1/Q(receipt['divide_by']));v=sp(receipt['normalized']);assert mv(fixed,w)==v
  reps[p]=w;norms[p]=v
 active=[0,1,4,5,7];variables=active+[8];blocks=[(2-q,w) for q in range(3) for w in itertools.combinations(active,q)];index={w:i for i,(_,w) in enumerate(blocks)};raw=[]
 for q in range(4):
  for w in itertools.combinations(active,q):
   p=3-q
   for j in range(14):
    col={}
    if p:col[(14*index[w]+j)*6+5]=Q(1)
    for s,z in enumerate(w):col=add(col,{(14*index[w[:s]+w[s+1:]]+j)*6+variables.index(z):Q((-1)**(p+s))})
    raw.append(col)
 assert len(raw)==364
 cert=D['actual-G3-plane-condition-fixed-rank-and-kernel.json'];kernels=[sp(v['original238_source_combination']) for v in cert['full_constant_kernel_columns']]
 assert kernels==list(map(sp,E['source_kernel_enlarged'])) and len(kernels)==218
 Efree=[{j:Q(1)} for j in t['active_F2_free42']]
 assert all(mv(Pi42,col)=={j:Q(1)} for j,col in enumerate(Efree))
 corrections=[];retained_z=0
 for n,c in enumerate(kernels):
  w=mv(I3,c);assert w==sp(E['source_original_raw364_combinations'][n]) and not mv(D30,w)
  V=mv(raw,w);correction={};zfirst={}
  for var in range(6):
   q={i//6:a for i,a in V.items() if i%6==var};pi=mv(Pi42,q)
   assert not mv(PiY[:42],pi)
   if var in (4,5):assert not pi
   for i,a in pi.items():zfirst[6*i+var]=a
   corrected=add(q,mv(Efree,pi),Q(-1));work=dict(corrected);coeff={}
   assert not mv(Pi42,corrected)
   for p in sorted(norms):
    x=work.get(p,Q(0))
    if x:coeff=add(coeff,reps[p],x);work=add(work,norms[p],-x)
   assert not work and mv(fixed,coeff)==corrected
   for i,a in coeff.items():
    if i>=56:correction[6*(i-56)+var]=a
  assert zfirst==sp(E['actual_nonzero_first_z_linear_columns'][n]);retained_z+=bool(zfirst)
  assert correction==sp(E['actual_enlarged_selected126_correction_linear_columns'][n]);corrections.append(correction)
 P=list(map(sp,E['actual_selected126_rawD31_Pi42_linear_columns']));NF=list(map(sp,E['PiY_L48_degree2_NF42_times20_columns']));pairs=list(map(tuple,E['degree2_active_pair_order20']));pi={p:i for i,p in enumerate(pairs)}
 envelopes=[];termcount=0;highest=0
 for corr in corrections:
  bounds=[Q(0) for _ in range(6)];den=1
  for k,x in corr.items():
   source,u=divmod(k,6)
   for i,y in P[source].items():
    row,v=divmod(i,6);pair=tuple(sorted((u,v)))
    if pair==(5,5):continue
    for target,z in NF[20*row+pi[pair]].items():
     term=abs(x*y*z);den=math.lcm(den,term.denominator);bounds[target]+=term;termcount+=1
  top=max([0]+[q.numerator*(den//q.denominator) for q in bounds]);envelopes.append([str(den),str(top)]);highest=max(highest,den.bit_length(),top.bit_length())
 assert envelopes==E['positive_bound_per_enlarged_common_D_U'] and termcount==E['signed_product_summand_count']
 assert highest==75==E['positive_output_bound_bits']
 maxinput=max([0]+[max(abs(c.numerator).bit_length(),c.denominator.bit_length()) for v in corrections+P+NF for c in v.values()])
 assert maxinput==105==E['max_factored_input_bits'] and max(highest,3*maxinput+2)==317==E['product_and_output_operational_bound_bits']
 assert sum(E['next_slot_breakdown'].values())==24127==E['next_whole_slot_bound']<60000
 assert E['next_raw_byte_bound']==9974816<32*1024*1024 and E['source_phase_owned_slots']==24434<150000
 assert E['limits']['guard_bits']==2048 and not E['signed_actual_quadratic_values_evaluated']
 proof={
  'corrected_span':'K0=[E56,D30sel126,Efree42]. For arbitrary q, q-Efree42*Pi42q is in span(first182). K0^-1 Efree42 has first182 coordinates zero, so A126 extraction is unchanged. The retained z projection was checked, not discarded.',
  'actual_formula':'The exact P(t)/B3 argument previously peered gives G3_2 c=-Pi42 D31sel A126(D31I3c). All218 corrected extraction coefficients are independently verified.',
  'later_unit_derivatives':'Plane-kernel source peer verifies full r/ell firstzero and projected zfirstzero. Later quotient/bottom derivatives multiply retained zfirst; all378 z*m normal forms kill these quadratic terms.',
  'bound':'Only absolute term products and positive envelopes were evaluated. No signed quadratic vectors or new rank computed. Product/intermediate bound317 and output75 are below the separately requested2048 guard; value job still needs root GO.',
  'scope':'New nonzero image would be a true quadratic initial relation and an upper-cover input only. Allzero would stop only this sufficient family; no full HF2/total length follows.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);path=out/'independent-enlarged-corrected-span-and-positive-budget.json'
 path.write_text(json.dumps({'all218_corrected182_span_and_A126_identities':True,'retained_nonzero_z_first_columns':retained_z,'verified_P_NF_maps_unchanged':True,
                             'positive_summands':termcount,'max_positive_output_bits':75,'operational_product_bits':317,'numeric_slot_bound':24127,
                             'raw_byte_bound':9974816,'source_proof':proof},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'enlarged218_corrected_source_formula_and_bound_accepted':True,'signed_quadratic_values_or_rank_computed':False,
         'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
         'evidence_sha256':{path.name:sha(path)},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'retained_z':retained_z,'positive_summands':termcount}))
if __name__=='__main__':main()
