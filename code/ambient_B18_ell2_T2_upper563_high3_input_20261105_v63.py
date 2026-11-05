#!/usr/bin/env python3
"""Positive source gate for direct OLD highest-third sensitivities only."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_UPPER563_HIGH3_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={};maximum=0
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 sup='reports/AMBIENT_B18_ELL2_T2_UPPER563_SUPPORT_PREFLIGHT_20261105_v63/';src=read(sup+'one563-coefficient-free-cubic-second-and-unit-carry-ancestor-inventory.json');assert src['case619']==563
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';B=[dict(v) for v in field(mug,'old640_source_columns')];ops=dict(field(mug,'old640_typed_closure_nodes'));P=field(H.PATHS['C0'],'honest_D0_relation_columns');gamma=src['Gamma_two'];old=[(j,c) for j,c in src['original768_constant_combination'] if j<640]
 terms=[];needed=set();opaque=set()
 for ch,row in enumerate(gamma):
  for j,c in old:
   for i,g in row:
    if i<48:
     k=B[j].get(i)
     if k is None or ops[k][0]=='c':continue
     kind='typed';needed.add(k)
    else:kind='literal_C0';k=P[j][i-48+1];opaque.add(k)
    terms.append({'channel':ch,'original_source':j,'source_row56':i,'root_kind':kind,'root':k,'coefficient_factors':[g,c]})
 closure=set()
 def visit(k):
  if k in closure:return
  closure.add(k);n=ops[k]
  for j in ([j for j,a,b in n[1]] if n[0]=='sum' else n[3] if n[0]=='product' else []):visit(j)
 for k in needed:visit(k)
 program=[[i,ops[i]] for i in sorted(closure)];atoms=sorted({n[1] for i,n in program if n[0]=='actual_C0_scalar_function'});allfirst=dict(field('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json','first640')['all79_atom_values_full10']);origins=[]
 for i in atoms:
  assert i in allfirst and not H.parse(allfirst[i]).get(0,F(0));origins.append([i,'0'])
 # The two directP8 origins are explicit cells of the ORIGINAL source array.
 X0=field(H.PATHS['origin'],'full_X0_56_by768');assert len(X0)==56
 literalproof=[]
 for t in terms:
  if t['root_kind']=='literal_C0':assert not F(X0[t['source_row56']][t['original_source']]);literalproof.append([t['root'],t['source_row56'],t['original_source'],'0'])
 del X0,B,P,allfirst
 zero=(1,0)
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());return a
 def ex(x):x=F(x);return mark((x.denominator,abs(x.numerator))) if x else zero
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=math.lcm(a[0],b[0]);return mark((d,a[1]*(d//a[0])+b[1]*(d//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else zero
 def add(a,b,c):
  z=dict(a)
  for k,x in b.items():z[k]=plus(z.get(k,zero),mul(c,x))
  return z
 dual={};centers={}
 for i,n in program:
  if n[0]=='actual_C0_scalar_function':center=zero;v={n[1]:ex(1)}
  elif n[0]=='c':center=ex(F(n[1],n[2]));v={}
  elif n[0]=='sum':
   center=zero;v={}
   for j,a,b in n[1]:c=ex(F(a,b));center=plus(center,mul(c,centers[j]));v=add(v,dual[j],c)
  else:
   assert n[0]=='product';center=ex(F(n[1],n[2]));v={}
   for j in n[3]:v=add({k:mul(x,centers[j]) for k,x in v.items()},dual[j],center);center=mul(center,centers[j])
  centers[i]=center;dual[i]=v
 targets=[{},{}]
 for t in terms:
  w=mul(ex(t['coefficient_factors'][0]),ex(t['coefficient_factors'][1]));v=dual[t['root']] if t['root_kind']=='typed' else {t['root']:ex(1)};targets[t['channel']]=add(targets[t['channel']],v,w)
 assert maximum<384
 data={'case619':563,'formal_scope':'ONLY directOLD GammaXc3 highestthirdsymbollinearresponse;allotheruppertermsretained','original_factor_terms':terms,'typed_source_program':program,'actual_scalar_leaf_origins':origins,'literalP8_point_origin_source_cells':literalproof,'positive_node_origin_envelopes':[[i,list(map(str,q))] for i,q in centers.items()],'positive_node_high3_sensitivity_envelopes':[[i,[[k,list(map(str,q))] for k,q in sorted(v.items())]] for i,v in dual.items()],'positive_two_weighted_high3_envelopes':[[[k,list(map(str,q))] for k,q in sorted(v.items())] for v in targets],
 'all_scalar_third_coefficients_remain_independent':True,'no_signed_gradient_or_actual_scalar_jet_value_computed':True,'positive_prior_committed_bits':maximum,'temporary_bound':3*maximum+1,'future_formal_only_scope':'Origin-dual rootfold only. No oldknownpart/native/b1/A2/M1/coherentcopy response values.',
 'no_true_higher_source_order_necessity_claim':True}
 inputcount=count(data);future=2*(inputcount+600)+6*len(program)+8*len(atoms)+4000;raw=1000000+future*(2*math.ceil(maximum*math.log10(2))+40)
 assert future<60000 and raw<32*1024*1024
 data.update({'future_whole_numeric_bound':future,'future_RAW_bound':raw,'metadata_reserve':600})
 ep=out/'direct-old-highestthird-origin-dual-source-and-positive-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'positive_source_gate_only':True,'typed_nodes':len(program),'typed_scalar_atoms':len(atoms),'literalP8_scalar_atoms':len(opaque),'positive_prior_bits':maximum,'temporary_bits':3*maximum+1,'future_numeric_bound':future,'future_RAW_bound':raw,'source_owned_bound':count([src,ops,program,dual,centers,data])+70000,'input_sha256':hs,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start,'formal_fold_or_upper_values_executed':False};assert report['source_owned_bound']<150000;H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','typed_nodes','typed_scalar_atoms','positive_prior_bits','future_numeric_bound','future_RAW_bound']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
