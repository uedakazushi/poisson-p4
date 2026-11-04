#!/usr/bin/env python3
"""Typed pair-only cofactor pure row and source-only11 column bounds."""
import argparse,hashlib,importlib.util,itertools,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_PAIR_ONLY_PURE_COVECTOR_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='certificate.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
 old=load('S61_NONSQUARE_PREPOWER_NORMALIZED_RAW8_ENTRY_20261104_v63','report.json');assert old['all_passed']and old['normalized_minor_is_zero']
 bottom=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6']
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63')
 path=ROOT/'code/s61_nonsquare_prePOWER_h22_kernel9_support_preflight_20261104_v63.py';pins[str(path.relative_to(ROOT))]=sha(path);sp=importlib.util.spec_from_file_location('bounds',path);b=importlib.util.module_from_spec(sp);sp.loader.exec_module(b)
 V=[]
 for n in q['scalar_DAG']['nodes']:
  op=n[0]
  if op=='z':x=None
  elif op=='one':x=b.ONE
  elif op in ['m','a','s']:x=b.mm(V[n[1]],V[n[2]])if op=='m'else None if op=='s'and n[1]==n[2]else b.aa(V[n[1]],V[n[2]])
  else:
   num,den=cp['matrix25x32'][n[1]][n[2]]if op=='checkpoint'else pd['Lframe_dual6x2'][n[1]][1]if op=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][n[1]]
   x=b.atom(num)
   if len(den)==1:assert den[0][0][1:]==[0,0];x=b.shift(x,-den[0][0][0])
   else:
    dd={tuple(m):c for m,c in den};assert set(dd)=={tuple(m)for m,c in b.BP};x=None if x is None else(x[0],x[1],1)
  V.append(x)
 rows=[2,4,5,8];compact=[1,3,4,7];cols=q['new22_target_original32_columns'];remaining=[1,2,3,4,6,7,8,9,14,17,18]
 A=[[b.atom(hc['raw_compact14x5'][r][j])for j in[0,4]]+[b.atom(bottom[i][0])]for i,r in enumerate(compact)]
 # Multiplication by known nonzeroQQ row constants changes no support/degree.
 nodes=[]
 def node(op,args=(),**kw):i=len(nodes);nodes.append(dict(id=i,op=op,args=list(args),**kw));return i
 atomroots=[]
 for i,r in enumerate(compact):
  H=[node('accepted_normalized_heavy_atom',compact_row=r,column=j,QQ_multiplier=hc['known_row_units_rho14'][r]['constant'])for j in[0,4]]
  H.append(node('accepted_bottom24_atom',row=i,column=0));atomroots.append(H)
 cofs=[];cofroots=[]
 for i in range(4):
  rr=[r for r in range(4)if r!=i];bound=None;terms=[]
  for p in itertools.permutations(range(3)):
   z=b.ONE
   for j in range(3):z=b.mm(z,A[rr[j]][p[j]])
   bound=b.aa(bound,z);parity=i+sum(p[j]>p[k]for j in range(3)for k in range(j+1,3));terms.append(dict(sign=-1 if parity%2 else 1,factors=[atomroots[rr[j]][p[j]]for j in range(3)]))
  cofs.append(bound);cofroots.append(node('signed_three_by_three_cofactor',sign=(-1)**i,terms=terms,bound=b.enc(bound)))
 profiles=[]
 for lab in remaining:
  vals=[V[q['remaining15x22_scalar_roots'][r][cols.index(lab)]]for r in rows];power=max(0 if x is None else x[2]for x in vals);total=None;terms=[]
  for i,(c,v)in enumerate(zip(cofs,vals)):
   if v is None:vv=None
   else:
    vv=b.numB(v,power-v[2]);vv=(vv[0],vv[1],0)
   z=b.mm(c,vv);total=b.aa(total,z);terms.append(dict(cofactor=cofroots[i],source_DAG_root=q['remaining15x22_scalar_roots'][rows[i]][cols.index(lab)],B_clear_to_power=power,source_bound=b.enc(vv),product_support_upper=b.size(z)))
  root=node('sum_four_source_products',terms=terms,bound=b.enc(total),known_B_column_clearing_power=power)
  profiles.append(dict(original32_target=lab,root=root,clearing_B_power=power,output_bound=b.enc(total),sum_source_entry_support=sum(b.size(v)for v in vals),largest_product_support_upper=max(t['product_support_upper']for t in terms),scalar_value_not_evaluated=True))
 profiles.sort(key=lambda x:(x['output_bound']['support_upper'],x['sum_source_entry_support'],x['original32_target']))
 recipe=dict(nodes=nodes,source_current15_rows=rows,source_compact14_rows=compact,anchor_original32=22,
  normalized_heavy_COLUMNS_ONLY=['divide original27 by L=aFC2','divide original31 by Lgamma'],
  cofactor_roots=cofroots,remaining11_profiles=profiles,
  universal_identity='kappa*[H0,H1,P22]=0 from signed maximal cofactors',
  six_other_zeros_source_proof='All90 two-minors of bottom4x6 zero; nonzero anchorentry gives domain cancellation, no pointwise inverse. kappa kills [22,11,13,25,15,19].',
  generic_kappa_nonzero_gate='Hpair rank2 and Saux nonzero force rank[H0,H1,P22]=3 overQQ(a,e,d)',
  chart_only_source_replacement='On D(kappa_component_for_current15row8), replace that row by kappa. Coefficient not declared FULLUunit.',
  no_new_scalar_values=True,all8_Q9_tails_retained=True,no_fullrank_or_bad_locus_claim=True)
 dump(out/'operator-and-support.json',recipe)
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 result=dict(id=ID,all_passed=True,frozen=True,source_only=True,typed_nodes=len(nodes),cofactor_support_upper=[b.size(c)for c in cofs],
  cheapest_profile=profiles[0],all11_output_support_upper=sum(x['output_bound']['support_upper']for x in profiles),
  coefficientheight_lifetime_gates_NOT_yet_evaluated=True,no_new_scalar_job_or_rank=True,input_sha256=pins,evidence_sha256={'operator-and-support.json':sha(out/'operator-and-support.json')},seconds=time.monotonic()-started)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;dump(out/'report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
