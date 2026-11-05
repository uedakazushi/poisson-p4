#!/usr/bin/env python3
"""Independent source-transport/NF check of full-kernel G3 quadratics."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_G3_QUADRATIC_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(a):return {int(i):Q(c) for i,c in a}
def add(a,b,c=Q(1)):
 d=dict(a)
 for i,v in b.items():
  q=d.get(i,Q(0))+c*v
  if q:d[i]=q
  else:d.pop(i,None)
 return d
def scale(a,c):return {i:c*v for i,v in a.items() if c*v}
def mv(cols,v):
 d={}
 for j,c in v.items():d=add(d,cols[j],c)
 return d
def main():
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 pp=ROOT/'reports'/A/'report.json';pr=json.loads(pp.read_text());assert pr['all_passed']
 inputs={'primary':{'path':str(pp.relative_to(ROOT)),'sha256':sha(pp)}}
 D={}
 for p,h in pr['input_sha256'].items():
  path=ROOT/p;assert sha(path)==h;inputs[p]={'path':p,'sha256':h}
  if path.name!='actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json':D[path.name]=json.loads(path.read_text())
 for stem in ['AMBIENT_B18_ELL2_T2_FULL_G3_LINEAR_KERNEL_20261105_v63_REVIEW','AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_GB_COVER_20261105_v63_REVIEW']:
  path=ROOT/'reports'/stem/'report.json';assert json.loads(path.read_text())['all_passed'];inputs[stem]={'path':str(path.relative_to(ROOT)),'sha256':sha(path)}
 ep=pp.parent/'actual176-G3quadratic-factored-source-and-known-NF-input.json';assert sha(ep)==pr['evidence_sha256'][ep.name]
 E=json.loads(ep.read_text());inputs['primary_evidence']={'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)}
 assert sha(ROOT/'code/ambient_B18_ell2_T2_G3_quadratic_source_20261105_v63.py')==pr['code_sha256']
 assert sha(ROOT/'research/notes'/f'{A}.md')==pr['note_sha256']
 prefix=next(p for p in pr['input_sha256'] if p.endswith('actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json'))
 F=fields(ROOT/prefix,['actual_active_I3_origin_columns','PiY3_from102_columns'])
 I3=list(map(sp,F['actual_active_I3_origin_columns']));PiY=list(map(sp,F['PiY3_from102_columns']))
 t=D['actual-active-unit-bases-and26-frame-transport.json'];pref=D['active-five-d3-source-split-and-budget.json']
 Pi42=list(map(sp,t['active_Pi42_columns']));fixed=list(map(sp,t['active_F2_fixed182_columns']));D30=list(map(sp,pref['active_d3_columns']))
 norms={};reps={}
 for receipt in t['active_F2_fixed182_receipts']:
  j=receipt['source_column'];p=receipt['pivot'];w={j:Q(1)}
  for k,c in receipt['subtract']:w=add(w,reps[k],-Q(c))
  w=scale(w,1/Q(receipt['divide_by']));v=sp(receipt['normalized'])
  assert mv(fixed,w)==v and min(v)==p and v[p]==1
  norms[p]=v;reps[p]=w
 active=[0,1,4,5,7];vars6=active+[8]
 b2=[(2-len(w),w) for q in range(3) for w in itertools.combinations(active,q)]
 lookup={w:i for i,(_,w) in enumerate(b2)}
 raw=[]
 for q in range(4):
  for w in itertools.combinations(active,q):
   p=3-q
   for j in range(14):
    col={}
    if p:col[(14*lookup[w]+j)*6+5]=Q(1)
    for s,z in enumerate(w):col=add(col,{(14*lookup[w[:s]+w[s+1:]]+j)*6+vars6.index(z):Q((-1)**(p+s))})
    raw.append(col)
 assert len(raw)==364 and all(not mv(D30,c) for c in I3)
 P=[]
 for j in t['active126_source_columns']:
  col={}
  for k,c in raw[j].items():
   row,var=divmod(k,6)
   col=add(col,{6*i+var:a*c for i,a in Pi42[row].items()})
  P.append(col)
 assert P==list(map(sp,E['actual_selected126_rawD31_Pi42_linear_columns']))
 fixedcert=D['actual-full-G3linear-active-fixed-rank-and-kernel.json']
 kernels=[sp(x['original238_source_combination']) for x in fixedcert['full_constant_kernel_columns']]
 assert kernels==list(map(sp,E['source_kernel176'])) and len(kernels)==176
 corrections=[]
 for n,c in enumerate(kernels):
  w=mv(I3,c);assert w==sp(E['source_original_raw364_combinations'][n]) and not mv(D30,w)
  linear=mv(raw,w);correction={}
  for var in range(6):
   v={i//6:a for i,a in linear.items() if i%6==var}
   assert not mv(Pi42,v)
   original=dict(v);coeff={}
   # Replay fixed normalized pivots, never choose new pivots.
   for p in sorted(norms):
    a=v.get(p,Q(0))
    if a:coeff=add(coeff,reps[p],a);v=add(v,norms[p],-a)
   assert not v and mv(fixed,coeff)==original
   for i,a in coeff.items():
    if i>=56:correction[6*(i-56)+var]=a
  assert correction==sp(E['actual176_selected126_correction_linear_columns'][n]);corrections.append(correction)
 lin=D['actual-complete-L48-and-degree2-rewrite-source-input.json'];nodes=lin['degree2_monomial_rewrite_nodes'];memo={}
 def nf(k):
  if k in memo:return memo[k]
  node=nodes[k]
  if node[0]=='standard':v={node[1]:Q(1)}
  elif node[0]=='zero_ell2':v={}
  else:
   v={}
   for j,c in node[2]:v=add(v,nf(j),Q(c))
  memo[k]=v;return v
 monidx={tuple(m):i for i,m in enumerate(lin['degree2_monomial_order44'])}
 pairs=[w for w in itertools.combinations_with_replacement(range(6),2) if w!=(5,5)]
 NF=[]
 for row in range(42):
  for u,v in pairs:
   m=[0]*9;m[vars6[u]]+=1;m[vars6[v]]+=1;mi=monidx[tuple(m)];col={}
   for target,a in PiY[row].items():col=add(col,nf(lin['all264_NF_root_nodes'][target][mi]),a)
   NF.append(col)
 assert NF==list(map(sp,E['PiY_L48_degree2_NF42_times20_columns']))
 pairidx={w:i for i,w in enumerate(pairs)};products=0;surviving=0
 for corr in corrections:
  for k,x in corr.items():
   source,u=divmod(k,6)
   for i,y in P[source].items():
    row,v=divmod(i,6);pair=tuple(sorted((u,v)))
    if pair==(5,5):continue
    products+=1
    z=NF[row*20+pairidx[pair]]
    assert not z,('nonzero_summand_requires_separate_signed_job',k,i)
    surviving+=len(z)
 assert surviving==E['signed_product_summand_count']==0
 assert not E['signed_actual_quadratic_values_evaluated']
 proof={
  'exact_projection':'B2=[E56,D3sel126,Efree42-E56K]. Subtracting E56K from last42 columns leaves last42 inverse rows P unchanged. P is last42 inverse[E56,D3sel126,Efree42].',
  'exact_source_clearing':'B3 only subtracts selectedD3 multiples. P*D3sel=0, so reducedG3=P*D3*Ecomp=P*D3*I3(0). D30*I3(0)=0 and rawD3 affine imply reducedG3=P*D31*I3(0) exactly.',
  'inverse_differentiation':'P1(q)=-Pi42*D31sel*A126(q), with q in fixed182 span. This is the only second-degree term because D31I3 has order1.',
  'final_minimal6_projection':'G3 constant and bottom30 identically zero. For fullG3_1*c=0, the derivative of every later target unit projection multiplies zero. Thus quadratic is PiY0*G3_2*c; no projected393 argument is used.',
  'scope':'All176 displayed quadratic classes lie in the known L48 plus ell² module. This family supplies no new class. Other source syzygies and higher relations remain uncomputed; actual HF2 and total length stay open.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 data={'all176_fixed182_span_and_selected126_coefficients':True,'all840_projection_normal_forms':True,
       'all_composable_nonell2_factor_terms_have_zero_normal_form':True,'factor_pairs_tested_before_NF':products,'nonzero_signed_summands':0,
       'source_proof':proof,'no_new_physical_implicit_or_Taylor_array_query':True}
 path=out/'independent-G3quadratic-source-transport-and-zero-NF.json';path.write_text(json.dumps(data,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'exact_internal_unit_Taylor2_formula_accepted':True,'full_kernel176_known_quadratic_classes_zero_mod_L48':True,
  'actual_HF2_or_total_length_claimed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{path.name:sha(path)},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'factor_pairs_tested':products,'nonzero_summands':0}))
if __name__=='__main__':main()
