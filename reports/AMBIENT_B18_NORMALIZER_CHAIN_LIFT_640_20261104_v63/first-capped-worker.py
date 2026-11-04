#!/usr/bin/env python3
"""Source-bound finite operator B with UP=QB; never evaluates old r/g germs."""
import argparse
import hashlib
import itertools
import json
import resource
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63'
PATHS = {
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'relative':'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'normalizer':'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json',
 'C0_peer':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
 'relative_peer':'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW/report.json',
 'Q_syzygy_peer':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63_REVIEW/report.json'}


def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()


def padd(a,b,c=F(1)):
 a=dict(a)
 for m,x in b.items():
  v=a.get(m,F(0))+c*x
  if v:a[m]=v
  else:a.pop(m,None)
 return a


def pmul(a,b):
 out={}
 for m,x in a.items():
  for n,y in b.items():
   k=tuple(sorted(m+n));v=out.get(k,F(0))+x*y
   if v:out[k]=v
   else:out.pop(k,None)
 return out


class Operators:
 def __init__(self):
  self.nodes=[];self.cache={};self.poly_cache={};self.poly_terms=0
  self.zero=self.const(0);self.one=self.const(1)
 def put(self,n):
  if n in self.cache:return self.cache[n]
  assert len(self.nodes)<50000,'typed50k_guard'
  i=len(self.nodes);self.nodes.append(n);self.cache[n]=i;return i
 def const(self,x):
  x=F(x);return self.put(('c',x.numerator,x.denominator))
 def ref(self,i):return self.put(('actual_C0_scalar_function',i))
 def add(self,*xs):
  terms={};c=F(0)
  for x in xs:
   n=self.nodes[x]
   if n[0]=='c':c+=F(n[1],n[2]);continue
   ts=n[1] if n[0]=='sum' else ((x,1,1),)
   for i,a,b in ts:
    v=terms.get(i,F(0))+F(a,b)
    if v:terms[i]=v
    else:terms.pop(i,None)
  if c:terms[self.one]=terms.get(self.one,F(0))+c
  if not terms:return self.zero
  if len(terms)==1:
   i,q=next(iter(terms.items()))
   if q==1:return i
   if i==self.one:return self.const(q)
  return self.put(('sum',tuple((i,q.numerator,q.denominator) for i,q in sorted(terms.items()))))
 def scale(self,x,c):
  c=F(c)
  if not c:return self.zero
  if c==1:return x
  return self.mul(self.const(c),x)
 def neg(self,x):return self.scale(x,-1)
 def mul(self,*xs):
  c=F(1);f=[]
  for x in xs:
   n=self.nodes[x]
   if n[0]=='c':c*=F(n[1],n[2])
   elif n[0]=='product':
    c*=F(n[1],n[2]);f.extend(n[3])
   else:f.append(x)
  if not c:return self.zero
  if not f:return self.const(c)
  f.sort()
  if len(f)==1 and c==1:return f[0]
  return self.put(('product',c.numerator,c.denominator,tuple(f)))
 def poly(self,i):
  if i in self.poly_cache:return self.poly_cache[i]
  n=self.nodes[i];o=n[0]
  if o=='c':p={():F(n[1],n[2])} if n[1] else {}
  elif o=='actual_C0_scalar_function':p={(i,):F(1)}
  elif o=='sum':
   p={}
   for j,a,b in n[1]:p=padd(p,self.poly(j),F(a,b))
  elif o=='product':
   p={():F(n[1],n[2])}
   for j in n[3]:p=pmul(p,self.poly(j))
  else:raise ValueError(n)
  assert len(p)<150000
  self.poly_terms+=len(p);assert self.poly_terms<150000,'universal_polynomial_support_guard'
  self.poly_cache[i]=p;return p
 def equal(self,a,b):return self.poly(a)==self.poly(b)


def main():
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True)
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert all(D[k]['all_passed'] for k in ['C0_peer','relative_peer','Q_syzygy_peer'])
 S,G,C,N=(D[k] for k in ['source','relative','C0','normalizer'])
 full=C['DAG_scalar_nodes'];base=S['DAG_scalar_nodes']
 assert full[:len(G['DAG_scalar_nodes'])]==G['DAG_scalar_nodes']
 assert full[:len(base)]==base
 records=[r for r in G['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if r['weight']==0]
 assert len(records)==68
 oldpairs=G['actual_Qn_R0_squared_generator_matrices']['0']
 P=C['honest_D0_relation_columns'];assert len(P)==640 and len(oldpairs)==612
 assert P[:28]==C['all28_C0_Koszul_kernel_columns'] and P[28:]==C['all612_Q0_lifted_columns']
 assert all(P[28+k][0]==p[0] for k,p in enumerate(oldpairs))
 (out/'progress.json').write_text(json.dumps({'phase':'source_and_size_preflight','actual_source_prefix_bound':True,
    'max_B_scalar_slots':48*640,'max_Hadamard_scalar_slots':68*12*4,
    'typed_nodes_guard':50000,'universal_poly_terms_guard':150000})+'\n')
 pool=Operators();Z=pool.zero;O=pool.one
 leaf={};aliases={};source_identity_records=[]
 def source(i):
  i=aliases.get(i,i)
  if i in leaf:return leaf[i]
  n=full[i];o=n[0]
  if o=='c':v=pool.const(F(n[1],n[2]))
  elif o=='a':v=pool.add(source(n[1]),source(n[2]))
  elif o=='n':v=pool.neg(source(n[1]))
  elif o=='m':v=pool.mul(source(n[1]),source(n[2]))
  elif o in ['r','g','i','v']:v=pool.ref(i)
  else:raise ValueError((i,n))
  leaf[i]=v;return v
 def opaque(i):
  leaf[i]=pool.ref(i);return leaf[i]
 ca,cb=S['families']['CA'],S['families']['CB']
 rows=[18,19,28,29,60,67,72,63,64,68,73,74]
 def linear_form(i,memo):
  if i in memo:return memo[i]
  n=base[i];o=n[0]
  if o=='c':a={-1:F(n[1],n[2])} if n[1] else {}
  elif o=='n':a={j:-v for j,v in linear_form(n[1],memo).items()}
  elif o=='a':a=padd(linear_form(n[1],memo),linear_form(n[2],memo))
  elif o=='m':
   x,y=linear_form(n[1],memo),linear_form(n[2],memo)
   if not x or not y:a={}
   elif set(x)<={-1}:a={j:x[-1]*v for j,v in y.items()}
   elif set(y)<={-1}:a={j:y[-1]*v for j,v in x.items()}
   else:a={i:F(1)}
  else:a={i:F(1)}
  memo[i]=a;return a
 lfm={}
 for family in [ca,cb]:
  for r,i in family['actual_zero_odd_I0_generators']:
   j=family['complete_normalized_physical_w_coefficients'][0][r][1]
   assert linear_form(i,lfm)==linear_form(j,lfm)
   aliases[j]=i
   source_identity_records.append({'kind':'exact_QQ_linear_source_alias','odd_generator':i,'physical_odd_coefficient':j})
 # The eight new C0 generators are actual regular functions; no function division is performed.
 cgens=C['actual_C0_basis8_scalar_nodes'];assert len(cgens)==8
 for i in cgens:opaque(i)
 EA=[opaque(r['node']) for r in C['all35_regular_coordinate_divisions'][0]['all7_regular_coordinate_divisions']]
 fA=ca['actual_zero_odd_I0_generators'][7][1]
 leaf[fA]=pool.add(source(cgens[7]),*(pool.mul(e,source(z)) for e,z in zip(EA,cgens[:7])))
 source_identity_records.append({'kind':'accepted_actual_CA_f63_tail','function':fA,'gA':cgens[7],
                                'coordinate_quotients':[r['node'] for r in C['all35_regular_coordinate_divisions'][0]['all7_regular_coordinate_divisions']]})
 for t,div in zip(C['all35_regular_coordinate_divisions'][1:],C['all4_known_regular_gA_divisions']):
  i=t['original_pure_CA_source_function']
  cs=[opaque(r['node']) for r in t['all7_regular_coordinate_divisions']]+[opaque(div['node'])]
  leaf[i]=pool.add(*(pool.mul(c,source(z)) for c,z in zip(cs,cgens)))
  source_identity_records.append({'kind':'accepted_actual_CA_zero_graph_membership','function':i,'coefficient_nodes':[r['node'] for r in t['all7_regular_coordinate_divisions']]+[div['node']]})
 EB=[opaque(i) for i in G['CB_f63_zero7_telescoping_coefficients']]
 bn=[opaque(i) for i in G['CB_scalar_odd_generators_normalized_coefficients'][:8]]
 fB=cb['actual_zero_odd_I0_generators'][7][1]
 leaf[fB]=pool.add(bn[7],*(pool.mul(e,b) for e,b in zip(EB,bn[:7])))
 source_identity_records.append({'kind':'accepted_actual_CB_f63_tail','function':fB,'gB':G['CB_g_restriction_node'],'coordinate_quotient_nodes':G['CB_f63_zero7_telescoping_coefficients']})
 memberships=dict(G['all87_actual_odd_membership_lifts'])
 hs={r:[opaque(i) if i else Z for i in memberships[r][:8]] for r in rows[8:]}
 for r in rows[8:]:
  b=cb['complete_normalized_physical_w_coefficients'][0][r][1]
  leaf[b]=pool.add(*(pool.mul(h,source(i)) for h,(_,i) in zip(hs[r],cb['actual_zero_odd_I0_generators'])))
  source_identity_records.append({'kind':'accepted_actual_CB_zero_graph_membership','physical_row':r,'function':b,'coefficient_nodes':memberships[r][:8]})
 # Opaque primitive even/odd normalizer coefficients remain the same original R0 functions.
 aA=[];aB=[];bA=[];bB=[]
 for r in rows:
  aa,ba=ca['complete_normalized_physical_w_coefficients'][0][r]
  ab,bb=cb['complete_normalized_physical_w_coefficients'][0][r]
  if aa not in leaf:opaque(aa)
  if ab not in leaf:opaque(ab)
  for i in [aliases.get(ba,ba),aliases.get(bb,bb)]:
   if i not in leaf:opaque(i)
  aA.append(source(aa));aB.append(source(ab));bA.append(source(ba));bB.append(source(bb))
 nuA=opaque(ca['actual_tracefree_nu']);nuB=opaque(cb['actual_tracefree_nu'])
 def va(*vs):return [pool.add(*(v[i] for v in vs)) for i in range(4)]
 def vs(v,c):return [pool.mul(x,c) for x in v]
 def vm(v,w):
  out=[Z]*4
  for i,x in enumerate(v):
   for j,y in enumerate(w):
    c=pool.mul(x,y)
    if i&j&1:c=pool.mul(c,nuA)
    if i&j&2:c=pool.mul(c,nuB)
    out[i^j]=pool.add(out[i^j],c)
  return out
 scalar=lambda x:[x,Z,Z,Z]
 # Raw centered odd coordinates X_i-aB_i, and their unipotent g_B replacement.
 YA=[[pool.add(a,pool.neg(b)),c,Z,Z] for a,b,c in zip(aA,aB,bA)]
 ZA=YA[:8]
 ZA[7]=va(YA[7],*(vs(YA[i],pool.neg(EB[i])) for i in range(7)))
 ZB=[[Z,Z,b,Z] for b in bn]
 L=[[O if i==j else Z for j in range(12)] for i in range(8)]
 for i in range(7):L[7][i]=pool.neg(EB[i])
 F12=[[pool.add(aa,pool.neg(ab)),ba,pool.neg(bb),Z] for aa,ab,ba,bb in zip(aA,aB,bA,bB)]
 hadamards=[];values=[];semantic=[]
 for r in records:
  lab=r.get('label','')
  if lab.startswith('regular_Koszul_') or lab.startswith('cross_'):
   ii=list(map(int,lab.split('_')[-2:]));i,j=ii
   H=[scalar(pool.add(pool.mul(bn[i],L[j][t]),pool.neg(pool.mul(bn[j],L[i][t])))) for t in range(12)]
   v=va(vs(ZA[j],bn[i]),vs(ZA[i],pool.neg(bn[j])))
   desc={'kind':'linear_Koszul','i':i,'j':j}
  elif lab.startswith('product_'):
   i,j=map(int,lab.split('_')[-2:])
   H=[va(vs(ZA[j],L[i][t]),vs(ZB[i],L[j][t])) for t in range(12)]
   v=va(vm(ZA[i],ZA[j]),scalar(pool.neg(pool.mul(nuB,bn[i],bn[j]))))
   desc={'kind':'quadratic_product','i':i,'j':j}
  else:
   r0=r['raw_row'];j=rows.index(r0)
   H=[scalar(hs[r0][t]) if t<8 else scalar(pool.neg(O) if t==j else Z) for t in range(12)]
   v=va(scalar(pool.add(aB[j],pool.neg(aA[j]))),[Z,pool.neg(bA[j]),Z,Z],
        *(vs(YA[t],hs[r0][t]) for t in range(8)))
   desc={'kind':'linear_graph','physical_row':r0}
  old=r['normalized_CA_trace_pair']
  actual=[source(old[0]),source(old[1]),Z,Z]
  assert all(pool.equal(x,y) for x,y in zip(v,actual)),('actual_polynomial_CA_pair',desc)
  q=[Z]*4
  for f,h in zip(F12,H):q=va(q,vm(f,h))
  assert all(pool.equal(x,y) for x,y in zip(v,q)),('universal_Hadamard',desc)
  hadamards.append(H);values.append(v);semantic.append(desc)
 (out/'phase1-Hadamard68-operator.json').write_text(json.dumps({'operator_nodes':pool.nodes,
    'Hadamard68_by12_C4':hadamards,'actual_CA_values68_C4':values,'semantics':semantic},separators=(',',':'))+'\n')
 (out/'progress.json').write_text(json.dumps({'phase':'all68_actual_CA_pairs_and_Hadamards_verified',
    'operator_nodes':len(pool.nodes),'seconds':time.monotonic()-start})+'\n')
 # Bind the exact existing Q48, not merely a new equivalent presentation.
 mp=dict(N['original_node_to_pruned_node']);nseed={}
 for i in list(leaf)+list(aliases):
  if i in mp:nseed[mp[i]]=source(i)
 nseed[0]=Z;nseed[1]=O
 def normal(i):
  if i in nseed:return nseed[i]
  n=N['DAG_scalar_nodes'][i];o=n[0]
  if o=='c':v=pool.const(F(n[1],n[2]))
  elif o=='a':v=pool.add(normal(n[1]),normal(n[2]))
  elif o=='n':v=pool.neg(normal(n[1]))
  elif o=='m':v=pool.mul(normal(n[1]),normal(n[2]))
  else:raise ValueError(('unbound_actual_normalizer_scalar',i,n))
  nseed[i]=v;return v
 qb=[]
 for j,f in enumerate(F12):
  for k in range(4):
   w=[Z]*4;w[k]=O;q=vm(f,w);actual=[normal(i) for i in N['normalizer_tensor_relation_columns4_by48'][4*j+k]]
   assert all(pool.equal(x,y) for x,y in zip(q,actual)),('Q48frame',j,k)
   qb.append(q)
 U=[[normal(i) for i in col] for col in N['actual_C0_free9_to_E4_columns']]
 assert len(U)==9
 oldmult=[scalar(O)]+[[Z,source(i),Z,Z] for _,i in ca['actual_zero_odd_I0_generators']]
 B=[{} for _ in range(28)];proofs=[];nonzero=0
 for k,(rec,H,v) in enumerate(zip(records,hadamards,values)):
  for m,mult in enumerate(oldmult):
   q=k*9+m;oldpair=oldpairs[q]
   expected=vm(v,mult)
   actual=[source(oldpair[0]),source(oldpair[1]),Z,Z]
   assert all(pool.equal(x,y) for x,y in zip(expected,actual)),('actual_Q0_pair_order',k,m)
   honest=[Z]*4
   for u,c in zip(U,P[28+q]):honest=va(honest,vs(u,source(c)))
   assert all(pool.equal(x,y) for x,y in zip(honest,actual)),('honest_C0_map',k,m)
   col={}
   for j,h in enumerate(H):
    for c,x in enumerate(vm(h,mult)):
     if x!=Z:col[4*j+c]=x
   rhs=[Z]*4
   for j,c in col.items():rhs=va(rhs,vs(qb[j],c))
   assert all(pool.equal(x,y) for x,y in zip(honest,rhs)),('UP_equals_QB',k,m)
   B.append(col);nonzero+=len(col)
   proofs.append({'D0_column':28+q,'CB_polynomial':k,'old_CA_multiplier':m,
                  'actual_Q0_pair':oldpair,'honest_C0_lift':P[28+q]})
 for k,p in enumerate(P[:28]):
  v=[Z]*4
  for u,c in zip(U,p):v=va(v,vs(u,source(c)))
  assert all(pool.equal(x,Z) for x in v),('intrinsic_C0_zero',k)
 assert len(B)==640 and nonzero<=30720 and len(pool.nodes)<150000
 evidence={'operator_scalar_nodes':pool.nodes,'actual_scalar_reference_artifact':'C0',
  'scalar_source_artifacts':{k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()},
  'scalar_references_keep_actual_original_r_g_functions_opaque':True,
  'accepted_actual_regular_function_identities':source_identity_records,
  'physical_target_rows12':rows,'coefficient_C_basis':['1','tauA','tauB','tauA*tauB'],
  'nu_A_B_operator_nodes':[nuA,nuB],
  'CB_polynomial68_source_records':records,'CB_polynomial68_semantics':semantic,
  'normalized_target_gradient8_by12':L,'CB_g_target_unipotent_coefficients7':EB,
  'actual_Hadamard68_by12_C4':hadamards,'actual_F12_C4':F12,
  'exact_existing_Q48_C4':qb,'exact_existing_U9_C4':U,
  'old_actual_CA_Q0_multipliers9_C4':oldmult,
  'old_multiplier_last_is_f63_not_new_gA':True,
  'B48_by640_sparse_columns':[[[i,v] for i,v in sorted(col.items())] for col in B],
  'all612_source_pair_and_honest_lift_bindings':proofs,
  'all68_actual_CA_evaluation_pairs_verified':True,'all68_Hadamard_identities_verified':True,
  'all48_existing_Q_frame_columns_verified':True,'all640_UP_equals_QB_verified':True,
  'all28_intrinsic_C0_columns_have_zero_lift':True,
  'typed_operator_nodes':len(pool.nodes),'B_nonzero_scalar_operator_entries':nonzero,
  'expanded_scalar_coefficients_evaluated':False,'new_division_or_source_jet_or_rank_search':False,
  'complete_low_comparison_ancestry':{'P':'actual9by640','U':'actual4by9','Q':'actual4by48','B':'actual48by640'},
  'K_or_specialized_kernel_or_depth_computed':False}
 data=json.dumps(evidence,separators=(',',':'))+'\n';assert len(data.encode())<32<<20
 ep=out/'actual-UP-QB-640-typed-operator.json';ep.write_text(data)
 report={'all_passed':True,'actual_chain_identity_UP_QB':True,'B_shape':[48,640],
  'B_nonzero_entries':nonzero,'operator_nodes':len(pool.nodes),'bytes':len(data.encode()),
  'no_new_source_jets_nonunit_divisions_rank_or_kernel_computation':True,
  'old_actual_regular_r_g_functions_kept':True,
  'input_hashes':evidence['scalar_source_artifacts'],
  'evidence':{'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)},
  'code_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start,
  'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__=='__main__':main()
