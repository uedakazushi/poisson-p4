#!/usr/bin/env python3
"""Independent fixed rewrite/rank/kernel receipts and one cubic source channel."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
CLOSURE='AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63'
CHANNEL='AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def sp(raw):return {int(i):Q(q) for i,q in raw if Q(q)}
def ax(a,b,c=Q(1)):
 a=dict(a)
 for i,q in b.items():
  v=a.get(i,Q(0))+c*q
  if v:a[i]=v
  else:a.pop(i,None)
 return a
def enc(a):return [[i,str(q)] for i,q in sorted(a.items())]
def combine(coeff,cols):
 a={}
 for i,q in coeff.items():a=ax(a,cols[i],q)
 return a
def fixed_receipts(columns,cert,prefix=None):
 """Check supplied operations; never select a new pivot."""
 ops=cert['fixed_elimination_receipts'];ops=ops if prefix is None else ops[:prefix]
 assert len(ops)==len(columns)
 bases=[];basis_source=[];kernel=[];zero_positions=[]
 for j,(col,op) in enumerate(zip(columns,ops)):
  assert op['column']==j
  v=dict(col);c={j:Q(1)}
  for p,a in op['subtract']:
   assert p<len(bases);a=Q(a);assert v.get(min(bases[p]),Q(0))==a
   v=ax(v,bases[p],-a);c=ax(c,basis_source[p],-a)
  if 'pivot_row' in op:
   r=op['pivot_row'];d=Q(op['divide']);assert v and min(v)==r and v[r]==d and d
   v={i:q/d for i,q in v.items()};c={i:q/d for i,q in c.items()}
   assert v==sp(op['normalized']);assert combine(c,columns)==v
   bases.append(v);basis_source.append(c)
  else:
   assert op.get('zero') and not v and not combine(c,columns)
   assert c[j]==1 and all(i==j or i in cert['selected_original_columns'][:len(bases)] for i in c)
   kernel.append(c);zero_positions.append(j)
 return bases,basis_source,kernel,zero_positions
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 hashes={}
 def load(p):hashes[p]=sha(p);return read(p)
 def get(st,f='report.json'):return load('reports/'+st+'/'+f)
 def bind(st,want=None):
  r=get(st);assert r['all_passed']
  if want:assert hashes['reports/'+st+'/report.json']==want
  for p,h in r.get('input_sha256',{}).items():assert sha(p)==h;hashes[p]=h
  for n,h in r.get('evidence_sha256',{}).items():p='reports/'+st+'/'+n;assert sha(p)==h;hashes[p]=h
  if 'code_sha256' in r:
   expected=ROOT/'code'/('ambient_B18_T_known_quadratic_cubic_closure_20261104_v63.py' if st==CLOSURE else 'ambient_B18_T_linear_syzygy_pure_r3_20261104_v63.py')
   if st in (CLOSURE,CHANNEL):p=str(expected.relative_to(ROOT));assert sha(p)==r['code_sha256'];hashes[p]=sha(p)
  return r
 RC=bind(CLOSURE);RS=bind(CHANNEL,'7bbfbe0916a3719ac7d6b97cb0a05d2f54854342c4002e4232f2b02181c774bb')
 fullrank=get('AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW');assert fullrank['all_passed']
 assert hashes['reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW/report.json']=='143af8f9218d2e43626f58db9abd34827dee9e419f6ec4a9626b5ca5d36caae1'
 for st,want in [('AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63_REVIEW','c933016430a1779bed947b0c81484b0aae23afbee23e4df15a386ea959bf0dd9'),('AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW','ca127c4dc0f3c12a4d4d80dfc1ca316f42dee8b90ecc74cb3e5e423752c906b7'),('AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW','998bcf929242837dee44de326f05950b1b295a8828f9cfc8785c082503804293')]:
  d=get(st);assert d['all_passed'] and hashes['reports/'+st+'/report.json']==want
 Low=get('AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63','actual-linear-quadratic-T-certificate.json')
 Lin=get('AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63','phase1-actual-five-target-linear-relations.json')
 Quad=get('AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63','phase2-certified-quadratic-input.json')
 Rank=get('AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63','actual-three-by620-quadratic-image-certificate.json')
 Dual=get('AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63','three-quadratic-duals.json')
 C=get(CLOSURE,'actual-known-quadratic-cubic-closure-certificate.json');W=get(CLOSURE,'linear-module-rewrite-certificate.json')
 S=get(CHANNEL,'actual-linear-syzygies-and-pure-r3-source-values.json')
 linear=list(map(sp,Low['complete_linear_rank_certificate']['basis_normalized_columns']));assert len(linear)==37
 assert W['linear_basis37']==Low['complete_linear_rank_certificate']['basis_normalized_columns']
 leading={min(v):i for i,v in enumerate(linear)};assert len(leading)==37
 assert all(v[min(v)]==1 for v in linear)
 m2=list(combinations_with_replacement(range(8),2));m3=list(combinations_with_replacement(range(8),3));im2={m:i for i,m in enumerate(m2)};im3={m:i for i,m in enumerate(m3)}
 assert len(m2)==36 and len(m3)==120 and W['monomials3_120']==list(map(list,m3))
 free1=[i for i in range(40) if i not in leading]
 standard=[120*t+j for t in range(5) for j,m in enumerate(m3) if all(8*t+k in free1 for k in m)]
 assert len(free1)==3 and len(standard)==5 and W['standard_global600_positions']==standard
 si={i:j for j,i in enumerate(standard)};steps={i:(step,sp(v)) for i,step,v in W['all600_rewrite_steps']};assert set(steps)==set(range(600))
 def order(i):t,j=divmod(i,120);m=m3[j];return (-t,)+tuple(m.count(k) for k in range(8))
 normal={}
 # Independent order-sorted dynamic verification, rather than producer DFS.
 for i in sorted(range(600),key=order):
  t,j=divmod(i,120);m=m3[j];step,want=steps[i]
  eligible=[k for k in m if 8*t+k in leading]
  if not eligible:assert step=={'standard':True};v={si[i]:Q(1)}
  else:
   k=step['leading_variable'];p=step['linear_relation_basis_index'];assert k in eligible and leading[8*t+k]==p
   rest=list(m);rest.remove(k);children=[];v={}
   for idx,q in linear[p].items():
    if idx==8*t+k:continue
    tt,var=divmod(idx,8);child=120*tt+im3[tuple(sorted(rest+[var]))]
    assert order(child)<order(i) and child in normal
    children.append([child,str(-q)]);v=ax(v,normal[child],-q)
   assert children==step['children']
  assert v==want;normal[i]=v
 def reduce(raw):return combine(raw,normal)
 true177=list(map(sp,Low['quadratic_certificate']['basis_normalized_columns']))
 rankone=sp(Rank['complete_fixed_rank_certificate']['basis_normalized_columns'][0]);reps=Dual['quotient_representative_rows']
 new={reps[i]:q for i,q in rankone.items()};assert new and C['actual_new_quad_canonical_lift_from_rank1_source']['canonical_quad_representative']==enc(new)
 assert [sum((q*new.get(i,Q(0)) for i,q in sp(row).items()),Q(0)) for row in Dual['projection3_by180_rows_sparse']]==[rankone.get(i,Q(0)) for i in range(3)]
 quads=true177+[new];assert C['quadratic_basis178']==list(map(enc,quads))
 reduced=[]
 for j,q in enumerate(quads):
  for k in range(8):
   raw={}
   for idx,a in q.items():
    t,l=divmod(idx,36);pos=120*t+im3[tuple(sorted(m2[l]+(k,)))];raw=ax(raw,{pos:a})
   assert C['cubic_column_labels'][8*j+k]==[j,k];reduced.append(reduce(raw))
 assert list(map(enc,reduced))==C['all1424_reduced_cubic_columns']
 cert=C['complete_reduced_rank_certificate'];cbasis,_,_,_=fixed_receipts(reduced,cert)
 assert len(cbasis)==cert['rank']==3 and cbasis==list(map(sp,cert['basis_normalized_columns']))
 assert all(combine(sp(c),cbasis)==v for c,v in zip(cert['all_column_factorization_in_normalized_basis'],reduced))
 # Given 3x3 selected minor, evaluated directly with six signed terms.
 rows=cert['selected_minor_rows'];selected=cert['selected_original_columns'];d=Q(0)
 from itertools import permutations
 for pm in permutations(range(3)):
  sign=(-1)**sum(pm[i]>pm[j] for i in range(3) for j in range(i+1,3));term=Q(sign)
  for i in range(3):term*=reduced[selected[pm[i]]].get(rows[i],Q(0))
  d+=term
 assert d==Q(cert['selected_minor_determinant']) and d
 overlaps=W['all1332_linear_overlap_remainders_retained_and_spanned'];assert len(overlaps)==37*36
 for p,l in enumerate(linear):
  for j,m in enumerate(m2):
   raw={}
   for idx,a in l.items():t,k=divmod(idx,8);raw=ax(raw,{120*t+im3[tuple(sorted(m+(k,)))]:a})
   v=reduce(raw);saved=overlaps[p*36+j]
   assert saved['linear_basis']==p and saved['quadratic_multiplier']==list(m)
   assert sp(saved['normal_form'])==v and combine(sp(saved['cubic_relation_span_coordinates']),cbasis)==v
 assert C['certified_cubic_rank']==3 and C['degree3_HF_upper_bound_only']==2 and C['actual_degree3_HF'] is None
 assert not C['m_cubed_T_zero_proved'] and C['total_T_length_proved'] is None and not W['linear_basis_assumed_Groebner']
 # Rebuild all296 actual linear multiples independently from the original
 # selected actual R1 columns, not from normalized linear basis surrogates.
 selected37=Low['complete_linear_rank_certificate']['selected_original_columns'];assert len(selected37)==37
 actualL=list(map(sp,Lin['actual40_by657_linear_relation_columns']));constraints=[]
 for j,idx in enumerate(selected37):
  for k in range(8):
   v={}
   for pos,q in actualL[idx].items():t,x=divmod(pos,8);v=ax(v,{36*t+im2[tuple(sorted((x,k)))]:q})
   constraints.append(v);lab=Quad['source_columns464'][8*j+k]
   assert lab['kind']=='variable_times_actual_linear_relation' and lab['variable']==k and lab['source_relation']==Lin['actual657_minimal_relation_source_labels'][idx]
 assert constraints==list(map(sp,Quad['all464_certified_quadratic_columns_sparse'][:296]))
 basis,_,kernels,zeropos=fixed_receipts(constraints,Low['quadratic_certificate'],296)
 assert len(basis)==177 and len(kernels)==119
 assert len(S['all119_linear_coefficient_syzygies'])==119
 for i,(v,j,saved) in enumerate(zip(kernels,zeropos,S['all119_linear_coefficient_syzygies'])):
  assert saved['fixed_zero_receipt_source296']==j and sp(saved['coefficients296_sparse'])==v
  assert not combine(v,constraints)
 # Reconstruct selected37 original768 combinations from the accepted full
 # source constant-rank receipts; no producer receipt_combos helper is used.
 Mu=get('AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63','actual-linear-cycle-muT-certificate.json')
 wcols=list(map(sp,Mu['actual84_by741_w0_columns_sparse']))
 _,_,wzero,_=fixed_receipts(wcols,Mu['complete_rank_certificate']);assert len(wzero)==657
 base768=list(map(sp,Mu['actual_constant_source_combinations']));source37=[combine(wzero[i],base768) for i in selected37]
 X=get('AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63','actual-projected-three-class-X2-operator-preflight.json')
 X0=get('AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63','actual-X0-and-unit-cancellation.json')['full_X0_56_by768']
 assert Lin['five_remaining_free_source_rows'][4]==55
 assert all(sp(v).get(55,Q(0))==0 for v in X['actual_pivot_X0_columns27_sparse56'])
 assert all(sp(v).get(55*8+7,Q(0))==0 for v in X['actual_pivot_X1_by27_sparse56_times8'])
 assert all(sp(v).get(39,Q(0))==0 for v in Lin['normalized84_unitcycle_Q1_columns_sparse'])
 assert all(sp(v).get(39,Q(0))==0 for v in X['Ystar1_40_by84_columns'])
 assert X['pi_free32_rows'][-1]==[[31,'1']] and X['free_kernel_rows32'][31]==55
 for c in source37:assert all(sum((a*Q(row[j]) for j,a in c.items()),Q(0))==0 for row in X0)
 C0=get('AMBIENT_B18_T0_HONEST_C0_20261104_v63','honest-C0-D0-D4-DAG.json')
 native=get('AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63','actual-honest-T116-936-Euler-operator.json')['actual_X_source']['source_X_sparse_columns']
 assert all(all(i<48 for i,_ in col) for col in native[640:])
 ownsource=get('AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW','third-class-source-certificate.json')
 sourcegate=get('AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW')
 assert hashes['reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW/third-class-source-certificate.json']==sourcegate['evidence_sha256']['third-class-source-certificate.json']
 vals={i:(list(map(Q,v)),sp(g)) for i,v,g in ownsource['all90_known_univariate_parts_and_symbolic_second_sensitivities']}
 vals[0]=([Q(0)]*3,{});assert C0['DAG_scalar_nodes'][0]==['c',0,1]
 h=[];allroots=set();h_receipts=[]
 for idx,c,saved in zip(selected37,source37,S['actual37_e5_r2_source_coefficients']):
  rp={}
  for j,a in c.items():
   if j<640:rp=ax(rp,{C0['honest_D0_relation_columns'][j][8]:a})
  allroots.update(rp);value=Q(0)
  for root,a in rp.items():assert root in vals and not vals[root][1];value+=a*vals[root][0][2]
  assert saved['minimal_cycle_index657']==idx and sp(saved['actual768_constant_source_combination'])==c
  assert sp(saved['literal_rawP8_last_root_coefficients'])==rp and Q(saved['actual_e5_r2_coefficient'])==value
  h.append(value);h_receipts.append({'source_index657':idx,'original768':enc(c),'P_roots':enc(rp),'rr':str(value)})
 assert len(allroots)==22 and 0 in allroots and len(allroots-{0})==21
 assert all(q.get(179,Q(0))==0 for q in quads)
 results=[]
 for v,saved in zip(kernels,S['all119_linear_coefficient_syzygies']):
  a=sum((q*h[j//8] for j,q in v.items() if j%8==7),Q(0));assert a==0 and Q(saved['actual_e5_r3_cubic_relation_value'])==a;results.append(str(a))
 assert S['nonzero_pure_r3_channel_syzygy_indices']==[]
 out1=ROOT/'reports'/(CLOSURE+'_REVIEW');out2=ROOT/'reports'/(CHANNEL+'_REVIEW')
 for out in (out1,out2):assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 def save(out,n,d):(out/n).write_text(json.dumps(d,indent=2)+'\n')
 save(out1,'independent-rewrite-and-closure-proof.json',{'all600_steps_coefficient_and_strict_order_checks':True,'all1332_linear_overlap_remainders_kept_and_spanned':True,'all1424_known_quad_times_variable_columns_verified':True,'fixed_rank3_minor':str(d),'standard_cubic_positions':standard,'known_closure_quotient_dimension':2,'actual_HF3_or_total_length_not_determined':True,'linear_basis_not_assumed_Groebner':True})
 save(out2,'independent-kernel-and-source-channel.json',{'constraint_rank':177,'kernel119_coefficients':list(map(enc,kernels)),'unique_nonpivot_positions':zeropos,'all119_actual_constraints_zero':True,'actual37_source_rr_values':h_receipts,'all119_pure_r3_values':results,'all21_nonzero_id_roots_from_own_independent90_source_cache':True,'unknown_second_source_sensitivities_for_these_roots':0,'first27_and84_all_r_zero_gates_checked':True,'all178_known_quadratics_pure_rr_zero':True,'only_this_channel_zero':True})
 paths=[str(Path(__file__).relative_to(ROOT)),'research/notes/'+CLOSURE+'_REVIEW.md','research/notes/'+CHANNEL+'_REVIEW.md']
 for p in paths:hashes[p]=sha(p)
 for out,status in [(out1,{'known_linear_rank':37,'known_quadratic_rank':178,'known_cubic_closure_rank':3,'standard_cubic_count':5,'remaining_candidate_classes':2,'actual_HF3_or_total_length_claimed':False}),(out2,{'actual_linear_coefficient_syzygy_kernel_dimension':119,'source37_rr_values_bound_to_same_actual_Schur':True,'nonzero_pure_r3_channel_values':0,'all_other_cubic_channels_or_actual_HF3_claimed':False})]:
  report={'id':out.name,'all_passed':True,**status,'no_new_source_jet_solve_point_prime_GB_or_rank_search':True,'input_sha256':dict(hashes),'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
  save(out,'report.json',report);hashes[str((out/'report.json').relative_to(ROOT))]=sha(str((out/'report.json').relative_to(ROOT)))
  print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')}),flush=True)
if __name__=='__main__':main()
