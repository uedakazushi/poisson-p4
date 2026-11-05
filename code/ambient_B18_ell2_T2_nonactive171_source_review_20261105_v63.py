#!/usr/bin/env python3
"""Independent fixed-source review; never evaluate the nine lower images."""
import hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_ELL2_T2_NONACTIVE171_SOURCE_PREFLIGHT_20261105_v63'
STEM=BASE+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(v):return {int(i):Q(x) for i,x in v if Q(x)}
def add(a,b,c=Q(1)):
 out=dict(a)
 for i,x in b.items():
  y=out.get(i,Q(0))+c*x
  if y:out[i]=y
  else:out.pop(i,None)
 return out
def image(cols,v):
 out={}
 for j,x in v.items():out=add(out,cols[j],x)
 return out
def count(x):
 if isinstance(x,int):return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except (ValueError,ZeroDivisionError):return 0
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return 0
def main():
 t=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 inputs={}
 def read(p,keys=None):
  path=ROOT/p;inputs[p]={'path':p,'sha256':sha(path)}
  return fields(path,keys) if keys else json.loads(path.read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed']
 assert r['code_sha256']==sha(ROOT/'code/ambient_B18_ell2_T2_nonactive171_source_20261105_v63.py')
 assert r['note_sha256']==sha(ROOT/'research/notes'/f'{BASE}.md')
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/exact171-decomposition-162-zero-9-lower-first-source-input.json'
 d=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 g=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json')
 first=read('reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json', ['actual_full_G3linear102_by418_columns','PiY3_from102_columns','PiBottom3_from56_columns'])
 G=list(map(sp,first['actual_full_G3linear102_by418_columns']));PY=list(map(sp,first['PiY3_from102_columns']))
 active=read('reports/AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63/actual-G2-linear-fixed-source-profile.json')
 G2=list(map(sp,active['actual_linear_G2_active102_by42_columns']))
 support=[j for j,v in enumerate(G2) if any(k%6 in (4,5) for k in v)]
 assert support==[0,1,29,31,36,38]
 assert [j for j,v in enumerate(PY[:42]) if v]==support and not any(PY[42:])
 projected=list(map(sp,g['actual_projected_G3linear_columns']));selected=g['selected_columns']
 for j in range(418):
  proj={}
  for k,x in G[j].items():
   row,var=divmod(k,9)
   for out,y in PY[row].items():proj=add(proj,{9*out+var:x*y})
  assert proj==projected[j]
  assert all(k%9 in ((2,3,6) if j>=238 else (0,1,4,5,7,8)) for k in proj)
 non=[j for j in selected if j>=238];assert non==d['projected_nonactive_selected9']==[238,239,267,280,281,309,322,323,351]
 mu=read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json', ['actual_selected_Y0_columns102','normalized_basis_rows','normalized_basis_original_selected_source_representations'])
 Y=list(map(sp,mu['actual_selected_Y0_columns102']));norm={p:sp(v) for p,v in mu['normalized_basis_rows']};rep={p:sp(v) for p,v in mu['normalized_basis_original_selected_source_representations']}
 # Only positive integer envelopes below; no lower-image signed products.
 maxbits=0
 def mark(b):
  nonlocal maxbits
  maxbits=max(maxbits,b[0].bit_length(),b[1].bit_length());return b
 def exact(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def plus(a,b):
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def times(a,b):return mark((a[0]*b[0],a[1]*b[1]))
 def cls(values,l1=False):
  v=list(values);D=math.lcm(*(x.denominator for x in v)) if v else 1
  u=[abs(x.numerator)*(D//x.denominator) for x in v]
  return mark((D,sum(u) if l1 else max([0]+u)))
 RE={p:cls(v.values(),True) for p,v in rep.items()};YE=cls(x for v in Y for x in v.values())
 records=d['all171_nonactive_projected_kernel_basis'];assert len(records)==171
 rest={x['source_index']:x for x in d['remaining9_sources']};struct=[];templates={};allcomb={}
 for rec in records:
  j=rec['source_index'];assert j>=238 and j not in non
  factor=sp(g['all418_original_basis_factors'][j]);assert all(selected[k]>=238 for k in factor)
  c=add({j:Q(1)},{selected[k]:x for k,x in factor.items()},Q(-1))
  assert c==sp(rec['original418_combination']) and not image(projected,c)
  allcomb[j]=c
  if j>=364:assert not factor;reason='G1zero2_or_F0zero3_fullfirstplane0'
  else:
   iz,local=divmod(j-238,42);assert all((k-238)//42==iz for k in c)
   v={(k-238)%42:x for k,x in c.items()}
   if local in templates:assert templates[local]==v
   else:templates[local]=v
   reason='G2zero1_column_outside_six_plane_sources' if local not in (31,36,38) else None
   if reason:assert not factor and local not in support
  assert reason==rec['structural_zero_reason']
  if reason:
   assert all(k%9 not in (7,8) for k in image(G,c));struct.append(j)
  else:
   rr=rest[j];assert sp(rr['original418_combination'])==c
   for var,key in [(7,'b1_r_selected99'),(8,'b1_ell_selected99')]:
    bound={}
    for k,x in c.items():
     for slot,y in G[k].items():
      if slot%9==var:bound[slot//9]=plus(bound.get(slot//9,(1,0)),times(exact(x),exact(y)))
    aE=(1,0)
    for p in sorted(norm):
     x=bound.pop(p,(1,0))
     if not x[1]:continue
     aE=plus(aE,times(x,RE[p]))
     for k,y in norm[p].items():
      if k!=p:bound[k]=plus(bound.get(k,(1,0)),times(x,exact(y)))
    times(aE,YE)
    expected={i//9:x for i,x in image(G,c).items() if i%9==var}
    assert expected==image(Y,sp(rr[key]))
    assert all(42<=i<93 and (i-42)%17 in (12,13,14,15,16) for i in expected)
 assert maxbits==r['point_b1_prior_bits']==261<=384
 assert struct==d['structural_zero162_source_indices'] and len(struct)==162 and len(templates)==39
 assert set(rest)=={269,274,276,311,316,318,353,358,360}
 # Exact cached firsts are bound by the independent17-bit source stage.
 hp='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_HSTAR_FIRST_20261105_v63'
 hr=read(hp+'/report.json');assert hr['all_passed'] and hr['prior_bits']==17
 h=read(hp+'/actual-cached-selected99-Hstar-r-ell-and-prior17.json')
 assert h['actual_Hstar_r_ell_columns99']==[d['Hstar1_r_selected99_columns'],d['Hstar1_ell_selected99_columns']]
 kp=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_KNOWN_INPUT_GATE_20261105_v63/pruned-known619-lower4-input-and-prior-resource-gate.json')
 assert first['PiBottom3_from56_columns']==kp['Pi3_source56_columns']
 Hr,He=[list(map(sp,M)) for M in h['actual_Hstar_r_ell_columns99']]
 lin=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_INPUT_GATE_20261105_v63/actual-complete-L48-and-degree2-rewrite-source-input.json')
 nf=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_LINEAR_DEGREE2_20261105_v63/actual432-overlap-rank-and-known-linear-cover-degree2.json')
 memo={i:sp(v) for i,v in nf['all264_linear_NF_columns']};mi={tuple(m):i for i,m in enumerate(lin['degree2_monomial_order44'])}
 for target,z,var in itertools.product(range(6),range(7),range(9)):
  mon=[0]*9;mon[z]+=1;mon[var]+=1
  assert not memo[lin['all264_NF_root_nodes'][target][mi[tuple(mon)]]]
 std={(t,tuple(m)):i for i,(t,m) in enumerate(lin['degree2_linear_standard_module_monomials'])}
 Ds=[]
 for chi in [(0,1,4),(11,0,270)]:
  for var in (7,8):
   mon=[0]*9;mon[7]+=1;mon[var]+=1
   Ds.append({std[(3+k,tuple(mon))]:Q(x) for k,x in enumerate(chi) if x})
 assert Ds==list(map(sp,d['lower_four_duals6']))
 QC=list(map(sp,d['Q_two_upper_quadratic_columns']))
 assert all(sum((a*q.get(i,Q(0)) for i,a in row.items()),Q(0))==0 for row in Ds for q in QC)
 cols=d['fixed_lower_dual4_minor_columns'];det=sum(((-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))*math.prod(Ds[i].get(cols[p[i]],Q(0)) for i in range(4)) for p in itertools.permutations(range(4))),Q(0))
 assert det==Q(d['fixed_lower_dual4_minor_determinant'])==-121
 # Freshly reconstruct positive envelopes only. All signed lower dots/rank
 # remain the future worker's authorized scope.
 env=[]
 for rr in d['remaining9_sources']:
  br,be=sp(rr['b1_r_selected99']),sp(rr['b1_ell_selected99'])
  for chi in [(0,1,4),(11,0,270)]:
   for mixed in (False,True):
    D,U=1,0
    for j in sorted(set(br)|set(be)):
     for k,c in enumerate(chi):
      terms=[abs(Q(c)*Hr[j].get(k,Q(0))*(be if mixed else br).get(j,Q(0)))]
      if mixed:terms.append(abs(Q(c)*He[j].get(k,Q(0))*br.get(j,Q(0))))
      for x in terms:
       DD=math.lcm(D,x.denominator);U=U*(DD//D)+x.numerator*(DD//x.denominator);D=DD
    env.append([str(D),str(U)])
 assert env==d['all36_lower_positive_D_U']
 DG=math.lcm(*(int(a) for a,b in env));UG=max(int(b)*(DG//int(a)) for a,b in env)
 assert str(DG)==d['common_integer_lower_D'] and str(UG)==d['common_integer_lower_abs_U']
 assert 5*UG.bit_length()+8==r['future_lower_rank_height_bits']==148
 assert d['source_prior_owned_bound']==2*d['all_input_numeric']+24000==74880
 assert d['source_owned_bound']==16309<74880<150000
 assert d['future_lower_value_owned_bound']==11408 and d['future_lower_value_RAW_bound']==2643762<32*1024*1024
 # Worker/custody is checked last, independently of arithmetic provenance.
 custody=read('work/a-zero-proof-v63/nonactive9-lower-prepared-custody-v2.json')
 wp=ROOT/custody['worker'];assert sha(wp)==custody['worker_sha256']
 assert custody['value_execution_started'] is False and custody['source_gate_sha256']==sha(ROOT/rp)
 text=wp.read_text();assert 'actual_lower4_by9_columns' in text and 'upper_components_not_zero_filled' in text
 assert 'temporary<=1153' in text and 'rank_bound<=384' in text
 assert 'source_actual_lower_values_evaluated' in text and 'common_integer_lower_D' in text
 inputs['literal_prepared_worker']={'path':custody['worker'],'sha256':sha(wp)}
 proof={
 'tensor_partition':'The exact transported resolution is active tensor K(z2,z3,z6). In F3 every nonactive block has PiY0 raw coefficient exactly linear, or identically zero: PiY0 vanishes on mixed F2 targets. This allorder source identity, together with the independently checked first array, proves162 sufficient zeros; no higher coefficient is assigned zero.',
 'actual_lower_schur':'For a projected-zero G3 constant combination c, after27 units the relation is (G3c,0). The99-unit elimination has b0=0, b1=B99 G3_1c and Hstar0=0. Its bottom quadratic is exactly -Hstar1*b1. Thus neither B99_1 nor unknown b2 contributes to that lower coefficient.',
 'gauge_scope':'The admissible comparison normalization is triangular on(top3,bottom3), so bottom coordinates are fixed. Unknown top quadratic coefficients and normalizer K*w2 contribute, in the normalized six-dimensional degree2 cover, only through the full upper image Q2. The four literal duals annihilate both Q columns. For162 with full plane first zero, all remaining unit derivative terms have a z factor and vanish by378 NF identities.',
 'finite_rank_scope':'Only four lower dual images are authorized. Their rank gives HF2 <= 6-rank; upper components, fullHF2, totalT2, torsion and nilpotence remain undecided.',
 'resource_scope':'Prior261 committed/784 temporary; integer4x9 rank148 is safely within384. Only known cached firsts and fixed basis memberships were used. No signed lower9 product, rank search, X2/X3/M1/physical/gauge query in this review.'}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 e=out/'independent171-source-lower9-frame-prior-and-worker-review.json'
 e.write_text(json.dumps({'all171_combinations_and_full102_plane_memberships_verified':True,'structural_zero_count':162,'remaining9':sorted(rest),'all99_Hstar_exact_match_independent_prior17':True,'all378_z_times_linear_NF_zero':True,'lower_four_Q_duals_and_minor':str(det),'all36_positive_envelopes_verified_without_signed_images':True,'prepared_worker_sha256':sha(wp),'proof':proof},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'source_gate_sha256':sha(ROOT/rp),'source_method_and_worker_custody_accepted':True,'actual_lower_values_or_rank_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{e.name:sha(e)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
 print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'worker_sha256':sha(wp),'no_signed_lower_products':True}))
if __name__=='__main__':main()
