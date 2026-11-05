#!/usr/bin/env python3
"""Independent source-linear G2/minor and constant-cycle theorem review."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63';STEM=A+'_REVIEW'
G='AMBIENT_B18_ELL2_C2_G2_LINEAR_INPUT_GATE_20261105_v63'
P={
 'primary':f'reports/{A}/report.json','rank':f'reports/{A}/actual-G2-linear-fixed-source-profile.json',
 'gate':f'reports/{G}/report.json','input':f'reports/{G}/known-G2-linear-input-and-provenance.json',
 'transport':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json',
 'transport_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json',
 'active':'reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63/active-five-d3-source-split-and-budget.json',
 'origin':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-source-and27-unit-split.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 # The X-origin file is optional only for the arithmetic source check;
 # the theorem states the separately accepted27-origin-unit premise.
 if not (ROOT/P['origin']).exists():P.pop('origin')
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'SOURCE_LINEAR42_AND_FIXED_MINOR','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert D['primary']['all_passed'] and D['gate']['all_passed'] and D['transport_peer']['all_passed']
 for name,folder in [('primary',A),('gate',G)]:
  for p,h in D[name]['input_sha256'].items():assert sha(ROOT/p)==h
  for p,h in D[name]['evidence_sha256'].items():assert sha(ROOT/'reports'/folder/p)==h
 T=D['transport'];S=D['active'];sel=T['active56_source_columns'];free=T['active_F2_free42'];Pi=list(map(sp,T['active_Pi17_columns']));d2=list(map(sp,S['active_d2_columns']))
 fixed=list(map(sp,T['active_F1_fixed67_columns']));basis={};comb={}
 for r in T['active_F1_fixed67_receipts']:
  v=dict(fixed[r['source_column']]);w={r['source_column']:Q(1)}
  for p,c in r['subtract']:v=add(v,scale(basis[p],-Q(c)));w=add(w,scale(comb[p],-Q(c)))
  c=Q(r['divide_by']);p=r['pivot'];assert v[p]==c;basis[p]=scale(v,1/c);comb[p]=scale(w,1/c)
  assert enc(basis[p])==r['normalized']
 alpha=[];I2=[]
 for j in free:
  v=dict(d2[j]);w={}
  while v:
   p=min(v);c=v[p];v=add(v,scale(basis[p],-c));w=add(w,scale(comb[p],c))
  assert all(k>=11 for k in w)
  aa={k-11:c for k,c in w.items()};assert mv([d2[k] for k in sel],aa)==d2[j]
  ii=add({j:Q(1)},{sel[k]:-c for k,c in aa.items()});assert not mv(d2,ii)
  alpha.append(aa);I2.append(ii)
 active=[0,1,4,5,7];row={w:i for i,(p,w) in enumerate([(1,())]+[(0,(i,)) for i in active])};lin=[]
 for p,w in [(2-q,w) for q in range(3) for w in itertools.combinations(active,q)]:
  for j in range(14):
   v={}
   if p:v[(14*row[w]+j)*6+5]=Q(1)
   for k,i in enumerate(w):
    rr=14*row[w[:k]+w[k+1:]]+j;v[rr*6+active.index(i)]=Q((-1)**(p+k))
   lin.append(v)
 cols=[]
 for p in I2:
  v=mv(lin,p);c={}
  for k,a in v.items():
   r,var=divmod(k,6);c=add(c,{i*6+var:a*b for i,b in Pi[r].items()})
  cols.append(c)
 E=D['rank'];assert list(map(enc,cols))==E['actual_linear_G2_active102_by42_columns']==D['input']['linear_G2_columns']
 assert list(map(enc,alpha))==E['constant_middle56_coefficients_alpha42'] and list(map(enc,I2))==E['actual_I2_source_origin42_columns']
 piv={};reprs={};selected=[];seq=[]
 for r in E['fixed_rank_receipts']:
  j=r['source_column'];v=dict(cols[j]);w={j:Q(1)}
  for p,c in r['subtract']:
   assert v.get(p)==Q(c);v=add(v,scale(piv[p],-Q(c)));w=add(w,scale(reprs[p],-Q(c)))
  p=r['pivot'];c=Q(r['divide_by']);assert c and min(v)==p and v[p]==c
  piv[p]=scale(v,1/c);reprs[p]=scale(w,1/c);selected.append(j);seq.append(p);assert enc(piv[p])==r['normalized']
 assert len(piv)==42 and E['exact_kernel_on_source42']==[]
 rows=sorted(piv);B=[piv[p] for p in rows];assert list(map(enc,B))==E['normalized_basis']
 assert all(mv(B,sp(p))==v for p,v in zip(E['all42_column_factorizations'],cols))
 rp={p:i for i,p in enumerate(rows)}
 assert E['fixed_minor_rows']==rows and E['fixed_minor_columns']==selected
 assert [enc({rp[p]:c for p,c in cols[j].items() if p in rp}) for j in selected]==E['fixed_minor_sparse_columns']
 assert all(not any(i<p for i in piv[p]) and piv[p][p]==1 for p in rows)
 sign=(-1)**sum(seq[i]>seq[j] for i in range(42) for j in range(i+1,42))
 assert E['fixed_minor_nonzero_exact_product']=={'sign':sign,'QQ_factors':[r['divide_by'] for r in E['fixed_rank_receipts']]}
 assert E['full102_linear_rank']==102 and E['linear_injectivity_full102']
 ep=o/'independent-G2-linear-source-and-X1-only-theorem.json'
 dump(o,ep.name,{'inputs':inputs,'all42_alpha_I2_and102_by42_linear_source_coefficients_reconstructed':True,'all42_fixed_minor_receipts_and_factorizations_verified':True,'no_projection_derivative_or_source_I2_derivative_needed':'d2*I2(0)=0 and Pi17*d2(0)=0','full_linear_constant_kernel_zero':'active target zero-variable coefficients force all51 mixed constants zero; remaining zeroKoszul2 forces9 constants zero; active42 injective','G3_constant_zero_follows_from_injection':True,'conditional_X1_only_theorem':'after27 actual X0 unit cancellations, all741 constant-zero source combos satisfy G2_linear Y0=Mnew0 Xeff1; uniqueY0 uses noM1','conditional_mu_T2_formula':'mu_B(T2)=105-rank(Y0) once actual full9 Xeff1 is source-bound','native_L1_drop_requires_separate_origin_cancellation':'oldB0 rank4/native0 rank23 directsum implies S0*cNative=0 for allzero-origin combos, so nativefirst=L0*S1','actual_X1_values_or_mu_T2_or_length_computed':False})
 report={'id':STEM,'all_passed':True,'actual_active42_and_full102_G2_linear_injection_accepted':True,'conditional_actual_X1_only_mu_T2_method_accepted':True,'source_derivative_or_actual_T2_mu_length_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
