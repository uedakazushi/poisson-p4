#!/usr/bin/env python3
"""Seven OLD-module socle-obstruction source functionals; NO new values.
All functional inputs are the already frozen37 actual e4*r² coefficients.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_TOLD_SOCLE_SEVEN_SOURCE_PREFLIGHT_20261105_v63';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 assert load('AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63_REVIEW')['all_passed']
 filt=load('AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63','actual-F2-source-products-and-m5-zero-proof.json')
 old=load('AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63');assert old['actual_T_length']==11 and old['actual_Tnil_length']==16
 src=load('AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63','phase1-actual-five-target-linear-relations.json');cert=load('AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63','actual-linear-quadratic-T-certificate.json')
 pure=load('AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63','actual-linear-syzygies-and-pure-r3-source-values.json');adj=load('AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_ADJOINT_20261104_v63','phase1-cubic-Pi2-and-source-adjoint.json')
 assert load('AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63_REVIEW')['all_passed']
 cache=load('AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63','phase3-composite-second-source-functional.json');assert load('AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW')['all_passed']
 formal=load('AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63','actual-minimal5-by825-Schur-operator.json')
 save('phase.json',dict(phase='all actual old source/cache and filtration peers pinned; no37-dot values'))
 c1=cert['complete_linear_rank_certificate'];R37=c1['selected_original_columns'];assert R37==pure['selected_actual_minimal_cycle_indices37']and len(R37)==37
 raw1=[Q.pairs(src['actual40_by657_linear_relation_columns'][i])for i in R37]
 def source_basis(raw,certificate,source_labels):
  norm=[Q.pairs(v)for v in certificate['basis_normalized_columns']];piv=certificate['selected_minor_rows'];T=[]
  for k,label in enumerate(source_labels):
   rec=certificate['fixed_elimination_receipts'][label];assert rec['column']==label
   z=raw[label]if len(raw)!=len(source_labels)else raw[k];t={label:F(1)}if len(raw)!=len(source_labels)else{k:F(1)}
   for j,c in rec['subtract']:
    z=Q.plus(z,Q.scaled(norm[j],-F(c)));t=Q.plus(t,Q.scaled(T[j],-F(c)))
   d=F(rec['divide']);assert d and Q.scaled(z,1/d)==norm[k]and rec['pivot_row']==piv[k]
   T.append(Q.scaled(t,1/d))
  return norm,T,piv
 b1,T1,p1=source_basis(raw1,c1,R37)
 def reduce_lift(p,basis,trans,piv,raw):
  z=p.copy();a={}
  # Frozen normalized columns are indexed by discovery, but their pivot
  # rows are the numerical leading rows.  Reduction follows pivot order.
  for j,row in sorted(enumerate(piv),key=lambda x:x[1]):
   c=z.get(row,F(0))
   if c:z=Q.plus(z,Q.scaled(basis[j],-c));a=Q.plus(a,Q.scaled(trans[j],c))
  assert Q.plus(z,*(Q.scaled(raw[j],c)for j,c in a.items()))==p
  return z,a
 lifts=[];other=[]
 for i in range(7):
  nf,a=reduce_lift({24+i:F(1)},b1,T1,p1,raw1);assert not nf,('zi_e3_not_in_linear_span',i,nf);lifts.append(a)
  nf,b=reduce_lift({32+i:F(1)},b1,T1,p1,raw1);assert set(nf)<={31,37};other.append(dict(zi=i,A_u=nf.get(31,F(0)),B_v=nf.get(37,F(0)),actual_linear_source_lift=b))
 mon2=list(itertools.combinations_with_replacement(range(8),2));ix={m:i for i,m in enumerate(mon2)};raw2=[]
 for col in raw1:
  for var in range(8):
   z={}
   for index,c in col.items():t,v=divmod(index,8);Q.put(z,36*t+ix[tuple(sorted((v,var)))],c)
   raw2.append(z)
 cq=cert['quadratic_certificate'];assert max(cq['selected_original_columns'])<296
 b2,T2,p2=source_basis(raw2,cq,cq['selected_original_columns']);assert len(b2)==177
 nf,a1=reduce_lift({4*36+ix[(5,5)]:F(1)},b2,T2,p2,raw2);assert not nf
 beta_r={j:a1[8*j+7]for j in range(37)if a1.get(8*j+7,F(0))}
 # Pi_cubic(g) is EXACT target4 r³ coefficient. Thus only these r source
 # coefficients and cached M2 target4 r² are relevant, not M3.
 assert adj['two_quotient_dual_rows600'][1]==[[599,'1']]
 assert pure['all178_known_quadratics_have_e5_r2_coefficient_zero']
 assert pure['first27_source55_r_linear_coefficients_all_zero']and pure['all84_unitcycle_source55_r_linear_coefficients_all_zero']
 assert pure['all37_source_values_from_PREEXISTING_108_root_cache']
 krows=pure['actual37_e5_r2_source_coefficients'];assert [x['minimal_cycle_index657']for x in krows]==R37
 # Source37 -> original768 and exact rawP8 coefficients have been frozen
 # independently. Retain these labels; do not replay source evaluations.
 c=F(11767,3652);functional37=[]
 for i in range(7):functional37.append(Q.plus(Q.scaled(lifts[i],-1),Q.scaled(lifts[5],-c*other[i]['A_u']),Q.scaled(beta_r,-c*other[i]['B_v'])))
 folded=[];needed=set();cachemap={i:(v,g)for i,v,g in cache['all108_composite_values_and_opaque_second_gradients']}
 assert cache['B11_monomial_order'][2]==[0,0,0,0,0,0,0,2]
 for f in functional37:
  z={}
  for j,a in f.items():
   for node,b in krows[j]['literal_rawP8_last_root_coefficients']:Q.put(z,node,a*F(b))
  z.pop(0,None);folded.append(z);needed.update(z)
 assert needed<=set(cachemap)and all(not cachemap[i][1]for i in needed)
 slots=sum(map(len,T1))+sum(map(len,T2))+sum(map(len,raw2))+sum(map(len,functional37))+sum(map(len,folded))+15000
 assert slots<150000
 sourcecount=sum(map(len,functional37));rootcount=sum(map(len,folded));bounds=dict(source_coefficient_slots=slots,seven_dot_multiplication_bound=min(sourcecount,rootcount),existing_root_count=len(needed),new_source_or_jet_coefficients=0,output_QQ_count=7)
 proof=dict(actual_old_minimal_target_sources=[22,23,29,30,55],actual_selected37_minimal_source_labels=R37,
  actual_selected37_original768_sources=[x['actual768_constant_source_combination']for x in krows],
  all7_zi_e3_constant_R37_lifts=lifts,all7_zi_e4_linear_NFs_and_source_lifts=other,
  all_exact_linear_lift_multiplybacks=True,z5_squared_e4_linear_coefficient_lift296=a1,z5_squared_e4_r_coefficients37=beta_r,exact_quadratic_lift_multiplyback=True,
  actual_scalar_input_field='PURE_R3.actual37_e5_r2_source_coefficients[*].actual_e5_r2_coefficient',
  scalar_input_scope='37 already evaluated actual M2(selected37) target4*rA²; no new evaluation or hiddenM3',
  alpha_formula='alpha_i=coefficient of oldg in zi*(e3r) = -lift(zi e3) dot k37, since r*q1=0, r*q2=g and oldm4=0',
  beta5_formula='beta5=coefficient of oldg in e4*z5² = -a1_r dot k37; M1R37*a1=e4*z5², and the higher M3*a1 is in oldm4=0; cubicdual is literal delta599',
  beta_i_formula='beta_i=A_i*alpha5+B_i*beta5 because zi e4=A_i*u+B_i*v+oldm² and z5*oldm²=0',
  q2_annihilation_proof='zi e4 lies in F2 by its linear NF with no e4*r term; r²F2=0 implies zi*q2=0. z5*q1 is a multiple of e4*z5² and z5*q2=qC, hence z5*oldm²=0 using oldm4=0.',
  seven_socle_obstruction_functionals_on_existing37=functional37,
  seven_equivalent_C0root_r_squared_functionals=folded,needed_cached_C0root_ids=sorted(needed),ordinary_opaque_second_sensitivities_allzero_in_existing108=True,
  same_S_and_r_tuple='old physical CA-F9 ell0, ordinary(rA,rA), B11 index2',
  abstract_source_t='u+c*v-alpha*q2, alpha from r*u=-c*q1+alpha*g is NOT evaluated; zi*q2=0 removes it from everyobstruction',
  sufficiency='ker(r:F2->F2) has dimension3 and basis q1,g,t. If m4Tnil=0, F3=m3Tnil has dimension3, is inside ker r, hence equals it; r-reduction then forces allzi to kill ker r. Any nonzero zi*t contradicts this, forcing HF55321. Seven zeros are only an inconclusive sufficient-test failure.',
  existing119_cache_scope='Uses the same actual37 source-value vector from the old119 cubic worker. No old source coefficient is recomputed. The7 contractions are new selected functionals, not declared equal to a saved119 syzygy value.',
  preflight_bounds=bounds,seven_obstruction_scalar_values_NOT_evaluated=True,no_917_or_new_source_jet_value=True)
 save('seven-old-socle-source-functionals-and-existing-cache-gate.json',proof)
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 rep=dict(id=ID,all_passed=True,frozen=True,all7_source_functionals_supported_by_existing37_cache=True,ordinary_opaque_second_sensitivities_zero=True,no_new_M2_M3_or_scalar_obstruction_evaluation=True,source_and_value_bounds=bounds,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20
 save('report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
