#!/usr/bin/env python3
"""Actual168 corrected d3 source combinations: an X2-only typed preflight.

Fixed saved elimination receipts only; no source coefficient evaluation.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_CORRECTED_D3_QUARTIC_PREFLIGHT_20261104_v63'
PATHS={
 'S':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'V':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'first640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'native128':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
 'A':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
 'm1':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
 'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
 'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
 'quad':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase2-certified-quadratic-input.json',
 'w0':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
 'resolution':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json',
 'boundary':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase1-corrected-d3-pure-r3-channel.json',
 'restricted':'reports/AMBIENT_B18_T_RESTRICTED118_QUARTIC_PREFLIGHT_20261104_v63/restricted118-quartic-source-operator-preflight.json',
 'restricted_values_peer':'reports/AMBIENT_B18_T_QUARTIC118_EXISTING_VALUES_20261104_v63_REVIEW/report.json',
 'Pi':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63_REVIEW/independent-cubic-witness71-source-and-channel.json'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(x):return {int(i):Q(c) for i,c in x if Q(c)}
def enc(x):return [[i,str(c)] for i,c in sorted(x.items()) if c]
def add(a,b,c=Q(1)):
 a=dict(a)
 for i,v in b.items():
  z=a.get(i,Q(0))+c*v
  if z:a[i]=z
  else:a.pop(i,None)
 return a
def put(a,k,c):
 z=a.get(k,Q(0))+c
 if z:a[k]=z
 else:a.pop(k,None)
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')

def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'corrected_d3_source_inputs_saved_before_fixed_receipt_replay','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert D['restricted_values_peer']['all_passed']
 S,X,A,L,W,R=[D[k] for k in ['S','X0','A','low','w0','resolution']]
 m2=list(itertools.combinations_with_replacement(range(8),2));i2={m:i for i,m in enumerate(m2)}
 source741=list(map(rec,S['source741_zero_constant_combinations']))
 source37=list(map(rec,A['actual_selected37_original768_constant_combinations']))
 star=[]
 for combo in S['second84_normalized_source741_combinations']:
  z={}
  for j,c in rec(combo).items():z=add(z,source741[j],c)
  star.append(z)
 first=D['first640']['source_columns640_rows56_sparse_dual9']+D['native128']['source_columns128_rows48_sparse_dual9']
 X1=[{8*r+k-1:Q(c) for r,p in col for k,c in p if k and Q(c)} for col in first]
 star1=[]
 for c in star:
  v={}
  for j,a in c.items():v=add(v,X1[j],a)
  star1.append(v)
 qnorm=list(map(rec,D['linear']['normalized84_unitcycle_Q1_columns_sparse']))
 qcols=list(map(rec,D['quad']['all464_certified_quadratic_columns_sparse']))
 assert len(qcols)==464
 # Reconstruct normalized177 basis SOURCE combinations from the saved
 # first296 elimination operations, without a fresh Gaussian/rank query.
 normsrc=[];basis=L['quadratic_certificate']['basis_normalized_columns']
 ops=L['quadratic_certificate']['fixed_elimination_receipts'][:296]
 for j,r in enumerate(ops):
  c={j:Q(1)}
  for p,a in r['subtract']:c=add(c,normsrc[p],-Q(a))
  if 'pivot_row' in r:
   c={i:a/Q(r['divide']) for i,a in c.items()};normsrc.append(c)
   z={}
   for i,a in c.items():z=add(z,qcols[i],a)
   assert z==rec(basis[len(normsrc)-1])
 assert len(normsrc)==177
 wb=list(map(rec,W['complete_rank_certificate']['basis_normalized_columns']));wl={min(v):j for j,v in enumerate(wb)}
 assert [enc(c) for c in wb]==S['actual_second_stage_U0_columns84']
 U=X['free_kernel_rows32'];pi=[{U[j]:c for j,c in rec(r).items()} for r in S['constant_free32_quotient_projection5_rows']]
 m0=[list(map(Q,r)) for r in D['V']['all24_by56_m0']]
 xp0=list(map(rec,D['m1']['actual_pivot_X0_columns27_sparse56']))
 receipts=[r for r in X['all768_fixed_elimination_receipts_and_constant_combinations'] if 'pivot' in r]
 pb={r['pivot']:rec(r['normalized']) for r in receipts};pindex={r['pivot']:j for j,r in enumerate(receipts)}
 actual37first=list(map(rec,A['actual_selected37_X1_sparse56_by8']))
 allcases=[];support=sum(map(len,star))+sum(map(len,star1))+sum(map(len,normsrc))+len(m0)*56
 for j,col in enumerate(R['differentials8_column_sparse'][2]):
  dlin={}
  for row,p in col:
   for mon,c in p:
    if sum(mon)==1:
     var=mon.index(1);assert var<7
     put(dlin,8*row+var,Q(c))
  b1={}
  for var in range(8):
   v={row:dlin[8*row+var] for row in range(84) if 8*row+var in dlin}
   while v:
    k=min(v);p=wl[k];c=v[k];put(b1,8*p+var,c);v=add(v,wb[p],-c)
  assert all(i%8!=7 for i in b1)
  q2={}
  for pv,c in b1.items():
   p,var=divmod(pv,8)
   for tv,a in qnorm[p].items():
    t,v=divmod(tv,8);put(q2,36*t+i2[tuple(sorted((var,v)))],-c*a)
  assert q2==qcols[296+j]
  assert all(not q2.get(36*t+35) for t in range(5))
  fac=rec(L['quadratic_certificate']['all_column_factorization_in_normalized_basis'][296+j]);a1={}
  for p,c in fac.items():a1=add(a1,normsrc[p],c)
  check={}
  for index,c in a1.items():check=add(check,qcols[index],c)
  assert check==q2
  c1=[{} for _ in range(8)];x2={}
  for index,c in a1.items():
   p,var=divmod(index,8);c1[var]=add(c1[var],source37[p],c)
   for rk,a in actual37first[p].items():
    row,k=divmod(rk,8);put(x2,36*row+i2[tuple(sorted((var,k)))],c*a)
  for index,c in b1.items():
   p,var=divmod(index,8);c1[var]=add(c1[var],star[p],c)
   for rk,a in star1[p].items():
    row,k=divmod(rk,8);put(x2,36*row+i2[tuple(sorted((var,k)))],c*a)
  # The full56 identity is stronger than checking the two formal kernels.
  free={36*i+k:x2.get(36*row+k,Q(0)) for i,row in enumerate(U) for k in range(36) if x2.get(36*row+k,Q(0))};a2={}
  for p in sorted(pb):
   for k in range(36):
    c=free.get(36*p+k,Q(0))
    if c:
     put(a2,36*pindex[p]+k,c);free=add(free,{36*i+k:a for i,a in pb[p].items()},-c)
  assert not free,('not_in_origin_pivot_span',j)
  full={}
  for index,c in a2.items():
   p,k=divmod(index,36)
   for row,a in xp0[p].items():put(full,36*row+k,c*a)
  assert full==x2,('full56_pivot_identity',j)
  for c in c1:
   for row in range(56):
    assert sum((Q(X['full_X0_56_by768'][row][s])*a for s,a in c.items() if s<640),Q(0))==0
    assert sum((Q(X['full_X0_56_by768'][row][s])*a for s,a in c.items() if s>=640),Q(0))==0
  support+=len(dlin)+len(b1)+len(q2)+len(a1)+len(a2)+len(x2)+sum(map(len,c1))
  assert support<150000
  allcases.append({'literal_d3_source':j,'actual_d3_linear84_times8':enc(dlin),
   'b1_U0_inverse_D3linear84_times8':enc(b1),'q_d3_2_sparse180':enc(q2),
   'a1_L37_lift_sparse296':enc(a1),'linear_original768_combinations_by8variables':list(map(enc,c1)),
   'known_x2_sparse56_times36':enc(x2),'actual_a2_sparse27_times36':enc(a2),
   'all56_x2_equals_Xp0_a2':True,'native_and_old_source_origin_parts_separately_zero':True})
 assert len(allcases)==168
 # Source-normalized extension to all known599 cubic coefficients.
 theta3=D['restricted']['Theta3_on_cubic600'];PiC=D['Pi']['canonical_Pi2_600_rows_independent'][0]
 assert rec(PiC).get(595)==1 and not rec(theta3['constant_part']).get(595) and not rec(theta3['canonical_qB_pure_cubic_coefficient_part']).get(595)
 assert not D['boundary']['nonzero_boundary_channel_indices']
 fields={'source_manifest':inputs,'all168':allcases,
  'actual_unitstar84_original768_combinations':list(map(enc,star)),
  'actual_unitstar84_first_sparse56_times8':list(map(enc,star1)),
  'actual_base':'QQ[[zA0,...,zA6,rA]]; lasttarget=e4=source55',
  'sign_convention':'corrected d3 residual has NEGATIVE the source functional below, modulo actual d3 times m',
  'typed_C3':'sum over variables c1[var]*var*actual_X2-Xp1*a2',
  'typed_beta2_rr':'r² coefficient of U0^-1*kappa_degree3(m0*C3)',
  'typed_Qproxy3':'pi*C3-Ystar1*U0^-1*kappa_degree3(m0*C3)',
  'typed_raw_G4':'sum c1[r]*rawP55,3-Xp55,2rr*a2rr-Ystar55,2rr*beta2rr',
  'typed_corrected168_quartic':'-typed_raw_G4+Theta4(typed_Qproxy3)',
  'Theta4_cubic600':{'constant_part':theta3['constant_part'],
   'times_b3':theta3['canonical_qB_pure_cubic_coefficient_part'],'times_c4':PiC},
  'b3_accepted':'-42215/3456',
  'c4_source_normalization':{'actual119_pivot':71,'divide':'-251099375/1051776',
   'formula':'restricted119_raw_G4[71]-Theta3(restricted119_Qproxy3[71]), divided by pivot',
   'canonical_qC_initial':'e4*zA5*rA²','Theta3_qC_zero':True},
  'm1_removed':'x0=x1=0 and full56 x2=Xp0*a2; z=x-Xp*a2 startsdegree3, so (mz)3=m0*C3',
  'Euler_boundary_control':'The corrected wedge plus kappa(mz) startsdegree2. Its exact d3 preimage startsm because d3_linear is injective on constant wedge³. Actual minimal d3 puredegree2/puredegree3 vanish, so d3*m has zero puredegree4. Theta4 kills their cubic initials because Theta2 annihilates d3_quadratic.',
  'no_nativeX3_or_commonCB4':'b1_r=0; native source55 identicallyzero; only c1_r old P3 is queried',
  'native_effective_S2_scope':'all168 c1[var] have separately zero native and old origins, checked full56; therefore accepted L2*S0 cancellation applies',
  'remaining_source_gate':'Contract Theta4 and C3, check requested X2 pairs and rawP3 labels against existing B11/common-source3 caches, keep unknown ordinary quotient coefficients symbolic until canceled. No evaluated168 values in this packet.',
  'whole_HF4_or_T_length_computed':False,'QQ_support_bound':support}
 ep=out/'actual168-corrected-d3-quartic-source-preflight.json';dump(out,ep.name,fields)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'all168_full56_source_pivot_identities_verified':True,
  'all168_old_and_native_constant_parts_separately_zero':True,'all168_b1_r_zero':True,
  'corrected_d3_quartic_X2_only_theorem_proved':True,'actual168_values_or_new_source_coefficients_computed':False,
  'actual_HF4_or_T_length_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},'QQ_support_bound':support,
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in out.iterdir())<32<<20
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)

if __name__=='__main__':main()
