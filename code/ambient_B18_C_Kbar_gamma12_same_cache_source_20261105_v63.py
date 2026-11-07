#!/usr/bin/env python3
"""Unpaired SOURCE: the two old filtered gamma dots use same37 cached t's."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
import ambient_B18_C_Kbar_r2_initial_source_20261105_v63 as I
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_GAMMA12_SAME_CACHE_SOURCE_20261105_v63'
PATHS={
 'initial':'reports/AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63/two-cached-alpha-reductions-and-projected-r2-seven-source-certificates.json',
 'initial_peer':'reports/AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63_REVIEW/report.json',
 'mu':'reports/AMBIENT_B18_C_KBAR_MU23_SOURCE_INVENTORY_20261105_v63/fourteen-z-actions-common37-M2-source-functionals-and-minimal-generator-criterion.json',
 'cache':'reports/AMBIENT_B18_C_KBAR_MU23_CACHE_SOURCE_20261105_v63/complete-old-Pg-second-cache-source-and-fourteen-action-positive-gate.json',
 'delta':'reports/AMBIENT_B18_C_KBAR_DELTA_SOURCE_INVENTORY_20261105_v63/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fs={'initial':['two_exact_linear_reductions_and_membership_cofactors','actual_alpha1_alpha2'],
 'mu':['normalized37_actual37_constant_lift','actual_selected37_sources'],
 'cache':['all37_Pg_source_weights']}
 d={k:S.selective(ROOT/p,fs[k],receipts)if k in fs else json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert d['initial_peer']['all_passed']
 T=[[I.parse(c)for c in row]for row in d['mu']['normalized37_actual37_constant_lift']];weights=[];rootweights=[]
 for record in d['initial']['two_exact_linear_reductions_and_membership_cofactors']:
  b=[Q(0)]*37
  for j,shift,c in record['exact_linear_source_cofactors']:
   assert shift==[0]*8;b[j]=I.addq(b[j],I.parse(c))
  w=[]
  for j in range(37):
   a=Q(0)
   for k in range(37):a=I.addq(a,I.mulq(b[k],T[k][j]))
   w.append(-a)
  weights.append(w);rw={}
  for a,row in zip(w,d['cache']['all37_Pg_source_weights']):
   for root,c in row:
    v=I.addq(rw.get(root,Q(0)),I.mulq(a,I.parse(c)))
    if v:rw[root]=v
    elif root in rw:del rw[root]
  rootweights.append(rw)
 evidence={'scope':'UNPAIRED SOURCE only; no gamma/C/F/Delta value. Reuses t37 alreadyidentifiedasrawPg,r²; doesnotchangeauthorizedmuworker.',
 'input_sha256':pins,'gamma1_gamma2_weights37':[[str(c)for c in row]for row in weights],
 'gamma1_gamma2_direct_given21_root_weights':[[[j,str(c)]for j,c in sorted(row.items())]for row in rootweights],
 'same_actual37_namespace':d['mu']['actual_selected37_sources'],
 'true_filtered_source':'rh1-alpha1*v=M1u1, rh2-rh3/4-alpha2*v=M1u2; module relationM1u=-M2u modoldm³. pi_q2 literal179, oldm³ tails do not enter gamma. Allcoefficients are fromtrue normalized→actual37 constantlifts, nothomogeneousM2.',
 'remaining_Delta_frontier':'C=-a3*gamma1-c3,F=gamma2-aS*gamma1-cS. Delta=lambdaU*F-lambdaS*(H+C), soafterthese2samecachelowerdots ONLYtrueG38cubiccoefficientHremains. NoH/primitiveM3cacheavailabilityclaim.',
 'signed_future_scope':'Any evaluation ofthese2dots requiresseparateROOTscope; currentpreparedmuworkercomputes37t+14actionsONLY andisunchanged.',
 'no_new_source_root_IFT_gauge_physical_jet':True,'max_pre_bits':I.MAXPRE,'selective_receipts':receipts}
 owned=S.count(d)+S.count(T)+S.count(weights)+S.count(evidence)+12000;assert owned<150000;evidence['whole_SOURCE_owned_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'two-unpaired-filtered-gamma-functionals-on-same-thirtyseven-cache-values.json';p.write_text(json.dumps(evidence,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_TWO_GAMMA_CACHE_ADJOINTS_UNPAIRED','input_sha256':pins,'code_sha256':sha(Path(__file__)),'evidence_sha256':{p.name:sha(p)},'whole_SOURCE_owned_upper':owned,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
