#!/usr/bin/env python3
"""Postverify provided finite105 corrections and exact C2 matrices; no solve."""
from fractions import Fraction as Q
from pathlib import Path
import collections,hashlib,json,resource,signal,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import add,scale,sp,enc,mv,mp,eye
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63';STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json','sourcevalues':f'reports/{A}/actual-source105-ell2-raw33-values.json',
 'C2':f'reports/{A}/actual-physical12-and-C2-right-exact-quotient.json',
 **{f'stage{i}':f'reports/{A}/stage{i}-source-residual-and-correction.json' for i in range(1,5)},
 'pref':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/whole-C2-source-and-budget-gate.json',
 'table':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json',
 'raw':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'raw_peer':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json',
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'J_peer':'reports/AMBIENT_B18_T_SOURCE4_ADJOINT9_20261104_v63_REVIEW/report.json',
 'oldpref':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_SOURCE_PREFLIGHT_20261104_v63/actual-Artin20-source-and-cost-preflight.json',
 'H20':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'oldsource':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-source105-Artin20-values.json',
 'oldactions':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-physical12-actions-and-C-special.json',
 'old_peer':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):
 (o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 assert sum(p.stat().st_size for p in o.iterdir())<32*1024*1024
class Algebra:
 def __init__(self,table):self.table=table;self.one={0:Q(1)};self.inverse_count=0
 def mul(self,a,b):
  out={}
  # Independent full table multiplication, with no degree shortcut.
  for i,c in a.items():
   for j,d in b.items():out=add(out,scale(self.table[i][j],c*d))
  return out
 def inv(self,p):
  c=p.get(0,Q(0));assert c
  h=dict(p);h.pop(0,None);h=scale(h,-1/c);v=t=dict(self.one)
  for _ in range(4):t=self.mul(t,h);v=add(v,t)
  assert not self.mul(t,h);v=scale(v,1/c);assert self.mul(p,v)==self.one;self.inverse_count+=1;return v
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'PROVIDED_FINITE_CORRECTIONS_AND_SOURCE_INPUTS_SAVED','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert all(D[k]['all_passed'] for k in ['primary','raw_peer','J_peer','old_peer'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 raw=D['raw'];basis=[tuple(m) for m in raw['actual33_basis']];assert len(basis)==33 and basis[0]==(0,0,0)
 actions=[[sp(p) for p in M] for M in raw['coordinate_actions9_sparse_columns']]
 table=[]
 for sector,a,e in basis:
  M=eye(33)
  for _ in range(a):M=mp(actions[7],M)
  for _ in range(e):M=mp(actions[8],M)
  if sector:M=mp(actions[sector-1],M)
  table.append(M)
 assert [[enc(p) for p in M] for M in table]==D['table']['table_columns_for_each_left_basis']
 degrees=[a+e+(sector!=0) for sector,a,e in basis];assert degrees==D['table']['basis_degrees'] and max(degrees)==4
 alg=Algebra(table);targets=[A0[0] for A0 in actions[:8]]+[add(actions[7][0],scale(actions[8][0],Q(1,5)))]
 assert list(map(enc,targets))==D['sourcevalues']['actual_CA_F9_targets']
 source=D['source'];nodes=source['DAG_scalar_nodes'];ids=D['joint']['all105_unknown_original_node_order'];eq=D['joint']['all105_equation_original_node_order']
 roots=D['pref']['source_requested164_roots'];assert len(roots)==len(set(roots))==164
 assert ids==D['sourcevalues']['joint105_unknown_order']==D['oldsource']['joint105_unknown_order']
 assert eq==D['sourcevalues']['joint105_equation_order']==D['oldsource']['joint105_equation_order']
 stops=set(ids);seen=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i)
  if i in stops:continue
  n=nodes[i];assert n[0] in ['a','m','n','i','c'],('UNBOUND_SOURCE',i,n)
  if n[0]!='c':todo.extend(n[1:])
 order=sorted(seen);assert order==D['pref']['source_closure_order']
 uses0=collections.Counter(roots)
 for i in order:
  if i not in stops and nodes[i][0]!='c':uses0.update(nodes[i][1:])
 peak=0
 def evaluate(phi):
  nonlocal peak
  uses=uses0.copy();v={};iv=[]
  for i in order:
   n=nodes[i];op=n[0]
   if i in phi:p=phi[i]
   elif op=='c':p={0:Q(n[1],n[2])} if n[1] else {}
   elif op=='a':p=add(v[n[1]],v[n[2]])
   elif op=='m':p=alg.mul(v[n[1]],v[n[2]])
   elif op=='n':p=scale(v[n[1]],-1)
   elif op=='i':p=alg.inv(v[n[1]]);iv.append([i,n[1],enc(v[n[1]]),enc(p)])
   else:raise AssertionError(n)
   v[i]=p;peak=max(peak,sum(map(len,v.values())));assert peak<150000
   if i not in phi:
    for j in n[1:] if op!='c' else []:
     uses[j]-=1
     if not uses[j]:del v[j]
  return {i:v[i] for i in roots},iv
 def residual(values):
  out=[values[i] for i in eq]
  for j,p in enumerate(targets):out[48+j]=add(out[48+j],scale(p,-1))
  return out
 jac=[{j:Q(c) for j,c in enumerate(row) if Q(c)} for row in D['joint']['fixed105_Jacobian']]
 phi={i:{} for i in ids};stages=[]
 for stage in range(1,5):
  values,iv=evaluate(phi);res=residual(values);given=D[f'stage{stage}'];corr=list(map(sp,given['correction105']))
  assert list(map(enc,res))==given['residual105']
  assert all(add(*(scale(corr[j],c) for j,c in row.items()))==scale(p,-1) for row,p in zip(jac,res))
  minimum=min((degrees[j] for p in res for j in p),default=5);assert minimum>=stage and minimum==given['minimum_residual_order']
  phi={i:add(phi[i],p) for i,p in zip(ids,corr)}
  stages.append({'stage':stage,'all105_J0_times_provided_correction_equals_minus_residual':True,'residual_minimum_order':minimum})
 final,iv=evaluate(phi);assert all(not p for p in residual(final))
 assert list(map(enc,[phi[i] for i in ids]))==D['sourcevalues']['all105_unknown_final_values']
 assert [[i,enc(final[i])] for i in roots]==D['sourcevalues']['all_requested164_root_values']
 assert iv==D['sourcevalues']['all6_unit_inverse_products_exact'] and len(iv)==6
 assert [enc(phi[i]) for i in range(2,11)]==D['sourcevalues']['physical_R0_readback']
 dump(o,'phase1-independent-provided105-and164-source-checks.json',{'all_four_provided_corrections_checked_without_solving':stages,
  'all105_final_equations_zero':True,'all164_actual_source_roots_and6_inverse_products_verified':True,'source_live_QQ_peak':peak})
 # Reconstruct both exact ell0 projection interfaces from monomial actions.
 H=D['H20'];oldA=[[sp(p) for p in A0] for A0 in H['multiplication8_column_sparse']];Q20=[]
 for sector,a,e in basis:
  p={} if e else {0:Q(1)}
  for _ in range(a):p=mv(oldA[7],p)
  if sector:p=mv(oldA[sector-1],p)
  Q20.append(p)
 assert [enc(p) for p in Q20]==D['C2']['raw33_to_H20_ell0_columns']
 assert all(mp(Q20,A0)==mp(B,Q20) for A0,B in zip(actions,oldA+[[{} for _ in range(20)]]))
 assert all(mv(Q20,phi[i])==sp(p) for i,p in zip(ids,D['oldsource']['all105_unknown_final_values']))
 oldroots={i:sp(p) for i,p in D['oldsource']['all_requested_original_source_root_values']}
 assert all(mv(Q20,final[i])==oldroots[i] for i in roots)
 branches=raw['canonical_branch_quotients6'];Qmap=[[sp(p) for p in b['quotient_map6_by33_columns']] for b in branches]
 b6=[tuple(m) for m in branches[0]['basis_plane_monomials']];assert b6==[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1)]
 idx={m:i for i,m in enumerate(basis)};lift=[{idx[0,a,e]:Q(1)} for a,e in b6]
 assert all(mv(Q0,p)=={i:Q(1)} for Q0 in Qmap for i,p in enumerate(lift))
 IA=[final[i] for row,i in H['source_conductor_generators_rows_and_ids']['CA']];IB=[final[i] for row,i in H['source_conductor_generators_rows_and_ids']['CB']]
 assert all(not mv(Qmap[0],p) for p in IA) and all(not mv(Qmap[1],p) for p in IB)
 assert all(not alg.mul(a,b) for a in IA for b in IB)
 nuA,nuB=[final[i] for i in D['oldpref']['nu_roots_CA_CB']]
 def block(a,b,c):return {**a,**{i+6:v for i,v in b.items()},**{i+12:v for i,v in c.items()}}
 physical=[];physical_data=[]
 for root in D['oldpref']['physical12_coefficient_roots']:
  aA,aB,bA,bB=[final[root[k]] for k in ['aA_node','aB_node','bA_node','bB_node']];delta=add(aA,scale(aB,-1))
  assert not mv(Qmap[0],bA) and not mv(Qmap[1],bB)
  assert all(not alg.mul(a,bB) for a in IA) and all(not alg.mul(b,bA) for b in IB)
  M=[]
  for p in lift:M.append(block(mv(Qmap[0],alg.mul(delta,p)),{},scale(alg.mul(bB,p),-1)))
  for p in lift:M.append(block({},mv(Qmap[1],alg.mul(delta,p)),alg.mul(bA,p)))
  for i in range(33):
   p={i:Q(1)};M.append(block(scale(mv(Qmap[0],alg.mul(alg.mul(nuB,bB),p)),-1),mv(Qmap[1],alg.mul(alg.mul(nuA,bA),p)),alg.mul(delta,p)))
  physical.append(M);physical_data.append({'physical_row':root['physical_row'],'roots':root,'delta':enc(delta),'bA':enc(bA),'bB':enc(bB)})
 assert physical_data==D['C2']['physical12_coefficients']
 assert [[enc(p) for p in M] for M in physical]==D['C2']['physical_actions12_sparse_columns']
 rawA=[[sp(p) for p in A0] for A0 in raw['raw45_coordinate_actions9_sparse_columns']]
 assert all(mp(A0,B)==mp(B,A0) for A0 in physical for B in rawA)
 assert all(mp(A0,B)==mp(B,A0) for i,A0 in enumerate(physical) for B in physical[i+1:])
 Q26=[{a:Q(1)} if not e else {} for a,e in b6]+[{a+3:Q(1)} if not e else {} for a,e in b6]+[{i+6:c for i,c in p.items()} for p in Q20]
 assert [enc(p) for p in Q26]==D['C2']['raw45_to_raw26_ell0_columns']
 oldF=[[sp(p) for p in M] for M in D['oldactions']['physical_actions12_column_sparse']]
 assert all(mp(Q26,A0)==mp(B,Q26) for A0,B in zip(physical,oldF))
 # Fixed rank31 receipts, quotient14 independence, and all induced actions.
 columns=[p for M in physical for p in M];assert [enc(p) for p in columns]==D['C2']['all540_relation_columns']
 piv={}
 for receipt in D['C2']['fixed_rank_receipts']:
  p=columns[receipt['source_column']]
  for k,c in receipt['subtract']:p=add(p,scale(piv[k],-Q(c)))
  k=receipt['pivot'];assert k not in piv and p.get(k)==Q(receipt['divide_by']) and all(i>=k for i in p)
  piv[k]=scale(p,1/Q(receipt['divide_by']));assert enc(piv[k])==receipt['normalized']
 assert len(piv)==31
 free=D['C2']['C2_basis_raw45_indices'];assert free==[i for i in range(45) if i not in piv] and len(free)==14
 projection=list(map(sp,D['C2']['all45_projection_columns']))
 assert all(projection[i]=={j:Q(1)} for j,i in enumerate(free)) and all(not mv(projection,p) for p in columns)
 Cactions=[[mv(projection,A0[j]) for j in free] for A0 in rawA]
 assert [[enc(p) for p in M] for M in Cactions]==D['C2']['C2_coordinate_actions9_sparse_columns']
 assert all(mp(projection,M)==mp(N,projection) for M,N in zip(rawA,Cactions))
 assert all(mp(M,N)==mp(N,M) for M in Cactions for N in Cactions)
 E=Cactions[8];assert all(not p for p in mp(E,E))
 epcols=[0,1,5,6,10,11];eprows=[3,4,8,9,12,13]
 assert [{i:p.get(r,Q(0)) for i,r in enumerate(eprows) if p.get(r,Q(0))} for p in [E[j] for j in epcols]]==eye(6)
 assert all(set(p)<=set(eprows) for p in E)
 ep=o/'independent-four105-corrections-physical12-and-C2-14-proof.json'
 dump(o,ep.name,{'inputs':inputs,'all_four_J0_provided_correction_checks':stages,'all105_final_equations_and6_unit_inverses_checked':True,
  'all164_roots_and_rawphysical9_readbacks_verified':True,'all105_and164_ell0_projection_checks':True,
  'all64_conductor_products_and12_branch_lift_annihilator_gates':True,'all108_scalar_intertwiners66_physical_commutators_and12_ell0_actions':True,
  'all540columns_and_fixed31_rank_receipts_verified':True,'actual_C2_dimension':14,'C2_basis_raw45_indices':free,
  'coordinate_order':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell'],
  'all9_source_bound_C2_actions14_verified':True,'ellrank6_fixed_unit_rows':eprows,'ellrank6_fixed_unit_columns':epcols,
  'C2_as_QQell_modell2_module':{'cyclic_length2':6,'cyclic_length1':2},
  'right_exact_scope':'Coker of the540 actual physical image columns gives C/ell²C; no kernel specialization used',
  'no_unknown_solve_inverse_Jacobian_rerun_T2_H0_or_nil_exponent_inference':True})
 note=ROOT/'research/notes'/f'{STEM}.md';report={'id':STEM,'all_passed':True,
  'provided_four105_source_corrections_and164_roots_independently_accepted':True,
  'actual_C2_dimension':14,'actual_physical_relation_rank':31,'actual_ell_action_rank':6,
  'C2_module_over_QQell_modell2':'A^6 + QQ^2','actual_T2_H0_nil_exponent_computed':False,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
