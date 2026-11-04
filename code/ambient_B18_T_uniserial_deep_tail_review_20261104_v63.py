#!/usr/bin/env python3
"""Independent reverse-order fixed divisions and uniserial-module proof gate."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
R=Path(__file__).resolve().parents[1];S='AMBIENT_B18_T_UNISERIAL_DEEP_TAIL_20261104_v63_REVIEW';P='AMBIENT_B18_T_UNISERIAL_DEEP_TAIL_20261104_v63';O=R/'reports'/S
FIXED={
 'AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63_REVIEW':'5bedafd8efeed3dee33bc43321a78d4caf5a83869ee033230c38b2ed070ab0bd',
 'AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63_REVIEW':'f830eaa73ce76134bcdd9bec9170bad1e3378550f86711ce320b68449ec4ed41',
 'AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW':'143af8f9218d2e43626f58db9abd34827dee9e419f6ec4a9626b5ca5d36caae1',
 'AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW':'ca127c4dc0f3c12a4d4d80dfc1ca316f42dee8b90ecc74cb3e5e423752c906b7',
 'AMBIENT_B18_CONDUCTOR_SQUARE_CLOSED348_20261104_v63_REVIEW':'b937f7f936645761081ffd4e5b8663d095806f19ce14a771b5ae0e6b204a4108'}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def save(n,x):(O/n).write_text(json.dumps(x,indent=2)+'\n')
def decode(g):return {(p,tuple(m)):F(c) for p,m,c in g if F(c)}
def add(p,key,c):
 v=p.get(key,F(0))+c
 if v:p[key]=v
 else:p.pop(key,None)
def enc(p):return [[i,list(m),str(c)] for (i,m),c in sorted(p.items())]
def lead(f):return max(f,key=lambda x:(-x[0],x[1]))
def mons(d):
 for xs in itertools.combinations_with_replacement(range(8),d):yield tuple(xs.count(i) for i in range(8))
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));O.mkdir(exist_ok=True);pins={}
 def load(st,name='report.json'):
  p=R/'reports'/st/name;pins[str(p.relative_to(R))]=sha(p);return json.loads(p.read_text())
 primary=load(P);assert pins[f'reports/{P}/report.json']=='a868358ec3e276e6a1577225cb63c67f2e24e258b6c81d9e91700b7d9646e498'
 for p,h in primary['input_sha256'].items():assert sha(R/p)==h
 for st,h in FIXED.items():assert load(st)['all_passed'] and pins[f'reports/{st}/report.json']==h
 data=load('AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63','known-initial40-and-length-upper17.json');g=list(map(decode,data['homogeneous_module_generators40']));assert len(g)==40
 assert data['actual_m10_T_zero_from_accepted_source_base'] and data['standard_degree_at_least3']==['e4*r^d']
 source=load(P,'seven-reductions-and-uniserial-tail-proof.json')
 lt=[lead(q) for q in g]
 for q in g:assert len({sum(m) for p,m in q})==1
 # Fixed source divisions in REVERSE generator order, independently of author's forward trace.
 def divide(f):
  todo=dict(f);rem={};trace=[];steps=0
  while todo:
   key=lead(todo);p,m=key;c=todo[key];divisor=None
   for j in range(39,-1,-1):
    pp,mm=lt[j]
    if p==pp and all(a>=b for a,b in zip(m,mm)):divisor=j;break
   if divisor is None:rem[key]=c;todo.pop(key);continue
   j=divisor;pm=tuple(a-b for a,b in zip(m,lt[j][1]));a=c/g[j][lt[j]]
   for (pp,mm),v in g[j].items():add(todo,(pp,tuple(x+y for x,y in zip(mm,pm))),-a*v)
   assert not todo or (-lead(todo)[0],lead(todo)[1])<(-p,m)
   trace.append([j,list(pm),str(a)]);steps+=1;assert steps<10000
  return rem,trace
 def combine(trace):
  back={}
  for j,m,c in trace:
   for (p,n),v in g[j].items():add(back,(p,tuple(a+b for a,b in zip(m,n))),F(c)*v)
  return back
 checks=[]
 for i in range(7):
  m=[0]*8;m[i]=1;m[7]=3;f={(4,tuple(m)):F(1)};rem,tr=divide(f);assert not rem and combine(tr)==f
  old=source['all7_zi_e4_r3_known_initial_zero_reductions'][i];assert old['variable_index']==i and combine(old['exact_known40_cofactors'])==f
  checks.append({'variable':i,'input':enc(f),'reverse_order_cofactors':tr,'literal_source_identity':True})
 standard={}
 for d in range(5):
  result=[]
  for p in range(5):
   for m in mons(d):
    if not any(p==pp and all(a>=b for a,b in zip(m,mm)) for pp,mm in lt):result.append([p,list(m)])
  standard[d]=result
 assert list(map(len,standard.values()))==[5,3,2,1,1]
 assert standard[3]==[[4,[0,0,0,0,0,0,0,3]]] and standard[4]==[[4,[0,0,0,0,0,0,0,4]]]
 for d in range(3,10):
  f={(4,(0,)*7+(d,)):F(1)};rem,tr=divide(f);assert rem==f and not tr
 # Logic: accepted initial containment, not equality, suffices for all conclusions.
 proof={'all7_literal_cofactor_identities_independent_reverse_division':checks,'standard_monomials_degrees0_4':standard,'pure_r3_to_r9_standard_in_known_upper_quotient_only':True,'M':'m^3*T','g':'rA^3*e4, last minimal5 target/sourceV56 row55','M_over_S_generated_by_g':'gr3(T) is at most kg, hence M=Sg+mM; Nakayama gives M=Sg, including zero case.','mM_equals_rM':'gr4(T) is at most k*rg, hence mM=rM+m^2M. Apply Nakayama to finite mM/rM.','M_over_k_r_generated_by_g':'M=kg+rM and finite length iterate; no rational integration or DVRell identification.','finite_cyclic_PID_classification':'M is k[[rA]]/(rA^s), s>=0; s<=7 since r^7M=m^10T=0.','all_S_submodules_linearly_ordered':'Every S-submodule is k[[rA]]-stable; all submodules of a cyclic truncated k[[rA]] module are r^jM.','hi_unique_and_in_r2':'zi*g has zero degree4 initial class by literal7identities, hence lies m^5T=m^2M=r^2M. Every zi commutes with r, so acts by unique hi in k[[r]]/(r^s).','actual_HF_form':'5,3,2 followed by s consecutive ones then zero, 0<=s<=7','actual_length_form':'10+s','no_known_upper_module_full_initial_equality_assumption':True,'no_actual_s_hi_H0_DVRell_exponent_or_nilcube_claim':True}
 save('independent-seven-reductions-and-module-proof.json',proof)
 report={'id':S,'all_passed':True,'all7_fixed_actual_known40_cofactor_identities':True,'independent_reverse_generator_order_division':True,'known40_standard_counts_degrees0_4':[5,3,2,1,1],'actual_tail':'m^3T is k[[rA]]/(rA^s), 0<=s<=7','zi_actions_in_rA_squared':True,'actual_HF_form':'5,3,2, then s consecutive ones, then zero','actual_total_length_form':'10+s','rA_not_DVRell':True,'actual_s_or_H0_nilcube_not_computed':True,'no_producer_import_rerun_source_jet_or_GB_search':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(R/'research/notes'/f'{S}.md'),'evidence_sha256':{'independent-seven-reductions-and-module-proof.json':sha(O/'independent-seven-reductions-and-module-proof.json')},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};save('report.json',report);print(json.dumps(report))
if __name__=='__main__':main()
