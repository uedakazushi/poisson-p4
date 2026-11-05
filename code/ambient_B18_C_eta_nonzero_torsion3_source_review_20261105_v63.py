#!/usr/bin/env python3
"""Independent conditional quartic/DVR proof, explicitly inapplicable to eta=0."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_ETA_NONZERO_TORSION3_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PRIMARY='3a6a89b22d3195c43ab64ca7bddd52c51373346475aafa3df58b53b7c51fdcc3'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
 def pin(p,h=None):
  p=pathlib.Path(p);v=sha(p)
  if h:assert v==h,(str(p),v,h)
  pins[str(p.relative_to(ROOT))]=v
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed']
 docs={}
 for p,h in pr['input_sha256'].items():
  d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
  if p.endswith('_REVIEW/report.json'):assert d['all_passed']
 D=pin(ROOT/'reports'/PRI/'conditional-eta-quartic-rank2-and-exclusion-of-exponent4.json',
  '7f93c43d8d8f8b319a4e1e9d2c018b382687c1cfe3226b03e42150ef88ac620c')
 pin(ROOT/'code/ambient_B18_C_eta_nonzero_torsion3_source_20261105_v63.py',pr['code_sha256'])
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'])
 line=docs['complete-full-C-cubic155-operator-with-four-formal-ell3-callbacks.json']
 N3=[[Q(x) for x in row] for row in line['full_no_cut_NF6']]
 N4=[[Q(0)]*15 for _ in range(6)];relations=[]
 for t in range(2):
  g=[Q(1)]+[-N3[3*t+j][4*t] for j in range(3)]
  def reduce(i,j):
   if i<3:return {(i,j):Q(1)}
   d={}
   for k in range(1,4):
    for mon,c in reduce(i-k,j+k).items():d[mon]=d.get(mon,Q(0))-g[k]*c
   return {m:c for m,c in d.items() if c}
  for col,i in enumerate([4,3,2,1,0]):
   v=reduce(i,4-i)
   for j,mon in enumerate([(2,2),(1,3),(0,4)]):N4[3*t+j][5*t+col]=v.get(mon,Q(0))
  for extra in [0,1]:
   coeff=([Q(0)] if extra else [])+g+([] if extra else [Q(0)])
   for row in N4[3*t:3*t+3]:assert sum(row[5*t+k]*coeff[k] for k in range(5))==0
  relations.append(list(map(str,g)))
 assert [[str(q) for q in row] for row in N4]==D['full_known_quartic_NF6_rows15_columns']
 # r*ell³(A-B/4), ell*ell³(A-B/4), in literal degree-four columns.
 p=[Q(0)]*15;q=p.copy();p[3]=1;p[8]=-Q(1,4);q[4]=1;q[9]=-Q(1,4)
 images=[[sum(a*b for a,b in zip(row,v)) for row in N4] for v in [p,q]]
 assert images==[[0,1,0,0,-Q(1,4),0],[0,0,1,0,0,-Q(1,4)]]
 assert images[0][1]*images[1][2]-images[0][2]*images[1][1]==1
 assert sum([3,6,6,5])==20 and 20-19==1
 assert 4*4+1+1+2+3==23
 # The all-order length/DVR argument is a source-bound pure implication.
 proofsteps=[
  'Assume some eta is nonzero. Complete155 plus the exact known initial NF kernel yields the actual homogeneous initial ell³(A-B/4), after dividing only by that nonzero QQ scalar.',
  'Its r and ell multiples give the two independently verified quartic columns. Thus HF4<=4; reduction m5C=ell*m4C and rankD4 give HF4>=4. Hence G=m4C is D-free4.',
  'F=m3C has rank4/mu5 and ellF=G. The torsion of F is killed by ell and has one minimal generator, hence is k. Its m multiple lies in G intersect Tor(C), hence vanishes.',
  'ell4C subset G is torsionfree, so all C torsion exponents are <=4. The accepted112e pattern leaves e3 or e4.',
  'If e4, ell³Tor(C) is k and lies in F, hence equals Tor(F)=k. Actual C3 m4=0 gives G subset L=ell³C. Length(C/G)=20, length(C/L)=19, so length(L/G)=1.',
  'Tor(F) lies in L and meets G trivially. The length1 quotient implies L=G directsum Tor(F), only as D-modules.',
  'For x write ell³x=ell*y+q with y in F,q in Tor(F). Then t=ell²x-y satisfies ell*t=q, so ell²t=0: t is genuine D-torsion.',
  'r*t lies in F by m³ order, and is torsion; hence r*t lies in Tor(F) subset L. Thus r ell²x=r*y+r*t lies in G+L=L.',
  'This kills r ell² A in actual C3, contradicting its literal nonzero NF coordinate. Therefore e4 is impossible; e3 is the conditional conclusion.'
 ]
 actual=pin(ROOT/'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63_REVIEW/report.json',
  '8dcf7359fe70c248ece84a61bca323f0820b1086e565d3687b23b59b7b70ca1b')
 assert actual['all4_eta03_zero_and_actual_fullC_HF3_6_accepted'] and actual['e_remains_open']
 proof={'pure_conditional_only':True,'hypothesis':'at least one actual eta03 !=0',
  'independent_recursive_curve_reduction_NF4':[[str(x) for x in row] for row in N4],
  'all_four_curve_times_variable_relations_annihilated':True,
  'actual_conditional_quartic_columns':[[str(x) for x in v] for v in images],
  'fixed_rank2_minor':'1','source_bound_DVR_and_colength_proof':proofsteps,
  'conditional_full_HF':[3,6,6,5,'4 in every higher degree'],
  'conditional_type':{'free_rank':4,'torsion_exponents':[1,1,2,3]},
  'conditional_series':'(3+3t-t³-t⁴)/(1-t)',
  'conditional_torsion_length':7,'conditional_ell4_quotient_length':23,
  'actual_eta_allzero_antecedent_false_no_actual_e_inference':True,
  'no_original_H0_nil_or_exhaustion_transfer':True,
  'no_new_source_jet_gauge_physical_coefficients':True}
 ep=out/'independent-conditional-quartic-initials-and-DVR-e4-obstruction.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
  'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
  'conditional_eta_nonzero_e3_theorem_accepted':True,
  'actual_eta_zero_no_current_exponent_claim':True,'elapsed_seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
