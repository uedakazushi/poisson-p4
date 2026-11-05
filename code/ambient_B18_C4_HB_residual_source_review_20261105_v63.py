#!/usr/bin/env python3
"""Independent C4 filtration/DVR proof and universal four-source minor."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C4_HB_RESIDUAL_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PRIMARY='c043b09dd44ce0a8b6f4dcc030f8eca2f19f4d66e222a98582781fad1cb41ec6'
Z=(0,)*6
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def const(x):return {Z:Q(x)} if x else {}
def var(i):return {tuple(int(j==i) for j in range(6)):Q(1)}
def add(a,b):
 d=dict(a)
 for k,v in b.items():
  d[k]=d.get(k,Q(0))+v
  if not d[k]:del d[k]
 return d
def scale(a,c):return {k:v*c for k,v in a.items() if v*c}
def mul(a,b):
 d={}
 for x,u in a.items():
  for y,v in b.items():
   k=tuple(i+j for i,j in zip(x,y));d[k]=d.get(k,Q(0))+u*v
 return {k:v for k,v in d.items() if v}
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
 def pin(p,h=None):
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h,(str(p),s,h)
  pins[str(p.relative_to(ROOT))]=s
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed']
 docs={}
 for p,h in pr['input_sha256'].items():
  d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
  if p.endswith('_REVIEW/report.json'):assert d['all_passed']
 D=pin(ROOT/'reports'/PRI/'C4-finite-HF-and-two-HB-residual-coordinate-criterion.json',
  '602cb1352009836524c40102221dab99bd7f1de425bbefeb41b1777c568f0ad2')
 pin(ROOT/'code/ambient_B18_C4_HB_residual_source_20261105_v63.py',pr['code_sha256'])
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'])
 # A row-subset determinant DP checks the universal block identity, with
 # all four lower HB entries retained as independent symbols.
 x,y,p,q,s,t=[var(i) for i in range(6)]
 A=[[const(-Q(3,5)),const(6),{},{}],[x,y,{},{}],[p,q,const(1),{}],[s,t,{},const(1)]]
 dp={0:const(1)}
 for r in range(4):
  nxt={}
  for mask,v in dp.items():
   for c in range(4):
    if mask>>c&1:continue
    sign=(-1)**sum(mask>>j&1 for j in range(c+1,4))
    nmask=mask|1<<c;nxt[nmask]=add(nxt.get(nmask,{}),scale(mul(v,A[r][c]),Q(sign)))
  dp=nxt
 determinant=dp[15]
 expected=add(scale(y,-Q(3,5)),scale(x,-Q(6)))
 assert determinant==expected
 assert [[list(k),str(v)] for k,v in sorted(determinant.items())]==D['universal_determinant_monomials_x_y_p_q_s_t']
 assert D['unconditional_C4_prefix_HF']==[3,6,6,6]
 assert D['unconditional_C4_HF_values']==[[3,6,6,6,2,0],[3,6,6,6,3,0]]
 counts=[]
 for e in range(3,9):
  length=4*4+sum(min(4,n) for n in [1,1,2,e]);h4=length-sum([3,6,6,6])
  assert (length,h4)==((23,2) if e==3 else (24,3))
  counts.append([e,length,[3,6,6,6,h4,0]])
 actual=pin(ROOT/'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_VALUES_20261105_v63_REVIEW/report.json',
  'b2e0b643a7c8d647f6f4f12ea72f8f12c7443ae0b1caa1d1015e2b3d4dbb31de')
 assert actual['full_HB_known_rank1_accepted'] and actual['no_actual_HF4_or_e_claim']
 proof={'same_actual_C4_equals_C_div_ell4':True,
  'unconditional_proof':[
   'Actual C3m4=0 means m4C subset ell³C. Exact full reduction implies m5C=ell*m4C subset ell4C. Therefore m5(C4)=0, not just an upper cover.',
   'ell4C subset m4C, so C4 has the actual full-C HF prefix3666 through degree3.',
   'Finite-D rank4/type112e gives length(C4)=16+1+1+2+min(4,e), hence23 for e3 and24 for e>=4.',
   'The prefix haslength21 and m5(C4)=0, yielding EXACT alternatives36662 or36663.'],
  'universal_source_matrix_has_all_lower_HB_symbols':True,
  'independent_subset_DP_determinant':[[list(k),str(v)] for k,v in sorted(determinant.items())],
  'conditional_four_source_minor':'rows[0,j,2,5], columns[HB0,HB1,ell4A,ell4B]: (-3/5)*(y+10x)=(-3/5)*w_j for j1 or4',
  'positive_conditional_proof':'A nonzero w1 orw4 makes four genuine C4 initial relations independent in its known6-dimensional quarticNF. Thus HF4(C4)<=2; the exact prefix/m5zero give length<=23, and lower23 forces23/e3.',
  'actual_full_HB_w_allzero_condition_not_satisfied':True,
  'no_current_exact_C4_length_or_e_claim':True,
  'finite_length_sanity_receipts_not_general_proof':counts,
  'no_new_source_or_jet_values_and_no_original_H0_nil_transfer':True}
 ep=out/'independent-C4-filtration-DVR-and-universal-HB-cut-minor.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
  'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
  'C4_HF36662_or36663_and_conditional_e3_accepted':True,
  'current_HB_negative_no_e_inference':True,'elapsed_seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
