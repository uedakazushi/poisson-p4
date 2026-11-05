#!/usr/bin/env python3
"""Independent primitive D-free complement, genuine T7 source and unrun gate."""
import ast,hashlib,json,math,pathlib,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import nc
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_GENUINE_T7_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
SOURCE='a55309b43bd58e43f9a8516fca219213e38881b6142e66f16a1b09c8505c4983'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def add(a,b):
 d=dict(a)
 for k,v in b.items():
  d[k]=d.get(k,Q(0))+v
  if not d[k]:del d[k]
 return d
def scale(a,c):return {k:v*c for k,v in a.items() if v*c}
def mult(a,b):
 d={}
 for x,u in a.items():
  for y,v in b.items():
   k=tuple(sorted(x+y));d[k]=d.get(k,Q(0))+u*v
 return {k:v for k,v in d.items() if v}
def c(q):return {():Q(q)} if q else {}
def v(j):return {(j,):Q(1)}
def mm(A,B):
 out=[]
 for row in A:
  nr=[]
  for j in range(len(B[0])):
   z={}
   for k,x in enumerate(row):z=add(z,mult(x,B[k][j]))
   nr.append(z)
  out.append(nr)
 return out
def decode4(a):
 d={}
 for exp,q in a:
  key=tuple(i for i,n in enumerate(exp) for _ in range(n));d[key]=Q(q)
 return d
def decode_aff(a):return {tuple(k):Q(q) for k,q in a}
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None):
  nonlocal RAW
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h,(str(p),s,h)
  k=str(p.relative_to(ROOT))
  if k not in pins:RAW+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text()) if p.suffix=='.json' else p.read_text()
 pr=pin(ROOT/'reports'/PRI/'report.json',SOURCE);docs={}
 for p,h in pr['input_sha256'].items():
  d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
  if p.endswith('_REVIEW/report.json'):assert d['all_passed']
 D=pin(ROOT/'reports'/PRI/'actual-filtered23-Dfree16-quotient7-and-genuine4-source.json',
   '2973604d746e6e98f45225ecb01591162d68e3fe3d2f178aad280dfc5bcee0c5')
 pin(ROOT/'code/ambient_B18_C_genuine_T7_source_20261105_v63.py',pr['code_sha256'])
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'])
 # Build literal graded bases independently; coordinates are (target,r,ell).
 basis=[]
 for target,r in [(0,0),(1,0),(0,1),(1,1)]:basis += [(target,r,e) for e in range(4)]
 tails=[(2,0,0),(2,1,0),(2,0,1),(0,2,0),(1,2,0),(0,2,1),(1,2,1)]
 basis+=tails
 assert basis==[(x['raw_target'],x['r_power'],x['ell_power']) for x in D['literal_filtered23_basis_F16_then_t7']]
 assert len(set(basis))==23 and Counter(r+e for t,r,e in basis)=={0:3,1:6,2:6,3:6,4:2}
 standards={0:{(t,0,0) for t in range(3)},1:{(t,r,1-r) for t in range(3) for r in range(2)},
  2:{(t,r,2-r) for t in range(2) for r in range(3)},
  3:{(t,r,3-r) for t in range(2) for r in range(3)},
  4:{(0,1,3),(1,1,3)}}
 assert all({m for m in basis if sum(m[1:])==n}==standards[n] for n in range(5))
 val=docs['all-full-tauAB12-values-and-fixed-known-rank-C4-certificates.json']
 cut=val['C4_rank4_certificate_if_positive'];assert cut['rows_ordered']==[0,3,2,5] and cut['universal_block_minor']=='-1331/15360'
 # The complement's grade4 standard rows1/4 are free because the known
 # actual relations have an invertible block on rows0/3 (cutrows2/5).
 freeimages=[(0,0,3),(1,0,3),(0,1,3),(1,1,3)]
 assert set(freeimages[:2])<=standards[3] and set(freeimages[2:])==standards[4]
 Fshift=[[int(j%4<3 and i==j+1) for j in range(16)] for i in range(16)]
 assert Fshift==D['F16_E_shift_matrix']
 assert D['actual_D_linear_quotient_q7_by23']==[[int(j==16+i) for j in range(23)] for i in range(7)]
 # Independent four coefficient operators: project -Pk3 to the two
 # surviving T7 cubic monomials using the actual curve NF3, not guesses.
 formal=docs['full-HB6-by2-affine-source-operator-and-old27-leaves.json']
 N3=[[Q(q) for q in row] for row in docs['complete-full-C-cubic155-operator-with-four-formal-ell3-callbacks.json']['full_no_cut_NF6']]
 operators=[]
 for label,side,k in [('alpha_A',0,1),('alpha_B',1,1),('beta_A',0,2),('beta_B',1,2)]:
  P=formal['P3_A' if side==0 else 'P3_B'][k];op={}
  for m in range(4):op=add(op,scale(decode_aff(P[m]),-N3[3*side][4*side+m]))
  assert all(key[0]=='cache' for key in op)
  saved=next(r for r in D['four_ell_coefficients_UNPAIRED'] if r['label']==label)
  assert op==decode_aff(saved['operator_old_C2_only']);operators.append((label,op))
 keys={k for label,op in operators for k in op};assert len(keys)==24 and sum(len(op) for label,op in operators)==60
 # Formal E7 and its COMPLETE kernel, with no numeric coefficient solve.
 E=[[{} for j in range(7)] for i in range(7)];E[2][0]=c(1)
 E[5][1]=v(0);E[6][1]=v(1);E[5][2]=v(2);E[6][2]=v(3);E[5][3]=c(1);E[6][4]=c(1)
 K=[[{} for j in range(4)] for i in range(7)]
 K[1][0]=c(1);K[3][0]=scale(v(0),-Q(1));K[4][0]=scale(v(1),-Q(1))
 K[2][1]=c(1);K[3][1]=scale(v(2),-Q(1));K[4][1]=scale(v(3),-Q(1));K[5][2]=c(1);K[6][3]=c(1)
 assert E==[[decode4(x) for x in row] for row in D['symbolic_E7_rows7_columns7']]
 assert K==[[decode4(x) for x in row] for row in D['symbolic_kernel7_by4']]
 assert mm(E,K)==[[{}]*4 for _ in range(7)]
 assert [[E[i][j] for j in [0,3,4]] for i in [2,5,6]]==[[c(i==j) for j in range(3)] for i in range(3)]
 assert [[K[i][j] for j in range(4)] for i in [1,2,5,6]]==[[c(i==j) for j in range(4)] for i in range(4)]
 E2=mm(E,E);assert E2[5][0]==v(2) and E2[6][0]==v(3)
 assert sum(bool(x) for row in E2 for x in row)==2 and mm(E2,E)==[[{}]*7 for _ in range(7)]
 # All nonzero multiplication paths are only constant*linear; no alpha
 # times beta product occurs before cancellation, so87bit applies globally.
 for A,B in [(E,K),(E,E),(E2,E)]:
  for row in A:
   for j in range(len(B[0])):
    terms=[mult(x,B[k][j]) for k,x in enumerate(row) if x and B[k][j]]
    assert len(terms)<=2 and all(len(key)<=1 for term in terms for key in term)
 known={tuple(k):Q(q) for k,q in D['known_old27_leaf_values_copied']};profiles=[];H=0
 for label,op in operators:
  U=0;den=1
  for key,w in op.items():
   q=known[key];U+=(abs(w.numerator)+w.denominator-1)//w.denominator*((abs(q.numerator)+q.denominator-1)//q.denominator)
   den=math.lcm(den,w.denominator*q.denominator)
  h=max((U*den).bit_length(),den.bit_length());H=max(H,h);profiles.append([label,str(U),str(den),h])
 assert profiles==D['four_future_cached_value_profiles'] and H==43 and 2*H+1==87
 gatepr=pin(ROOT/'reports/AMBIENT_B18_C_GENUINE_T7_VALUE_GATE_20261105_v63/report.json',
  '0bea531427072ef101882666b20aee3470d68a10defdda3796081c4a0bac8213')
 for p,h in gatepr['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 gate=pin(ROOT/'reports/AMBIENT_B18_C_GENUINE_T7_VALUE_GATE_20261105_v63/four-cached-values-code-and-whole-positive-custody.json',
  '155aa4c0d15279bbc904a530fb3725cb7362e97d48c9cf8e333c7caccaa8f00f')
 worker=pin(ROOT/'code/ambient_B18_C_genuine_T7_cached_values_20261105_v63.py',
  'fe20a0bf5e8709aaee8a78baa498c632f5b99140fd1e961b5e578d89dfbb7673')
 c2=pin(ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
  'c29b3794da3410228b7620102ec1906af68d5a12df3222d0ffa31c6ca8f439a6')
 delta={r['physical_row']:dict(r['delta']) for r in c2['physical12_coefficients']}
 assert all(Q(delta[k[1]].get(k[2],'0'))==known[k] for k in keys)
 fs={n.name:ast.get_source_segment(worker,n) for n in ast.parse(worker).body if isinstance(n,ast.FunctionDef)}
 for name,expr in [('add','a+b'),('times','a*b')]:assert fs[name].index('check(')<fs[name].index('return note('+expr+')')
 assert "for h in[SOURCE_SHA,sha(ROOT/GATE_REPORT),sha(Path(__file__))]" in fs['main']
 assert "assert c or e" in fs['main'] and 'assert mm(E,K)' in fs['main']
 assert "MAX_PRECHECK<=gate['whole_before_operation_bound']" in fs['main']
 future=nc(pr)+nc(D)+nc(gatepr)+nc(gate)+5000+nc(c2)+12000
 assert future<25000<30000 and RAW<1<<20
 assert not (ROOT/'reports/AMBIENT_B18_C_GENUINE_T7_CACHED_VALUES_20261105_v63').exists()
 proof={'literal23_filtered_basis_and_complementary_grade4_rows_proved':True,
  'primitive_D_free_proof':[
   'The four E³ images ofA,B,rA,rB are independent: two distinct grade3 units ell³A/B, and the two remaining grade4 units r ell³A/B after actual cutminor.',
   'Their count is the proved rankE³=4. Since ell³Tor C=0, multiplicationell³ identifies C/Tor C with ell³C. Thus these four source classes project to a basis moduloell in C/Tor.',
   'D-Nakayama and rank4 imply the fullDspanF mapsisomorphically to C/Tor. Therefore F is primitive D-free4 and C=F directsum Tor C asD-modules, with this actualchosenF.',
   'F is NOT R0-stable: r*(rA)=r²A has a nonzero T7 coordinate. No naive quotient r/z action is permitted.'
  ],
  'coefficient_source_order_proof':[
   'Allm4C4 is inF16 because its literal two gradedbasis vectors are r ell³A/B andm5C4=0.',
   'Every z-containing cubic lies in m2*(z target), andz target is alreadyinm2C4 by the actualfirstframe. Thus it lies inm4 andvanishesunderq. PlaneHcubics lie inm4 by theactualHquadratics.',
   'The actualrelations r ellH+P1³+terms>=4 andell²H+P2³+terms>=4 give E(rH)=-qP1³ and E(ellH)=-qP2³.',
   'Actualcurve cubics give q(r³A)=-a1q(r²ellA),q(r³B)=-b1q(r²ellB), with allremaining curve/carryterms inF16. Hence the fouroperators are exactly -N3 rows0/3 applied to Pk³. No λ,μ orδ4 coefficient survives.'
  ],
  'four_operators_old24_C2_leaves_60terms_checked':True,
  'formal_rank3_unit_minor_and_all28_kernel_products_and_freeI4':True,
  'complete_kernel_proof':'E*x=0 forces x0=0,x3=-alphaA*x1-betaA*x2,x4=-alphaB*x1-betaB*x2; x1,x2,x5,x6 arbitrary. This proves completeness, not only four independent cycles.',
  'genuine_target_proof':'Using actual D-split C=F+Tor, q induces kerE23/imE³ onto kerE7; kernelF16[ell]=ell³F16=imE³. This is a genuineC[ell] coordinate interface, not an identification with finitekerE23 alone.',
  'no_R0_transport_or_actual_formal_lifts_or_connecting_values_constructed':True,
  'unrun_worker_code_guard_and_whole_gate_PASS':True,'source_H43_and_all_operations_bound87':True,
  'whole_future_numeric_with5000_peer':future,'raw_unique_loaded_bytes':RAW,
  'four_response_values_not_evaluated':True,'parent_conditional_value_readiness':True,
  'hard_value':30000,'static_value':25000,'hard_global':150000,'hard_bits':2048,
  'hard_seconds':60,'hard_memory':2<<30,'hard_RAW':32<<20}
 ep=out/'independent-literal-Dfree16-T7-genuine4-source-and-unrun-code-gate.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'D_only_primitive_split_genuine_T7_source_and_code_PASS':True,
  'four_responses_R0_actions_and_connecting_values_uncomputed':True,'elapsed_seconds':time.monotonic()-start}
 assert nc(report)<5000
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
