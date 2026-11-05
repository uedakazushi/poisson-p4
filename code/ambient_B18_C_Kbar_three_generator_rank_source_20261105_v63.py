#!/usr/bin/env python3
"""Pure SOURCE: unknown deep coefficients, unit3-generator minor and r-rank4."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_THREE_GENERATOR_RANK_SOURCE_20261105_v63'
PATHS={
 'initial':'reports/AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63/two-cached-alpha-reductions-and-projected-r2-seven-source-certificates.json',
 'initial_peer':'reports/AMBIENT_B18_C_KBAR_R2_INITIAL_SOURCE_20261105_v63_REVIEW/report.json',
 'seven':'reports/AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63/actual-seven-barbasis-typed56-cycle-operator-and-full11-unit-complement.json',
 'seven_peer':'reports/AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63_REVIEW/report.json',
 'R3_peer':'reports/AMBIENT_B18_C_CONNECTING_R3_FACTORIZATION_SOURCE_20261105_v63_REVIEW/report.json',
 'old_exact':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json',
 'old_exact_peer':'reports/AMBIENT_B18_T_COMBINED619_QUARTIC_WITNESS_20261105_v63_REVIEW/report.json',
 'four_upper':'reports/AMBIENT_B18_C_K_R0_FOUR_GENERATOR_UPPER_SOURCE_20261105_v63/pure-four-R0-generator-upper-bound-and-exact-rU-recovery.json',
 'finite_K_peer':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
}
MAXPRE=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def h(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def pre(n):
 global MAXPRE
 MAXPRE=max(MAXPRE,n);assert n<=1153
def qa(a,b):pre(2*max(h(a),h(b))+1);q=a+b;assert h(q)<=384;return q
def qm(a,b):pre(h(a)+h(b));q=a*b;assert h(q)<=384;return q
def c(q):return{():Q(q)}if q else{}
def var(label):return{(label,):Q(1)}
def add(a,b):
 out=a.copy()
 for k,v in b.items():
  out[k]=qa(out.get(k,Q(0)),v)
  if not out[k]:del out[k]
 return out
def neg(a):return{k:-v for k,v in a.items()}
def mul(a,b):
 out={}
 for k,x in a.items():
  for l,y in b.items():
   m=tuple(sorted(k+l));out[m]=qa(out.get(m,Q(0)),qm(x,y))
 return{k:q for k,q in out.items()if q}
def sum_poly(xs):
 p={}
 for x in xs:p=add(p,x)
 return p
def mm(A,B):return[[sum_poly(mul(x,y)for x,y in zip(row,col))for col in zip(*B)]for row in A]
def enc(p):return[[list(k),str(v)]for k,v in sorted(p.items())]
def encM(A):return[[enc(p)for p in row]for row in A]
def det(A):
 # Subset DP gives a full polynomial certificate without assigning any
 # unknown actual deep coefficient or introducing a symbolic denominator.
 n=len(A);state={0:c(1)}
 for row in A:
  nxt={}
  for mask,v in state.items():
   for j in range(n):
    if mask>>j&1:continue
    term=mul(v,row[j])
    if sum(mask>>k&1 for k in range(j+1,n))%2:term=neg(term)
    key=mask|1<<j;nxt[key]=add(nxt.get(key,{}),term)
  state=nxt
 return state[(1<<n)-1]
def main():
 start=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()}
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 for k,v in d.items():
  if k.endswith('_peer'):assert v['all_passed']
 initial=d['initial'];assert initial['projected_rank']==1
 matrix=initial['actual_projected_r2_matrix2_by7'];lamU,lamS=map(Q,matrix[0][-2:]);assert lamS
 cc=Q(11767,3652);aa=qa(lamU,cc)
 # Fixed honest basis(v,q1,q2,g,w,U,S). B,C,D,E,F,G,H are TRUE but
 # uncomputed actual deeper coordinates; they stay algebraically free.
 R=[[{}for _ in range(7)]for _ in range(7)]
 R[1][0]=c(1);R[3][2]=c(1)
 R[1][4]=c(-cc);R[3][4]=var('H')
 for row,p in enumerate([c(aa),var('B'),var('C'),var('D'),c(1),{},{}]):R[row][5]=p
 for row,p in enumerate([c(lamS),var('E'),var('F'),var('G'),{},{},{}]):R[row][6]=p
 R2=mm(R,R);R3=mm(R,R2);assert all(not p for row in R3 for p in row)
 assert R2[1][5]==c(lamU)and R2[1][6]==c(lamS)
 assert R2[3][5]==add(var('H'),var('C'))and R2[3][6]==var('F')
 assert all(not p for i,row in enumerate(R2)for j,p in enumerate(row)if i not in[1,3]or j not in[5,6])
 rankrows=[1,3,4,0];rankcols=[0,2,5,6]
 minor=[[R[i][j]for j in rankcols]for i in rankrows];assert det(minor)==c(lamS)
 # Every column belongs to span(rv,rq2,rU,rS), exactly.
 C4=[[{}for _ in range(7)]for _ in range(4)]
 C4[0][0]=c(1);C4[1][2]=c(1);C4[0][4]=c(-cc);C4[1][4]=var('H');C4[2][5]=c(1);C4[3][6]=c(1)
 selected=[[row[j]for j in rankcols]for row in R];assert mm(selected,C4)==R
 e=lambda j:[c(int(i==j))for i in range(7)]
 columns=[e(5),e(6),e(2),[row[5]for row in R],[row[6]for row in R],e(3),[row[6]for row in R2]]
 prolonged=[list(row)for row in zip(*columns)];prolonged_det=det(prolonged)
 assert prolonged_det==c(qm(lamS,lamS))
 remaining_det=det([[R2[1][5],R2[1][6]],[R2[3][5],R2[3][6]]])
 assert remaining_det==add(mul(c(lamU),var('F')),neg(mul(c(lamS),add(var('H'),var('C')))))
 data={'scope':'Pure SOURCE actual original zero-degree K/ellK: three-generator upper bound and exact r-rank, with all unknown deeper coordinates formal. No new actual action query.',
  'honest_basis':['v=z5h4','q1=z5rh4','q2=r²h4','g=r³h4','w=rh3','U','S'],
  'known_initial_scalars':{'lambdaU':str(lamU),'lambdaS':str(lamS),'c':str(cc),'A=lambdaU+c':str(aa)},
  'seven_TRUE_uncomputed_deeper_coordinates':['B','C','D','E','F','G','H'],
  'formal_model_R7':encM(R),'formal_R_squared':encM(R2),'all49_R_cubed_zero_polynomial_identities':True,
  'actual_source_of_model':'ITold4=span(v,q1,q2,g), rv=q1,rq2=g,rq1=rg0. G38 gives rw=-c*q1+H*g since oldm³=QQ*g. Exact M7 relations put rU in w+ITold andrS inITold. Accepted initial reductions give the v coefficients A andlambdaS. Hence every unspecified actual coefficient is exactly one of B..H; no omitted formal action or guessed zero.',
  'rank4_witness_rows':rankrows,'rank4_witness_columns':rankcols,'rank4_witness_determinant':str(lamS),
  'all49_exact_rank4_factorizations':True,'factor4_by7':encM(C4),'actual_rank_r_Kbar':4,
  'three_candidate_actual_bar_generators':['U','S','q2=r²h4'],
  'all_seven_prolonged_columns':['U','S','q2','rU','rS','rq2','r²S'],
  'formal_prolonged7_matrix':encM(prolonged),'formal_prolonged7_determinant':enc(prolonged_det),
  'unit_minor_conclusion':'lambdaS!=0 is an actual QQunit; det=lambdaS² independently of allB..H. Seven prolongedcolumns areQQbasis, soU,S,q2 generateKbar overQQ[[r]] andthereforeRbar.',
  'upstairs_Nakayama':'K is finiteD andthusfiniteR0. Any actualK lifts ofU,S,q2 generateK byellNakayama. Hence mu_R0K<=3. FullR0lifts are not expanded.',
  'lower_bound_two':'U,S have independent constant-head images inTold/mTold; mKbar mapsintomTold, so mu_R0K>=2. Exact minimality2or3 remainsopen.',
  'r_Jordan_possibilities':[[3,3,1],[3,2,2]],
  'Jordan_proof':'dim7,rankr4,r³0 andrankr²>=1 imply exactly3blocks and thetwo displayedpartitions. This is r, notellDVR.',
  'remaining_one_scalar_obstruction':enc(remaining_det),
  'remaining_obstruction_scope':'rankr²2 iff lambdaU*F-lambdaS*(H+C)!=0 (partition331); otherwise rankr²1 (partition322). This scalar remainsUNCOMPUTED; alloldm³corrections retained.',
  'old_four_generator_result_retained_history':True,'no_H0_length_ell_torsion_nilcube_or_index_inference':True,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'maximum_QQ_preop_bits':MAXPRE,'all_passed':True}
 import re
 def count(x):
  if isinstance(x,(int,float,Q)):return 1
  if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
  if isinstance(x,(list,tuple)):return sum(map(count,x))
  if isinstance(x,dict):return count(list(x.keys()))+count(list(x.values()))
  return 0
 owned=count(d)+count(data)+count(R)+count(R2)+count(R3)+count(prolonged)+15000;assert owned<150000;data['whole_SOURCE_owned_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'formal-deep-action-unit-three-generator-minor-and-exact-rank-four.json';p.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'PURE_SOURCE_MU_R0_K_AT_MOST_THREE_RANK_R_FOUR_DEEP_SCALAR_FORMAL',
  'mu_R0_K_upper':3,'mu_R0_K_lower':2,'actual_rank_r_Kbar':4,'r_Jordan_possibilities':[[3,3,1],[3,2,2]],
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{p.name:sha(p)},'whole_SOURCE_owned_upper':owned,'maximum_QQ_preop_bits':MAXPRE,'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
