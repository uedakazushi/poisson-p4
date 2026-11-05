#!/usr/bin/env python3
"""Independent formal fixed-matrix proof of the ACTUAL M7-only splitting."""
import hashlib,itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1];PRI='AMBIENT_B18_C_CONNECTING_KERNEL3_SPLIT_20261105_v63';ID=PRI+'_REVIEW'
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def c(q):return {():Q(q)} if q else {}
def v(k):return {(k,):Q(1)}
def add(a,b):
 d=dict(a)
 for k,q in b.items():
  d[k]=d.get(k,Q(0))+q
  if not d[k]:del d[k]
 return d
def sc(p,q):return {k:x*q for k,x in p.items() if q*x}
def mul(a,b):
 d={}
 for x,q in a.items():
  for y,r in b.items():
   k=tuple(sorted(x+y));assert len(k)<=4
   d[k]=d.get(k,Q(0))+q*r
 return {k:q for k,q in d.items() if q}
def total(xs):
 p={}
 for x in xs:p=add(p,x)
 return p
def mm(A,B):return [[total(mul(x,B[k][j]) for k,x in enumerate(row)) for j in range(len(B[0]))] for row in A]
def det(A):
 p={};n=len(A)
 for perm in itertools.permutations(range(n)):
  sign=(-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n));term=c(sign)
  for i,j in enumerate(perm):term=mul(term,A[i][j])
  p=add(p,term)
 return p
def decode(A):return [[total(c(Q(q)) if key=='QQ' else sc(v(key),Q(q)) for key,q in x) for x in row] for row in A]
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None):
  nonlocal RAW
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h
  k=str(p.relative_to(ROOT))
  if k not in pins:RAW+=p.stat().st_size
  pins[k]=s;return json.loads(p.read_text())
 pr=pin(ROOT/'reports'/PRI/'report.json');docs={}
 for p,h in pr['input_sha256'].items():
  d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
  if p.endswith('_REVIEW/report.json'):assert d['all_passed']
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 code=ROOT/'code/ambient_B18_C_connecting_kernel3_split_corollary_20261105_v63.py';assert sha(code)==pr['code_sha256'];pins[str(code.relative_to(ROOT))]=sha(code)
 note=ROOT/'research/notes'/f'{PRI}.md';assert sha(note)==pr['note_sha256'];pins[str(note.relative_to(ROOT))]=sha(note)
 values=docs['phase2-full-coherent-head23-and-low-two-values.json'];assert values['simple_target_rank2_witness_nonzero']
 low=[[Q(x) for x in row] for row in values['fivehead_low_rH_ellH_rows']];chi=list(map(Q,values['actual_cyclic_top_chi']))
 d,a0,a4=chi[4],low[0][1],low[4][1];p1,b1=map(Q,values['head23_canonical_J2_k2_socle_coordinates'])
 assert chi[:4]==[0,0,0,0] and all(low[j]==[0,0] for j in [1,2,3])
 assert d and a0 and b1
 Sdata=docs['honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json']
 assert Sdata['source7_basis']==['h0','h1','h2','h3','h4','r*h3','r*h4']
 Rm=[[c(Q(x)) for x in row] for row in Sdata['actual_r_action7']]
 Rg=[[c(Q(x)) for x in row] for row in docs['genuine-four-R0-action-transport-and-cached-HB-operator-identities.json']['canonical_r_Jordan_action']]
 assert Rg==[[c(i==1 and j==0) for j in range(4)] for i in range(4)]
 # Fixed known values and all8remaining highcoordinates kept indeterminate.
 M=[[{} for j in range(7)] for i in range(4)]
 M[0][4]=c(d);M[1][6]=c(d);M[2][0]=c(a0);M[2][4]=c(a4);M[1][1]=c(p1);M[3][1]=c(b1)
 for j in [0,2,3,4]:M[1][j]=v('p'+str(j));M[3][j]=v('b'+str(j))
 assert M==decode(D['actual_edge4_by7_with_eight_formal_coords'])
 owned=count([pr,D,docs,M,Rm,Rg,pins])+15000;assert owned<50000 and RAW+1000000<32<<20
 save(out/'positive-before-fixed-source-identities.json',{'source_only':True,'eight_unknowns_formal':True,'fixed4det_only_24permutations':True,'max_polynomial_degree':4,'whole_numeric_upper':owned,'hard_bits':2048})
 assert mm(M,Rm)==mm(Rg,M)
 S=[[c(i==j) for j in [4,6,0,1]] for i in range(7)]
 assert S==decode(D['S_selector7_by4']) and mm(Rm,S)==mm(S,Rg)
 frame=mm(M,S);assert frame==decode(D['delta_S_matrix4'])
 determinant=det(frame);assert determinant==c(d*d*a0*b1)==c(Q(11979,2260608800))
 # Constructkernel independently by subtracting theactual section's known
 # head1 andrh4 columns fromh3 andfromh2-h3/4. No unknown inversion.
 K=[[{} for j in range(3)] for i in range(7)]
 K[3][0]=c(1);K[1][0]=sc(v('b3'),-1/b1)
 K[6][0]=sc(add(v('p3'),sc(v('b3'),-p1/b1)),-1/d)
 K[5][1]=c(1);K[2][2]=c(1);K[3][2]=c(Q(-1,4))
 beta=add(v('b2'),sc(v('b3'),Q(-1,4)));p=add(v('p2'),sc(v('p3'),Q(-1,4)))
 K[1][2]=sc(beta,-1/b1);K[6][2]=sc(add(p,sc(beta,-p1/b1)),-1/d)
 assert K==decode(D['kernel_M_basis7_by3_with_four_formal_coords'])
 assert mm(M,K)==[[{} for j in range(3)] for i in range(4)]
 Rk=[[c(i==1 and j==0) for j in range(3)] for i in range(3)]
 assert mm(Rm,K)==mm(K,Rk) and det([K[i] for i in [2,3,5]])==c(1)
 assert Sdata['fixed_I_Told_source_basis']['first']['source_V56_row']==55
 assert len(Sdata['fixed_I_Told_source_basis']['other_three'])==3
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'all28_actual_edge_equivariance_identities':True,'all28_literal_S_stability_identities':True,
  'fixed_S_determinant_independent_of8unknowns':str(d*d*a0*b1),
  'all12_deltaK_and21_rK_products_and_fixed_unit3minor':True,
  'actual_M7_sequence_R0_split':True,'actual_S_type':'J2(r) plus k²','actual_kernel_type':'J2(r) plus k',
  'split_proof':'S=<h4,rh4,h0,h1> is stableinM7 anddelta|S hasnonzeroQQconstantdeterminantd²a0b1. ThusitsR0-linear inverse isasectionandM7=S directsumkerdelta. KnownKcolumnsare3independentcyclesandkerdimension3,soarecomplete.',
  'five_known_actual_Kbar_classes':['z5h4','z5rh4','r²h4','r³h4','rh3'],
  'Kbar_basis_scope':'ITold4 isalreadytheactual4-dimensionalkernel submodule. Sincehead3lowtwo0,delta(h3)isgenuinesocle andrdelta(h3)=0;rh3 isactualKbar anditsnonzeroM7image provesindependencefromITold. Addingthetwoliftsu,s withtheirfourTRUEbutuncomputed parameters gives7barclasses. This doesnotserializefullformalK-lifts.',
  'no_lifting_of_quotient_split':'r²h4 isNONZEROinTold. Theliteral4-spanS is stableONLYinM7,so noTold→Gsplit orKbar-extension split orR0typeKbar follows.',
  'no_original_H0_length_exponents_or_nilcube_claim':True,'owned_upper':owned}
 ep=out/'independent-actual-M7-splitting-and-kernel-source-identities.json';save(ep,proof)
 rep={'id':ID,'all_passed':True,'frozen':True,'pure_source':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'actual_M7_only_R0_split_and_kernel_J2k_PASS':True,'Told_or_Kbar_extension_split_claimed':False,'seconds':time.monotonic()-started}
 save(out/'report.json',rep);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'seconds':rep['seconds']}))
if __name__=='__main__':main()
