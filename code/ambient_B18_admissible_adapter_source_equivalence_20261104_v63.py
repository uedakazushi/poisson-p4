#!/usr/bin/env python3
"""Universal admissible adapter source-unit identity and a jet nonlift test."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ADMISSIBLE_ADAPTER_SOURCE_EQUIVALENCE_20261104_v63'
PATHS={
 'regular8':'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63_REVIEW/report.json',
 'actual_interface':'reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63_REVIEW/report.json',
 'actual_native_operator':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def zero():return {}
def one():return {():Q(1)}
def var(*s):return {(tuple(s),):Q(1)}
def add(a,b,sign=Q(1)):
 z=dict(a)
 for m,c in b.items():
  t=z.get(m,Q(0))+sign*c
  if t:z[m]=t
  else:z.pop(m,None)
 return z
def mul(a,b):
 z={}
 for m,c in a.items():
  for n,d in b.items():
   k=tuple(sorted(m+n));v=z.get(k,Q(0))+c*d
   if v:z[k]=v
   else:z.pop(k,None)
 return z
def tr(A):return list(map(list,zip(*A)))
def mm(A,B):
 return [[sum_poly(mul(a,b) for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def sum_poly(ps):
 z={}
 for p in ps:z=add(z,p)
 return z
def counts(A):return sum(len(p) for row in A for p in row)
def K(G):
 M=[[zero() for _ in range(28)] for _ in range(8)]
 for p,(i,j) in enumerate(itertools.combinations(range(8),2)):
  M[j][p]=G[i];M[i][p]={m:-c for m,c in G[j].items()}
 return M
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'universal_source_identity_inputs_saved','inputs':inputs})
 assert all(json.loads((ROOT/p).read_text())['all_passed'] for p in PATHS.values())
 A=[[var('A',i,j) for j in range(8)] for i in range(8)]
 B=[[var('B',i,j) for j in range(4)] for i in range(8)]
 D=[[var('D',i,j) for j in range(4)] for i in range(4)]
 H=[[var('H',i,j) for j in range(4)] for i in range(28)]
 G=[var('G',i) for i in range(8)];KG=K(G);AG=[r[0] for r in mm(A,[[g] for g in G])];KAG=K(AG)
 assert not any(sum_poly(mul(G[i],KG[i][p]) for i in range(8)) for p in range(28))
 Ct=mm(KG,H);E=[A[i]+B[i] for i in range(8)]+[tr(Ct)[i]+D[i] for i in range(4)]
 S=[KG[i]+[zero() for _ in range(4)] for i in range(8)]
 S.extend([[zero() for _ in range(28)]+[one() if i==j else zero() for j in range(4)] for i in range(4)])
 Sp=[KAG[i]+[zero() for _ in range(4)] for i in range(8)]
 Sp.extend([[zero() for _ in range(28)]+[one() if i==j else zero() for j in range(4)] for i in range(4)])
 pairs=list(itertools.combinations(range(8),2))
 wedge=[[add(mul(A[i][p],A[j][q]),mul(A[j][p],A[i][q]),-1) for i,j in pairs] for p,q in pairs]
 bottom=mm(tr(B),KAG);T=[wedge[i]+H[i] for i in range(28)]+[bottom[i]+tr(D)[i] for i in range(4)]
 left=mm(tr(E),Sp);right=mm(S,T)
 assert left==right
 natural=mm(tr(A),KAG);assert natural==mm(KG,wedge)
 live=sum(map(counts,[A,B,D,H,KG,KAG,Ct,E,S,Sp,wedge,bottom,T,left,right,natural]));assert live<150000
 # A formal quadratic cover, not an actual Poisson/source point: R=k[[x,s]],
 # C=R[t]/(t²-s). Over R/m³, C'=x²(0,t) passes CG=0, but is not the
 # reduction of an exact syzygy. The fixed functional proves nonliftability.
 base=list(itertools.chain.from_iterable([(a,d-a) for a in range(d+1)] for d in range(3)))
 mon=[(a,b,t) for a,b in base for t in range(2)]
 def toy_mul(u,v):
  x,s,t=(u[i]+v[i] for i in range(3));s+=t//2;t%=2
  return (x,s,t) if x+s<3 else None
 def psi(v):return v[0].get((1,1,0),0)+v[1].get((2,0,1),0)
 images=[]
 for m in mon:
  first=toy_mul((0,0,1),m);second=toy_mul((1,0,0),m)
  v=[{first:-1} if first is not None else {},{second:1} if second is not None else {}]
  assert psi(v)==0;images.append({'h_monomial':m,'separating_functional':0})
 target=[{},{(2,0,1):1}];assert psi(target)==1 and toy_mul((2,0,1),(0,0,1)) is None
 ev={'source_manifest':inputs,'universal_polynomial_method':'independent sparse monomial dictionaries, all A8x8/B8x4/D4x4/H28x4/G8 indeterminates',
  'Koszul_sign':'K_G[:,i<j]=G_i*e_j-G_j*e_i','all28_Gt_KG_zero':True,
  'all224_naturality_entries_zero':'A^t K(AG)=K(G) wedge²(A^t)',
  'all384_source_equivalence_entries_exact':True,
  'formal_parameterization':'E=[[A,B],[(K_G H)^t,D]], Gprime=A G; A,D nearI, B near0; Eunit',
  'source32_matrix':'[[wedge²(A^t),H],[B^t K(AG),D^t]]',
  'source32_unit_reason':'At the original C9 closed point G=0,A=D=I, so this is upper triangular with identity diagonal, regardless H0',
  'actual_source_application':'Nprime=E N gives Nprime^t S(Gprime)=N^t S(G) T32; oldphysicalHadamard640 unchanged. Res4 gives diag(I640,Res4(T32)) sourceunit before subsequent ell/artin basechange.',
  'regularity_location':'G8 is regular in original C9, NOT assumed regular after ell or in C44; matrix identity and unit survive basechange',
  'coherent_invariance_scope':'complete presentation/coker and intrinsic initial Hilbert function; individual fixed619 coordinate values need coherent lower-frame transport',
  'finite_jet_admissibility':'actual admissible jets are jets of the displayed allorder parameterization; a truncated row annihilator condition alone is insufficient',
  'toy_nonlift_countercondition':{'ring':'C=k[[x,s]][t]/(t²-s), baseideal(x,s), G=(x,t) regular in C',
   'candidate_bottom_row_mod_base_m3':'(0,x²t)','candidate_times_G_mod_m3':0,
   'exact_syzygy_form':'(-t*h,x*h)',
   'separating_functional':'coefficient(first even,xs)+coefficient(second odd,x²)',
   'all12_monomial_h_images_annihilated':images,'candidate_functional':1,
   'proof':'h=a+b*t gives firsteven=-b*s, secondodd=x*b; the xs and x² coefficients cancel on every exact syzygy'},
  'arbitrary_N2_zero_or_free_fibre_jet_directions_not_declared_admissible':True,
  'no_actual_adapter_coefficient_value_source_jet_or_HF4_computed':True,'universal_live_QQ_terms':live}
 ep=out/'universal-admissible-source-unit-and-jet-countercondition.json';dump(out,ep.name,ev)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'universal384_source_matrix_entries_verified':True,
  'actual_admissible_adapter_source_unit_lemma_proved':True,'truncated_admissibility_sufficiency_rejected_by_formal_countercondition':True,
  'arbitrary_N2_zero_or_actual_N2_values_or_HF4_claimed':False,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
