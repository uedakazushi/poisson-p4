#!/usr/bin/env python3
"""Independent exponent-ring/socle/source-input proof, no scalar DAG outputs."""
import itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1];PRI='AMBIENT_B18_C_CONNECTING_B63_SAFE_OUTPUT_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PRE=0
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def calc(a,b,m=False):
 global PRE
 if a and b:
  na,da=abs(a.numerator).bit_length(),a.denominator.bit_length();nb,db=abs(b.numerator).bit_length(),b.denominator.bit_length()
  h=max(na+nb,da+db)if m else max(na+db,nb+da)+1;PRE=max(PRE,h);assert h<=2048
 return a*b if m else a+b
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None,parse=True):
  nonlocal raw
  p=pathlib.Path(p);k=str(p.relative_to(ROOT));s=sha(p)
  if h:assert h==s,(k,h,s)
  if k not in pins:raw+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text())if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json')
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 pin(ROOT/'code/ambient_B18_C_connecting_B63_safe_output_source_20261105_v63.py',pr['code_sha256'],False)
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 def src(s):return ROOT/next(p for p in pins if p.endswith(s))
 joint=src('joint-degree3-certificate.json');q20=src('phase1-common105-and-114-source-equations.json')
 r_s=field(joint,'monomials');assert r_s==[[d-j,j]for d in range(4)for j in range(d+1)]
 U=field(joint,'all105_unknown_original_node_order');F=field(joint,'all105_equation_original_node_order')
 jets=dict(field(joint,'all105_unknown_jets'));old=dict(field(q20,'all105_cached_composite_values'));b27=field(q20,'base27_monomials')
 assert set(U)==set(jets)==set(old)and len(U)==len(F)==105 and field(joint,'CA_F9_targets')==['0']*7+['r','s']
 ad=src('one-selected-delta63-fourth-origin-adjoint-and-unrun-B9-source-closure.json')
 assert field(ad,'actual_g_original_root_pair')==[20144,40922]and not field(ad,'origin_gradient_g_wrt105_sparse')and not field(ad,'fixed_rho_equation_order_sparse')
 def exp(m):return tuple(m.count(i)for i in range(9))
 source=D['basis63_monomials'];E=list(map(exp,source));index={e:i for i,e in enumerate(E)}
 expected=set()
 for d in range(5):
  for m in itertools.combinations_with_replacement(range(9),d):
   e=exp(m);nz=sum(e[:7]);plane=sum(e[7:])
   if nz==0 or(nz==1 and plane<=2)or(nz==2 and not plane and all(not e[i]for i in[1,2,3,6])):expected.add(e)
 assert set(E)==expected and len(E)==63 and E[0]==(0,)*9
 for e in E:
  for i,x in enumerate(e):
   if x:a=list(e);a[i]-=1;assert tuple(a)in index
 table=D['product63_indices_minus1_zero'];assert len(table)==63
 for i,e in enumerate(E):
  for j,f in enumerate(E):assert table[i][j]==index.get(tuple(a+b for a,b in zip(e,f)),-1)
 # A direct all63^3 associativity check, independent of variable-commutator
 # fixtures, together with downward closure certifies the monomial quotient.
 for i,j,k in itertools.product(range(63),repeat=3):
  a,b=table[i][j],table[j][k]
  assert(-1 if a<0 else table[a][k])==(-1 if b<0 else table[i][b])
 vars=[index[exp([i])]for i in range(9)]
 socle=[j for j in range(63)if all(table[i][j]<0 for i in vars)]
 assert socle==D['socle32_indices']and len(socle)==32
 known=[j for j in range(63)if j not in set(socle)];assert known==D['known31_indices']and len(known)==31
 assert all(table[i][j]<0 for i,j in itertools.product(socle,repeat=2))
 assert sum(not any(e[:7])for e in E)==15 and sum(sum(e[:7])==1 for e in E)==42 and sum(sum(e[:7])==2 for e in E)==6
 assert sum(not any(E[i][:7])for i in socle)==5 and sum(sum(E[i][:7])==1 for i in socle)==21
 for record in D['trusted_quotient_homomorphisms']:
  B=list(map(exp,record['target_monomials']));bi={e:i for i,e in enumerate(B)};P=[bi.get(e,-1)for e in E]
  assert P==record['projection63']and all(P[i]<0 for i in socle)
  for i,j in itertools.product(range(63),repeat=2):
   k=table[i][j];left=-1 if k<0 else P[k]
   right=bi.get(tuple(a+b for a,b in zip(E[i],E[j])),-1)if P[i]>=0 and P[j]>=0 else -1
   assert left==right
 assert len(D['trusted_quotient_homomorphisms'][0]['target_monomials'])==10 and D['trusted_quotient_homomorphisms'][1]['target_monomials']==b27
 owned=count([pr,D,r_s,U,F,jets,old,b27,E,table,pins])+25000;assert owned<150000 and raw+1000000<32<<20
 save(out/'positive-before-input-conversion-only.json',{'whole_SOURCE_upper':owned,'only_provided_cache_inputs_converted':True,'no_scalar_DAG_output_composed':True,'hard2048':True})
 lower=dict(D['provided105_lower_extension63']);overlaps=0;terms=0
 for node in U:
  a=list(map(Q,jets[node]));assert len(a)==10
  v={}
  # Coefficient extraction first fixes the desired ell exponent; each source
  # s-power supplies binom(b,j)/5^j under s=r+ell/5.
  for degree in range(4):
   for ell in range(degree+1):
    v0=Q(0)
    for (ra,sb),q in zip(r_s,a):
     if ra+sb==degree and ell<=sb:v0=calc(v0,calc(q,Q(math.comb(sb,ell),5**ell),True))
    v[exp([7]*(degree-ell)+[8]*ell)]=v0
  for m,q in zip(b27,old[node]):
   e=exp(m);q=Q(q)
   if not any(e[:7]):assert v[e]==q;overlaps+=1
   else:v[e]=q
  assert set(v)=={E[i]for i in known}
  actual=[str(v.get(e,Q(0)))for e in E];assert actual==lower[node]
 assert overlaps==630
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'all3969_products_and250047_full_associativity_checks':True,'monomial_complement_ideal_downward_closure_verified':True,
  'basis63_parts':[15,42,6],'known31_unknown32':[31,32],'socle_parts_plane4_zi_plane2_activezz':[5,21,6],
  'all288_m_times_unknown_and1024_socle_squared_zero':True,'both3969_trusted_quotient_homomorphisms_checked':True,
  'all630_provided_input_overlaps_and105_lowerextensions_verified':True,'same_original105_and_r_s_to_r_ell_ordinary_frame':True,
  'entire_safe_delta_proof':'All unprovided input coefficients lie in Soc(B63), not an arbitrary truncation ideal. For two extensions Phi,L with the same known31, every regular rational original DAG f satisfies f(Phi)-f(L)=Df(Phi0)*(Phi-L) because m*Soc=Soc^2=0. Units stay units at the same origin. The independently proved ALL105 gradient zero for20144-40922 therefore gives equality of the ENTIRE B63 delta63 output, without computing it.',
  'other_output_scope':'For bA,bB,nu and other nonzero-gradient functions, ONLY plane10/B27 output projections are certified: both quotient maps kill every unknown32 input slot. Full B63 values of these functions are not adopted.',
  'no_105_implicit_solve_or_new_physical_jet_or_delta_output_or_edge_value':True,
  'lower_skeleton_satisfies_full105_equations_NOT_claimed':True,'head30_coherence_and_future_worker_height_custody_GO_still_separate':True,
  'source_owned_upper':owned,'input_conversion_prebits':PRE,'RAW':raw}
 ep=out/'independent-exponent63-socle-safe-output-and-all-cache-input-overlaps.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_only':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'entire_delta63_lower_cache_sufficiency_proved':True,'actual_scalar_output_evaluated':False,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'prebits':PRE,'seconds':report['seconds']}))
if __name__=='__main__':main()
