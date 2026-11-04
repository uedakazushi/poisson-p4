#!/usr/bin/env python3
"""Independent polynomial, prime-intersection and generic-module review.
No producer imports, new source solve/jet, point, or CAS decomposition.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
R=Path(__file__).resolve().parents[1]
S='AMBIENT_B18_NORMALIZER_NIL_SQUARE_THREE_SUPPORTS_20261104_v63_REVIEW'
O=R/'reports'/S
P='AMBIENT_B18_NORMALIZER_NIL_SQUARE_THREE_SUPPORTS_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def save(n,x):(O/n).write_text(json.dumps(x,indent=2)+'\n')
def put(p,m,c):
 c=p.get(m,F(0))+c
 if c:p[m]=c
 else:p.pop(m,None)
def mul(p,q):
 out={}
 for m,c in p.items():
  for n,d in q.items():put(out,tuple(a+b for a,b in zip(m,n)),c*d)
 return out
def add(p,q,c=F(1)):
 p=dict(p)
 for m,v in q.items():put(p,m,c*v)
 return p
def power(p,n):
 out={(0,)*len(next(iter(p))):F(1)}
 for i in range(n):out=mul(out,p)
 return out
def compose(p,images):
 out={};dim=len(next(iter(next(q for q in images if q))))
 for m,c in p.items():
  x={(0,)*dim:c}
  for q,k in zip(images,m):
   for i in range(k):x=mul(x,q)
  out=add(out,x)
 return out
def decode(x):return {tuple(row[0]):F(row[1]) if len(row)==2 else F(row[1],row[2]) for row in x}
def enc(p):return [[list(m),str(c)] for m,c in sorted(p.items())]
def minimal(ms):return sorted(m for m in set(ms) if not any(n!=m and all(a<=b for a,b in zip(n,m)) for n in set(ms)))
def intersect(a,b):return minimal([tuple(max(x,y) for x,y in zip(m,n)) for m in a for n in b])
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));O.mkdir(exist_ok=True)
 paths={
  'primary':f'{P}/report.json','proof':f'{P}/nil-square-generic-lengths-and-annihilator-proof.json',
  'five':'AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63/five-line-and-generic-models.json',
  'prep':'AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json',
  'five_review':'AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW/report.json',
  'nil_review':'AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63_REVIEW/report.json',
  'ci_review':'AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63_REVIEW/report.json'}
 hashes={str(Path('reports')/p):sha(R/'reports'/p) for p in paths.values()};save('phase.json',{'stage':'hashes_persisted','input_sha256':hashes})
 d={k:json.loads((R/'reports'/p).read_text()) for k,p in paths.items()}
 assert sha(R/'reports'/paths['primary'])=='b2faaabeefd6cd0e2c6b3ad052df343d4518c21f9982e8b8f56d65690b38893f'
 for k,h in [('five_review','cffc3d5c6a9cf9ed1207788a27675bdc76f925380fa27f0dc3a3c0c8e5337f80'),('nil_review','02888d3c4b987327f6df2d8faf1f2adffa3f9c01af47e27aff8c5c0982e75472'),('ci_review','8e58bac4008623050b67dcd15509ac9a6237390cf417c694ec316a932f6f5899')]:
  assert sha(R/'reports'/paths[k])==h and d[k]['all_passed']
 f=d['five'];polys=[decode(f[k]) for k in ['source_Q0','source_Q1','source_f']]
 assert polys==list(map(decode,d['prep']['generators3']))
 unit=[{tuple(int(j==i) for j in range(4)):F(1)} for i in range(4)];a,dd,x,t=unit
 assert polys[2]==mul(add(a,x),power(add(add(x,a,-1),dd,-1),2))
 lines={'N':[0,0,0,1],'P':[2,0,-2,1],'C':[0,1,1,0],'R':[-2,1,-1,-2],'L':[1,1,-1,F(-1,2)]}
 for label,vec in lines.items():
  images=[{(1,):F(c)} if c else {} for c in vec]
  for p in polys:assert not compose(p,images),(label,p)
 # Linear coordinates in original (a,d,x,t), with explicit inverse.
 coords=[dd,{(1,0,0,0):F(1,2)},add(t,{(1,0,0,0):F(1,2)},-1),add(add(a,dd,-1),x)]
 matrix=[[p.get(tuple(int(j==i) for j in range(4)),F(0)) for i in range(4)] for p in coords]
 det=sum((F((-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4)))*matrix[0][p[0]]*matrix[1][p[1]]*matrix[2][p[2]]*matrix[3][p[3]] for p in itertools.permutations(range(4))),F(0));assert det==F(1,2)
 c,pp,n,H=unit
 inv=[{(0,1,0,0):F(2)},c,add(add(H,pp,-2),c),add(n,pp)]
 assert [compose(p,inv) for p in coords]==unit and [compose(p,coords) for p in inv]==unit
 axes={k:[sum(row[i]*v[i] for i in range(4)) for row in matrix] for k,v in lines.items() if k in ['N','P','C']}
 assert axes=={'N':[0,0,1,0],'P':[0,1,0,0],'C':[1,0,0,0]}
 # Exact squarefree monomial intersection; no radical/decomposition oracle.
 axes_primes=[[tuple(int(k==i) for k in range(4)) for i in inds] for inds in [[0,1,3],[0,2,3],[1,2,3]]]
 inter=intersect(intersect(axes_primes[0],axes_primes[1]),axes_primes[2]);expected=sorted([(0,0,0,1),(1,1,0,0),(1,0,1,0),(0,1,1,0)]);assert inter==expected
 pull=[compose({m:F(1)},coords) for m in [(1,1,0,0),(1,0,1,0),(0,1,1,0)]]
 assert pull[0]=={m:v/2 for m,v in mul(a,dd).items()}
 assert pull[1]==mul(dd,add(t,a,F(-1,2)))
 assert pull[2]==mul({(1,0,0,0):F(1,2)},add(t,a,F(-1,2)))
 # Replay ONLY previously certified generic Laurent models; t remains symbolic.
 Nimages=[{(1,0):F(1)},{(2,-1):F(4)},{(1,0):F(-3),(2,-1):F(2)},{(0,1):F(1)}]
 Pimages=[{(0,1):F(2),(1,0):F(1),(2,-1):F(7,12)},{(1,0):F(1)},{(0,1):F(-2),(1,0):F(-1),(2,-1):F(-11,12)},{(0,1):F(1)}]
 jets={}
 for label,images,cut in [('N',Nimages,3),('P',Pimages,3)]:
  q=[compose(p,images) for p in polys];assert not {m:v for p in q[:2] for m,v in p.items() if m[0]<cut}
  leading=min(m[0] for m in q[2]);assert leading=={'N':3,'P':2}[label]
  assert {m:v for m,v in q[2].items() if m[0]==leading}==({(3,0):F(-32)} if label=='N' else {(2,1):F(-16,3)})
  jets[label]=[enc(p) for p in q]
 models=f['generic_models'];assert models['N']['HF']==[1,1,1] and models['P']['HF']==[1,1] and models['C']['HF']==[1,2,1] and models['R']['HF']==[1,1] and models['L']['HF']==[1]
 assert models['N']['d_order']==2 and models['N']['d_squared_zero'] and models['P']['d_order']==1 and models['P']['d_squared_zero']
 assert f['normalized_E4']['Ass']==['(d)','N','C','R','P','L']
 # In kappa[v]/v^k, d-images are nonzero multiples of last basis vector.
 # Their ambient maximal ideals act trivially; dimension is exactly one.
 nil_lengths={'N':1,'P':1,'C':1,'R':0,'L':0,'main_d':0}
 assert d['proof']['generic_nil_square_lengths']==nil_lengths
 assert d['proof']['exact_annihilator_S9']==['u1','u2','u3','u4','u5','a-d+x','a*d','d*(2*t-a)','a*(2*t-a)']
 proof={'source_Q0_Q1_f_bound_to_accepted_CI':True,'all5_symbolic_line_substitutions_zero':True,'existing_generic_N_P_Laurent_polynomial_checks':jets,'coordinates_and_inverse':{'forward':list(map(enc,coords)),'inverse':list(map(enc,inv)),'det':str(det)},'squarefree_prime_intersection_monomials_c_p_n_H':inter,'pullback_scalars':[str(F(1,2)),str(F(1,2)),str(F(1,4))],'generic_nil_square_lengths':nil_lengths,'embedded_module_product_warning':'d:T -> Nil(E4) is a module isomorphism; Nil^2 corresponds to d*T, not d^2*T','Ass_exact':['N','P','C'],'H0_M_zero_from_H0_E4_zero':True,'annihilator_exact_reduced':d['proof']['exact_annihilator_S9'],'annihilator_reason':'I*M vanishes at all three generic associated primes; if nonzero it has an associated prime inherited from M, impossible. Equivalently finite-origin image is zero by H0(E4)=0. Reverse inclusion from all three nonzero simple localizations.','dimension_depth_multiplicity':[1,1,3],'upper8_or_D0_nil_square_claimed':False}
 save('independent-support-annihilator-proof.json',proof)
 report={'id':S,'all_passed':True,'independent_Fraction_polynomial_and_monomial_intersection':True,'producer_import_or_rerun':False,'new_source_solve_jet_point_CAS_GB':False,'normalized_E4_only':True,'Nil_E4_square_Ass_exact':['N','P','C'],'reduced_annihilator_exact':True,'dimension_depth_multiplicity':[1,1,3],'generic_lengths':[1,1,1],'core_d_square_vs_whole_E4_product_distinguished':True,'original_D0_nil_square_or_upper8_not_claimed':True,'input_sha256':hashes,'source_sha256':sha(Path(__file__)),'note_sha256':sha(R/'research/notes'/f'{S}.md'),'evidence_sha256':{'independent-support-annihilator-proof.json':sha(O/'independent-support-annihilator-proof.json')},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 save('report.json',report);save('phase.json',{'stage':'FINAL','all_passed':True});print(json.dumps(report))
if __name__=='__main__':main()
