#!/usr/bin/env python3
"""Tiny associative D-algebra: listed module and barcube facts do not improve index7."""
import hashlib,itertools,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_BAR_CUBE_ZERO_INDEX7_COUNTERMODEL_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
names=['1','x','x2','y']+[f'u{i}'for i in range(1,10)]+['v1','v2']
limits={3:11,13:1,14:1}
def term(i,k=0):return {}if i in limits and k>=limits[i]else {(i,k):1}
def basisprod(i,j):
 if i==0:return term(j)
 if j==0:return term(i)
 table={(1,1):(2,0),(1,2):(3,1),(1,3):(3,3),(2,2):(3,4),(2,3):(3,6),(3,3):(3,8)}
 key=tuple(sorted((i,j)))
 return term(*table[key])if key in table else {}
def mul(p,q):
 z={}
 for (i,a),c in p.items():
  for (j,b),d in q.items():
   for (h,k),v in basisprod(i,j).items():
    for key,w in term(h,k+a+b).items():z[key]=z.get(key,0)+c*d*v*w
 return {k:c for k,c in z.items()if c}
def moduleprod(A,B):
 out={}
 for i,a in A.items():
  for j,b in B.items():
   for (h,k),c in mul(term(i,a),term(j,b)).items():out[h]=min(out.get(h,10**6),k)
 return out
def subset(A,B):return all(i in B and k>=B[i]for i,k in A.items())
def labelmod(A):return {names[i]:f'ell^{k}'for i,k in sorted(A.items())}
def encoded(p):return [[names[i],k,c]for (i,k),c in sorted(p.items())]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 actual='reports/AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63_REVIEW/report.json';actualreport=json.loads((ROOT/actual).read_text());assert actualreport['all_passed']and actualreport['genuine_Nbar_cube_zero']and not actualreport['original_cube_vanishing_proved']
 assoc=hashlib.sha256();assoc_count=0
 for i,j,k in itertools.product(range(15),repeat=3):
  a=mul(mul(term(i),term(j)),term(k));b=mul(term(i),mul(term(j),term(k)));assert a==b
  assoc.update(json.dumps([i,j,k,encoded(a)],separators=(',',':')).encode()+b'\n');assoc_count+=1
 assert all(mul(term(0),term(i))==term(i)==mul(term(i),term(0))for i in range(15))
 assert all(basisprod(i,j)==basisprod(j,i)for i in range(15)for j in range(15))
 # Bilinearity descends through every D-torsion relation.
 assert all(mul(term(i,e),term(j))=={}for i,e in limits.items()for j in range(15))
 N={i:0 for i in range(1,15)};H={i:0 for i in limits};K={i:0 for i in[3,4,5,6,7,13,14]}
 assert subset(moduleprod(K,N),K)
 powers=[N]
 for _ in range(6):powers.append(moduleprod(powers[-1],N))
 assert powers[2]=={3:1}and powers[3]=={3:4}and powers[4]=={3:7}and powers[5]=={3:10}and powers[6]=={}
 ellH={i:k+1 for i,k in H.items()if i not in limits or k+1<limits[i]};assert powers[2]==ellH
 K2=moduleprod(K,K);assert K2=={3:8}and moduleprod(N,K2)=={}
 assert all(k>=1 for k in powers[2].values()) # Nbar^3=0; no scalar source values involved.
 xp=term(1)
 for _ in range(5):xp=mul(xp,term(1))
 assert xp=={(3,10):1}and mul(xp,term(1))=={}
 freeN=[i for i in N if i not in limits];torN=[i for i in N if i in limits];freeK=[i for i in K if i not in limits];torK=[i for i in K if i in limits]
 assert (len(freeN),len(torN),len(freeK),len(torK))==(11,3,4,3)
 assert set(limits).issubset(K)and len(set(freeN)-set(freeK))==7
 data={'base':'D=QQ[[ell]]','D_module_basis15':names,'torsion_basis_annihilator_powers':{names[i]:e for i,e in limits.items()},
 'all225_D_bilinear_basis_products':[[i,j,encoded(basisprod(i,j))]for i in range(15)for j in range(15)],'all3375_associativity_checks':assoc_count,'associativity_stream_sha256':assoc.hexdigest(),'unit_commutativity_and_torsion_well_defined':True,
 'four_generator_core_presentation':'D[x,y]/(x^3-ell*y, x*y-ell^3*y, y^2-ell^8*y, ell^11*y); constructed D-module basis1,x,x²,y proves no additional y relation. Added u1..u9 are D-free square-zero and annihilate N; v1,v2 are k-square-zero and annihilate N.',
 'core_derived_products':{'x2_y':'ell^6*y','x2_x2':'ell^4*y','x3_y':'ell^9*y'},
 'exact_nilradical_source':'N=(x,x²,y,u1..u9,v1,v2), N^7=0 and algebra/N=D is reduced, hence it is exactly the nilradical.',
 'exact_N_powers1_through7':list(map(labelmod,powers)),'x6_nonzero_witness':encoded(xp),'x7_zero':True,'exact_nilpotency_index':7,
 'H0_D_torsion':labelmod(H),'N3_equals_ell_H0':True,'Nbar_cube_zero':True,'Nbar_dimension_over_QQ':14,
 'K_ideal':labelmod(K),'K_square':labelmod(K2),'N_times_K_square_zero':True,'D_rank_N':11,'torsion_minimal_generators_N':3,'mu_D_N':14,'D_rank_K':4,'torsion_minimal_generators_K':3,'mu_D_K':7,'N_over_K_D_free_rank':7,
 'scope':'Only the listed nil identities and D-module rank/minimal-generator facts are modeled. This is NOT a counterexample satisfying all actual S=QQ[[z7,r,ell]] source matrices, exact sequences, Loewy data, or normalization maps. It proves those listed facts alone do not improve the original index[3,7].',
 'actual_original_nilindex_changed':False,'actual_source_gauge_jet_values_evaluated':False}
 ep=out/'associative-table-and-index7-Dmodule-countermodel.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 result={'id':STEM,'all_passed':True,'status':'PROVED_tiny_associative_countermodel_for_listed_Dmodule_facts_only','input_sha256':{actual:sha(ROOT/actual)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'associativity_checks':assoc_count,'nilpotency_index_countermodel':7,'actual_nilindex_improvement_claimed':False,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));signal.alarm(0)
if __name__=='__main__':main()
