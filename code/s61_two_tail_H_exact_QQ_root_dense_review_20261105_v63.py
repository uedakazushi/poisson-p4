"""Independent reverse-order dense three-variable product of exact H factors."""
from pathlib import Path
import hashlib,json,math,time
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();pins={};start=time.monotonic()
out=R/'reports/S61_TWO_TAIL_H_EXACT_QQ_ROOT_DENSE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;return json.loads(q.read_text())
prefix='reports/S61_TWO_TAIL_H_EXACT_QQ_ENTRY_20261105_v63/'
r=load(prefix+'report.json','835f6a6404958bf947efa5e1809f49bbe95872ce2033fd89bdab0d574c640a5b');assert r['all_passed'] and r['primitive_positive_factor_count']==2
c=load(prefix+'canonical-factor-candidates.json',r['evidence_sha256']['canonical-factor-candidates.json'])
hp='reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/factor-sparse.txt';pins[hp]=sha(R/hp);assert pins[hp]=='241488fafa9b45e5353f0b707f75902bc5bc32e85ee093b9f2e0a5f07c7d571c'
(out/'phase.json').write_text(json.dumps(dict(SOURCE_prior=True,method='reverse factors, dense7x10x5 tensor, no native API or packed factoraddress reader',committed_output_bound357=True,identity_hard_bits=384,owned_upper=50000,RAW_cap=4194304))+'\n')
def blank():return [[[0 for d in range(5)]for e in range(10)]for a in range(7)]
target=blank()
for line in (R/hp).read_text().splitlines()[1:]:
 a,e,d,v=map(int,line.split());assert target[a][e][d]==0;target[a][e][d]=v
factors=[]
for f in sorted(c['factors'],key=lambda f:f['total_degree']):
 assert f['multiplicity']==1
 terms=[(tuple(k),int(v))for k,v in f['terms']];assert math.gcd(*(v for k,v in terms))==1
 assert all(0<=a<=6 and 0<=e<=9 and 0<=d<=4 and abs(v).bit_length()<=357 for (a,e,d),v in terms)
 factors.append((f['total_degree'],terms))
assert [n for n,t in factors]==[4,7] and [len(t)for n,t in factors]==[11,29]
current=blank();current[0][0][0]=1;pairs=0;bits=0;stagebounds=[]
for degree,terms in factors:
 ab=max(abs(v).bit_length()for block in current for row in block for v in row);bb=max(abs(v).bit_length()for k,v in terms);bound=ab+bb+18;assert bound<=384;stagebounds.append(bound)
 z=blank()
 for (a,e,d),v in reversed(terms):
  for i in range(7-a):
   for j in range(10-e):
    for k in range(5-d):
     x=current[i][j][k]
     if not x:continue
     t=x*v;pairs+=1;assert abs(t).bit_length()<=384;z[i+a][j+e][k+d]+=t
     bits=max(bits,abs(t).bit_length(),abs(z[i+a][j+e][k+d]).bit_length());assert bits<=bound
 current=z
assert current==target
poly={n:{k:v for k,v in t}for n,t in factors};S=poly[4]
assert {k:v for k,v in S.items()if k[0]==2}=={(2,0,2):125}
assert sum(v for (a,e,d),v in S.items()if a==1 and d==0)==-4000
# Exact discriminant in a: sparse coefficient polynomials in e,d.
co=[{(e,d):v for (a,e,d),v in S.items()if a==i}for i in range(3)]
def prod(x,y):
 z={}
 for (e,d),v in x.items():
  for (f,g),w in y.items():z[e+f,d+g]=z.get((e+f,d+g),0)+v*w
 return z
D=prod(co[1],co[1])
for k,v in prod(co[0],co[2]).items():D[k]=D.get(k,0)-4*v
D={k:v for k,v in D.items()if v};assert max(e for e,d in D)==3
assert {d:v for (e,d),v in D.items()if e==3}=={3:54000}
ep=out/'independent-dense-all350-product-and-quadratic-discriminant.json';ep.write_text(json.dumps(dict(dense_coefficients350=current,reverse_factor_total_degrees=[4,7],stage_bounds=stagebounds,discriminant_S4_e_d=[[e,d,v]for (e,d),v in sorted(D.items())],S4_discriminant_e3_leading='54000*d^3',all350_coefficients_match=True),indent=2)+'\n')
r=dict(id=out.name,all_passed=True,input_sha256=pins,independent_exact_reverse_dense_product_all350=True,H_QQ_reducible_proved=True,factor_terms=[11,29],factor_total_degrees=[4,7],factor_variable_degrees=[[2,3,2],[4,6,2]],
 S4_geometrically_irreducible_proved=True,S4_absolute_proof='As polynomial in a itsLC125d² and a-linear coefficient atd0=-4000 give geometric contentunit. Discriminant has odd e-degree3 withLC54000d³, hence cannot be a square in C(d)(e); Gauss gives absolute irreducibility.',
 R7_QQ_or_geometric_irreducibility_independently_proved=False,API_completeness_not_used_for_decomposition_or_S4_absolute_proof=True,
 integer_products=pairs,observed_identity_bits=bits,hard_identity_bits=384,owned_upper=50000,native_internal_growth_remains_opaque=True,
 code_sha256=sha(Path(__file__)),evidence_sha256={ep.name:sha(ep),'phase.json':sha(out/'phase.json')},seconds=time.monotonic()-start)
p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');print(sha(p),pairs,bits)
