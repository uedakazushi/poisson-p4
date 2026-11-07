"""Independent full coefficient check of base-Fp d1 factor pair; no char0 factor."""
from pathlib import Path
import hashlib,json,time
R=Path(__file__).resolve().parents[1];P=32003;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out=R/'reports/S61_TWO_TAIL_H_D1_MODFACTOR_ROOT_REVIEW_20261105_v63';out.mkdir(exist_ok=False);start=time.monotonic();pins={}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;return json.loads(q.read_text())
r=load('reports/S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_ENTRY_20261105_v63/report.json','3028597f392c2c16d3a1aa51065b59cb5702401df3ed2f93f045742dc2ad4fa6')
assert r['all_passed'] and r['fixed_image_factor_found_only'] and not r['geometric_H_irreducible_PRIMARY']
p='reports/S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_ENTRY_20261105_v63/complete-six-root-certificate.json'
c=load(p,r['evidence_sha256']['complete-six-root-certificate.json']);m=c['matched_polynomial_factor_pair'];assert m['mask']==3 and m['LC_allocation']==0
hp=R/'reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/factor-sparse.txt';pins[str(hp.relative_to(R))]=sha(hp);assert pins[str(hp.relative_to(R))]=='241488fafa9b45e5353f0b707f75902bc5bc32e85ee093b9f2e0a5f07c7d571c'
phase=dict(SOURCE_prior=True,method='independent sparse bivariate product directly from original94-term ZZinput; no Hensel/seed/producerreader',same_d1_prime32003=True,field_bits=30,existing_ROM_bits=40,owned_upper=24000,RAW_cap=4194304,no_char0_factor_claim=True)
(out/'phase.json').write_text(json.dumps(phase)+'\n')
assert all(P%d for d in range(2,179))
G={};Q={}
for dest,key,degree in [(G,'G',2),(Q,'Q',4)]:
 rows=m[key];assert len(rows)==10*(degree+1)
 for i,(x,y)in enumerate(rows):
  assert y==0 and 0<=x<P
  if x:dest[i//10,i%10]=x
assert max(a for a,e in G)==2 and max(a for a,e in Q)==4
F={};alpha_inv=pow(17578125,-1,P)
for line in hp.read_text().splitlines()[1:]:
 a,e,d,v=map(int,line.split());key=(a,e);F[key]=(F.get(key,0)+int(v)%P)%P
F={key:(v*alpha_inv)%P for key,v in F.items()if v}
product={};count=0;bits=0
for (a,e),x in reversed(list(Q.items())):
 for (b,f),y in reversed(list(G.items())):
  z=x*y;bits=max(bits,z.bit_length());assert bits<=30;count+=1;k=(a+b,e+f);product[k]=(product.get(k,0)+z)%P
product={k:v for k,v in product.items()if v};assert product==F
assert all(product.get((a,e),0)==F.get((a,e),0)for a in range(7)for e in range(19))
ev=dict(G_sparse=[[a,e,x]for (a,e),x in sorted(G.items())],Q_sparse=[[a,e,x]for (a,e),x in sorted(Q.items())],all133_coefficient_differences_zero=True,original94_term_reduction_independent=True)
p=out/'independent-sparse-full-product.json';p.write_text(json.dumps(ev,indent=2)+'\n')
d=dict(id=out.name,all_passed=True,independent_fixed_bivariate_image_factor_identity=True,input_sha256=pins,
 factor_a_degrees=[2,4],factor_terms=[len(G),len(Q)],all133_coefficients_match_original_H_d1_div_alpha=True,
 no_H_char0_factor_reducibility_irreducibility_or_geometric_component_claim=True,field_products=count,max_field_temp_bits=bits,existing_ROM_bits=40,owned_upper=24000,
 code_sha256=sha(Path(__file__)),evidence_sha256={p.name:sha(p),'phase.json':sha(out/'phase.json')},seconds=time.monotonic()-start)
p=out/'report.json';p.write_text(json.dumps(d,indent=2)+'\n');print(sha(p),count,d['factor_terms'])
