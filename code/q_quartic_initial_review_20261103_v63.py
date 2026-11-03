#!/usr/bin/env python3
"""Complete cubic-kernel coverage and actual raw quartic-plane carry audit."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha

def add_scaled(out,p,c):
 for m,v in p.items():
  z=out.get(m,F(0))+c*v
  if z:out[m]=z
  else:out.pop(m,None)
def unpack(p):return{tuple(m):F(c)for m,c in p}
def modular_rank(rows,prime):
 basis={}
 for original in rows:
  row={m:int(c.numerator)*pow(int(c.denominator),-1,prime)%prime for m,c in original.items()};row={m:c for m,c in row.items()if c}
  while row:
   pivot=min(row)
   if pivot not in basis:
    inv=pow(row[pivot],-1,prime);basis[pivot]={m:c*inv%prime for m,c in row.items()};break
   factor=row[pivot]
   for m,c in basis[pivot].items():
    v=(row.get(m,0)-factor*c)%prime
    if v:row[m]=v
    else:row.pop(m,None)
 return len(basis)

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--mixed-replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_ORDINARY_QUARTIC_INITIAL_20261103_v63';mainreport=read(folder/'report.json');mixed=read(folder/'mixed-syzygies/report.json');fresh=read(args.mixed_replay/'report.json');assert{k:v for k,v in mixed.items()if k not in('seconds','peak_rss_kb','certificate_sha256')}=={k:v for k,v in fresh.items()if k not in('seconds','peak_rss_kb','certificate_sha256')}
 paths=[Path(__file__),ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for reportpath in[folder/'report.json',folder/'mixed-syzygies/report.json',folder/'coverage.json']:
  report=read(reportpath);assert report['all_passed'];paths.append(reportpath)
  for file,h in report['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 tailpath=folder/'quartic-tails.json.gz';assert sha(tailpath)==mainreport['quartic_tails_sha256'];tails=read(tailpath)
 certpath=folder/'mixed-syzygies/certificate.json.gz';assert sha(certpath)==mixed['certificate_sha256']
 decoded=gzip.decompress(certpath.read_bytes());assert decoded==gzip.decompress((args.mixed_replay/'certificate.json.gz').read_bytes());certificate=json.loads(decoded);paths.append(certpath)
 for name in['integer-blocks.in','kernels.txt']:
  p=folder/f'mixed-syzygies/{name}.gz';assert gzip.decompress(p.read_bytes())==(args.mixed_replay/name).read_bytes();paths.append(p)
 q=read(ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz');c=read(ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63/cubic-initials.json.gz');planes=read(ROOT/'reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63/certificate.json.gz')
 Q=[{(i,j):F(v)for i,j,v in row}for row in q['all350_quadratic_rows']];C=list(map(unpack,c['raw350_cubic_terms']));G={row['source']:row for row in tails['all293_constant_quadratic_cancellation_rows']}
 # Replay the invertible constant row replacement independently.  Selected
 # original A57 rows plus e_source+span(A57) for all other293 are triangular.
 normalized={};selected=[];cofactors={}
 for operation in c['exact_constant_row_elimination']:
  source=operation['source'];cf={source:F(1)}
  for pivot,coefficient in operation['reductions']:add_scaled(cf,normalized[tuple(pivot)],-F(coefficient))
  if'new_pivot'in operation:
   normalized[tuple(operation['new_pivot'])]={i:v/F(operation['normalizer'])for i,v in cf.items()};selected.append(source)
  else:cofactors[source]=cf
 assert len(selected)==57 and len(cofactors)==293 and selected==c['quadratic_linear_multiples_recipe']['quadratic_original_rows']
 oldc={row['source']:unpack(row['terms'])for row in c['cubic_initial_relations']};nonzero=[];zero=[]
 for source,row in G.items():
  assert cofactors[source]=={int(i):F(v)for i,v in row['constant_full350_cofactors']}
  assert cofactors[source][source]==1 and set(cofactors[source])-{source}<=set(selected)
  row['C']=unpack(row['cubic_terms']);assert row['C']==oldc[source]
  (nonzero if row['C']else zero).append(source)
 assert len(nonzero)==135 and len(zero)==158
 syzygies=certificate['all219_exact_syzygies'];assert len(syzygies)==219;rankrows=defaultdict(list)
 for syz in syzygies:
  cubic={};actual=defaultdict(dict);linearcomb={}
  for term in syz['row_recipe']:
   if term[0]=='W':
    _,source,variable,value=term;value=F(value);p={tuple(sorted(m+(variable,))):v for m,v in Q[source].items()};add_scaled(cubic,p,value);add_scaled(actual[source],{(variable,):F(1)},value);key=('W',source,variable)
   else:
    _,source,value=term;value=F(value);add_scaled(cubic,G[source]['C'],value)
    for raw,coefficient in cofactors[source].items():add_scaled(actual[raw],{():F(1)},value*coefficient)
    key=('G',source)
   linearcomb[key]=linearcomb.get(key,F(0))+value
  assert not cubic
  saved={i:unpack(p)for i,p in syz['actual350_cofactors']};assert saved=={i:p for i,p in actual.items()if p}
  rankrows[syz['weight']].append({m:v for m,v in linearcomb.items()if v})
 prime=1000003;independent=sum(modular_rank(rows,prime)for rows in rankrows.values());assert independent==219
 rankreport=read(ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63/exact-rank/report.json');assert sum(row['QQ_cubic_initial_rank']for row in rankreport['block_ranks'])==1911
 assert 1995+135-1911==219 and 1995+293-1911==377
 # Restrict the ACTUAL original350 raw D_i, not the saved reduced G quartics,
 # then form each true cofactor's degree4 contribution independently.
 for branch in planes['branches']:
  names=q['variables'];eta=names.index('z22');free=sorted([names.index(n)for n in branch['free_variables']]+[eta]);R,*v=ring(','.join(names[i]for i in free),s.QQ);variables=dict(zip(free,v));graph=[]
  for encoded in branch['all35_graph']:
   p=R.zero
   for mon,coefficient in encoded:
    assert len(mon)<=1;value=F(coefficient);p+=s.QQ(value.numerator,value.denominator)*(variables[eta]if not mon else variables[mon[0]])
   graph.append(p)
  cache={():R.one}
  def restrict(poly):
   result=R.zero
   for mon,coefficient in poly:
    mon=tuple(mon)
    if mon not in cache:
     value=R.one
     for i in mon:
      value*=graph[i]
      if not value:break
     cache[mon]=value
    if cache[mon]:
     value=F(coefficient);result+=s.QQ(value.numerator,value.denominator)*cache[mon]
   return result
  qr=[restrict([[m,v]for m,v in row.items()])for row in Q];assert not any(qr)
  cr=[restrict(p)for p in c['raw350_cubic_terms']];dr=[restrict(p)for p in tails['all350_quartic_terms']]
  for source in G:
   assert sum((s.QQ(v.numerator,v.denominator)*cr[i]for i,v in cofactors[source].items()),R.zero)==0
  for source in zero:assert sum((s.QQ(v.numerator,v.denominator)*dr[i]for i,v in cofactors[source].items()),R.zero)==0
  for syz in syzygies:
   result=R.zero
   for source,poly in syz['actual350_cofactors']:
    for mon,coefficient in poly:
     value=F(coefficient);term=dr[source]if not mon else graph[mon[0]]*cr[source];assert len(mon)<=1
     result+=s.QQ(value.numerator,value.denominator)*term
   assert not result
 paths += [ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63_REVIEW.json',ROOT/'reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63_REVIEW.json',ROOT/'research/notes/Q_ORDINARY_QUARTIC_INITIAL_20261103_v63.md',ROOT/'research/notes/Q_ORDINARY_QUARTIC_INITIAL_20261103_v63_REVIEW.md']
 report=dict(id='Q_ORDINARY_QUARTIC_INITIAL_20261103_v63_REVIEW',all_passed=True,
  full202second_actual_DAG_and_physical4jet_receipt_audited_hash_pinned_not_repeated=True,
  bounded_mixed_producer_decoded_certificate_and_exact_native_inputs_outputs_match=True,
  invertible_original350_to_selected57_plus293_row_replacement_independently_verified=True,
  all219_global_cubic_kernel_and_actual350_cofactor_recipes_independently_verified=True,
  kernel_basis_QQ_independence_lower219_mod_prime=prime,old_exact_QQ_rank1911_proves_complete_kernel219_plus_zero158=True,
  independent_RAW350_quartic_plane_restrictions_and_all377_actual_carries_zero=True,
  full_true_initial_degree4_spanning_formula_accepted=True,
  entire_true_initial_through4_ideal_contains_both15planes_and_any12linear_section_dim_at_least3=True,
  full_initial_or_actual_local_dimension_lower15_or_actual15branch_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS COMPLETE377 kernel and actual raw4 cofactor restrictions; true initials through4 retain15planes only')

if __name__=='__main__':main()
