#!/usr/bin/env python3
"""Independent carry replay, exhaustive15-axis check, and exact-rank audit."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def sub(row,other,c):
 for m,v in other.items():
  x=row.get(m,Fraction(0))-c*v
  if x:row[m]=x
  else:row.pop(m,None)

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--cubic-replay',type=Path,required=True)
 ap.add_argument('--QQ-replay',type=Path,required=True);ap.add_argument('--mod-replay',type=Path,required=True)
 ap.add_argument('--minimal-replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 directory=ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63';manifest=load(directory/'manifest.json');paths=[Path(__file__),directory/'manifest.json']
 for name,h in manifest['source_sha256'].items():p=ROOT/name;assert sha(p)==h,name;paths.append(p)
 for name,h in manifest['evidence_sha256'].items():p=directory/name;assert sha(p)==h,name;paths.append(p)
 fresh=load(args.cubic_replay/'report.json');original=load(directory/'report.json')
 ignore={'seconds','peak_rss_kb'}
 assert {k:v for k,v in fresh.items()if k not in ignore}=={k:v for k,v in original.items()if k not in ignore}
 QQ=load(directory/'exact-rank/report.json');freshQQ=load(args.QQ_replay/'report.json')
 for k in ['block_ranks','quadratic_linear_multiple_rank_QQ','full_cubic_initial_rank_QQ','new_cubic_rank_mod_quadratic_multiples_QQ','ordinary_Hilbert_values_through3_QQ','integer_matrix_sha256']:
  assert QQ[k]==freshQQ[k]
 mod=load(directory/'modular-rank/report.json');assert mod['block_ranks']==load(args.mod_replay/'report.json')['block_ranks']
 assert QQ['quadratic_linear_multiple_rank_QQ']==1894 and QQ['full_cubic_initial_rank_QQ']==1911
 assert QQ['ordinary_Hilbert_values_through3_QQ']==[1,35,573,5859]
 with gzip.open(directory/'cubic-initials.json.gz','rt')as f:d=json.load(f)
 oldpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz'
 with gzip.open(oldpath,'rt')as f:old=json.load(f)
 Q=[{(i,j):Fraction(c)for i,j,c in row}for row in old['all350_quadratic_rows']]
 C=[{tuple(m):Fraction(c)for m,c in row}for row in d['raw350_cubic_terms']]
 expected={r['source']:{tuple(m):Fraction(c)for m,c in r['terms']}for r in d['cubic_initial_relations']}
 basis={};carries={};seen=set()
 for op in d['exact_constant_row_elimination']:
  i=op['source'];assert i not in seen;seen.add(i);q,c=Q[i].copy(),C[i].copy()
  for mon,v in op['reductions']:
   mon=tuple(mon);v=Fraction(v);assert q[mon]==v and min(q)==mon
   sub(q,basis[mon],v);sub(c,carries[mon],v)
  if'new_pivot'in op:
   m=tuple(op['new_pivot']);v=Fraction(op['normalizer']);assert q[m]==v and v.numerator%32003 and v.denominator%32003
   basis[m]={a:b/v for a,b in q.items()};carries[m]={a:b/v for a,b in c.items()}
  else:assert not q and c==expected[i]
 assert seen==set(range(350))and len(basis)==57 and len(expected)==293 and sum(bool(v)for v in expected.values())==135
 active={d['variables'].index(n)for n in d['plane_active_variables']}
 factor={tuple(sorted([d['variables'].index('z8')]*2+[d['variables'].index(n)])):Fraction(v)for n,v in [('z14',12),('z17',1),('z3',6)]}
 for i,row in expected.items():
  restricted={m:c for m,c in row.items()if set(m)<=active};v=Fraction(3,16)if i==234 else Fraction(-1,16)if i==349 else 0
  assert restricted=={m:v*c for m,c in factor.items()}if v else not restricted
 assert all(not(set(m)<=active)for row in Q for m in row)
 # Exhaust all axis subsets of the15 variables not forbidden by square
 # monomials. No quadratic-graph maximality claim is taken on trust.
 forbidden={i for row in Q for (i,j),c in row.items()if i==j};allowed=[i for i in range(35)if i not in forbidden];assert len(allowed)==15
 edges={tuple(sorted([allowed.index(i),allowed.index(j)]))for row in Q for(i,j),c in row.items()if i in allowed and j in allowed}
 masks=[(1<<i)|(1<<j)for i,j in edges];maximum=0;maxima=[]
 for m in range(1<<15):
  size=m.bit_count()
  if size<maximum or any(m&e==e for e in masks):continue
  if size>maximum:maximum=size;maxima=[]
  maxima.append({allowed[i]for i in range(15)if m>>i&1})
 assert maximum==13 and maxima==[active]
 minimum=ROOT/'reports/Q_ORDINARY_CUBIC_MINIMAL_GENERATORS_20261103_v63/report.json';mn=load(minimum);mf=load(args.minimal_replay/'report.json')
 assert {k:v for k,v in mn.items()if k!='seconds'}=={k:v for k,v in mf.items()if k!='seconds'}
 for file,h in mn['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 paths += [oldpath,minimum,ROOT/'research/notes/Q_ORDINARY_CUBIC_INITIAL_20261103_v63_REVIEW.md',ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63_REVIEW.json']
 result=dict(id='Q_ORDINARY_CUBIC_INITIAL_20261103_v63_REVIEW',all_passed=True,
  fresh_all350_exact_QQ_degree3_DAG_and_physical_Jacobi_replay_matches=True,
  fresh_exact_integer_QQ_ranks_and_modular_blocks_match=True,all350_saved_constant_carries_independently_replayed=True,
  true_degree3_initial_space_exhaustion_proof_accepted=True,quadratic_multiple_rank=1894,true_cubic_initial_rank=1911,
  new_minimal_cubic_generators=17,ordinary_Hilbert_through3=[1,35,573,5859],
  all32768_allowed_coordinate_axis_subsets_independently_checked=True,unique_quadratic13_axis_plane_verified=True,
  whole_false13plane_restriction='z8^2*(12*z14+z17+6*z3)',maximum_coordinate_plane_in_J_through3_dimension=12,
  explicit17_generators_mod1995_multiples_independently_replayed=True,
  full_initial_cone_Krull_dimension_or_actual_local_dimension_claimed=False,completed_local_ideal_generated_by74_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS true cubic initial component and exact restricted divisor; no full-cone claim')

if __name__=='__main__':main()
