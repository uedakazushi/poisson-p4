#!/usr/bin/env python3
"""A finite four-row MODULE reduction of actual D at generic B18.

Universal polarized Cayley-Hamilton is checked in two generic quadratic
algebras. All target word indices <=3 are frozen; algebraic coefficient
series remain specified by the prior rank2 source circuit, not expanded.
"""
import argparse,gzip,hashlib,itertools,json,math,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def partitions(items):
 if not items:yield [];return
 first,*rest=items
 for p in partitions(rest):
  yield [(first,),*p]
  for i in range(len(p)):yield [*p[:i],tuple(sorted((first,*p[i]))),*p[i+1:]]
def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic();tA,nA,tB,nB=s.symbols('tauA nuA tauB nuB');zero=(s.S.Zero,)*4;one=(s.S.One,s.S.Zero,s.S.One,s.S.Zero)
 def add(x,y):return tuple(s.expand(a+b)for a,b in zip(x,y))
 def scale(c,x):return tuple(s.expand(c*a)for a in x)
 def mul(x,y):
  a,b,c,d=x;e,f,g,h=y
  return tuple(map(s.expand,(a*e+nA*b*f,a*f+b*e+tA*b*f,c*g+nB*d*h,c*h+d*g+tB*d*h)))
 def prod(indices):
  ans=one
  for i in indices:ans=mul(ans,gens[i])
  return ans
 def trace(x):return s.expand(2*x[0]+tA*x[1]+2*x[2]+tB*x[3])
 gens=[tuple(s.symbols(f'aA{i} bA{i} aB{i} bB{i}'))for i in range(4)]
 traces={tuple(idx):trace(prod(idx))for k in range(1,5)for idx in itertools.combinations(range(4),k)}
 def mixed_sigma(indices):
  ans=0
  for blocks in partitions(list(indices)):
   term=(-1)**(len(indices)-len(blocks))
   for block in blocks:term*=math.factorial(len(block)-1)*traces[tuple(block)]
   ans+=term
  return s.expand(ans)
 rhs=zero;coefficients=[]
 for k in range(1,5):
  for indices in itertools.combinations(range(4),k):
   co=(-1)**(k+1)*math.factorial(4-k)*mixed_sigma(indices);rest=[i for i in range(4)if i not in indices];rhs=add(rhs,scale(co,prod(rest)));coefficients.append(dict(indices=list(indices),coefficient=str(co),remaining_product=rest))
 assert all(s.expand(a-b)==0 for a,b in zip(scale(24,prod(range(4))),rhs))
 # Independent universal rank2 CH and a degree3->words<=1 reduction.
 aa,bb,cc,dd,ee,ff,t,n=s.symbols('a b c d e f tau nu')
 def m2(x,y):return(s.expand(x[0]*y[0]+n*x[1]*y[1]),s.expand(x[0]*y[1]+x[1]*y[0]+t*x[1]*y[1]))
 def tr2(x):return s.expand(2*x[0]+t*x[1])
 def add2(x,y):return tuple(s.expand(a+b)for a,b in zip(x,y))
 def sc2(c,x):return tuple(s.expand(c*a)for a in x)
 b1,b2,b3=(aa,bb),(cc,dd),(ee,ff);t1,t2,t3=map(tr2,[b1,b2,b3]);c12=s.expand(t1*t2-tr2(m2(b1,b2)));c13=s.expand(t1*t3-tr2(m2(b1,b3)));c23=s.expand(t2*t3-tr2(m2(b2,b3)))
 pair=add2(add2(sc2(t1,b2),sc2(t2,b1)),(-c12,0));assert sc2(2,m2(b1,b2))==pair
 triple=add2(add2(sc2(2*tr2(m2(b1,b2)),b3),sc2(t1*t3,b2)),sc2(t2*t3,b1));triple=add2(triple,(-t1*c23-t2*c13,0));assert sc2(4,m2(m2(b1,b2),b3))==triple
 words=[[]]+[list(w)for k in range(1,4)for w in itertools.combinations_with_replacement(range(97),k)];counts=[sum(len(w)==k for w in words)for k in range(4)];assert counts==[1,97,4753,156849]and len(words)==161700
 archive=out/'paired-words-degree0-through3.json.gz';archive.write_bytes(gzip.compress(json.dumps(words,separators=(',',':')).encode(),mtime=0))
 dependency=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json';old=json.loads(dependency.read_text());assert old['all_passed']
 for f,h in old['input_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h,f
 recipe=dict(base='The SAME R=C[[r28]] from the frozen common-target rank2 reduction',normalizer_bases=['1,pA','1,pB'],quadratics=['pA²=tauA*pA+nuA','pB²=tauB*pB+nuB'],target_generator_pairs='b_j=(a_Aj,b_Aj,a_Bj,b_Bj) in B_A directsum B_B, 0<=j<97',coefficient_representation='Algebraic formal series specified by exact finite source+unitIFT+quadratic Weierstrass; not expanded.',M_generator_count=196,M_columns='(1,0) and(0,1), plus(z_j^A,0),(0,z_j^B); 4 rows in B free4 basis',U_generator_count=161700,U_columns='Products of paired b_j indexed by frozen word archive, including diagonal1',module_exact_sequence='0 -> U -> A_A directsum A_B -> actualD ->0',last_map='(a,b) -> lift(a)-lift(b) modulo sum of the TWO ACTUAL image kernels',ring_structure_warning='U is not an ideal of normalizer directproduct; use paired TARGET actions, not arbitrary normalizer multiplication',target_action_matrices='blockdiag([[a_Aj,nuA*b_Aj],[b_Aj,a_Aj+tauA*b_Aj]],[[a_Bj,nuB*b_Bj],[b_Bj,a_Bj+tauB*b_Bj]])',unit_class='class of (1,0) in A_A directsum A_B; paired target actions make this the full cyclic T-module actualD',rank4_product_reduction='24*b1*b2*b3*b4 = sum_nonemptyS (-1)^(|S|+1)*(4-|S|)!*mixed_sigma(S)*product_remaining_b',mixed_sigma='sum_setpartitions pi (-1)^(|S|-|pi|)*product_blocks (|block|-1)!*Tr(product_block b)',rank2_pair_reduction='2*b1*b2=Tr(b1)*b2+Tr(b2)*b1-(Tr(b1)Tr(b2)-Tr(b1b2))*1',actual_image_kernels_not_computed=True,full_primary_structure_not_claimed=True)
 (out/'finite-module-recipe.json').write_text(json.dumps(recipe,indent=2)+'\n')
 files=[Path(__file__),dependency,ROOT/'code/ambient_B18_quadratic_image_reduction_20261103_v63.py',ROOT/'research/notes/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63.md'];report=dict(id='AMBIENT_B18_FINITE_UNION_MODULE_20261103_v63',all_passed=True,universal_two_quadratic_algebras_used=True,rank4_CH_all4_components_symbolically_zero=True,rank2_CH_pair_and_triple_symbolically_verified=True,paired_word_counts_by_degree=counts,paired_union_U_finite_generator_bound=161700,ambient_submodule_M_generator_bound=196,normalizer_ambient_module_rank=4,exact_actualD_module_subquotient=True,finite_target_actions_recover_actual_ring=True,algebraic_coefficients_expanded=False,actualD_ring_type_or_Ass_or_depth_claimed=False,word_archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in files},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('B18 ACTUAL D FINITE FOUR-ROW MODULE REDUCTION PASS',report['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
