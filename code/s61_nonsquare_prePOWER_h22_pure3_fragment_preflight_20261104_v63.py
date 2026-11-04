#!/usr/bin/env python3
"""Exact integer l1/common-denominator MAJORANTS; no parameter products."""
import argparse,hashlib,itertools,json,math,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_H22_PURE3_FRAGMENT_PREFLIGHT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);paths=[]
 def load(stem,file='certificate.json'):
  p=ROOT/'reports'/stem/file;paths.append(p);return read(p)
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');h=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json')
 support=load('S61_NONSQUARE_PREPOWER_H22_KERNEL9_SUPPORT_PREFLIGHT_20261104_v63');nodes=q['scalar_DAG']['nodes'];cols=q['new22_target_original32_columns'][:17];rows=[1,3,6];chosen=[11,13,25];pivot=cols.index(22)
 roots=[q['remaining15x22_scalar_roots'][i][cols.index(j)]for i in rows for j in chosen+[22]];need=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i);n=nodes[i]
  if n[0]in['a','s','m']:todo+=n[1:]
 count=0;maxbits=0
 def atom(p,unit=Q(1)):
  nonlocal count,maxbits
  c=[Q(v)/unit for m,v in p];D=math.lcm(*(v.denominator for v in c))if c else 1;N=sum(abs(v.numerator)*(D//v.denominator)for v in c);count+=len(c);maxbits=max(maxbits,max((max(abs(v.numerator).bit_length(),v.denominator.bit_length())for v in c),default=0));return D,N
 def add(a,b):
  D=math.lcm(a[0],b[0]);return D,a[1]*(D//a[0])+b[1]*(D//b[0])
 def mul(a,b):return a[0]*b[0],a[1]*b[1]
 L={};denominator_atoms=[]
 for i in sorted(need):
  n=nodes[i]
  if n[0]=='z':v=(1,0)
  elif n[0]=='one':v=(1,1)
  elif n[0]in['checkpoint','accepted_pair_second','new_unit_atom']:
   a,d=cp['matrix25x32'][n[1]][n[2]]if n[0]=='checkpoint'else pd['Lframe_dual6x2'][n[1]][1]if n[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][n[1]]
   assert len(d)==1 and d[0][0][1:]==[0,0]and Q(d[0][1]);v=atom(a,Q(d[0][1]));denominator_atoms.append(dict(DAG_root=i,a_power=d[0][0][0],QQ_unit=d[0][1]))
  elif n[0]=='m':v=mul(L[n[1]],L[n[2]])
  else:assert n[0]in['a','s'];v=add(L[n[1]],L[n[2]])
  L[i]=v
 small={j:atom(p)for j,p in zip(h['original_target_labels'],h['source_h9'])};M=[]
 for i in rows:M.append([add(mul(small[22],L[q['remaining15x22_scalar_roots'][i][cols.index(j)]]),mul(L[q['remaining15x22_scalar_roots'][i][pivot]],small[j]))for j in chosen])
 result=(1,0)
 for perm in itertools.permutations(range(3)):
  x=(1,1)
  for i,j in enumerate(perm):x=mul(x,M[i][j])
  result=add(result,x)
 D,N=result;assert N.bit_length()<=292 and D.bit_length()<=89
 candidate=support['cheapest12_pure3_fixed_minor_candidates'][0];assert candidate['targets']==chosen and candidate['support_upper']==10080
 raw=10080*(len(str(N))+80)+4630*(max(len(str(x[1]))for r in M for x in r)+80)
 data=dict(source_rows15=rows,actual135=[18,39,68],fixed_target_original32_columns=chosen,h_pivot_original32=22,
  actual_entry_recipe='h22*P_i[j]-P_i[22]*Q9[j]; thethree Hzero rows use no H adjugate and need no rho normalization.',
  needed_scalar_DAG_roots=sorted(need),needed_scalar_count=len(need),all_monomial_a_QQ_denominator_atoms=denominator_atoms,
  atom_terms_counting_repeats=count,max_atom_QQ_coefficient_bits=maxbits,
  nine_entry_integer_common_denominator_and_l1_majorants=[[[str(D),str(N)]for D,N in r]for r in M],
  determinant_common_QQ_denominator=str(D),determinant_integer_coefficient_abs_upper=str(N),
  determinant_support_upper=10080,selected9_entry_support_upper=4630,
  determinant_integer_coefficient_bits=N.bit_length(),determinant_common_denominator_bits=D.bit_length(),raw_JSON_bytes_upper=raw,
  no_parameter_polynomial_coefficients_products_or_minor_evaluated=True,
  next_job='ONE fixed actualnine-entry/3det computation; ifzero STOPnoalternativecolumnsearch. Ifnonzero onlythispure3blockisproved, notfullkernel9.',
  chart_complement='FULLU intersect D(existingDelta5*h22) is the chosen9-kernel chart; this3blockitselfusesonlyFULLUandh22. No complementremoved.')
 (out/'certificate.json').write_text(json.dumps(data,indent=2)+'\n');paths += [Path(__file__),ROOT/f'research/notes/{ID}.md'];report=dict(id=ID,all_passed=True,source_only=True,frozen=True,
  scalar_nodes=len(need),literal_atom_terms=count,determinant_support_upper=10080,determinant_numerator_bits=N.bit_length(),determinant_QQ_denominator_bits=D.bit_length(),raw_JSON_bytes_upper=raw,
  no_new_coefficient_value_minor_rank_point_prime_GB=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},certificate_sha256=sha(out/'certificate.json'),seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
