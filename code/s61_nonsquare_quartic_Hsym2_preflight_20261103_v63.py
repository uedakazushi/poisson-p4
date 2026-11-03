#!/usr/bin/env python3
"""Bounded lazy924-coefficient preflight of actual33 quartic rows.

Only structural support and formal-atom identities are tested. Core parameter
polynomials are not expanded and no rank or pointwise nonzero claim is made.
"""
import argparse,hashlib,itertools,json,time
from pathlib import Path
from fractions import Fraction
import sympy as s
from s61_nonsquare_lazy_quartic_oracle_20261103_v63 import Oracle,ROOT
ID='S61_NONSQUARE_QUARTIC_HSYM2_PREFLIGHT_20261103_v63'

def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);path=ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63/oracle.json';before=hashlib.sha256(path.read_bytes()).hexdigest();o,data=Oracle.load(path);oldnodes=len(o.nodes);rows=data['genuine_quartic_original_rows'];roots=dict(zip(rows,data['genuine_quartic_roots']));light=[0,1,2,3,4,6,7];pairs=list(itertools.combinations_with_replacement(light,2));monomials=[]
 for i,j in pairs:
  m=[0]*8;m[5]=1;m[i]+=1;m[j]+=1;monomials.append(tuple(m))
 raw=[[o.coefficient(roots[r],m)for m in monomials]for r in rows];H2=tuple(2 if i==5 else 0 for i in range(8));b={r:o.coefficient(roots[r],H2)for r in rows};hp=ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.json';h=json.loads(hp.read_text());K,R=o.K,o.R;a,e,d=o.a,o.e,o.d;parse=lambda text:K.from_expr(s.sympify(text));u,v=map(parse,h['row_combination_coefficients']);b77,b131=map(parse,h['coefficient_formulas']);B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2;unit=K((50*a*B*C)**2);assert u*b77+v*b131==unit
 # No symbolic core coefficient expansion is performed in these transformations.
 Rraw=dict(zip(rows,raw));chi=[o.lin([(u,x),(v,y)])for x,y in zip(Rraw[77],Rraw[131])];psi=[o.lin([(-b131,x),(b77,y)])for x,y in zip(Rraw[77],Rraw[131])];newlabels=[['psi',77,131]]+[r for r in rows if r not in[77,131]];new=[psi]+[[o.lin([(1,x),(-1/unit,o.mul(b[r],y))])for x,y in zip(Rraw[r],chi)]for r in rows if r not in[77,131]];assert len(new)==32
 # Formal polynomial normalization in independent atom symbols is a SAFE
 # identity detector. It does not identify a nonzero form with a nonzero core
 # function; nontrivial rational-function relations among atoms are ignored.
 cache={0:{}};totalterms=0;limit=250000
 class Bound(Exception):pass
 def atom_form(atom):
  c=o.atoms[atom]
  if c.numer.is_ground and c.denom.is_ground:return {():Fraction(str(c))}if c else{}
  return {(atom,):Fraction(1)}
 def product(p,q):
  if len(p)*len(q)>200000:raise Bound()
  result={}
  for x,c in p.items():
   for y,d in q.items():
    z=tuple(sorted(x+y));result[z]=result.get(z,Fraction())+c*d
  result={m:c for m,c in result.items()if c}
  if len(result)>2000:raise Bound()
  return result
 def normal(node):
  nonlocal totalterms
  if node in cache:return cache[node]
  if totalterms>=limit:raise Bound()
  nd=o.nodes[node];assert o.weights[node]==0
  if nd[0]=='c':ans=atom_form(nd[1])
  elif nd[0]=='m':ans=product(normal(nd[1]),normal(nd[2]))
  else:
   assert nd[0]=='l';ans={}
   for c,n in nd[1]:
    for m,vv in product(atom_form(c),normal(n)).items():ans[m]=ans.get(m,Fraction())+vv
   ans={m:c for m,c in ans.items()if c}
   if len(ans)>2000:raise Bound()
  totalterms+=len(ans);cache[node]=ans;return ans
 nfs={};refused=[]
 for node in sorted(set(x for row in raw+new for x in row)):
  try:nfs[node]=normal(node)
  except Bound:refused.append(node)
 def supports(matrix):
  exact_zero=[];upper=[]
  for row in matrix:
   zero=[j for j,n in enumerate(row)if not n or (n in nfs and not nfs[n])];exact_zero.append(zero);upper.append([j for j in range(28)if j not in zero])
  return exact_zero,upper
 z0,s0=supports(raw);z1,s1=supports(new)
 def proportional(matrix):
  fingerprints={};relations=[]
  for i,row in enumerate(matrix):
   if not all(n in nfs for n in row):continue
   terms=[(j,m,c)for j,n in enumerate(row)for m,c in sorted(nfs[n].items())]
   if not terms:continue
   scale=terms[0][2];key=tuple((j,m,c/scale)for j,m,c in terms)
   if key in fingerprints:
    old,oldscale=fingerprints[key];relations.append({'row':i,'proportional_to':old,'QQfactor':str(scale/oldscale)})
   else:fingerprints[key]=(i,scale)
  return relations
 evidence={'source_rows':rows,'light_pair_columns':[list(p)for p in pairs],'H_times_light_pair_monomials':[list(m)for m in monomials],'raw33_by28_coefficient_nodes':raw,'raw_H2_coefficient_nodes':[[r,b[r]]for r in rows],'chi_Hsym2_coefficient_nodes':chi,'H2_eliminated32_labels':newlabels,'H2_eliminated32_by28_coefficient_nodes':new,'raw_exact_zero_columns':z0,'raw_possible_support_columns':s0,'H2_eliminated_exact_zero_columns':z1,'H2_eliminated_possible_support_columns':s1,'raw_direct_QQ_proportional_rows':proportional(raw),'H2_eliminated_direct_QQ_proportional_rows':proportional(new),'uncompleted_formal_atom_normalizations':refused,'formal_atom_cache_terms':totalterms,'formal_atom_zero_is_exact_QQ_identity_but_nonzero_not_a_nonvanishing_certificate':True,'arithmetic':o.encoded()}
 text=json.dumps(evidence,separators=(',',':'))+'\n';assert len(o.nodes)<=200000 and len(text)<32*1024*1024;(out/'preflight.json').write_text(text);assert hashlib.sha256(path.read_bytes()).hexdigest()==before
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();paths=[Path(__file__),ROOT/'code/s61_nonsquare_lazy_quartic_oracle_20261103_v63.py',ROOT/'research/notes'/(ID+'.md'),path,hp]
 report={'id':ID,'all_bounded_symbolic_checks_passed':True,'Hsym2_queries':924,'required_additional_H2_queries':33,'core_parameter_polynomials_not_expanded':True,'actual_H2_unit_source_transform_applied':True,'source_oracle_unchanged':True,'original_DAG_nodes':oldnodes,'new_DAG_nodes':len(o.nodes),'artifact_bytes':len(text),'raw_exact_zero_coefficient_counts':[len(v)for v in z0],'raw_possible_support_counts':[len(v)for v in s0],'H2_eliminated_possible_support_counts':[len(v)for v in s1],'raw_direct_QQ_proportional_rows':proportional(raw),'H2_eliminated_direct_QQ_proportional_rows':proportional(new),'formal_atom_normalization_refusals':len(refused),'no_nonzero_coefficient_or_rank_claim':True,'no_field_point_GB_or_large_matrix_expansion':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'preflight_sha256':sha(out/'preflight.json'),'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_bounded_symbolic_checks_passed','new_DAG_nodes','artifact_bytes','raw_possible_support_counts','H2_eliminated_possible_support_counts','formal_atom_normalization_refusals','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
