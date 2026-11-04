#!/usr/bin/env python3
"""Read-only individual common-monomial test of frozen high-weight identities.

No reconstruction, matrix, GB, point, field, normal graph or pivot calculation.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,collections,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_BIANCHI_MONOMIAL_DIVISION_PREFLIGHT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b):
 z=dict(a)
 for m,c in b.items():z[m]=z.get(m,Q())+c
 return{m:c for m,c in z.items()if c}
def mul(a,b):
 z={}
 for m,c in a.items():
  for n,d in b.items():k=tuple(x+y for x,y in zip(m,n));z[k]=z.get(k,Q())+c*d
 return{k:c for k,c in z.items()if c}
def decode(ts):
 z={}
 for ex,c in ts:
  assert len(ex)==3 and min(ex)>=0;z=add(z,{tuple(ex):Q(c)})
 return z
def encode(p):return[[list(m),str(c)]for m,c in sorted(p.items())]
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(500000000,500000000));signal.alarm(30);out.mkdir(parents=True,exist_ok=False)
 parent=ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63';op_path=parent/'operator.json';pr_path=parent/'report.json';pr=json.loads(pr_path.read_text());assert pr['all_full_raw_map_and_literal_identity_checks_passed']and sha(op_path)==pr['operator_sha256']
 for rel,want in pr['input_sha256'].items():assert sha(ROOT/rel)==want,rel
 data=json.loads(op_path.read_text());fixed=data['local_fixed_q_index'];columns=data['local_other_q_order'];weights=data['local_coordinate_weights'];eqw=data['equation_weights'];assert len(weights)==45 and len(eqw)==135
 chart_path=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json';chart=json.loads(chart_path.read_text());zero_coords=[i for i,w in enumerate(weights)if not w];assert zero_coords==[19,23]
 for i in zero_coords:
  assert all(Q(chart['T34'][i][j])==0 for j in chart['normal34'])
  assert all(Q(c)==0 for c in chart['U34'][i])
 F=[];source_terms=0;removed_source_terms=0
 for i,poly in enumerate(data['full_homogeneous135_equations']):
  f={}
  for mon,encoded in poly:
   assert len(mon)==2;core=decode(encoded)
   local=tuple(sorted(44 if v==45 else columns.index(v)for v in mon if v!=fixed))
   assert sum(weights[v]for v in local)==eqw[i]
   source_terms+=bool(core)
   if any(v in zero_coords for v in local):removed_source_terms+=bool(core);continue
   f[local]=add(f.get(local,{}),core)
  F.append({m:c for m,c in f.items()if c})
 assert len(F)==135
 checked_cofactor_terms=0;scans=[];candidates=[]
 for index,rec in enumerate(data['local_weight_pieces']):
  coeff={};discarded=0
  for g,var,encoded in rec['coefficients']:
   assert 0<=g<135 and(-1<=var<45)
   p=decode(encoded);mon=()if var==-1 else(var,)
   assert eqw[g]+sum(weights[v]for v in mon)==rec['weight'];checked_cofactor_terms+=bool(p)
   if var in zero_coords:discarded+=bool(p);continue
   if p:coeff.setdefault(g,{})[mon]=add(coeff.get(g,{}).get(mon,{}),p)
  coeff={g:{m:c for m,c in f.items()if c}for g,f in coeff.items()};coeff={g:f for g,f in coeff.items()if f}
  if rec['weight']<=4:continue
  allmons=[m for f in coeff.values()for m in f]
  common=collections.Counter(allmons[0])if allmons else collections.Counter()
  for mon in allmons[1:]:common&=collections.Counter(mon)
  factor=tuple(sorted(v for v,n in common.items()for _ in range(n)));fweight=sum(weights[v]for v in factor)
  item={'source_weight_piece_index':index,'operator_row':rec['operator_row'],'spatial_degree':rec['spatial_degree'],'spatial_monomial':rec['spatial_monomial'],'total_weight':rec['weight'],'retained_generator_columns':sorted(coeff),'cofactor_term_count':len(allmons),'discarded_only_verified_zero_nuisance_terms':discarded,'common_monomial_local_indices':list(factor),'common_monomial_weight':fweight,'zero_vector_after_nuisance_substitution':not allmons,'new_total_weight_if_divided':rec['weight']-fweight if factor else None}
  scans.append(item)
  if factor and fweight>0 and rec['weight']-fweight<=4:
   divided={};check={}
   for g,poly in coeff.items():
    divided[g]={}
    for mon,p in poly.items():
     cc=collections.Counter(mon);cc.subtract(common);assert min(cc.values(),default=0)>=0;mm=tuple(sorted(v for v,n in cc.items()for _ in range(n)));assert eqw[g]+sum(weights[v]for v in mm)==rec['weight']-fweight
     divided[g][mm]=p
     for fm,fc in F[g].items():check[tuple(sorted(mm+fm))]=add(check.get(tuple(sorted(mm+fm)),{}),mul(p,fc))
   check={m:c for m,c in check.items()if c};assert not check
   candidates.append({**item,'all135_columns_retained_during_common_factor_test':True,'divided_coefficients':[[g,list(m),encode(c)]for g,p in sorted(divided.items())for m,c in sorted(p.items())],'literal_135_polynomial_identity_after_division_verified':True})
 evidence={'source_operator_sha256':sha(op_path),'verified_zero_coordinates':zero_coords,'all135_F_local_weights_checked':True,'source_terms_checked':source_terms,'source_terms_removed_only_by_zero_nuisance':removed_source_terms,'all_nonzero_cofactor_term_weights_checked':checked_cofactor_terms,'higher_weight_piece_count':len(scans),'scans':scans,'candidates':candidates,'all135_generator_columns_used_for_every_common_factor_test':True,'NO_higher_equation_column_dropped':True,'NO_linear_combination_gcd_or_saturation_search':True}
 (out/'preflight.json').write_text(json.dumps(evidence,indent=2)+'\n');paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),op_path,pr_path,chart_path,ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63_REVIEW/report.json'];report={'id':ID,'all_read_only_source_weight_and_monomial_checks_passed':True,'source135_local_homogeneity_verified':True,'all_piece_cofactor_homogeneity_verified':True,'only_literal_zero_graph_nuisance_19_23_substituted':True,'all135_equation_columns_retained':True,'scanned_individual_higher_weight_identities':len(scans),'positive_common_monomial_rows':sum(bool(x['common_monomial_local_indices'])for x in scans),'divisions_reaching_weight_at_most4':len(candidates),'literal_divided_identities_checked':len(candidates),'tG6_not_divided_unless_all_coefficients_divisible':True,'individual_monomial_division_scope_only':True,'no_normal_graph_pivot_rank_matrix_GB_field_point_or_regeneration':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'preflight_sha256':sha(out/'preflight.json'),'seconds':time.monotonic()-start,'frozen':True};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
