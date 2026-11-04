#!/usr/bin/env python3
"""Bounded full80/140 Jacobi map and homogeneous135 Bianchi operator data.

No quartic Macaulay matrix, point, field or Groebner basis is computed.
"""
import argparse,collections,hashlib,itertools,json,resource,signal,time
from pathlib import Path
from fractions import Fraction
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61,X
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_binary_uniform_cubic_20261101 import full_raw_jacobi_maps
from elliptic_weight1_dominance import projective_affine
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(33);out.mkdir(parents=True,exist_ok=False);R,a,e,d=ring('a,e,d',s.QQ);zero=R.zero;one=R.one
 def add(target,key,value):
  if value:
   z=target.get(key,zero)+value
   if z:target[key]=z
   else:target.pop(key,None)
 def encode(p):return [[list(ex),str(c)]for ex,c in p.items()]
 model=S61();chartpath=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json';chart=json.loads(chartpath.read_text());oldpath=ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json';old=json.loads(oldpath.read_text())['linear'];extra_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';extra=json.loads(extra_path.read_text());ids_path=ROOT/'legacy/v02/data/symplectic_bianchi_identities.json';ids=json.loads(ids_path.read_text());beta={int(k[1:])-1:s.sympify(v)for k,v in chart['beta'].items()};betap={i:R.from_expr(c)for i,c in beta.items()}
 zero_nuisance=[i for i,w in enumerate(chart['local_coordinate_weights'])if w==0]
 assert zero_nuisance==[19,23]
 for i in zero_nuisance:
  assert all(R.from_expr(s.sympify(chart['T34'][i][j]))==zero for j in chart['normal34'])
  assert all(R.from_expr(s.sympify(c))==zero for c in chart['U34'][i])
 print('START_RAW80_140',flush=True);raw,i4,i5,lower=full_raw_jacobi_maps(ROOT,model);assert[len(r)for r in raw]==[80,140]and lower==60;print('RAW_MAPS_EXACT',time.monotonic()-start,flush=True)
 old_to_new={old:i for i,old in enumerate(model.reduced['old_generator_indices'])};eliminated={r['eliminated_old_generator']:r for r in ids['generator_redundancy_identities']}
 def old_generator(j):
  if j in old_to_new:return{old_to_new[j]:one}
  rec=eliminated[j];out={}
  for k,n,den in rec['constant_coefficients_on_retained_old_generators']:add(out,old_to_new[k],R(s.Rational(n,den)))
  for k,u,n,den in rec['linear_coefficients_on_weight4_old_generators']:add(out,old_to_new[k],R(s.Rational(n,den))*betap[u])
  return out
 triples=list(itertools.combinations(range(4),3));Jparts=[];raw_encoded=[]
 for entries,oldids in zip(raw,[i4,i5]):
  T={tr:{}for tr in triples};encoded=[]
  for tr,mon,coords in entries:
   coeff={}
   for k,c in coords.items():
    cc=R(s.Rational(c))
    for g,v in old_generator(oldids[k]).items():add(coeff,g,cc*v)
   for g,v in coeff.items():add(T[tr],(g,mon),v)
   encoded.append({'triple':list(tr),'spatial_monomial':list(mon),'old_basis_coordinates':[[k,str(c)]for k,c in sorted(coords.items())],'reduced134_coordinates':[[g,encode(v)]for g,v in sorted(coeff.items())]})
  Jparts.append(T);raw_encoded.append(encoded)
 two,_=specialize(model,beta,True);q=s.symbols('q0:45');values=[s.S.Zero]*61
 for i,w in enumerate(model.reduced['parameter_weights']):values[i]=beta[i]if w==1 else q[two.index(i)]
 y=s.symbols('y0:4');phys=projective_affine(model.reconstruct(values),X,y).applyfunc(s.expand);P={}
 for i,j in itertools.combinations(range(4),2):
  comp={}
  for mon,c in s.Poly(phys[i,j],*y).terms():
   for qm,coef in s.Poly(c,*q).terms():
    if not coef:continue
    assert sum(qm)<=1;var=qm.index(1)if sum(qm)else 45;add(comp,(var,tuple(mon)),R.from_expr(coef))
  P[i,j]=comp
 print('FULL_HOMOG_P_COEFFICIENTS',sum(len(p)for p in P.values()),time.monotonic()-start,flush=True)
 def pair(i,h):
  if i==h:return{}
  return P[i,h]if i<h else{k:-v for k,v in P[h,i].items()}
 def trivector(indices):
  if len(set(indices))<3:return 0,()
  return(-1)**sum(indices[i]>indices[j]for i in range(3)for j in range(i+1,3)),tuple(sorted(indices))
 T={tr:dict(Jparts[0][tr])for tr in triples}
 for tr in triples:
  for key,v in Jparts[1][tr].items():add(T[tr],key,v)
 operator={}
 def term(mon,g,var,value):
  if mon not in operator:operator[mon]={}
  add(operator[mon],(g,var),value)
 for i in range(4):
  other=tuple(k for k in range(4)if k!=i)
  for h in range(4):
   for (var,pm),pc in pair(i,h).items():
    for (g,tm),tc in T[other].items():
     if tm[h]:
      mon=tuple(pm[k]+tm[k]-int(k==h)for k in range(4));term(mon,g,var,((-1)**i)*tm[h]*pc*tc)
 for i,k in itertools.combinations(range(4),2):
  other=tuple(v for v in range(4)if v not in[i,k])
  for h in range(4):
   sign,tr=trivector((h,*other))
   if not sign:continue
   for (var,pm),pc in P[i,k].items():
    if not pm[h]:continue
    for (g,tm),tc in T[tr].items():
     mon=tuple(pm[j]+tm[j]-int(j==h)for j in range(4));term(mon,g,var,((-1)**(i+k))*sign*pm[h]*pc*tc)
 operator={m:v for m,v in operator.items()if v};assert all(sum(m)in[4,5,6]for m in operator)
 # Full q/t homogeneous equations, including G6. These are the SAME source135.
 _,F=specialize(model,beta,True);Fh=[]
 for poly in F:
  h={}
  for mon,c in poly.items():assert len(mon)<=2;add(h,tuple(sorted(mon+(45,)*(2-len(mon)))),R.from_expr(c))
  Fh.append(h)
 assert len(Fh)==135
 def verify(row):
  check={}
  for(g,var),c in row.items():
   for mon,v in Fh[g].items():add(check,tuple(sorted((var,*mon))),c*v)
  assert not check
 records=[]
 for D in[4,5,6]:
  count=0
  for mon in itertools.product(range(D+1),repeat=4):
   if sum(mon)!=D:continue
   row=operator.get(mon,{});verify(row);count+=1;records.append({'spatial_degree':D,'spatial_monomial':list(mon),'coefficients':[[g,var,encode(c)]for(g,var),c in sorted(row.items())]})
  assert count=={4:35,5:56,6:84}[D];print('FULL_IDENTITIES_VERIFIED_DEG',D,time.monotonic()-start,flush=True)
 # The q-part of each degree4 row is precisely the old canonical35, not a replacement.
 for monindices,stored in zip(extra['affine_quartic_monomial_order'],extra['canonical_35_Bianchi_vectors']):
  mon=tuple(monindices.count(i)for i in range(4));actual={(g,var):c for(g,var),c in operator.get(mon,{}).items()if var<45};wanted={}
  for index,c in stored:j,i=divmod(index,70);add(wanted,(i,j),R(s.Rational(c)))
  assert actual==wanted
 extra_row={(134,45):one}
 for index,c in extra['selected_non_Bianchi_syzygy']:
  j,i=divmod(index,70);add(extra_row,(i,j),-R(s.Rational(c)))
 verify(extra_row)
 checkpoint={'raw_Jacobi_maps':raw_encoded,'homogeneous_bivector_q_t_coefficients':[[i,j,[[v,list(m),encode(c)]for(v,m),c in sorted(row.items())]]for(i,j),row in P.items()],'full_homogeneous135_equations':[[[list(m),encode(c)]for m,c in sorted(f.items())]for f in Fh],'full175_Schouten_rows':records,'extra_tG6_row':[[g,var,encode(c)]for(g,var),c in sorted(extra_row.items())],'all175_and_extra_literal_identities_verified':True,'verified_zero_weight_nuisance_coordinates':zero_nuisance}
 checkpoint_text=json.dumps(checkpoint,separators=(',',':'))+'\n';assert len(checkpoint_text)<32*1024*1024;(out/'source-verified-checkpoint.json').write_text(checkpoint_text)
 fixed,cols=old['fixed'],old['cols'];cw=chart['local_coordinate_weights'];rw=chart['equation_weights'];pieces=[];counts=collections.Counter();quartic_rows=[i for i,w in enumerate(rw)if w==4 and i not in chart['selected34_rows']];assert len(quartic_rows)==33
 for idx,rec in enumerate(records+[{'spatial_degree':'extra_tG6','spatial_monomial':[],'coefficients':[[g,var,encode(c)]for(g,var),c in sorted(extra_row.items())]}]):
  groups={}
  row=operator.get(tuple(rec['spatial_monomial']),{})if idx<len(records)else extra_row
  for(g,var),c in row.items():
   local=-1 if var==fixed else(44 if var==45 else cols.index(var));sigma=rw[g]+(0 if local==-1 else cw[local]);target=groups.setdefault(sigma,{})
   add(target,(g,local),c)
  for weight,group in sorted(groups.items()):
   counts[(str(rec['spatial_degree']),weight)]+=1;pieces.append({'operator_row':idx,'spatial_degree':rec['spatial_degree'],'spatial_monomial':rec['spatial_monomial'],'weight':weight,'coefficients':[[g,v,encode(c)]for(g,v),c in sorted(group.items())]})
 block=[];blocklabels=[];nuisance_terms_discarded=0
 for i,rec in enumerate(pieces):
  if rec['weight']!=4:continue
  row={}
  for g,var,c in rec['coefficients']:
   if g in quartic_rows:
    if var==-1:row[g]=c
    else:assert var in zero_nuisance;nuisance_terms_discarded+=1
  block.append([row.get(g,[])for g in quartic_rows]);blocklabels.append(i)
 # Preserve compact source/formula maps, raw135 and all175 operator identities.
 evidence={'original_source_count':135,'local_fixed_q_index':fixed,'local_other_q_order':cols,'local_coordinate_weights':cw,'equation_weights':rw,'raw_J3_J4_verified_counts':[80,140],'raw_lower_zero_coefficients':lower,'old154_weight4_generator_indices':i4,'old154_weight5_generator_indices':i5,'raw_Jacobi_maps':raw_encoded,'homogeneous_bivector_q_t_coefficients':[[i,j,[[v,list(m),encode(c)]for(v,m),c in sorted(row.items())]]for(i,j),row in P.items()],'full_homogeneous135_equations':[[[list(m),encode(c)]for m,c in sorted(f.items())]for f in Fh],'full175_Schouten_rows':records,'extra_tG6_row':[[g,var,encode(c)]for(g,var),c in sorted(extra_row.items())],'local_weight_pieces':pieces,'local_weight_piece_counts':[[D,w,n]for(D,w),n in sorted(counts.items())],'genuine_quartic_original_rows':quartic_rows,'weight4_quartic_constant_block_row_piece_indices':blocklabels,'weight4_quartic_constant_block':block,'zero_weight_nuisance_terms_dropped_after_exact_graph_zero_proof':nuisance_terms_discarded,'matrix_entry_coefficient_core_order':['a','e','d'],'homogeneous_variable_order':'q0,...,q44,t (index45); dehom q_fixed=1','all175_full_identities_and_extra_tG6_literally_verified':True,'all_degree4_q_parts_equal_stored35canonical':True,'NO_normal_graph_transport_or_319_deletion_performed':True}
 text=json.dumps(evidence,separators=(',',':'))+'\n';assert len(text)<32*1024*1024;(out/'operator.json').write_text(text)
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),chartpath,oldpath,extra_path,ids_path,ROOT/'code/s61_binary_uniform_cubic_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'legacy/v02/code/reduce_symplectic_bianchi.py',ROOT/'legacy/v02/code/build_symplectic_slice.py',ROOT/'legacy/v02/code/build_universal_slice.py',ROOT/'legacy/v02/data/symplectic_slice_certificate.json',ROOT/'legacy/v02/data/symplectic_slice_equations.json',ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json']
 report={'id':ID,'all_full_raw_map_and_literal_identity_checks_passed':True,'raw_Jacobi_coefficients_verified':[80,140],'full_Schouten_degree_rows':[35,56,84],'all175_literal135_syzygies_verified':True,'homogenized_t_G6_identity_verified':True,'all_degree4_q_parts_match_old35canonical':True,'actual_local_weight_piece_counts':evidence['local_weight_piece_counts'],'weight4_quartic_constant_block_shape':[len(block),33],'weight4_block_nonzero_rows':sum(any(c for c in row)for row in block),'operator_bytes':len(text),'verified_zero_nuisance_coordinates':zero_nuisance,'discarded_zero_nuisance_multiplier_terms':nuisance_terms_discarded,'source_checkpoint_sha256':sha(out/'source-verified-checkpoint.json'),'no_rank_pivot_GB_field_point_or_quartic_matrix_job':True,'no80_or_new_unit_or_fullquartic_claim':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'operator_sha256':sha(out/'operator.json'),'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
