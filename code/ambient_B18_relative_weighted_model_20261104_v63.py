#!/usr/bin/env python3
"""Frozen-source weight binding for the relative B18 full-intersection model."""
import argparse,collections,gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/report.json';gr=read(gp);assert gr['all_exact_physical_tangent_and_stabilizer_weight_checks_passed'];gc=read(gp.parent/'certificate.json')
 for rel,h in gr['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 assert sha(gp.parent/'certificate.json')==gr['certificate_sha256']
 rp=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json';r=read(rp);assert r['all_passed']and r['both_normalizers_finite_free_rank']==2
 tp=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63/certificate.json.gz';tc=read(tp)
 op=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz';oc=read(op)
 mp=ROOT/'reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63/report.json';mr=read(mp);assert mr['all_passed']and mr['whole_target_kernel_99_relations']
 mcp=mp.parent/'relations99.json.gz';assert sha(mcp)==mr['evidence_sha256']['relations99.json.gz'];mc=read(mcp)
 nilp=ROOT/'reports/AMBIENT_B18_NILCOTANGENT_BASIS_20261103_v63/report.json';assert read(nilp)['all_passed']
 weights=gc['physical150_torus_weights'];O=s.Matrix(oc['physical150_by19_scalar_then_GL']);T=s.Matrix(tc['common150_by39_normalized_tangent']);C39=tc['common39_physical_coordinate_rows']
 Rrows=tc['common29_even_target_rows'];assert Rrows==r['common28_linear_target_rows']+[r['common29th_target_row']]
 def weight(v):
  ws={weights[i]for i,c in enumerate(v)if c};assert len(ws)==1;return ws.pop()
 rbasis=list(O[Rrows,:].T.rref()[1]);assert len(rbasis)==19
 groupcuts=[Rrows[j]for j in rbasis];basecuts=[j for j in Rrows if j not in groupcuts];det=s.factor(O[groupcuts,:].det());assert det
 assert collections.Counter(weights[j]for j in groupcuts)=={0:15,1:4}
 assert collections.Counter(weights[j]for j in basecuts)=={0:9,1:1}
 positivebase=[j for j in basecuts if weights[j]==1];assert len(positivebase)==1
 params=s.symbols('p k0:9');P,K=params[0],params[1:];records={}
 for name in['CA','CB']:
  fam=r['families'][name];old=tc['families'][name]
  J=s.Matrix(old['full150_by28_physical_even_linear_J']);E=s.Matrix(old['full150_by29_even_cotangent_embedding']);B=s.Matrix(old['full150_by10_actual_odd_cotangent_embedding']);W=s.Matrix(old['physical_pure_kernel_second_W'])
  jw=[weight(J[:,j])for j in range(28)];ew=[weight(E[:,j])for j in range(29)];bw=[weight(B[:,j])for j in range(10)]
  assert collections.Counter(jw)=={0:23,1:5};assert collections.Counter(ew)=={0:24,1:5};assert collections.Counter(bw)=={0:8,1:2};assert weight(W)==0
  assert E[Rrows,:]==s.eye(29)
  M=s.Matrix(fam['source_kernel_Lie_matrix']);D=s.diag(0,0,-1,0,0);assert D*M-M*D==s.zeros(5)
  assert s.Rational(fam['source_velocity'][3])==0
  extras=fam['source_fiber_basis_indices'];kw=[int(j==3)for j in extras];assert kw==[0,0,1,0,0,0,0,0,0]
  assert jw[19:]==kw
  counts=collections.Counter();nonzero=collections.Counter();maxordinary=0;maxp=0
  for i,text in enumerate(fam['core_full150_coordinates']):
   poly=s.Poly(s.sympify(text),*params)
   if poly.is_zero:continue
   nonzero[weights[i]]+=1
   for mon,c in poly.terms():
    assert mon[3]==weights[i] # params=(p,k0,k1,k2,k3,...), positive k2 only.
    counts[weights[i]]+=1;maxordinary=max(maxordinary,sum(mon));maxp=max(maxp,mon[0])
   assert weights[i]in[0,1]
  records[name]={'physical_J28_weights':jw,'actual_even_R29_weights':ew,'actual_odd_minimal10_weights':bw,'kernel_square_W_weight':0,'kernel_Lie_matrix_commutes_with_fixed_torus':True,'kernel_source_velocity_has_no_positive_fiber_component':True,
   'full_source_core10_parameter_order':list(map(str,params)),'full_source_core10_weights':[0,*kw],
   'all150_exact_frozen_core_polynomials_homogeneous':True,'core_nonzero_physical_coordinate_counts_by_weight':dict(nonzero),'core_terms_by_weight':dict(counts),'core_max_ordinary_parameter_degree':maxordinary,'core_max_p_degree':maxp,
   'source_normalizer29_weights':{'0':24,'1':5},'source_normalizer_after19_common_cuts_weights':{'0':9,'1':1},'tracefree_rank2_generator_weight':0}
 # Canonical ambient126 is a torus representation. Its old normalized basis
 # consists of weight eigenvectors; no fresh physical tensor is generated.
 canonical=s.Matrix(r['canonical_reconstruction_from126']);rawfree=r['common_target126_raw_coordinate_indices'];assert canonical.shape==(150,126)
 aw=[weight(canonical[:,j])for j in range(126)];ac=collections.Counter(aw);assert ac=={0:36,1:60,2:24,3:6}
 tcounters=collections.Counter(weights[j]for j in C39);assert tcounters=={0:32,1:7}
 graphcounts={i:ac[i]-tcounters[i]for i in range(4)};assert graphcounts=={0:4,1:53,2:24,3:6}
 pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2));physicalkeys=[(pair,mon)for pair in pairs for mon in mons]
 cubicraw=[i for i,w in enumerate(weights)if w==3];assert len(cubicraw)==6
 # Compatibility of the old abstract99 presentation with the same weight
 # counts. This does NOT identify its CA/CB relative coordinate changes.
 names=[*mc['base_R_core_variables'],*mc['odd_target_variables']];vars=s.symbols(' '.join(names));vmap=dict(zip(names,vars));assignment={n:int(n=='w'or n in['Y7','Y8'])for n in names}
 relationcounts=collections.Counter()
 for text in mc['all99']:
  poly=s.Poly(s.sympify(text,locals=vmap),*vars);ww={sum(mon[j]*assignment[n]for j,n in enumerate(names))for mon,c in poly.terms()}
  assert len(ww)==1;relationcounts[ww.pop()]+=1
 assert relationcounts=={0:64,1:32,2:3}
 evidence={'actual_four_term_B18_point_only':True,'common_physical_group_cuts19':groupcuts,'common_group_cut_minor':str(det),'remaining_actual_even_base10_raw_coordinates':basecuts,'positive_even_base_raw_coordinate':positivebase[0],
 'group_cut_weights':{'0':15,'1':4},'remaining_base_weights':{'0':9,'1':1},'component_records':records,
 'canonical126_weights':dict(ac),'common_tangent39_weights':dict(tcounters),'eliminated87_ambient_graph_weights':graphcounts,
 'possible_cubic_discrepancy_raw_coordinates':[{'raw':i,'pair':list(physicalkeys[i][0]),'quadratic_monomial':list(physicalkeys[i][1])}for i in cubicraw],
 'abstract99_weight_compatibility_counts':dict(relationcounts),'abstract99_coordinate_alignment_NOT_assumed':True,
 'exact_structural_model':'D_perp = completion of B[Y0,Y1,Y2]/K, B=C[[X0,...,X16]], K homogeneous in Y and generated in degrees0..3',
 'minimal_component_kernel_bound_positive_degree':2,'whole_D_bound_positive_degree':3,'whole_D_number_of_generators_upper_bound':285,'graph_discrepancies_retained_all87':True,'possible_cubic_graph_discrepancies':6,
 'positive_degree_monomial_counts':[1,3,6,10],'total_coefficient_families':20,
 'fixed_zero_weight_D':{'presentation':'C[[X17]]/K0, K0=(K modulo positive3)','cotangent_dimension':17,'linear_visible_nil_cotangent_dimension':1,'explicit_linear_class':'u63','raw_u63_not_itself_claimed_nilpotent':True,'nil_lift_exists_from_weight0_projection_of_accepted_actual_nil_class':True},
 'uncomputed':'all actual zero17 coefficient functions, relative CA/CB graph differences, full nil ideal, socle, depth19/20, degree2 improvement',
 'no_GB_elimination_new_point_primary_solver_source_regeneration':True}
 (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),gp,gp.parent/'certificate.json',rp,tp,op,mp,mcp,nilp]
 report={'id':ID,'all_frozen_source_weight_and_sparsity_checks_passed':True,'normalizers_transverse10_weights':{'0':9,'1':1},'actual_common_even_base10_weights':{'0':9,'1':1},'actual_odd10_weights':{'0':8,'1':2},'component_kernel_positive_degree_bound':2,'whole_D_positive_degree_bound':3,'whole_D_possible_cubic_graph_discrepancies':6,'fixed_zero_D_embedding_dimension':17,'fixed_zero_D_linear_nil_cotangent_dimension':1,'actual_coefficient_functions_or_full_D_nil_depth_computed':False,'formal_orbit_IFT_truncation_not_assumed':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())
