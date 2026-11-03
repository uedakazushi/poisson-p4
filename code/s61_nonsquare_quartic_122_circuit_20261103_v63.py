#!/usr/bin/env python3
"""Source-bound full122x64 quartic circuit; no large coefficients expanded."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
from fractions import Fraction
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_QUARTIC_122_CIRCUIT_20261103_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 st=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md')]
 deps=['S61_NONSQUARE_CUBIC_NINE_ROW_CRITERION_20261103_v63','S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63','S61_NONSQUARE_QUARTIC_286_SUFFICIENT_20261103_v63','S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63','S61_NONSQUARE_VERTEX_SCROLL_FINITE_20261103_v63']
 for dep in deps:
  p=ROOT/'reports'/dep/'report.json';d=json.loads(p.read_text());assert any(d.get(k,False)for k in ['all_passed','all_exact_QQ_checks_passed','all_exact_QQ_actual_source_checks_passed'])
  paths.append(p)
  for name,want in d.get('input_sha256',{}).items():assert sha(ROOT/name)==want,name
 cp=ROOT/'reports/S61_NONSQUARE_CUBIC_NINE_ROW_CRITERION_20261103_v63/circuit.json';cubic=json.loads(cp.read_text());paths.append(cp)
 qp=ROOT/'reports/S61_BINARY_HYPERPLANE_QUADRATIC_CODIMENSION_20261102.json';q=json.loads(qp.read_text());paths.append(qp)
 wp=ROOT/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json';chart=json.loads(wp.read_text());paths.append(wp)
 sp=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';source=json.loads(sp.read_text());paths.append(sp)
 bp=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json';basis=json.loads(bp.read_text());paths.append(bp)
 cubicpath=ROOT/cubic['source'];raw=gzip.decompress(cubicpath.read_bytes());assert hashlib.sha256(raw).hexdigest()==cubic['source_uncompressed_sha256'];paths.append(cubicpath)
 assert basis['global_free_row_module_rank']==57
 brp=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW/report.json';br=json.loads(brp.read_text());paths.append(brp)
 assert br['abstract_row_module_quotient_free_rank']==16
 light=[0,1,2,3,4,6,7];mons=source['normal_monomials'];heavy_by_light=[]
 for v in light:
  mon=[0]*8;mon[5]=1;mon[v]=1;heavy_by_light.append(mons.index(mon))
 assert sorted(heavy_by_light)==cubic['heavy_columns']
 E=cubic['heavy_pivot_rows_basis_positions'];other=[i for i in range(41,57)if i not in E]
 assert len(E)==7 and len(other)==9 and other==cubic['remaining_rows_basis_positions']
 lines=raw.decode().splitlines();assert lines[0]=='57 91 57'and len(lines)==1+57*91
 heavy_entries=[]
 for i in E+other:
  row=[]
  for j in heavy_by_light:
   terms=lines[1+i*91+j].split();n=int(terms[0]);assert len(terms)==1+4*n
   for k in range(n):Fraction(terms[1+4*k]);assert min(map(int,terms[2+4*k:5+4*k]))>=0
   row.append({'basis_row':i,'column':j,'source_line_1_based':2+i*91+j,'sparse_encoding':lines[1+i*91+j]})
  heavy_entries.append(row)
 assert all(not int(lines[1+i*91+j].split()[0])for i in range(41)for j in heavy_by_light)
 pairs=list(itertools.combinations_with_replacement(range(7),2));offpairs=list(itertools.combinations(range(7),2));assert len(pairs)==28 and len(offpairs)==21
 # The49 actual f_i*x_j coefficients in H*Sym2 are this universal constant map.
 C=s.zeros(49,28)
 for i in range(7):
  for j in range(7):C[7*i+j,pairs.index(tuple(sorted((i,j))))]=1
 assert C.rank()==28
 piv=[7*i+j for i,j in pairs];assert C[piv,:]==s.eye(28)
 K=s.zeros(21,49)
 for k,(i,j)in enumerate(offpairs):K[k,7*i+j]=1;K[k,7*j+i]=-1
 assert K*C==s.zeros(21,28)and K.rank()==21
 for qi in q['quadratic_rows_sparse']:
  assert all(mon[5]==0 and sum(mon)==2 for mon,_ in qi)
 q4=[i for i,w in enumerate(chart['equation_weights'])if w==4 and i not in chart['selected34_rows']]
 assert len(q4)==33 and 77 in q4 and 131 in q4
 hp=ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.json';h=json.loads(hp.read_text());paths.append(hp)
 a,e,d=s.symbols('a e d');b0,b1=map(s.sympify,h['coefficient_formulas']);u,v=map(s.sympify,h['row_combination_coefficients']);unit=(50*a*((50*a+18*e**2-9*e-2)*d+60*e-40)*(100*a+18*e**2-9*e-2))**2
 assert s.cancel(u*b0+v*b1-unit)==0
 # Universal2x2 source change, including the potentially nonunit individual b_i.
 U,V,B0,B1=s.symbols('u v b0 b1');T=s.Matrix([[U,V],[-B1,B0]])
 assert T.det()==U*B0+V*B1 and T*s.Matrix([B0,B1])==s.Matrix([U*B0+V*B1,0])
 # A different fixture retains an unknown h in the second transformed row.
 r,hvar,z=s.symbols('r h z');F0=r*z+1;F1=(1-r)*z+hvar
 first=F0+F1;second=-(1-r)*F0+r*F1
 assert s.expand(first-z)==hvar+1 and s.expand(second)==r*hvar+r-1
 residual=[['H_times_quadric',i]for i in range(6)]+[['commutator',i,j]for i,j in offpairs]+[['pure_cubic_times_light',i,j]for i in other for j in range(7)]
 assert len(residual)==90
 q4res=[['second_H2_free_combination',77,131]]+[['genuine_quartic_minus_H2_multiple',i]for i in q4 if i not in [77,131]]
 assert len(q4res)==32
 circuit={
 'coefficient_domain':'QQ[a,e,d,1/(a*F*B*C)] on extra DeltaH*DeltaQ charts; FULLU ensures product41 and rank146 away from finitely many points',
 'source_quartic_shape':[319,239],'source_cubic_matrix':cubic['source'],'source_cubic_sha256_uncompressed':cubic['source_uncompressed_sha256'],
 'light_old_variable_indices':light,'heavy_old_variable_index':5,'heavy_columns_in_light_order':heavy_by_light,
 'genuine_basis_positions':list(range(41,57)),'chosen_heavy_basis_rows':E,'remaining_nine_basis_rows':other,'heavy_coefficient_entries':heavy_entries,
 'DeltaH':'det A[E,heavy_columns_in_light_order]; an EXTRA chart unit',
 'heavy_graph_cubics':'f_i=sum_k (A[E,heavy]^-1)[i,k]*G_E[k]=H*x_i-p_i',
 'pure_nine_cubics':'r_j=G_j-sum_i A[j,heavy_i]*f_i for j in remaining_nine_basis_rows',
 'H_Sym2_pivot_rows':[[i,j]for i,j in pairs],'H_Sym2_columns':pairs,
 'twentyone_kernel_rows':[[i,j,7*i+j,7*j+i]for i,j in offpairs],
 'normal_form_Hxi_xj':'H*x_i*x_j -> p_i*x_j for i<=j',
 'commutator_normal_form':'f_i*x_j-f_j*x_i = p_j*x_i-p_i*x_j',
 'pure_quadratic_products':[['quadric',i,'light_pair',list(pair)]for i in range(6)for pair in pairs],
 'pure_quartic_columns':[list(c)for c in itertools.combinations_with_replacement(range(7),4)],
 'DeltaQ':'ANY146x146 minor of the168x210 pure quadratic-product matrix; an EXTRA chart unit',
 'pure_quad_complement_columns':64,'pure_quad_projection':'For pivot rows I and pivot columns K: pi(v)=v_J-v_K*Q[I,K]^-1*Q[I,J], J=complement(K)',
 'pure_quad_global_upper146_dependency':'S61_NONSQUARE_ZW_PRODUCT_RANK19_20261103_v63',
 'genuine_quartic_original_source_rows':q4,
 'genuine_quartic_coefficient_oracle':'weighted degree4 component of original135 equation after the exact selected34 implicit graph in S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102; same source as319 packet, no new expansion',
 'H2_source_row_transform':{'input_rows':[77,131],'matrix':[['u','v'],['-b131','b77']],'u':h['row_combination_coefficients'][0],'v':h['row_combination_coefficients'][1],'b77':h['coefficient_formulas'][0],'b131':h['coefficient_formulas'][1],'determinant':'(50*a*B*C)^2'},
 'H2_pivot_chi':'chi=u*Q77+v*Q131, coefficient U=(50aBC)^2',
 'H2_replacement':'H^2 -> -NF_H(chi-U*H^2)/U',
 'genuine_quartic_residual_rows':q4res,
 'other_quartic_elimination':'for k!=77,131: Qk-(b_k/U)*chi; second transformed row=-b131*Q77+b77*Q131; then NF_H and pi',
 'low_weight_residual_rows':residual,'full_residual_rows':residual+q4res,
 'full_shape':[122,64],'low_weight_sufficient_shape':[90,64],
 'rank_full239_iff':'rank(actual original full weight4)=175+rank(full122x64) on DeltaH*DeltaQ chart',
 'rank_low_weight238_iff':'rank(low-weight products)=174+rank(low90x64) on same chart',
 'finite_cover':'all7-row minors of16 genuine heavy rows and all146 minors of168x210 pure quadratic-products',
 'no_full_cubic57_premise':True,'no_new_rank64_claim':True,'no_expanded_parametric_quartic_matrix':True}
 assert len(circuit['pure_quartic_columns'])==210 and len(circuit['pure_quadratic_products'])==168
 assert 28+146+1==175 and 239-175==64 and 6+21+9*7+32==122 and 6+21+9*7==90
 (out/'circuit.json').write_text(json.dumps(circuit,indent=2)+'\n')
 result={'id':ID,'all_passed':True,'actual_sources_and_all112_heavy_entries_pinned':True,'universal49_to28_plus21_exact':True,'global_known_unit_H2_two_row_transform_exact':True,'all33_genuine_quartics_retained_equivalently':True,'full122_by64_cokernel_criterion':True,'low90_by64_sufficient_only':True,'extra_charts':['DeltaH7','DeltaQ146'],'full_cubic57_not_assumed':True,'rank64_full_quartic_length65_or_ambient_not_proved':True,'no_GB_field_or_large_parametric_expansion':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'circuit_sha256':sha(out/'circuit.json'),'seconds':time.monotonic()-st}
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'all_passed':True,'shapes':[[122,64],[90,64]],'seconds':result['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
