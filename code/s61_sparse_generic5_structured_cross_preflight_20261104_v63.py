#!/usr/bin/env python3
"""Read-only support/source preflight for a structured cross inverse.

Does NOT expand a new cross block or perform a field solve/rank test.
"""
import argparse,gzip,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_STRUCTURED_CROSS_PREFLIGHT_20261104_v63'
def read(p):
 b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pd(p):return {tuple(m):Q(c)for m,c in p}
def add(a,b,c=Q(1)):
 z=dict(a)
 for m,v in b.items():
  z[m]=z.get(m,Q(0))+c*v
  if not z[m]:del z[m]
 return z
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
 idir=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63';fdir=ROOT/'reports/S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63';ndir=ROOT/'reports/S61_SPARSE_GENERIC5_NORMAL42_EXACT_PRESENTATION_20261104_v63'
 known=read(idir/'known-five-input.json.gz');raw=read(idir/'pulled134-input.json.gz');full=read(fdir/'full134-direct-circuit.json.gz');normal=full['normal_q_rows42'];idx={j:i for i,j in enumerate(normal)};piv=full['pivot_q_rows'];assert piv==[3,6,13]
 center=list(map(pd,known['q_center']));d=list(map(pd,known['Q52_minus_center']));q51=list(map(pd,known['Q51']));q52=list(map(pd,known['Q52']))
 assert all(add(q52[i],center[i],Q(-1))==d[i]and add(q51[i],center[i],Q(-1))=={m:-v/2 for m,v in d[i].items()}for i in range(45))
 dp=[d[i]for i in piv];a0={(1,0,0,1,1):Q(1,144)};assert dp==[a0,{},{}]
 support={}
 for name,v in [('d',d),('Q51',q51),('Q52',q52)]:
  mat=[]
  for eq in raw['equations']:
   row=[set()for _ in normal]
   for mon,e,c in eq:
    if len(mon)==1 and name!='d'and mon[0]in idx:row[idx[mon[0]]].add(tuple(e))
    if len(mon)!=2:continue
    i,j=mon
    if i in idx:row[idx[i]].update(tuple(x+y for x,y in zip(e,t))for t in v[j])
    if j in idx:row[idx[j]].update(tuple(x+y for x,y in zip(e,t))for t in v[i])
   mat.append(row)
  support[name]=dict(input_terms=sum(map(len,v)),input_nonzero=sum(bool(z)for z in v),support_upper_total=sum(len(c)for row in mat for c in row),support_upper_max=max(len(c)for row in mat for c in row),entry_slots=134*42,structurally_possible_nonzero_entries=sum(bool(c)for row in mat for c in row),parameter_degree_upper=max((sum(m)for row in mat for c in row for m in c),default=0),not_exact_nonzero_or_term_counts=True)
 qterms=[(mon,e,Q(c))for eq in raw['equations']for mon,e,c in eq if len(mon)==2]
 denq=math.lcm(*(c.denominator for _,_,c in qterms));dend=math.lcm(*(c.denominator for z in d for c in z.values()));qmax=max(abs(c*denq).numerator for _,_,c in qterms);dmax=max(abs(c*dend).numerator for z in d for c in z.values())
 bits=(2*len(qterms)*qmax*dmax).bit_length();denbits=(denq*dend).bit_length();weights=(1,2,3,5,7);Tdeg=max(sum(w*(x+y)for w,x,y in zip(weights,e,m))for mon,e,c in qterms for j in mon for m in d[j])
 assert support['d']['support_upper_total']==13506 and bits==89 and denbits==44 and Tdeg==20
 source=dict(coefficient_field='QQ(H,L,s0,q0,t0)',q_pivot_rows=piv,d_definition='Q52-q_center=(2/3)*(Q52-Q51)',d_pivot=[[[list(m),str(c)]for m,c in z.items()]for z in dp],E_column0='144*d/(H*q0*t0)',proof='Accepted W=P*[NT0,NT1,d] and E=W*adj(Wpivot)/Gamma give E*d_pivot=d. Since d_pivot=(H*q0*t0/144,0,0), E[:,0]=d/a0 exactly.',actual_45_coordinate_difference_identities_checked=True,no_E_or_W_fraction_expansion=True,all_inverses_parameter_only=True)
 (out/'first-column-source-contraction.json').write_text(json.dumps(source,indent=2)+'\n')
 operator=read(ndir/'typed-normal42-operator.json.gz');xr=operator['actual_X130x126'];node=operator['nodes'][xr];assert node['op']=='matrix_from_entries'and len(node['inputs'])==130*126
 recipe=dict(actual_X_source=str((ndir/'typed-normal42-operator.json.gz').relative_to(ROOT)),actual_X_matrix_root=xr,first42_block_entry_refs=[[node['inputs'][126*i+j]for j in range(42)]for i in range(130)],first42_block_formula='X0=C_d_res/a0; C_d_res=rawPolar_d-lambda*selectedPolar_d, withsameactual five-plane selector; G6 row included from accepted X.',rawPolar_d_formula='For F=sum Q_ij q_i q_j+..., C_d[k]=sum_j Q_kj*d_j with diagonal doubled; normalq indices are current42 std_nonpivot.',remaining84_columns='X[:,42:126], exactexistingunexpandedatoms',block_split_not_executed='Z=X0, B=Xremaining. LZ=(ZtZ)^-1Zt; Y=(I-ZLZ)B; LY=(YtY)^-1Yt. Leftinverse y -> (LZ*(y-B*LY*y),LY*y).',block_injectivity='Accepted rankX126 makes Z injective42 and Y injective84; formallyrealK proves both Gram determinants nonzero.',not_a_computed_leftinverse=True,not_a_total_cost_reduction_proof=True,largest_unresolved_coefficient_gates=['actual five-plane quotient scalar coefficients useall692953 terms unless kept external','42- and84-block scalar inverse growth notbounded','two toric plane directions remaincoupled; no commonweightdecomposition supplied'])
 (out/'source-bound-42-plus84-recipe.json').write_text(json.dumps(recipe,indent=2)+'\n')
 costs=dict(actual_K5=support,raw_d_polar_coefficient_numerator_bits_upper=bits,raw_d_polar_common_denominator_bits_upper=denbits,actual_plane_functional_input_terms=692953,actual_plane_functional_exceeds150000=True,whole_K5_47column_input_terms_upper=692953+13506,formal_curve_distinct=dict(map='(T,T²,T³,T⁵,T⁷)',old_S5_actual_collapsed_terms=55594,new_raw_d_polar_Tsupport_upper=13506,total_upper=55594+13506,raw_d_polar_Tdegree_upper=Tdeg,det47_Tdegree_upper=5*111+42*Tdeg,not_a_generic_parameter_coefficient_bound=True,no_formal_rank_job_recommended_merely_to_reprove_known_rank=True),no_new_coefficients_or_rank_or_inverse_computed=True,no_known_plane_or_normal_unknown_discarded=True)
 (out/'support-and-inverse-frontier.json').write_text(json.dumps(costs,indent=2)+'\n')
 paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',idir/'known-five-input.json.gz',idir/'pulled134-input.json.gz',fdir/'full134-direct-circuit.json.gz',ndir/'report.json',ndir/'typed-normal42-operator.json.gz',ROOT/'reports/S61_SPARSE_GENERIC5_CROSS_RANK_COVARIANCE_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_FIVE_UNIVARIATE_20261104_v63/native-certificate.json.gz',ROOT/'reports/S61_FOUR_DIMENSIONAL_INVERSE_FORM_BIANCHI_20261104_v63_REVIEW/report.json']
 report=dict(id=ID,all_passed=True,scope='ActualsourcefirstEcolumncontraction+supportpreflight; NOstructuredinversecoefficients/rankjob',actual_E_first_column_simplification_proved=True,known_45_coordinate_difference_identities_checked=True,raw_d_polar_terms_upper13506=True,full_plane_selector_692953_terms_not_hidden=True,no_K5_field_Gaussian_Gram_expansion_point_prime_GB_or_formal_curve_job=True,source_actual_normal42_adapter_bound=True,structured_inverse_still_open=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)
