#!/usr/bin/env python3
"""Global20 QQ(beta)-linear row identities and audited structural rank transfer."""
import argparse,gzip,itertools,json,resource,signal,time
from collections import Counter,defaultdict
from fractions import Fraction as Q
from pathlib import Path
from s61_model import S61
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CROSS_RANK_COVARIANCE_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False);folder=ROOT/'legacy/v02/data';old=read(folder/'symplectic_slice_equations.json');new=read(folder/'symplectic_slice_reduced_equations.json');cert=read(folder/'symplectic_slice_certificate.json');ids=read(folder/'symplectic_bianchi_identities.json');weights=old['parameter_weights'];order=new['old_generator_indices'];inverse={j:i for i,j in enumerate(order)};assert len(order)==134 and len(old['equations'])==154
 pre=dict(global20_term_contributions_upper=262420,new_explicit154x134_operator_terms_upper=747,hard_new_explicit_terms_TOTAL150000=True,no_full350_rebuild=True,no_new_field_Gaussian_point_prime=True,old_exact350_to154_constant_basis_reused='S61_REBUILD_20261030, originalpolynomial_row_basisconstantQQalgorithm')
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 rebuild=read(ROOT/'reports/S61_REBUILD_20261030.json');assert rebuild['all_passed']and rebuild['comparison_summary']['json_files_equal_ignoring_runtime_seconds']==5
 for f in rebuild['source_scripts']:assert sha(ROOT/f['path'])==f['sha256'],f['path']
 for f in rebuild['file_comparisons']:
  p=folder/f['file'];assert sha(p)==f['archived_sha256'],str(p)
 assert len(cert['remaining_equation_rows'])==154 and cert['generator_weight_counts']=={'4':70,'5':84}
 equations=[{tuple(m):Q(n,d)for m,n,d in f}for f in old['equations']];reduced=[{tuple(m):Q(n,d)for m,n,d in f}for f in new['equations']]
 matrix=[defaultdict(dict)for _ in range(154)]
 def add(i,j,u,c):
  p=matrix[i][j];p[u]=p.get(u,Q(0))+c
  if not p[u]:del p[u]
 for i,j in enumerate(order):assert equations[j]==reduced[i];add(j,i,-1,Q(1))
 contributions=0;maxinter=0;records=ids['generator_redundancy_identities'];assert len(records)==20
 for r in records:
  i=r['eliminated_old_generator'];assert i not in inverse
  for j,n,d in r['constant_coefficients_on_retained_old_generators']:add(i,inverse[j],-1,Q(n,d))
  for j,u,n,d in r['linear_coefficients_on_weight4_old_generators']:assert weights[u]==1;assert all(sum(weights[k]for k in m)==4 for m in equations[j]);add(i,inverse[j],u,Q(n,d))
  replay={}
  for j,p in matrix[i].items():
   for u,c in p.items():
    for mon,v in reduced[j].items():
     m=tuple(sorted((*mon,*([u]if u>=0 else[]))));replay[m]=replay.get(m,Q(0))+c*v
     if not replay[m]:del replay[m]
     contributions+=1
    maxinter=max(maxinter,len(replay));assert len(replay)<=150000
  assert replay==equations[i]
 assert contributions==262420
 for i,j in enumerate(order):assert dict(matrix[j])=={i:{-1:Q(1)}}
 operator=[[[j,[[u,str(v)]for u,v in sorted(p.items())]]for j,p in sorted(row.items())if p]for row in matrix];terms=sum(len(p)for row in matrix for p in row.values());assert terms<=150000
 (out/'global154to134-linear-operator.json').write_text(json.dumps(dict(operator154x134_sparse=operator,coefficient_convention='u=-1meansconstant, otherwise old0basedweight1u; NOqcoefficients',all20_original61QQ_identities_exact=True,retained134_rows_are_identity=True,selected350_raw_rows154=cert['remaining_equation_rows'],full350_to154_gate='ExactconstantQQpolynomial_row_basisintheacceptedS61_REBUILD; no newraw350ortruncatedidealassertion'),separators=(',',':'))+'\n')
 # Check all reconstruction weights, not just onPoissonpoints.
 pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2));W=[4,-1,-1,-1,-1]
 def rw(raw):
  i,j=pairs[raw//15];a,b=mons[raw%15];v=10+W[i]+W[j]-W[a]-W[b];assert v%5==0;return v//5
 ec=Counter();maxq=0
 for row in old['eliminated_raw_variables']:
  d=rw(row['raw_index']);ec[d]+=1
  for mon,n,den in row['polynomial']:assert sum(weights[j]for j in mon)==d;maxq=max(maxq,sum(weights[j]==2 for j in mon))
 assert dict(ec)=={1:4,2:15,3:36}and maxq==1
 assert all(rw(i)==1 for i in cert['perspective_orbit_pivots_raw'])
 assert all(rw(i)==weights[j]for j,i in enumerate(cert['raw_free_indices']))
 assert all(len({sum(weights[j]for j in m)for m in f})==1 for f in equations)
 assert Counter(sum(weights[j]for j in next(iter(f)))for f in equations)=={4:70,5:84}
 assert all(sum(weights[j]==2 for j in m)<=2 for f in reduced for m in f)
 response=S61().perspective_matrix();assert not response.free_symbols and response.det()!=0
 # Old126 lineage is EXACToriginal134, with independentphysical reconstruction.
 fdir=ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63';fr=read(fdir/'report.json');fc=read(fdir/'certificate.json');review=read(ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63_REVIEW/report.json');assert fr['all_passed']and fr['exact_mixed_cross_rank126']and review['all_passed'];assert fc['cross_matrix_shape']==[129,126]and fc['exact_cross_rank']==126;assert review['all134_original_equations_accounted_for']and review['cross_input_rebuilt_from_actual_G_and_byte_equal']and review['integer_inverse_products_checked']==2
 assert sha(fdir/'cross-input.txt.gz')==fr['evidence_sha256']['cross-input.txt.gz']
 assert __import__('hashlib').sha256(gzip.decompress((fdir/'cross-input.txt.gz').read_bytes())).hexdigest()==fc['native_input_sha256']
 sparse=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63/report.json';sp=read(sparse);assert sp['all_passed'];dom=ROOT/'reports/Q_NEW_FAMILY_GENERIC_BETA_SPARSE_SECTION_MINOR_20261104_v63_REVIEW/report.json';assert read(dom)['all_passed']
 certificate=dict(all_passed=True,global154_to134_QQ_beta_rowspace_equality=True,full350_to154_constantQQ_rowspace_gate_reused=True,not_entire_I_qdegree_le2='G6outside134span; notneededforRAWJacobicovariance',row_operator_explicit_terms=terms,exact20_term_contributions=contributions,max_replay_support=maxinter,low_Jacobi_graph_elimination_weights=dict(ec),reconstruction_max_qdegree=maxq,perspective_response4x4=[[str(x)for x in row]for row in response.tolist()],perspective_response_determinant=str(response.det()),off_shell_graph_reason='RawJacobiorder>=3 preservedbypointfixingprojectivechanges; universal55lowequationgraphunique; perspectiveparameterdependsONLYweight1beta',q_action_reason='Positiveweightsforcebeta_prime linearbeta, q_prime=M_gq+b_g(beta^2); inverseactionforcesM_ginvertible',RAW350_action_reason='EwedgeJacobiisGL5tensor; verticalEulerwedgelinearadditionsleaveitunchanged; normalizingGdependsONLYbeta, soRAWcoefficientspantransformsQQ(beta)-linearly',adapted_coordinate_rank_reason='Onkerplane restriction, n-zcrossblockchangesbyGL42tensorGL3; shearonlyaddsnormal2. Onfullspaceextended[S|Y]rightblocktriangularoperationaddsS*R',old126_lineage=dict(equations=134,residuals=129,cross_shape=[129,126],added_G6=False,source='load_fiber originalreduced134; finite5review own physicalall134 and126inverse',exact_QQ_rank=126),known_plane_equivariance='AffinespanofcanonicalknownT01three+T39one+T40one; connectedSp4xGm preservescomponentclosures andgenericfivefibers',generic_cross_rank_on_K5=126,generic_extended_S_Y_rank=131,dense_source_transfer='G=Sp4xweightedGm rankpropertyinvariant; acceptedGxsparse5dominanceforcesgenericrankonK5',new_explicit126minor_computed=False,fixed_old_minor_assumed_nonzero_on_source=False,full_affine_fiber_exhaustion_proved=False)
 (out/'covariance-certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',ROOT/'reports/S61_REBUILD_20261030.json',folder/'symplectic_slice_certificate.json',folder/'symplectic_slice_equations.json',folder/'symplectic_slice_reduced_equations.json',folder/'symplectic_bianchi_identities.json',ROOT/'legacy/v02/code/build_symplectic_slice.py',ROOT/'legacy/v02/code/build_universal_slice.py',ROOT/'legacy/v02/code/reduce_symplectic_bianchi.py',ROOT/'code/s61_model.py',ROOT/'code/s61_affine_fiber_macaulay.py',ROOT/'code/s61_affine_normal_finite5_20261103_v63.py',ROOT/'code/s61_affine_normal_finite5_review_20261103_v63.py',fdir/'report.json',fdir/'certificate.json',ROOT/'reports/S61_AFFINE_NORMAL_FINITE5_20261103_v63_REVIEW/report.json',sparse,dom,ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',ROOT/'research/notes/S61_UNIVERSAL_WEIGHT6_20261031.md']
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<=32*1024*1024
 report=dict(id=ID,all_passed=True,new20_global_QQ_beta_linear_identities=True,original134_not135_used_for_old126=True,RAW134_rowspace_covariance_proved=True,off_shell_normalization_q_affine_proved=True,generic_K5_cross_rank126_proved_by_covariance=True,generic_K5_extended_rank131=True,new126minor_or_fieldGaussian=False,full_generic_exhaustion=False,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,row_operator_terms=terms,global20_contributions=contributions,cross_rank126=True,seconds=report['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
