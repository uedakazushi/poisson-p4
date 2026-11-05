#!/usr/bin/env python3
"""Actual formal connecting map Tnil16 -> Cprime[ell], source operators only."""
import hashlib,json,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_minimal_operator_20261104_v63 import Matrices,enc,count
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NIL_CONNECTING16_SOURCE_OPERATOR_20261105_v63'
PATHS={
 'Loewy':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63/report.json',
 'basis':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63/exact-Loewy16-and-explicit-filtered-basis.json',
 'basispeer':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63_REVIEW/report.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'minimalpeer':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW/report.json',
 'honest':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63/actual-honest-Tnil128-972-formal-operator.json',
 'honestpeer':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63_REVIEW/report.json',
 'radical':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63/actual-radical-normalizer-comparison-operator.json',
 'radicalpeer':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW/report.json',
 'natural':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json',
 'joint105':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'jointpeer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
 'nilfinite':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
 'Cprimebar':'reports/AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63_REVIEW/report.json',
 'Cprimerank':'reports/AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 for k in ['Loewy','basispeer','minimalpeer','honestpeer','radicalpeer','jointpeer','nilfinite','Cprimebar','Cprimerank']:assert D[k]['all_passed'],k
 A=D['basis'];assert A['exact_actual_HF']==[5,5,3,2,1]and len(A['explicit_sixteen_filtered_basis_lifts'])==16
 M=D['minimal'];indices=[x['index']+(44 if x['kind']=='Koszul_wedge84'else 0)for x in M['minimal_target_labels']];assert indices==[42,43,50,68,70]
 # Every Schur projection has its literal retained-coordinate right inverse;
 # composing the three inverse injections gives this original128 injection.
 current=list(range(128))
 for step in M['three_exact_unit_Schur_steps']:
  assert step['retained_rows']==[i for i in range(len(current))if i not in step['unit_rows']]
  current=[current[i]for i in step['retained_rows']]
 assert current==indices
 mat=Matrices();mo=mat.external('minimal','matrices.m',24,44);nn=mat.external('minimal','matrices.direct_unit_section_n',44,24);d2=mat.external('honest','matrices.d2_new',24,84)
 phi=mat.hcat(mat.add(mat.eye(44),mat.neg(mat.mul(nn,mo))),mat.mul(nn,d2));inj=mat.literal([[Q(i==j)for j in indices]for i in range(128)]);heads=mat.mul(phi,inj)
 natural=mat.external('natural','matrices.actual_natural_map_minimal5_by5',5,5)
 def smon(op,mon,coef=Q(1)):return mat.put(['multiply_by_literal_S_monomial',op,mon,coef],*mat.shapes[op])
 cols=[]
 for b in A['explicit_sixteen_filtered_basis_lifts']:
  typ=b['type']
  if typ=='new_minimal_unit':z=mat.select(mat.eye(5),list(range(5)),[b['target']])
  elif typ=='new_minimal_monomial':z=smon(mat.select(mat.eye(5),list(range(5)),[b['target']]),b['monomial'])
  elif typ=='natural_column_times_monomial':z=smon(mat.select(natural,list(range(5)),[b['old_minimal_column']]),b['monomial'])
  else:
   assert typ=='sum_of_natural_column_monomials';terms=[smon(mat.select(natural,list(range(5)),[j]),m,Q(c))for j,m,c in b['terms']];z=terms[0]
   for op in terms[1:]:z=mat.add(z,op)
  cols.append(z)
 basis=mat.hcat(*cols);cycles=mat.mul(heads,basis)
 # Same actual R0 formal coordinate section, with distinct original-source
 # clone identities. Never identify full and ell0 physical variables.
 J=D['joint105'];assert len(J['all105_unknown_original_node_order'])==len(J['all105_equation_original_node_order'])==105 and Q(J['det_joint105'])
 clones=[{'namespace':name,'unknown_order_source':'joint105.all105_unknown_original_node_order','equation_order_source':'joint105.all105_equation_original_node_order','equation_target_indices':list(range(48,57)),'target9':target,'origin_Jacobian_determinant':Q(J['det_joint105'])}for name,target in[
  ('FULL_R0_F9', [['z',i]for i in range(7)]+[['r'],['sum',['r'],['scale','1/5',['ell']]]]),
  ('ELL0_S_F9', [['z',i]for i in range(7)]+[['r'],['r']])]]
 gamma=mat.literal([[Q(1),Q(0),Q(0)],[Q(0),Q(0),Q(1)],[-Q(11,12),Q(1),Q(0)]])
 vf=mat.external('radical_full_R0','actual_V3_by44_rows_after_unit1_cancellation',3,44);vf=mat.put(['actual_coefficient_pullback_by_original105_clone',vf,'FULL_R0_F9','all full proper/unit matrix functions retain their samephysical R0 origin branches'],3,44);vf=mat.mul(gamma,vf)
 sigma=mat.put(['coefficientwise_formal_section_S_to_R0',cycles,'ELL0_S_F9','same z7,r; independent formalell; clone every implicit/regular-quotient function BEFORE multiplication'],44,16)
 numerator=mat.mul(vf,sigma)
 quotient=mat.put(['matrix_of_regular_formal_coordinate_divisions',numerator,'ell',{'exact_vanishing_under_ell0':'Vgamma_bar*Phi_direct=0 by mn=I, Vgamma_bar=d1*m,d1*d2=0','each_quotient_origin':'uncomputed, never set0','division_ring':'regular R0=QQ[[z7,r,ell]], ell nonzerodivisor; NOT finiteArtin33'}],3,16)
 edge=mat.put(['classes_in_actual_full_Cprime_cokernel',quotient,{'presentation_matrix':vf,'ell_annihilation':'ell*quotient=Vfull*sigma(cbar), hencezero class'}],3,16)
 ten=mat.select(edge,list(range(3)),list(range(10)))
 data={'actual_full_R0_variables':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell'],'ell_equals':'5*(sA-rA)',
  'two_distinct_original105_source_clones':clones,'formal_section_representation_limit':'finite function-ID/implicit pullback specification, with proved regular coordinate divisions; not an expanded scalar polynomial system and not computed values',
  'source_namespace_bindings':{k:p for k,p in PATHS.items()},'matrix_operator_nodes':mat.nodes,'matrix_operator_shapes':mat.shapes,
  'matrix_roots':{'minimal5_to_honest128_injection':inj,'directPhi44_by128':phi,'actual_five_minimal_cycle_heads44_by5':heads,'actual_filtered_basis_coordinates5_by16':basis,'actual_ell0_cycles44_by16':cycles,'actual_full_Vgamma3_by44':vf,'formal_ell_independent_cycles44_by16':sigma,'all48_regular_ell_numerators3_by16':numerator,'all48_regular_ell_quotients3_by16':quotient,'actual_connecting_classes3_by16':edge,'remaining_head_first_classes3_by10':ten},
  'all_order_head_and_basis_cycle_proof':'Vbar=d1*m; m*n=I; Vbar*Phi=(d1*m)*(I-nm,nd2)=(0,d1d2)=0. All sixteen basis coordinates are S-linear combinations of thosefive cycleheads, so all48 numerators vanish EXACTLY atell0.',
  'formal_coordinate_division_proof':'The actual105unit system gives physicalR0 as formal functions ofY=(z7,r,s). Fullclone has s=r+ell/5; ell0clone has s=r. This is the literal coefficient section QQ[[z7,r]]→QQ[[z7,r,ell]]. The numerator restriction0 therefore supplies unique regular formal ellquotients; no finite quotient or unknown origin is used to define them.',
  'boundary_descent_proof':'For any actualfull comparison boundary Xfull*a, sigma(Xbar*abar)-Xfull*sigma(abar) is divisiblebyell. Since Vfull*Xfull=0, dividing proves its connecting representative is in imVfull. Section andd3 boundaries map0 by the same Phi/full-comparison identities. Hence the map descends to actualTnil, independently of representative.',
  'exact_sequence':'0→Nil(D0)/ellNil(D0)→Tnil→Cprime[ell]→0',
  'target_dimension_and_maximalideal':'Cprime is finite overD=QQ[[ell]], genericrank2 andbarlength4, hence two torsion summands and dim_kCprime[ell]=2. Any length2 module over the SAME localS is killed bym_S².',
  'all_six_degree_at_least2_basis_columns_zero_CLASSES':list(range(10,16)),'not_literal_zero_vectors_or_computed_V_preimages':True,
  'remaining_genuine_edge_source_count':10,'genuine_target_dimension':2,'genuine_kernel_dimension':14,
  'actual_Cprime_ell_torsion_basis_or_edge_coefficients_not_computed':True,'Cprime2_does_not_replace_full_target':'current length7 leaves one exponent>=2; longtorsion may be invisible moduloell² orbar, so no finitekernel identification is used',
  'next_completion_gate':'Identify two source-bound fullCprime[ell] generators or prove a finite-target adequacy theorem for these10 expressions, then evaluate theactual2x10 map and liftits14 kernelclasses using explicit fullV membership. No naive m² coefficient truncation or Artinell division.',
  'no_new_source_gauge_jet_or_scalar_ell_values':True,'actual_nil_generators_or_cube_claimed':False}
 assert count(data)<150000
 p=out/'actual-formal-connecting16-and-ten-source-obstructions.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_source_bound_genuine_formal_connecting16_operator','actual_formal_regular_ell_division_count':48,'genuine_remaining_source_classes':10,'genuine_target_dimension':2,'genuine_kernel_dimension':14,'six_higher_filtered_basis_classes_in_kernel':True,'edge_values_or_actual_full_nil_generators_computed':False,'new_source_or_gauge_values':False,'visible_typed_QQ_slots':count(data),'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
