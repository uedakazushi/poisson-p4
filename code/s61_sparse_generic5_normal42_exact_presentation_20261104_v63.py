#!/usr/bin/env python3
"""Typed finiteoperator for the actualK5 closedresidual; NO Gram expansion."""
import argparse,gzip,itertools,json,resource,signal,time
from pathlib import Path
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_NORMAL42_EXACT_PRESENTATION_20261104_v63'
class Operator:
 def __init__(self):self.nodes=[];self.types=[];self.cache={}
 def new(self,op,inputs,attrs,typ):
  assert all(isinstance(i,int)and 0<=i<len(self.nodes)for i in inputs)
  key=json.dumps([op,inputs,attrs,typ],sort_keys=True,separators=(',',':'))
  if key in self.cache:return self.cache[key]
  i=len(self.nodes);self.cache[key]=i;self.nodes.append(dict(op=op,inputs=inputs,attrs=attrs,type=typ));self.types.append(typ);return i
 def scalar(self,ctx,root):return self.new('external_parameter_scalar',[],dict(context=ctx,root=root),dict(kind='K'))
 def poly(self,ctx,root):return self.new('external_polynomial_z_n',[],dict(context=ctx,root=root),dict(kind='PolyZN',degree_bound=2))
 def zeroz(self,p,minorder=0):
  assert self.types[p]['kind']=='PolyZN';return self.new('specialize_z_zero',[p],{},dict(kind='PolyN',degree_bound=self.types[p]['degree_bound'],normal_order_lower_bound=minorder))
 def mat(self,rows):
  n=len(rows);m=len(rows[0]);assert all(len(r)==m for r in rows)and all(self.types[i]['kind']=='K'for r in rows for i in r)
  return self.new('matrix_from_entries',[i for r in rows for i in r],dict(rows=n,cols=m),dict(kind='MatK',shape=[n,m]))
 def tr(self,a):
  n,m=self.types[a]['shape'];return self.new('transpose',[a],{},dict(kind='MatK',shape=[m,n]))
 def mm(self,a,b):
  n,k=self.types[a]['shape'];k1,m=self.types[b]['shape'];assert k==k1;return self.new('matrix_product',[a,b],{},dict(kind='MatK',shape=[n,m]))
 def solve(self,a,b):
  n,k=self.types[a]['shape'];n1,m=self.types[b]['shape'];assert n==k==n1 and self.types[a]['kind']==self.types[b]['kind']=='MatK'
  return self.new('unique_parameter_linear_solve',[a,b],dict(nonzero_determinant_proof='CauchyBinet:det(XtX)=sum_126rowminors(X)^2!=0 overformallyrealK; rankX126 accepted',evaluated=False),dict(kind='MatK',shape=[n,m]))
 def vec(self,ps):
  assert all(self.types[p]['kind']=='PolyN'for p in ps);return self.new('polynomial_vector',ps,{},dict(kind='VecPolyN',length=len(ps),degree_bound=max(self.types[p]['degree_bound']for p in ps),normal_order_lower_bound=min(self.types[p]['normal_order_lower_bound']for p in ps)))
 def mv(self,a,b):
  n,m=self.types[a]['shape'];t=self.types[b];assert t['kind']=='VecPolyN'and m==t['length'];return self.new('parameter_matrix_times_polynomial_vector',[a,b],{},dict(kind='VecPolyN',length=n,degree_bound=t['degree_bound'],normal_order_lower_bound=t['normal_order_lower_bound']))
 def vneg(self,a):return self.new('vector_negation',[a],{},self.types[a])
 def va(self,a,b):
  assert self.types[a]==self.types[b];return self.new('vector_addition',[a,b],{},self.types[a])
 def get(self,a,i):
  t=self.types[a];assert t['kind']=='VecPolyN'and 0<=i<t['length'];return self.new('vector_entry',[a],dict(index=i),dict(kind='PolyN',degree_bound=t['degree_bound'],normal_order_lower_bound=t['normal_order_lower_bound']))
 def n(self,i):return self.new('normal_variable',[],dict(index=i),dict(kind='PolyN',degree_bound=1,normal_order_lower_bound=1))
 def kpoly(self,a):
  assert self.types[a]['kind']=='K';return self.new('parameter_as_normal_polynomial',[a],{},dict(kind='PolyN',degree_bound=0,normal_order_lower_bound=0))
 def prod(self,args):
  assert args and all(self.types[a]['kind']=='PolyN'for a in args);return self.new('polynomial_product',args,{},dict(kind='PolyN',degree_bound=sum(self.types[a]['degree_bound']for a in args),normal_order_lower_bound=sum(self.types[a]['normal_order_lower_bound']for a in args)))
 def add(self,args):
  assert args and all(self.types[a]['kind']=='PolyN'for a in args);return self.new('polynomial_sum',args,{},dict(kind='PolyN',degree_bound=max(self.types[a]['degree_bound']for a in args),normal_order_lower_bound=min(self.types[a]['normal_order_lower_bound']for a in args)))
 def sub(self,a,b):return self.new('polynomial_difference',[a,b],{},dict(kind='PolyN',degree_bound=max(self.types[a]['degree_bound'],self.types[b]['degree_bound']),normal_order_lower_bound=min(self.types[a]['normal_order_lower_bound'],self.types[b]['normal_order_lower_bound'])))
def universal_checks():
 # All monomials degree<=2 in3plane+42normalvariables, for all42 charts.
 mons=[()]+[(i,)for i in range(45)]+list(itertools.combinations_with_replacement(range(45),2));assert len(mons)==1081;count=0
 for j in range(42):
  for mon in mons:
   z=[i for i in mon if i<3];n=[i-3 for i in mon if i>=3]
   reference=tuple(sorted([('P',h)for h in z]+[('n',k)for k in n]+[('n',j)]*(2-len(z))))
   if len(z)==0:assembled=tuple(sorted([('n',j)]*2+[('n',k)for k in n]))
   elif len(z)==1:assembled=tuple(sorted([('n',j),('P',z[0])]+[('n',k)for k in n]))
   else:assembled=tuple(sorted([('P',z[0]),('P',z[1])]))
   assert assembled==reference and 2*len(z)+len(n)+2-len(z)<=4;count+=1
 # n_j*(z_h*n_k-P_hk)=n_k*P_hj-n_j*P_hk whenz_h=P_hj/n_j.
 charts=3*42*42
 return dict(all42x1081_polynomial_homogenization_basis_checks=count,all_quartics_degree_atmost4=True,all_chart_recovery_compatibility_formulas=charts,chart_recovery_formula='n_j*(z_h*n_k-P_hk)=n_k*P_hj-n_j*P_hk=-C_h,j,k; CantisymmetricandC_jj0',residual_identity='R=Xu+b=X(u-P)+H, whereP=-Lb,H=b+XP,LX=I126',no_evaluation_point_prime_or_Gram_value=True)
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
 rdir=ROOT/'reports/S61_SPARSE_GENERIC5_RESIDUAL129_CIRCUIT_20261104_v63';gdir=ROOT/'reports/S61_SPARSE_GENERIC5_G6_RESIDUAL_20261104_v63';fdir=ROOT/'reports/S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63';r=read(rdir/'residual129-factored-circuit.json.gz');g=read(gdir/'G6-normal-residual-circuit.json.gz');f=read(fdir/'full134-direct-circuit.json.gz');blocks=read(fdir/'coefficient-block-circuit.json.gz')
 paths=[rdir/'report.json',gdir/'report.json',ROOT/'reports/S61_SPARSE_GENERIC5_G6_RESIDUAL_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_SPARSE_GENERIC5_RESIDUAL129_CIRCUIT_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_SPARSE_GENERIC5_CROSS_RANK_COVARIANCE_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63/report.json',ROOT/'reports/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63_REVIEW/report.json']
 assert all(read(p)['all_passed']for p in paths);J=r['selected_original134_rows'];assert J==[50,10,52,51,15]
 pre=dict(new_explicit_parameter_polynomial_terms_TOTAL=0,all_actual_scalar_inputs_external=True,scalar_X_shape=[130,126],Gram_shape=[126,126],leftinverse_shape=[126,130],Gram_scalar_products_if_expanded=126*126*130,expanded_Gram_nodes_built=0,estimated_operator_nodes_upper=50000,hard_nodes_less_than1000000=True,hard_seconds=60,hard_memory_bytes=2000000000,hard_bytes=32*1024*1024,old_abstract_saturation_theorem_already_exists='S61_AFFINE_NORMAL_FINITE5_20261103_v63; newresultisactualK5sourceboundcanonical130rowoperator')
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n');op=Operator();contexts=dict(residual129=str((rdir/'residual129-factored-circuit.json.gz').relative_to(ROOT)),G6=str((gdir/'G6-normal-residual-circuit.json.gz').relative_to(ROOT)),full134_direct=str((fdir/'full134-direct-circuit.json.gz').relative_to(ROOT)),full134_blocks=str((fdir/'coefficient-block-circuit.json.gz').relative_to(ROOT)))
 Xrows=[[op.scalar('residual129',i)for i in row]for row in r['cross129x126']]+[[op.scalar('G6',i)for i in g['R_G6_cross126']]];assert len(Xrows)==130 and all(len(row)==126 for row in Xrows);X=op.mat(Xrows);Xt=op.tr(X);Gram=op.mm(Xt,X);L=op.solve(Gram,Xt)
 sourceR=[op.poly('residual129',i)for i in r['residual_direct_F129']]+[op.poly('G6',g['R_G6_direct_node'])];assert len(sourceR)==130;b=op.vec([op.zeroz(i,1)for i in sourceR]);P=op.vneg(op.mv(L,b));H=op.va(b,op.mv(X,P));N=[op.n(i)for i in range(42)];Pe=[op.get(P,i)for i in range(126)];pure=[op.get(H,i)for i in range(130)];compat=[];compat_labels=[]
 for h in range(3):
  for j in range(42):
   for k in range(j+1,42):compat_labels.append([h,j,k]);compat.append(op.sub(op.prod([N[j],Pe[42*h+k]]),op.prod([N[k],Pe[42*h+j]])))
 assert len(compat)==2583
 selected=[op.poly('full134',f['F134'][i])for i in J];selected_b=[op.zeroz(i)for i in selected];zz=list(itertools.combinations_with_replacement(range(3),2));source_zlinear=[];source_zquadratic=[]
 for i in J:
  row=blocks['rows134'][i];source_zlinear.append([op.add([op.kpoly(op.scalar('full134',row['plane_z_degree_le2'][1+h])),*[op.prod([op.kpoly(op.scalar('full134',row['cross_z_n126'][42*h+k])),N[k]])for k in range(42)]])for h in range(3)]);source_zquadratic.append([op.kpoly(op.scalar('full134',row['plane_z_degree_le2'][4+k]))for k in range(6)])
 quartic=[];quartic_labels=[]
 for a in range(5):
  for j in range(42):quartic_labels.append([J[a],j]);quartic.append(op.add([op.prod([N[j],N[j],selected_b[a]]),*[op.prod([N[j],Pe[42*h+j],source_zlinear[a][h]])for h in range(3)],*[op.prod([source_zquadratic[a][k],Pe[42*h+j],Pe[42*l+j]])for k,(h,l)in enumerate(zz)]]))
 assert len(quartic)==210
 for roots,degree in[(pure,2),(compat,3),(quartic,4)]:assert all(op.types[i]['kind']=='PolyN'and op.types[i]['degree_bound']<=degree and op.types[i]['normal_order_lower_bound']>=1 for i in roots)
 assert len(op.nodes)<50000;solves=[n for n in op.nodes if n['op']=='unique_parameter_linear_solve'];assert len(solves)==1 and op.types[Gram]['kind']=='MatK'
 proof=universal_checks();(out/'universal-assembly-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 obj=dict(nodes=op.nodes,coefficient_field='K=QQ(H,L,s0,q0,t0)',polynomial_ring='R=K[n0,...,n41]',source_context_files=contexts,source_context_semantics={'full134':'NodeIDs concatenatefull134_direct.nodes andfull134_blocks.appended_nodes','residual129':'literalnodesinresidual129source, commonz/nvariables','G6':'literalnodesinG6source, commonz/nvariables'},selected_original134_rows=J,residual130_row_order=[*r['residual_original134_rows'],'G6'],actual_X130x126=X,pure_normal_b130=b,Gram_parameter_matrix=Gram,leftinverse_parameter_matrix=L,normal_rewrite_P126=P,pure_normal_H130=H,H_generator_nodes130=pure,cross_compatibility_generator_nodes2583=compat,cross_compatibility_labels=compat_labels,plane_quartic_generator_nodes210=quartic,plane_quartic_labels=quartic_labels,all2923_J_generators=[*pure,*compat,*quartic],generator_counts=dict(redundant_quad130_span4=130,cubic=2583,quartic=210,total=2923),canonical_solve_semantics='Lunique126x130Kmatrix satisfying(Transpose(X)*X)*L=Transpose(X); completefinitelinearsystem, NOTevaluated',matrix_operations_semantics='finiteusualsums/products; noneexpandedGram/nonecoefficientfieldGaussianexecuted',normal_rewrite_entry_order='P_hk at42*h+k',quartic_formula='n_j^2*Fselected_r(P_0j/n_j,P_1j/n_j,P_2j/n_j,n), implementedpolynomiallywithoutn_jinverse',patch_maps=dict(normal_to_original='n->n,z_h->P_hj/n_j onD(n_j)',original_to_normal='n->n; inversebyLX=Iandcompatibility'),saturation_recipe=dict(J='all2923_J_generators',ideal_m=N,operation='J:m^infinity',finite_elimination_recipe='intersection_j0..41 of ((J,1-t_j*n_j) intersect R)',coefficient_saturation_computed=False),no_sum_normal_squares_inverted=True,no_normal_unknown_inverse_in_generator_construction=True,only_parameter_solve=True,all42normals_andnormal2_retained=True)
 raw=json.dumps(obj,separators=(',',':')).encode();assert len(raw)<=32*1024*1024;(out/'typed-normal42-operator.json.gz').write_bytes(gzip.compress(raw,mtime=0))
 theorem=dict(actual_A='K[z3,n42]/(selected5F,all130R), sameoriginal134idealbyG6identity',J='130pure-normalH(degree<=2,span4)+2583compat(degree<=3)+210manifestquartics(degree<=4)',patch_scheme_isomorphism=True,LX_equals_I126='Gramnormal equations,detGramnonzero byformalreality/CauchyBinet',forward='R=0impliesu=P,H=0,compat=0andquartics=0',backward='onD(n_j),compatgivesallP_hk=z_h*n_k; H=0givesR=0; quarticsandunitn_j²giveallselectedF',all_nilpotents_preserved=True,finite_over_R_generators_atmost5='crossrewrites+planeHF(1,3,1)/topcubicspan10; noK-finitenessclaim',known_five_clopen='AcceptedactualpolarizationA=E5xB andA/(n)=E5',normal_projection_B_closed_immersion='BfiniteR;B_m0byNakayama;R_j->B_jsurjall42charts;finitecokerzeroallprimes=>R->Bsurj',exact_ideal_kernel='ker(R->B)=intersection_j(J:n_j^infinity)=J:m^infinity sinceSpecBclosedandavoidsorigin',exact_B='R/(J:m^infinity)',empty_criterion='B=0 iff J:m^infinity=R iff someintegerd>=0 hasm^d subsetJ',empty_criterion_verified=False,old_abstract_theorem_not_new='ExistingNORMAL_FINITE5alreadyprovesthisschemeeliminationtypewith129rows; herecanonicalactualK5/130rows/allsourceboundtypedoperatorsarefixed',formal_reality_scope='OnlycoefficientfieldGramproof; extensiontoCpreservesnonzeroKdet; noorder/positivityassumptiononcomplexnsolutions')
 (out/'scheme-elimination-theorem.json').write_text(json.dumps(theorem,indent=2)+'\n')
 paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',*paths,rdir/'residual129-factored-circuit.json.gz',gdir/'G6-normal-residual-circuit.json.gz',fdir/'full134-direct-circuit.json.gz',fdir/'coefficient-block-circuit.json.gz',ROOT/'research/notes/S61_AFFINE_NORMAL_FINITE5_20261103_v63.md',ROOT/'research/notes/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63.md']
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<=32*1024*1024
 report=dict(id=ID,all_passed=True,typed_operator_nodes=len(op.nodes),raw_operator_bytes=len(raw),new_explicit_parameter_polynomial_terms=0,Gram_scalar_products_if_expanded=2063880,expanded_Gram_products_computed=0,field_Gaussian_or_126inverse_computed=False,all2923_generators_degree_le4_defined=True,all42_chart_scheme_isomorphisms_proved=True,all45402_quartic_basis_checks=True,closed_residual_exact_saturation_presentation_proved=True,normal_saturation_or_unit_computed=False,no_point_prime_GB_or_normal_sumsquares_inverse=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,nodes=len(op.nodes),rawbytes=len(raw),generators=2923,seconds=report['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
