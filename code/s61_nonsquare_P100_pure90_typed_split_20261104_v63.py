#!/usr/bin/env python3
"""Typed global8-target split from external five-minor witness; coefficient terms0."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
class Typed:
 def __init__(self):self.nodes=[]
 def add(self,op,shape,**data):
  n={'id':len(self.nodes),'op':op,'shape':list(shape),**data};self.nodes.append(n);return n['id']
 def shape(self,i):return self.nodes[i]['shape']
 def eye(self,n):return self.add('identity_matrix',(n,n))
 def mul(self,a,b):
  aa,bb=self.shape(a),self.shape(b);assert aa[1]==bb[0]
  return self.add('matrix_product',(aa[0],bb[1]),args=[a,b])
 def scale(self,s,a):assert self.shape(s)==[1,1];return self.add('central_scalar_multiple',self.shape(a),args=[s,a])
 def summ(self,args):
  assert args and all(self.shape(i)==self.shape(args[0])for i in args)
  return self.add('matrix_sum',self.shape(args[0]),args=args)
 def neg(self,a):return self.add('negate',self.shape(a),args=[a])
 def sub(self,a,b):assert self.shape(a)==self.shape(b);return self.summ([a,self.neg(b)])
 def sel(self,a,cols):return self.add('select_columns',(self.shape(a)[0],len(cols)),args=[a],columns=cols)
 def adj(self,a):
  n,m=self.shape(a);assert n==m==8
  return self.add('polynomial_adjugate_8x8',(8,8),args=[a],definition='Signed7x7cofactor determinants; no scalar expansion or matrix inverse.')
 def invunit(self,a,proof):assert self.shape(a)==[1,1];return self.add('certified_unit_inverse',(1,1),args=[a],unit_proof=proof)
 def rows(self,a,rows):return self.add('select_rows',(len(rows),self.shape(a)[1]),args=[a],rows=rows)
 def transpose(self,a):s=self.shape(a);return self.add('transpose',(s[1],s[0]),args=[a])
 def hcat(self,args):assert len({self.shape(i)[0]for i in args})==1;return self.add('horizontal_block', (self.shape(args[0])[0],sum(self.shape(i)[1]for i in args)),args=args)
 def vcat(self,args):assert len({self.shape(i)[1]for i in args})==1;return self.add('vertical_block',(sum(self.shape(i)[0]for i in args),self.shape(args[0])[1]),args=args)
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 cb=base/'S61_NONSQUARE_P100_PURE15_COVER_20261104_v63';cr=read(cb/'report.json');cov=read(cb/'native.json');cbind=read(cb/'source-binding.json');assert cr['complete']and cr['whole_P100_quartic_pointwise_lower149']
 assert sha(cb/'native.json')==cr['evidence_sha256']['native.json']and cov['literal_sum_Ai_actualMinor_i_equals_E_times_gcd_verified']
 peer=base/'S61_NONSQUARE_P100_PURE15_COVER_20261104_v63_REVIEW/report.json';assert sha(peer)=='0e830d25d8631a7e38da157d15e7be181dac4dc80bc6c282381cc9e64f901f1a'
 rb=base/'S61_NONSQUARE_P100_PURE8_GENERIC_RANK_20261104_v63';rd=read(rb/'native.json');rlabels=read(rb/'source-labels.json')
 qb=base/'S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63';gb=base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63';kb=base/'S61_NONSQUARE_P100_K2_PURE7_20261104_v63'
 kcb=base/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63';kp=kcb/'pure-cubic84-2.json'
 rowpaths=[qb/'pure98.json']+[gb/f'pure98-row{i}.json'for i in range(7)]+[kb/f'pure98-row{i}.json'for i in range(7)]
 qp=read(rowpaths[0]);labels98=qp['target122_columns'];assert len(labels98)==98
 for p in rowpaths[1:]:d=read(p);assert d['target122_columns']==labels98 and len(d['values98'])==98
 hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';h=read(hp)
 fd=base/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json'
 lp=base/'S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63/native.json';lr=read(lp);assert lr['all23_source_identities_both_source_frames_verified']
 catp=base/'S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json';cat=read(catp);catalogue=[cat['source_rows'][i]for i in cat['retained_row_indices']];assert len(catalogue)==319
 assert catalogue[174:286]==h['all112_cubic_source319_definitions']
 q5H=[i for i,v in enumerate(catalogue)if v[0]=='quadric'and v[1]==5 and v[2]==[0,0,0,0,0,1,0,0]];assert len(q5H)==1
 C=rd['selected_columns'];T=[i for i in range(98)if i not in C];assert len(C)==8 and len(T)==90 and len(set(C+T))==98
 assert cbind['fixed_old_selected_target98_columns']==C and cbind['fixed_original122_target_labels']==[labels98[i]for i in C]
 basis=[v['basis_index']for v in cov['actual_minor_polynomial_Bezout_contributors']];assert basis==[0,7,9,10,15]
 Epow=cov['common_actual_minor_unit_powers5'];assert Epow==[21,93,6,0,1]
 def ref(path,pointer):return {'path':str(path.relative_to(ROOT)),'sha256':sha(path),'json_pointer':pointer}
 # No input polynomial or RF coefficient is copied into the typed operator.
 t=Typed();M=t.add('external_row_matrix',(15,98),row_refs=[ref(p,'/values98')for p in rowpaths],coordinate_context='P100/actualNF117/H2-and-H1-cleared/pure98')
 A=t.sel(M,C);B=t.sel(M,T);I8=t.eye(8);I15=t.eye(15);I90=t.eye(90)
 E=t.add('external_polynomial_scalar',(1,1),reference=ref(cb/'native.json','/common_actual_minor_unit_polynomial'),unit_monomial_powers5=Epow)
 Ei=t.invunit(E,{'ring':'T100','known_factor_names':['z','z+20','z+60','203z+2700','z+12'],'powers':Epow,'proof':ref(cb/'native.json','/literal_sum_Ai_actualMinor_i_equals_E_times_gcd_verified')})
 minorops=[];terms=[];sourceOld=rd['selected_rows'];scalesOld=cbind['old_selected_source_binding']['selected_rowLCD_powers']
 for n,b in enumerate(basis):
  rows=sourceOld[:];powers=[x[:]for x in scalesOld]
  if b==0:minorref=ref(cb/'native.json','/old_minor_strip/full');candidate=None
  else:
   cc=cov['candidates'][b-1];light=cc['K2_light_position'];replace=cc['replaced_old_selected_row_position'];rows[replace]=8+light;powers[replace]=cbind['new_sources_in_order'][light]['known_unit_rowLCD_powers5'];minorref=ref(cb/'native.json',f'/candidates/{b-1}/split/full');candidate={'K2_light_position':light,'replace_position':replace}
  assert len(set(rows))==8
  P=t.add('weighted_source_selector',(8,15),selected_rows=rows,row_unit_powers5=powers,unit_factor_names=['z','z+20','z+60','203z+2700','z+12'])
  S=t.mul(P,A);Adj=t.adj(S);D=t.add('external_polynomial_scalar',(1,1),reference=minorref,certified_equal_to_determinant_of=S,source_minor_proof=ref(cb/'source-binding.json','/minor_formula'))
  ai=t.add('external_polynomial_scalar',(1,1),reference=ref(cb/'native.json',f'/actual_minor_polynomial_Bezout_contributors/{n}/coefficient'),basis_index=b)
  terms.append(t.scale(ai,t.mul(Adj,P)));minorops.append({'basis_index':b,'candidate':candidate,'weighted_P':P,'S':S,'Adj':Adj,'D':D,'Ai':ai,'source_rows_in_order':rows,'LCDpowers_in_order':powers})
 K=t.scale(Ei,t.summ(terms));KA=t.mul(K,A);KB=t.mul(K,B);Psource=t.sub(I15,t.mul(A,K));R=t.mul(Psource,B)
 Z90x8=t.add('zero_matrix',(90,8));Z8x90=t.add('zero_matrix',(8,90))
 U=t.add('block_matrix',(98,98),row_partition=[8,90],column_partition=[8,90],blocks=[[I8,t.neg(KB)],[Z90x8,I90]],acts_after_column_order=C+T)
 Ui=t.add('block_matrix',(98,98),row_partition=[8,90],column_partition=[8,90],blocks=[[I8,KB],[Z90x8,I90]],acts_after_column_order=C+T)
 Mperm=t.hcat([A,B]);MU=t.mul(Mperm,U);KMU=t.mul(K,MU)
 V=t.add('arbitrary_regular_row',(1,98),context='Anyactualpure98source; neverassume its remaining90tail iszero.')
 V90=t.sub(t.sel(V,T),t.mul(t.sel(V,C),KB))
 # Bind the FULL319 source context; all otherwise unresolved pure tails are lazy.
 lightframe=gb/'joint-cubic-frame.json';assert sha(lightframe)==read(gb/'report.json')['evidence_sha256']['joint-cubic-frame.json']
 PhysicalFull=t.add('external_lazy_full_weight4_source_matrix',(319,239),catalogue_ref=ref(catp,'/retained_row_indices'),source_catalogue_rows=catalogue,frame='AcceptedphysicalTn normalized quadrics/cubics and Qhat truequartics; knownunit row scalings only.',coefficient_queries_performed=0)
 Full=t.add('accepted_weighted4_coordinate_transport',(319,239),args=[PhysicalFull],old_physical_monomials_ref=ref(catp,'/weight4_columns'),new_joint_monomials_ref=ref(hp,'/parent_joint_target239_monomials'),physical_light_inverse_ref=ref(lightframe,'/physical_light_inverse7x7'),rule='Pure210 usesSym4(B); H*Sym2light28 usesSym2(B); physicalH=v5 equalsjointH(var7),H2unchanged. No expansion.')
 NF=t.add('apply_external_linear_operator',(319,122),args=[Full],reference=ref(fd,'/full_target_projection239_to122'),strict_Hdegree=True)
 F0=t.sel(NF,h['target122_Hdegree_classes']['0']);F1=t.sel(NF,h['target122_Hdegree_classes']['1']);F2=t.sel(NF,h['target122_Hdegree_classes']['2'])
 assert h['target122_Hdegree_classes']['0']==labels98
 U2=t.add('external_certified_unit_scalar',(1,1),reference=ref(hp,'/H2_known_unit_root'),formula='(50*a*B*C)^2 restrictedtoP100',known_P100_unit_proof='a,B,C areproducts/ratios ofthefive knowncurveunits.')
 U2i=t.invunit(U2,{'proof':ref(hp,'/source_pair_determinant_equals_knownunit_exact_QQ_checked'),'formula':'(50aBC)^2; onlyknownunits'})
 chi0=t.add('lazy_polynomial_then_NF_pure_part',(1,98),polynomial_ref=ref(hp,'/chi_full_joint_weight4_polynomial_root'),polynomial_nodes_ref=ref(hp,'/chi_polynomial_nodes'),NF_ref=ref(fd,'/full_target_projection239_to122'),target122_labels=labels98)
 chi1=t.add('external_lazy_scalar_root_row',(1,23),reference=ref(hp,'/chi_H1_post117_scalar_roots'))
 PiOrig=t.sel(t.rows(NF,list(range(174,286))),labels98)
 delta=t.add('certified_regular_unit_formula',(1,1),formula='-64000*C^2 restrictedtoP100',known_unit_factor_context='fiveP100linearunits')
 Ddelta=t.add('delta7_source_preconditioner',(112,112),delta=delta,rule='For eachphysical light j: row6*7+j becomesdelta*(row6*7+j-row1*7+j); otherrowsunchanged.')
 PiActual=t.mul(Ddelta,PiOrig);L=t.add('external_regular_RF_matrix',(112,23),reference=ref(lp,'/right_inverse112x23'),source_frame='actualdelta7; allunitdenominatorsverifiedbyacceptedH1peer')
 Section=t.mul(t.transpose(L),PiActual)
 f2u=t.scale(U2i,F2);H2pure=t.sub(F0,t.mul(f2u,chi0));H2h1=t.sub(F1,t.mul(f2u,chi1));Gamma=t.sub(H2pure,t.mul(H2h1,Section))
 Gamma90=t.sub(t.sel(Gamma,T),t.mul(t.sel(Gamma,C),KB));Presentation90=t.vcat([Gamma90,R])
 # These are exact source embeddings; polynomial/RF coefficients remain external.
 embeddings=[{'kind':'q5H_after_H1_clear','original319_source_index':q5H[0],'producer_cache_row0_and_QPAIR57_proof':ref(qb/'source-binding.json','/literal_q5H23_cache_matches_producer')}]
 embeddings +=[{'kind':'literalG18_product','original319_source_index':174+j,'physical_light':v}for j,v in enumerate([0,1,2,3,4,6,7])]
 embeddings +=[{'kind':'correctedK2_product','physical_light':v,'original319_source_indices':[174+7*s+j for s in range(16)],'original16_cofactor_ref':ref(kp,'/old16_source_cofactors'),'actual_delta7_cofactor_ref':ref(kp,'/actual_delta7_source_cofactors')}for j,v in enumerate([0,1,2,3,4,6,7])]
 assert len(embeddings)==15
 J=t.add('certified_sparse_source_embedding',(15,319),row_recipes=embeddings,coefficient_expansion_performed=False)
 JGamma=t.mul(J,Gamma)
 # Shape/provenance assertions, plus universal algebra rewrites; no scalar arithmetic.
 assert t.shape(K)==[8,15]and t.shape(KA)==[8,8]and t.shape(R)==[15,90]and t.shape(V90)==[1,90]
 assert t.shape(Presentation90)==[334,90]and t.shape(Section)==[23,98]and t.shape(Gamma)==[319,98]
 for v in minorops:
  assert t.nodes[v['S']]['args']==[v['weighted_P'],A]and t.nodes[v['Adj']]['args']==[v['S']]
  assert t.nodes[v['D']]['certified_equal_to_determinant_of']==v['S']
 inverses=[n for n in t.nodes if n['op']=='certified_unit_inverse'];assert {n['args'][0]for n in inverses}=={E,U2}
 proof={'KA_equals_I8':{'steps':['Expand K*A and associate each Adj(Si)*Pi*A.','Pi*A=Si bydefinition, includingallrowLCDunits.','Adj(Si)*Si=det(Si)*I8 byuniversalcofactoridentity.','det(Si)=Di bythefixed source/minorcertificate.','sum Ai*Di=E bytheexternal literalpolynomialidentity.','E^-1*E=1 usingonlythefiveknowncurveunits.'],'external_identity':ref(cb/'native.json','/literal_sum_Ai_actualMinor_i_equals_E_times_gcd_verified'),'left_node':KA,'right_node':I8},
  'target_GL98':{'node':U,'inverse':Ui,'determinant':1,'proof':'Blockunitriangular; the two offdiagonal blocks cancel inbothproducts.'},
  'unit_rows':{'left_node':KMU,'right_expression':'[I8,0_8x90]','proof':'K[A,B]U=[KA,KB-KA*KB]=[I8,0].'},
  'selected_residual':{'node':R,'expression':'(I15-A*K)*B','proof':'M U=[A,R]; K R=0 and(I-AK)^2=I-AK. Sourcequotienthasprojectiverank7; all15rowsretained,nobasischosen.'},
  'general_row':{'node':V90,'expression':'V_T-V_C*K*B','all_puretails_retained':True},
  'source_embedding':{'left_node':M,'right_node':JGamma,'proof':'q5actualproducerROW0=QPAIR57withoriginq5H, literalG18products, andtheacceptedold16K2cofactors commutewiththelinearNF/H2/H1targetprojection. All15cachedsourceidentitieswereindependentlyverified.'},
  'full_source_same_coker':{'node':Presentation90,'proof':'KnownNF117/H2/H1sectionsgivecoker(original319x239)=coker(Gamma). M=J*Gamma bythe15exactsourceembeddings. AppendingMchangesnoimage. ApplytargetGL, useK M U=[I8,0], anddeleteonlytheunit8targetdirections. RemainingallGamma90plusall15R presentsamecoker.','source_count334_is_redundant_no_freebasis_claim':True},
  'no_new_bound':True,'full239_criterion':'Originalcokerzero overT100 iff the displayed334x90 residualpresentation haszero cokernel. No residualrowor rankwas evaluated.'}
 context={'ring':'T100=QQ[z,1/(z(z+20)(z+60)(203z+2700)(z+12))]','pure98_target122_labels':labels98,'selected_local98_columns8':C,'selected_original122_labels8':[labels98[i]for i in C],'remaining_local98_columns90':T,'remaining_original122_labels90':[labels98[i]for i in T],
          'H1_source_frame':'UsePiActual=delta7*PiOrig withLactual; doNOTmixLactualandoriginalpuretails.','general_full_source_formula':'Gamma=f0-(f2/U2)*chi0-(f1-(f2/U2)*chi1)*(Lactual^T*PiActual). Allchi0/PiActual/rawf0 staylazyandareNOTzero.','source_embeddings15':embeddings}
 record={'coefficient_ring':context['ring'],'explicit_new_polynomial_or_RF_coefficient_terms':0,'external_Bezout_atoms_remain_external':True,'nodes':t.nodes,'exports':{'M15x98':M,'A15x8':A,'B15x90':B,'K8x15':K,'KA8x8':KA,'KB8x90':KB,'target_GL98':U,'target_GL98_inverse':Ui,'residual15x90':R,'general_V90':V90,'source_embedding15x319':J,'full_Gamma319x98':Gamma,'full_Gamma319x90':Gamma90,'full_redundant334x90':Presentation90},'five_minors':minorops,'context':context,'proof':proof}
 (out/'operator.json').write_text(json.dumps(record,indent=2)+'\n')
 paths=[Path(__file__),cb/'report.json',cb/'native.json',cb/'source-binding.json',peer,rb/'report.json',rb/'native.json',rb/'source-labels.json',qb/'source-binding.json',gb/'report.json',lightframe,kb/'report.json',kp,hp,fd,lp,catp]+rowpaths
 report={'id':ID,'all_shape_and_source_bindings_passed':True,'new_explicit_coefficient_terms':0,'typed_nodes':len(t.nodes),'global_KA_identity_from_accepted_cover':True,'target_GL98_det1':True,'actual_pure90_samecoker_constructor':True,'retained_residual15_no_source_freebasis_assumed':True,'all_full319_puretails_retained_lazy':True,'no_polynomial_NF_rank_query_or_residual90_row_evaluation':True,'uniform149_cover_peer_pinned':sha(peer),'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'operator_sha256':sha(out/'operator.json'),'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert len(t.nodes)<200000 and sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
