#!/usr/bin/env python3
"""Ten pureR0 ordinary coefficient functionals: literalIAIB and660NF source."""
from pathlib import Path
from fractions import Fraction as Q
import itertools,json,math,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_NATIVE_NONZZ10_DESCENT_SOURCE_20261105_v63'
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):
  (o/n).write_text(json.dumps(v,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in o.iterdir())<4000000
 src='reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json';S=fld(src,['actual_physical_rows9','source_F9_ids','YB_in_YA_linear','A0','actual_generator_normalization','source_conductor_generators_rows_and_ids','gA_plane3','gB_own_plane3','ell_exact'])
 assert S['actual_physical_rows9']==[27,108,109,112,113,117,118,119,61] and S['ell_exact']=='5*(s_CA-r_CA)'
 ca=[7148,7199,7454,7505,7688,7943,8045];cb=[27897,27948,28203,28254,28437,28692,28794]
 assert S['source_F9_ids']['CA']==ca+[17493,19106] and S['source_F9_ids']['CB']==cb+[38273,39886]
 assert [j for row,j in S['source_conductor_generators_rows_and_ids']['CA']]==ca+[7790] and [j for row,j in S['source_conductor_generators_rows_and_ids']['CB']]==cb+[28539]
 full='reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json';T=fld(full,['source2_monomials','both57_source2']);mon2=[];u=[]
 for slots in T['source2_monomials']:
  m=[0]*9
  for i in slots:m[i]+=1
  mon2.append(tuple(m))
 for i in range(9):m=[0]*9;m[i]=1;u.append(tuple(m))
 zero=(0,)*9;indices=[mon2.index(zero)]+[mon2.index(m) for m in u];rawf=[]
 for side,f63,coords in [('CA',7790,ca+[17493,19106]),('CB',28539,cb+[38273,39886])]:
  family=T['both57_source2'][side];assert family['F9_nodes']==coords and family['f63']==f63 and family['all57_equations_all55_coefficients_verified']
  v=family['pure_function2jets'][str(f63)];assert all(Q(v[k])==0 for k in indices);rawf.append({'family':side,'raw_f63_node':f63,'own_F9_constant_and_first_indices':indices,'all_constant_firsts_exact0':True})
 # Literal B coordinate linearization in the commonACTUAL CA F9 chart.
 B1=[[Q(x) for x in row] for row in S['YB_in_YA_linear']];assert len(B1)==9 and all(len(row)==9 for row in B1)
 assert [row[:7] for row in S['YB_in_YA_linear'][:7]]==S['A0'] and all(B1[i][7]==B1[i][8]==0 for i in range(7))
 # No coordinate inverse or jet solve: centeredCBcompositionpreserveszero0/1ofrawf63.
 IA=[{u[i]:Q(1)} for i in range(7)]+[{}];IB=[{u[j]:B1[i][j] for j in range(7) if B1[i][j]} for i in range(7)]+[{}]
 c0='reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json';old=fld(c0,['actual_C0_basis8_scalar_nodes']);assert old['actual_C0_basis8_scalar_nodes']==ca+[505938]
 # Keep thecurrentpointcomparisonframe literal, withoutusingitinthelambda arithmetic.
 mg='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';M=fld(mg,['Mnew0_columns'])['Mnew0_columns'];fr='reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json';assert fld(fr,['actual_Mnew0_columns'])['actual_Mnew0_columns']==M
 pnf='reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json';D=fld(pnf,['monomials660_degree_le4','all660_actual_raw33_normal_forms','actual33_basis_monomial_indices']);mons=list(map(tuple,D['monomials660_degree_le4']));NF=list(map(sp,D['all660_actual_raw33_normal_forms']));basisidx=D['actual33_basis_monomial_indices'];assert len(mons)==len(NF)==660 and len(basisidx)==33
 raw='reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json';rawb=fld(raw,['basis74','free_column_indices']);basis=[]
 for k in rawb['free_column_indices']:
  sector,r,e=rawb['basis74'][k];m=[0]*9;m[7]=r;m[8]=e
  if sector:m[sector-1]=1
  basis.append(tuple(m))
 assert basis==[mons[j] for j in basisidx] and len(set(basis))==33
 premises=[('reports/AMBIENT_B18_ELL2_RAW45_TRUNCATION_PREFLIGHT_20261105_v63_REVIEW/report.json',None),('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json',None),('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63_REVIEW/report.json','d07787a307570c50a903ad4532dd94d58e0c2a73c72e62ebcc7c78b213cc58e9'),('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/report.json','6428e2178c41cfd7cfcf271402c4d13898879414a0b5129fd409e3c97300d460')]
 rawpeer=None
 for p,h in premises:
  q=rd(p);assert q['all_passed']
  if h:assert pins[p]==h
  if 'RAW33_QUOTIENT_' in p:rawpeer=q;assert q['actual_product_m5_zero']
 expected_mons=[]
 for degree in range(5):
  for slots in itertools.combinations_with_replacement(range(9),degree):
   m=[0]*9
   for i in slots:m[i]+=1
   if m[8]<2:expected_mons.append(tuple(m))
 assert mons==expected_mons
 pairs=[(0,7),(0,8),(1,7),(1,8),(4,7),(4,8),(5,7),(5,8),(7,7),(7,8)];targets=[]
 for i,j in pairs:m=[0]*9;m[i]+=1;m[j]+=1;targets.append(tuple(m))
 lambdas=[{i:Q(1) for i,m in enumerate(basis) if m==t} for t in targets];assert all(len(v)==1 for v in lambdas)
 # SourcecountsandpositiveenvelopesFIRST;noκ/G8component readbackisloaded.
 prior=1
 for lam in lambdas:
  for nf in NF:
   terms=[(x,nf[i]) for i,x in lam.items() if i in nf]
   if terms:
    den=math.lcm(*(x.denominator*y.denominator for x,y in terms));U=sum(abs(x.numerator*y.numerator)*(den//(x.denominator*y.denominator)) for x,y in terms);prior=max(prior,den.bit_length(),U.bit_length())
 assert prior<=384
 sourcebase=count([S,T,mon2,indices,rawf,B1,IA,IB,old,M,D,mons,NF,basisidx,rawb,basis,pairs,targets,lambdas,pins]);owner=sourcebase+30000;assert owner<150000
 rawinputs=sum((ROOT/p).stat().st_size for p in pins);assert rawinputs+4000000<33554432
 put('phase-positive-before-source-checks.json',{'phase':'native10functional_SOURCE_only_before640initial_and6600NFchecks','positive_lambda_NF_prior_bits':prior,'positive_literal_product_bound':'eachCAlinearformisunitcoordinateor0,so64initialproductsarecopiesofliteralBlinearcoefficientsornothing;nosum/inverse','source_owned_prior':owner,'whole_RAW_prior':rawinputs+4000000,'no_new_unitIFT_functionresponse':True})
 products=[]
 for i,a in enumerate(IA):
  for j,b in enumerate(IB):
   init={}
   for m,x in a.items():
    for n,y in b.items():
     k=tuple(m[t]+n[t] for t in range(9));assert k not in init;init[k]=x*y
   assert all(sum(m)==2 and sum(m[:7])==2 and m[7]==m[8]==0 for m in init)
   assert all(not init.get(t,Q(0)) for t in targets)
   products.append({'literal_CA_generator_index':i,'literal_CB_generator_index':j,'source_nodes':[S['source_conductor_generators_rows_and_ids']['CA'][i][1],S['source_conductor_generators_rows_and_ids']['CB'][j][1]],'homogeneous_degree2_initial':[[list(m),str(x)] for m,x in sorted(init.items())],'all10_nonzz_coefficients_exact0':True})
 values=[];bits=0
 for m,nf in zip(mons,NF):
  out=[]
  for t0,(target,lam) in enumerate(zip(targets,lambdas)):
   q=Q(0)
   for i,x in lam.items():q+=x*nf.get(i,Q(0));bits=max(bits,abs(q.numerator).bit_length(),q.denominator.bit_length());assert bits<=prior
   assert q==(1 if m==target else 0)
   if q:out.append([t0,str(q)])
  values.append(out)
 assert all(NF[j]=={i:Q(1)} for i,j in enumerate(basisidx))
 data={'physical_R0_variable_order':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell'],'ell_exact':'5*(sA-rA)','same_actual_physical_rows9':S['actual_physical_rows9'],'same_CA_CB_F9_nodes':S['source_F9_ids'],'original_raw_IA8_IB8_generator_ids':S['source_conductor_generators_rows_and_ids'],'raw_f63_zero_constant_first_source':rawf,'actual_common_zB_linear_forms7':[[[list(m),str(x)] for m,x in sorted(v.items())] for v in IB[:7]],'all64_literal_product_degree2_initials':products,'pair_order10':pairs,'ordinary_target_monomials10':targets,'raw33_basis_monomials':basis,'lambda10_by33_sparse_rows':[[[i,str(x)] for i,x in row.items()] for row in lambdas],'lambda10_by33_dense_rows':[['1' if i in row else '0' for i in range(33)] for row in lambdas],'all660_monomial_inputs':mons,'all660_lambda_values10_sparse':values,'all6600_source_NF_identities_including_nonhomogeneouscarry':True,'ideal_annihilation_proof':'ALL16rawgeneratorshavezeroorigin. CAlinearformsareziand0; CBlinearformsareA0ziand0 (last2commonplaneentries0);rawf63own0/1zeroisstableundercenteredcommonpullback. ThereforeALL64productquadraticinitialsarezzonly. Anyformalmultipliercontributestodegree2onlythroughitsconstant;nonzzdegree2functionalsannihilateI. Theyalsoannihilateell²andm³. ThisusesALLliteralproductgenerators,notjustanormalizedheuristic.','exact_descent_proof':'ActualR/(I,ell²)hasm5=0andraw33basisbyacceptedsourcecertificate. lambdaon33literalbasismonomialsannihilatesactualideal;ALL660savedNFcolumnsverifythatlambdaNF(m)=ordinarycoefficient(m). NonhomogeneoushighercarryisretainedbythefullNFcomparison.','OLD_Q9_and_current_M0_bindings':'OLD_C0basis8=actualCAz7+[gA505938],notrawf637790;physicalF9rowsandellchartareunchanged. CurrentM0isexactlytheacceptedactivepointtransport;noM1/comparisongaugeisassignedzero.','future_unitIFT_scope':'lambda maylaterbeappliedcomponentwisetoSAMEphysicalR0canonicalunitIFTκ/G8imagesinthefree4basisafteraseparatebranch/table/overlap/resourcegate. Noκ/G8/ν/functionresponseisreadhere.','honest_frontier':'Theoremdoesnotsupplyfiniteimagesofgenuinehonestg/hquotientsanddoesnotauthorizeb/gdivisionintheArtinring. TheyremainopaqueuntilwholeADcancellationoraseparatesourcecertificate.','plane3_4_scope_separate':'Noextensionofthetensecondfunctionaltohighermixedcoefficients. Pureplane3/4requiresitsseparateP9quotient/sourcegate.','positive_prior_bits':prior,'actual_check_bits':bits,'no_native_OLD_honest_or_upper_functionresponse_values':True}
 ep='actual-ten-nonzz-functionals-and-all660-source-NF-descent.json';put(ep,data)
 report={'id':STEM,'all_passed':True,'SOURCE_descent_only_no_functionresponses':True,'literal_product_generators_checked64':64,'lambda_shape':[10,33],'all660_NF_columns_all10functionals_checked':True,'nonhomogeneous_carry_retained':True,'positive_prior_bits':prior,'actual_check_bits':bits,'source_owned_prior':owner,'whole_RAW_prior':rawinputs+4000000,'no_kappa_G8_honest_newcoeff_or_upper_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'seconds':time.monotonic()-start};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owner,'positive':prior,'lambda_indices':[next(iter(v)) for v in lambdas]}))
if __name__=='__main__':main()
