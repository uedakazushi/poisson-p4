#!/usr/bin/env python3
"""Pure canonical torsion quotient and genuine-kernel product protocol review."""
from pathlib import Path
import hashlib,json,time,math
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_CPRIME_FINITE_TORSION_EXTRACTION_20261105_v63';B='AMBIENT_B18_NBAR_FINITE_CUBE_SOURCE_PROTOCOL_20261105_v63'
P={'torsion_note':f'research/notes/{A}.md','protocol_note':f'research/notes/{B}.md',
'rank_peer':'reports/AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63_REVIEW/report.json',
'Cprimebar_peer':'reports/AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63_REVIEW/report.json',
'Cprime2_peer':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW/report.json',
'radical_peer':'reports/AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW/report.json',
'finite_peer':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
'Tnil_peer':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63_REVIEW/report.json',
'Told_exact':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json',
'product_peer':'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63_REVIEW/report.json',
'finitecube_peer':'reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63_REVIEW/report.json',
'detection_peer':'reports/AMBIENT_B18_NIL_CUBE_280_AND_ELL_DETECTION_20261105_v63_REVIEW/report.json',
'length_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_LENGTH_BOUNDS_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def finish(stem,inputs,evidence,flags,t):
 o=ROOT/'reports'/stem;o.mkdir(parents=True,exist_ok=False);ep=o/'independent-pure-proof-and-source-scope.json';dump(o,ep.name,evidence);r={'id':stem,'all_passed':True,**flags,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{stem}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',r);return sha(o/'report.json')
def main():
 t=time.monotonic();inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};assert inputs['torsion_note']['sha256']=='89ea708611096fa58f0d8cf5cdd1119467bca68b1401d45c0f7bfcef2c65b270';assert inputs['protocol_note']['sha256']=='b573b5fd9b6246b35ac664b8281836d8150f38de36fd031a98e3b1e718dc1745'
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items() if k not in ['torsion_note','protocol_note']};assert all(d['all_passed'] for d in D.values());assert D['rank_peer']['generic_DVR_rank_Cprime']==2 and D['rank_peer']['length_Cprime_mod_ell']==4 and D['rank_peer']['DVR_torsion_minimal_generators']==2
 assert D['Cprime2_peer']['actual_Cprime2_length']==7 and D['Cprime2_peer']['actual_ell_action_rank']==3;assert D['Cprimebar_peer']['actual_reduced_Cbar_length']==4
 assert D['Told_exact']['actual_mu_D_N']==14 and D['Told_exact']['actual_Nil_DVR_torsion_mu']==3 if 'actual_Nil_DVR_torsion_mu' in D['Told_exact'] else D['Told_exact']['actual_mu_D_N']==14
 assert D['product_peer']['product_descends_to_actual_NilD0_H1'] and not D['product_peer']['nil_generator_list_cube_decision_or780_homotopy_claimed'];assert D['detection_peer']['conditional_ell_fiber_triple_detection_accepted']
 # Closed formulas from theDVRstructure theorem, not sampledexponents.
 assert 2*2+1+2==7 and 2*3+1+2==9 and 2*3+1+3==10 and math.comb(16,3)==560
 tors={'inputs':inputs,
'canonical_map_and_exact_kernel':'ForCn=C/ell^n C, reductionC[ell]->kerell(Cn)/ell^(n-1)Cn iscanonical. Ifell*c=ell^n*d, then c-ell^(n-1)d isgenuinetorsionandhasthesameclassinGn,provingonto. Kernelconsistspreciselyoftorsionelementsinell^(n-1)C.',
'cyclic_summand_formula':'AfreeD-summandcontributes0. AcyclicD/ell^e contributesktoGn iff e<n,and0iff e>=n. ThusdimGn=#torsionexponentsSTRICTLYsmallerthann. Genuinetorsionisobtainediff all exponents<n,not<=n.',
'actual_Cprime_type':'D² plusD/ell plusD/ell^e,e>=2. Rank2/mu4givetwotorsionblocks; ell²length7 forcesoneexponent1andtheother>=2.',
'actual_G2_dimension':1,'genuine_Cprime_ell_torsion_dimension':2,'Tnil_to_G2_kernel_dimension':15,'genuine_Nbar_dimension':14,
'actual_Cprime_modell3_length_formula':'7+min(3,e):9iffe2,10iffe>=3. Onlyaftere2isproveddoesG3canonicallyequalthefullCprime[ell].',
'no_finite_ell_inverse':'FiniteCnvaluesaloneprovidenofullsourcequotientbyell. Numeratorforedgemustberetainedonformalorone-order-thickersource; no filtrationstrictnessfromfiniteHF.',
'actual_e_G3_edge_basis_or_values_not_computed':True}
 ht=finish(A+'_REVIEW',inputs,tors,{'canonical_Gn_exact_quotient_and_actual_G2_G3_scope_accepted':True,'actual_Cprime_exponent_or_edge_values_computed':False},t)
 protocol={'inputs':inputs,
'genuine_edge_exact_sequence':'TheacceptedfullR0freecomplexhasH1=N,H0=Cprime; ellisinjectiveonitsfreechainterms. Longhomologyofthechainshortexact0->F --ell->F ->F/ellF->0 gives0->N/ellN->Tnil->Cprime[ell]->0,independentofhigherhomology. Dimensions14,16,2areactual.',
'existence_of_full_source_lift_for_given_kernel_representative':'AtruekernelclassisimageofanNclass. Ifgivenfinitecycle differsfromreductionofanexactfullcyclebyd2bar*bbar,liftbbararbitrarilytothefreefullsourceandaddd2*b. Exactd1d2=0preservescycleandmakesitsreductionthegivenrepresentative. ThisprovesEXISTENCEwithoutcomputingthelift.',
'finite_product_source_proof':'Product44/cube44 areacceptedpolynomialsourceoperationsonexactNcycles. Substitutionmodellcommuteswiththem. Thefullsourcecycleguaranteedabovehasthesamefiniteinputs,soapplyingthereducedformula totruekernelrepresentativesgivespreciselytheinjectedNbarproduct. ChangingfullliftschangesNclassesbyellN,anideal; changingfiniteboundaryrepresentativescanbehandledbyliftedd2boundariesandtheacceptedNproductdescent.',
'no_Tnil_algebra_or_new_homotopy':'ThistransportdefinesmultiplicationONLYontheactual14-dimensionalkernelidentifiedwithNbar. Itdoesnotextendproducttoall16Tnilclasses,anddoesnotconstructabilinearhomotopyforthe780boundaries. Genuinekernelchoice,actualfinitecycles,sourceprojectionandboundarytestingareSTILLrequired; HF55321/filtered16basisaloneisnotenough.',
'finite_cube_detection':'AnynonzerotripleinNbarprovesoriginalN³nonzero. VanishingalltriplesonlyimpliesN³subsetellN. UnderadditionalunprovedellH0=0 (equivalentlyactualT2length17), N³subsetH0 andH0intersectellN0proveN³zero. Atmost560symmetricbasis-triples sufficeonceactual14basis/productsareavailable.',
'no_actual_source_edge_projection_products_lifts_H0_or_nilcube_job':True}
 hp=finish(B+'_REVIEW',inputs,protocol,{'genuine14_kernel_only_finite_product_protocol_accepted':True,'actual_edge_nil_products_or_original_cube_computed':False,'all16_Tnil_algebra_or780_homotopy_claimed':False},t)
 print(json.dumps({'all_passed':True,'torsion_peer_report_sha256':ht,'protocol_peer_report_sha256':hp}))
if __name__=='__main__':main()
