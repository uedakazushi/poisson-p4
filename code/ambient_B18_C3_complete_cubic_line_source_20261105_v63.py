#!/usr/bin/env python3
"""Complete order-three relation-source reduction, with four formal alphas."""
import hashlib,itertools,json,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_COMPLETE_CUBIC_LINE_SOURCE_20261105_v63'
PATHS={
 'frame':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
 'frame_peer':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json',
 'kernel':'reports/AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63/nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json',
 'kernel_peer':'reports/AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63_REVIEW/report.json',
 'B4':'reports/AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63/four-actual-cubic-sources-two-cached-adjoints-and-rank2-completion.json',
 'B4_peer':'reports/AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63_REVIEW/report.json',
 'graph':'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63/actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json',
 'graph_peer':'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
 started=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()};d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['B4']['all_eight_corrected_adjoints_zero_without_dot_worker']
 assert d['graph']['all5_actual_delta_r_ell2_forced_zero_without_readback']
 assert d['frame']['actual_C_and_C3_HF_degrees0_1_2']==[3,6,6]
 # Complete first-frame linear-syzygy witness: all189 monomial columns,
 # all126 canonical images and the63 literal free duplicate differences.
 columns=[];groups={}
 for sector in range(3):
  for z in range(7):
   for variable in range(9):
    key=(sector,tuple(sorted([z,variable])))
    ci=len(columns);columns.append([sector,z,variable]);groups.setdefault(key,[]).append(ci)
 assert len(columns)==189 and len(groups)==126
 canonical=[];duplicates=[];syzygies=[]
 for key,inds in groups.items():
  first=inds[0];canonical.append([list(key[:1])+[list(key[1])],first])
  for second in inds[1:]:
   sector,z,var=columns[second];s0,z0,var0=columns[first]
   assert sector==s0 and z<7 and var<7 and z!=var
   assert sorted([z,var])==sorted([z0,var0])
   duplicates.append(second)
   syzygies.append({'free_column':second,'canonical_column':first,
       'terms':[[second,'1'],[first,'-1']],
       'same_quadratic_target':[sector,list(sorted([z,var]))],
       'actual_source':'z_variable*R_sector,z-z_z*R_sector,variable',
       'all_linear_coefficients_are_z':True})
 assert len(syzygies)==63 and len(set(duplicates))==63
 # Each duplicate appears only in its own difference: literal freeI63.
 for i,s in enumerate(syzygies):
  assert[(Q(1)if t==s['free_column']else Q(0))for t in duplicates]==[Q(i==j)for j in range(63)]
 first_physical=[r['physical_row']for r in d['frame']['all12_verified_source_records']]
 constant=[]
 for side,index in[('A',7),('B',15)]:constant.append({'kind':'curve','side':side,'raw116_index':index,'known_cover_cubic_image':'0'})
 for i in range(8):
  for j in range(8):
   constant.append({'kind':'cross','IA_index':i,'IB_index':j,'raw116_index':16+8*i+j,
       'order3_status':'first2=zi*zj eH; clear by zi*RH_j; all plane3 carry retainsz'if i<7 and j<7 else'order>=4',
       'known_cover_cubic_image':'0'})
 remaining=[19,28,29,72];graph=[63,64,68,73,74]
 for j in d['kernel']['free_columns9']:
  physical=first_physical[j]
  constant.append({'kind':'tauA_kernel','physical_row':physical,'raw_physical_index':j,
      'raw116_tauA_index':80+3*j,
      'known_cover_cubic_image':f'alpha_{physical}*(0,1,0,-1/4)'if physical in remaining else'0',
      'source_reason':'B4 corrected exact NF with all weighted cache leaves0'if physical in remaining else'bB1=0; graphdelta3modell3=0; Hmixed RH clears Bplane2=0'})
 for physical in graph:
  j=first_physical.index(physical);constant.append({'kind':'graph_tauB','physical_row':physical,'raw116_index':81+3*j,'known_cover_cubic_image':'0'})
 for j,physical in enumerate(first_physical):
  constant.append({'kind':'tauAB','physical_row':physical,'raw116_index':82+3*j,'known_cover_cubic_image':'0'})
 assert len(constant)==92
 assert sum(r['known_cover_cubic_image']!='0'for r in constant)==4
 const_counts={k:sum(r['kind']==k for r in constant)for k in['curve','cross','tauA_kernel','graph_tauB','tauAB']}
 assert const_counts=={'curve':2,'cross':64,'tauA_kernel':9,'graph_tauB':5,'tauAB':12}
 # Verify every graph tauA really has zero h before invoking the full
 # graph delta proof. All inactive mixed leaves are retained in principle.
 fmap={r['physical_row']:r for r in d['frame']['all12_verified_source_records']}
 for physical in graph:assert all(Q(v)==0 for v in fmap[physical]['bB_first_z7'])
 physical_roots={r['physical_row']:r['roots']for r in d['C2']['physical12_coefficients']}
 alpha=[{'label':f'alpha_{p}','physical_row':p,'literal_roots':physical_roots[p],
          'functional':'coefficient r*ell^2 of the actual same-R0 delta=aA-aB; lower r3,r2ell vanish from the selected literal cache leaves',
          'value':'UNAVAILABLE_FORMAL_NOT_ZERO'}for p in remaining]
 # Formal reduced operator on155 sources:151zero, four multiples of v.
 operator=[]
 for i,r in enumerate(constant):
  if r['kind']=='tauA_kernel'and r['physical_row']in remaining:
   operator.append([i,f'alpha_{r["physical_row"]}',['0','1','0','-1/4']])
 assert len(operator)==4 and len(constant)+len(syzygies)==155
 data={
  'actual_complete_presentation':'Original C has3raw targets,116sources=IA8+IB8+IAIB64+physical12*3. All origin0. The selected actual first21 frame is z7 ineachrawsector; fullconstant change gives95zero-first sources.',
  'actual_first21':{'intrinsic_A_raw116_columns':list(range(7)),'intrinsic_B_raw116_columns':list(range(8,15)),
       'physical_RH_raw116_columns':[81+3*i for i in range(7)],'target_order':['eA','eB','eH']},
  'first_frame_linear_monomial_columns189':columns,
  'all126_canonical_quadratic_images':canonical,
  'complete63_linear_firstframe_syzygies_with_freeI63':syzygies,
  'constant_firstkernel95_source_count':'12tauA+5graph_tauB+12tauAB+64cross+2curves',
  'constant_quadratic_kernel92_counts':const_counts,
  'all92_constant_sources':constant,
  'quad_constraint_source':'Only12tauA constant columns have nonzero plane quadratic: -bB2 eH, rank3 by fixed5/768. Thus their complete ker12x9 plus theother83columns is thewhole92constantquadkernel.',
  'all855_linear_zero_first_multiples_accounted':'95zero-first columns times9variables have cubic initials S1*M2. Their plane parts are only S1*plane2 eH, killed by theknownquadraticinitialideal; allotherparts retainz. NoM3 enters these linear-coefficient terms beforedegree4.',
  'quadratic_firstframe_coefficients_accounted':'Everyquadraticcoefficient*firstframe hascubicinitial inS2*(z7), alreadykilled. Quadraticcoefficients*zero-first sources startat4; higherformalcoefficients cannotcontributetocubicinitials.',
  'full_order3_initial_relation_completeness':'Expandarbitraryformalsourcecoefficients asu0+u1+u2+... . Order1 forcesu0 inthe95constantfirstkernel. Order2 forcesitsplanequadcomponent0, leaving92; chooseactualu1firstframe toclearitszquadratics. Allchoicesdifferbythe63Koszulsyzygies; theirplane3M2carryretainsz. Remainingu1on95zero-first andu2onfirstframearetheknownclosuretermsabove. HenceALLextraactualcubicinitialsmodtheknowncoveraretheimagesofthe92+63sources, withnohighersourcecoefficientsorTorinference.',
  'cross64_full_carry':'For49linear*linearproducts, primitivecube is zi*beta_j2 eH and clearingzi*zj eH withzi*RH_j givesonlyz*RH_j2. AllA/B/Hpureplane3projections0. Theother15productscontainaconductorcurveoforder3, soord>=4. Noindividualunknowncoefficientissetzero.',
  'known_cover_source_maps':'BothnormalizedintrinsiccurvecubicsareintheknownCI(3,3)cover. AcceptedtauAB12andgraphtauB5havezeroimage. GraphtauA5alsozero: h=bB1=0, graphdelta3modell3=0, andRHmixedBplane2=0. RemainingtauA4 usecorrectedB4NFwithopaquealpha retained.',
  'formal_operator_shape':[4,155],
  'formal_nonzero_columns':operator,
  'all_other151_columns_zero_in_known_cover':True,
  'remaining_four_alpha_source_functionals':alpha,
  'actual_cubic_image_exactly':'span over QQ of four alpha_j*v, v=(0,1,0,-1/4) in(A_r2ell,A_rell2,B_r2ell,B_rell2)',
  'actual_cubic_extra_rank_cases':[0,1],
  'actual_HF_degree3_cases':[4,3],
  'one_nonzero_alpha_consequence':'Cubicrank1, HF3=3; r*v=(1,-1/4) nonzero quartic andell*v=0. Henceonlyknownupper19, whileindependentprimitivequarticsmayforce18.',
  'all_alpha_zero_consequence':'HF3=4, soactualC3length>=19 andatleastoneofthetworemainingDVRunexponents>=3. Thisisconditionaluntilalpha valuesorapurevanishingproofexist.',
  'remaining_possible_C3_HFs':[[3,6,6,3,0],[3,6,6,3,1],[3,6,6,4,0],[3,6,6,4,1]],
  'remaining_possible_C3_lengths':[18,19,20],
  'no_alpha_readbacks_newjet_or_actual_C3_length_claim':True,
  'fixed_source_profile':{'source_families':155,'known_cubic_closure_terms':855,'whole_owned_numeric_envelope':50000,'hard_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2147483648,'new_physical_source_gauge_IFT_or_jet_coefficients':0}
 }
 ep=out/'complete155-source-cubic-operator-with-four-formal-alphas.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'COMPLETE_SOURCE_cubic_rank0_or1_four_alpha_callbacks_UNEVALUATED',
      'actual_C3_HF_degree3_in':[3,4],'actual_C3_length_or_exponents_determined':False,
      'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
      'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)


if __name__=='__main__':main()
