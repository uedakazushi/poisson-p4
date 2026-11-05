#!/usr/bin/env python3
"""Complete full-C cubic source: four formal pure-ell3 callbacks on one line."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_FULL_CUBIC_LINE_SOURCE_20261105_v63'
PATHS={
 'complete':'reports/AMBIENT_B18_C3_COMPLETE_CUBIC_LINE_SOURCE_20261105_v63/complete155-source-cubic-operator-with-four-formal-alphas.json',
 'complete_peer':'reports/AMBIENT_B18_C3_COMPLETE_CUBIC_LINE_SOURCE_20261105_v63_REVIEW/report.json',
 'frame':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
 'frame_peer':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'alpha':'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63/eight-cached-even-root-values-and-four-complete-alpha-readbacks.json',
 'alpha_peer':'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'reduction':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63/full-C-ell-reduction-and-conditional-graph03-exponent3-proof.json',
 'reduction_peer':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63_REVIEW/report.json',
 'graph03':'reports/AMBIENT_B18_C_FULL_GRAPH03_FORCED_ZERO_20261105_v63/report.json',
 'graph03_peer':'reports/AMBIENT_B18_C_FULL_GRAPH03_FORCED_ZERO_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()}
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['graph03']['all_five_graph03_coefficients_zero_source_proved']
 complete=d['complete'];assert complete['constant_firstkernel95_source_count']==95
 assert complete['constant_quadratic_kernel92_counts']=={'curve':2,'cross':64,'tauA_kernel':9,'graph_tauB':5,'tauAB':12}
 assert len(complete['all92_constant_sources'])==92 and len(complete['complete63_linear_firstframe_syzygies_with_freeI63'])==63
 b33=d['raw33']['actual33_basis'];byrow={r['physical_row']:r for r in d['C2']['physical12_coefficients']}
 first={r['physical_row']:r for r in d['frame']['all12_verified_source_records']}
 alphas={r['physical_row']:r for r in d['alpha']['all4_actual_alpha_delta_B6']}
 nf=[[Q(s)for s in row]for row in d['reduction']['full_no_cut_cubic_NF6_rows12_columns']]
 pure=[nf[k][3]-nf[k][7]/4 for k in range(6)];assert pure==[0,0,1,0,0,Q(-1,4)]
 active_zeros=[];four=[]
 for old in complete['remaining_four_alpha_source_functionals']:
  row=old['physical_row'];vals=dict(byrow[row]['delta']);ix=list(first).index(row)
  assert list(map(Q,first[row]['bA_first_z7']))==[Q(k==ix)for k in range(7)]
  assert list(map(Q,first[row]['bB_first_z7']))==[Q(-1,4)*Q(k==ix)for k in range(7)]
  for z in[0,4,5]:
   for mon in[[z+1,1,0],[z+1,0,1]]:
    bi=b33.index(mon);v=Q(vals.get(bi,'0'));assert not v;active_zeros.append([row,mon,'0'])
  assert all(Q(vals.get(b33.index(mon),'0'))==0 for mon in[[0,3,0],[0,2,1]])
  assert Q(alphas[row]['delta_r_ell2'])==0
  four.append({'physical_row':row,'literal_roots':old['literal_roots'],
   'new_functional':'coefficient ell^3 of SAME actual delta=aA-aB at r=0; value remains FORMAL',
   'exact_actual_source':'F_j*(tauA-tauB/4) with evenpart removed by literal R0 unit target',
   'full_cubic_NF':'eta_'+str(row)+'*(0,0,1,0,0,-1/4)',
   'delta3_r3_r2ell_rell2':'all0 by reviewed C2 and alpha4 values; only ell3 remainsformal'})
 assert len(active_zeros)==24 and len(four)==4
 columns=[]
 for old,callback in zip(complete['formal_nonzero_columns'],four):
  assert old[0]in[66,67,68,69]
  columns.append([old[0],'eta_'+str(callback['physical_row']),list(map(str,pure))])
 categories=[
  {'kind':'curve','count':2,'full_no_cut_proof':'The two ACTUAL intrinsic curve cubic initials are precisely the relations defining NF6; their ell3 coefficients are retained in NF6, not discarded.'},
  {'kind':'cross','count':64,'full_no_cut_proof':'49zi*zjH initials clear with zi*RHj; every resulting cubicA/B carry retainszi and hasplane0. Other15 products start>=4. Replacing any remainingz bybeta2 only raisesdegree to>=4.'},
  {'kind':'tauA_kernel_graph','count':5,'full_no_cut_proof':'All four plane cubic delta coefficients now SOURCEzero, including NEWgraphell3. Hmixed clears haveRH second-planeB0; other z clears retainz.'},
  {'kind':'tauA_kernel_free','count':4,'full_no_cut_proof':'Literal q=Fj(tauA-tauB/4). A3=delta3. RB evaluation of B2=-delta2/4 adds +sum(delta2_zi*x*x*beta_i2)/4. Onlyi0,4,5 havebeta2; all24active mixed delta2 leaves are literally0. Hmixed RHclears haveB2_plane0. Thus B3=-delta3/4 withNOell3cut.'},
  {'kind':'graph_tauB','count':5,'full_no_cut_proof':'All graphdelta plane3 anddelta2mixed vanish; RBbeta carry0 by the reviewed fullgraph proof.'},
  {'kind':'tauAB','count':12,'full_no_cut_proof':'PrimitiveA/B products nu*b after own A/B graph evaluation start>=4. Full plane cubicH is a linear multiple of ACTUALH quadratics P0/P1/P2; subtracting this changesA/B only at4. This includesell3H=ell*P2 INITIAL, without adding a quotient cut.'},
  {'kind':'firstframe_syzygy','count':63,'full_no_cut_proof':'The three literalz7 Koszul blocks havez coefficients. Their M2carry has a z factor; any RB z→-beta2 substitution raisesit to>=4. Hence plane3A/B0, independent ofell3.'},
 ]
 assert sum(x['count']for x in categories)==155
 data={
  'scope':'COMPLETE degree3 initial space of FULL comparison C, NOT C/ell³C. Same actual R0, original116presentation and complete155source coordinates.',
  'no_cut_NF6_coordinates':['A*r²ell','A*r ell²','A*ell³','B*r²ell','B*r ell²','B*ell³'],
  'full_no_cut_NF6':[[str(v)for v in row]for row in nf],
  'whole_source_counts':{'first21':21,'zero_first_constant_kernel':95,'quadratic_constant_kernel':92,'linear_firstframe_syzygies':63,'complete_extra_cubic_columns':155,
   'zero_first_times_linear_coefficients':855,'firstframe_times_quadratic_coefficients':945},
  'complete_no_cut_source_categories':categories,
  'literal24_active_mixed_delta2_zero_checks':active_zeros,
  'four_FORMAL_pureell3_callbacks':four,
  'complete_formal_operator_shape':[6,155],'only_possible_nonzero_columns':columns,'all_other151_columns_zero':True,
  'full_formal_coefficient_completeness':'Expand arbitrary source coefficientsu0+u1+u2+... . Asinthe reviewed155theorem, order1 leaves95; order2 leaves92 after actual firstframe clearing; allchoices differ by63syzygies. Remaining855zero-firstlinear multiples and945firstframequadratic multiples are known closure. This bookkeeping is before any ell3cut; adding neither arbitrarysourceu3 nor a Tor assumption is needed atdegree3.',
  'full_C_HF_degree3':'5 if at least one eta_j!=0;6 if all4eta_j=0',
  'all_ell3_callbacks_UNREAD':True,
  'comparison_DVR_exponent_scope':'This decides full HF3 and mu_D Tor(m³C), not automatically the remaining torsion exponent e>=3. No exact exponent, original H0 or nilcube inference.',
  'no_new_physical_IFT_gauge_or_jet':True,
  'fixed_profile':{'degree':3,'numeric_envelope':50000,'hard_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2147483648,'hard_RAW_bytes':33554432},
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'complete-full-C-cubic155-operator-with-four-formal-ell3-callbacks.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_FULL_C_CUBIC_LINE_FOUR_PURE_ELL3_CALLBACKS_UNREAD',
  'complete_full_C_cubic_operator_proved':True,'actual_eta_values_evaluated':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
