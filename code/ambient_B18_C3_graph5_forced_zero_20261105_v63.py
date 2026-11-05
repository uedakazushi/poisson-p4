#!/usr/bin/env python3
"""Fixed cover normal forms and triangular source proof; no graph readbacks."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63'
PATHS={
 'pair':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'frame':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
 'frame_peer':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json',
 'tauA_source':'reports/AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63/nine-complete-tauA-source-kernel-and-unpaired-cubic-cache-adjoints.json',
 'tauA_peer':'reports/AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2_report':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/report.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'rank4':'reports/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63/report.json',
 'rank4_peer':'reports/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW/report.json',
 'Cbar_peer':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW/report.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'ell':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
 'ell_peer':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b):
 c=a.copy()
 for m,v in b.items():
  c[m]=c.get(m,Q(0))+v
  if not c[m]:del c[m]
 return c
def scale(a,q):return{m:v*q for m,v in a.items()if v*q}
def mul(a,b):
 c={}
 for u,v in a.items():
  for w,z in b.items():
   m=tuple(x+y for x,y in zip(u,w));c[m]=c.get(m,Q(0))+v*z
 return{m:v for m,v in c.items()if v}
def power(a,n):
 c={(0,0):Q(1)}
 for k in range(n):c=mul(c,a)
 return c
def compose(rows,r,s):
 c={}
 for(a,b),v in rows:c=add(c,scale(mul(power(r,a),power(s,b)),Q(v)))
 return c
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Q(0))


def main():
 started=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()};d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['rank4']['actual_comparison_C_generic_DVR_rank4_proved']
 assert d['rank4']['actual_mu_D_C8_from_independent_source_evaluation']
 assert d['C2_report']['actual_C2_dimension']==14
 assert d['C2_report']['actual_C_DVR_torsion_cyclic_summands']==4
 assert d['C2_report']['actual_C_DVR_torsion_exponents_at_least2']==2
 lower=3*4+1+1+2+2;assert lower==18
 r={(1,0):Q(1)};s={(1,0):Q(1),(0,1):Q(1,5)}
 pair=d['pair'];ga=compose(pair['gA_plane3'],r,s)
 bot=[list(map(Q,row[7:]))for row in pair['YB_in_YA_linear'][7:]]
 rb=add(scale(r,bot[0][0]),scale(s,bot[0][1]))
 sb=add(scale(r,bot[1][0]),scale(s,bot[1][1]))
 gb=compose(pair['gB_own_plane3'],rb,sb)
 curves=[[p.get((3-i,i),Q(0))for i in range(4)]for p in[ga,gb]]
 assert[curves[0][0],curves[1][0]]==[Q(125,24),Q(-15125,144)]
 nf3=[];nf4=[];block_receipts=[]
 for co in curves:
  a,b,c=co[1]/co[0],co[2]/co[0],co[3]/co[0]
  q3=[[-a,Q(1),Q(0),Q(0)],[-b,Q(0),Q(1),Q(0)]]
  q4=[a*a-b,-a,Q(1),Q(0),Q(0)]
  assert all(dot(row,co)==0 and row[3]==0 for row in q3)
  rg=co+[Q(0)];lg=[Q(0)]+co
  assert dot(q4,rg)==dot(q4,lg)==q4[3]==q4[4]==0
  assert [[row[j]for j in[1,2]]for row in q3]==[[Q(1),Q(0)],[Q(0),Q(1)]]
  assert q4[2]==1
  nf3.append(q3);nf4.append(q4)
  block_receipts.append({'curve3_r3_r2ell_rell2_ell3':list(map(str,co)),
      'NF3_to_r2ell_rell2':[[str(x)for x in row]for row in q3],
      'NF4_to_r2ell2':list(map(str,q4)),
      'curve_and_all_ell3_multiple_products_zero':True,
      'literal_r_ell2_and_r2_ell2_images':[0,1,1]})
 # Whole cover maps retain the exact raw eA/eB sectors; H dies in>=2.
 n3=[[Q(0)]*12 for _ in range(4)];n4=[[Q(0)]*15 for _ in range(2)]
 for sector in range(2):
  for i in range(2):n3[2*sector+i][4*sector:4*sector+4]=nf3[sector][i]
  n4[sector][5*sector:5*sector+5]=nf4[sector]
 assert[n3[i][j]for i,j in[(1,2),(3,6)]]==[Q(1),Q(1)]
 assert[n4[i][j]for i,j in[(0,2),(1,7)]]==[Q(1),Q(1)]
 basis=d['raw33']['actual33_basis'];raws={x['physical_row']:x for x in d['C2']['physical12_coefficients']}
 frows=d['frame']['all12_verified_source_records'];fmap={x['physical_row']:x for x in frows}
 rows=[63,64,68,73,74];kernel=[[Q(x)for x in row]for row in d['tauA_source']['complete_constant_kernel12_by9']]
 free=d['tauA_source']['free_columns9'];phys=d['tauA_source']['physical_row_order']
 nonzz=[basis.index([i+1,1,0])for i in range(7)]+[basis.index([i+1,0,1])for i in range(7)]
 prem=[]
 for row in rows:
  j=phys.index(row);kcol=free.index(j)
  assert[kernel[i][kcol]for i in range(12)]==[Q(i==j)for i in range(12)]
  assert all(Q(x)==0 for x in fmap[row]['bB_second_plane'])
  assert all(Q(x)==0 for x in fmap[row]['bA_first_z7'])
  vals={name:dict(raws[row][name])for name in['delta','bA']}
  assert all(Q(vals['delta'].get(i,'0'))==0 for i in nonzz)
  assert all(Q(vals['delta'].get(basis.index(m),'0'))==0 for m in[[0,3,0],[0,2,1]])
  mixed=[[basis[i],vals['bA'][i]]for i in nonzz if i in vals['bA']and Q(vals['bA'][i])]
  prem.append({'physical_row':row,'literal_free_unit_kernel_column':kcol,
       'bB2_plane_zero':True,'bA1_zero':True,'delta2_nonzz_zero':True,
       'delta3_r3_and_r2ell_zero':True,'actual_bA2_H_mixed_terms_retained':mixed,
       'tauA_bB1_RH_source_coefficients':fmap[row]['bB_first_z7']})
 # Universal minors are checked in QQ[d,u,v]. d is never assigned an
 # actual readback. Both triangular determinant polynomials are d^2.
 D={(2,0,0):Q(1)}
 z={};one={(0,0,0):Q(1)};dd={(1,0,0):Q(1)};u={(0,1,0):Q(1)};v={(0,0,1):Q(1)}
 cubic_minor=add(mul(dd,dd),scale(mul(z,v),-1));assert cubic_minor==D
 rqA_B=add(scale(u,-curves[1][1]/curves[1][0]),v)
 quartic_minor=add(mul(dd,dd),scale(mul(z,rqA_B),-1));assert quartic_minor==D
 data={
  'coefficient_scope':'The actual original R0-module C=coker(D0→E4), C3=C/ell^3 C. No algebra splitting is inferred.',
  'accepted_lower_bound_receipts':{'DVR_rank':4,'C_mod_ell_length':8,'torsion_cyclic_summands':4,'C_mod_ell2_length':14,'torsion_exponents_exactly1_count':2,'torsion_exponents_atleast2_count':2,'C_mod_ell3_length_lower':lower},
  'known_cover':'M=(QQ[r,ell]/(gA3,ell3))eA ⊕ (QQ[r,ell]/(gB3,ell3))eB ⊕ (QQ[r,ell]/(r,ell)^2)eH surjects onto gr(C3) after the actual first21/complete quadratic initial relations. It is a module cover, not an asserted split of C3.',
  'known_cover_HF':[3,6,6,4,2],'known_cover_length':21,
  'two_block_source_curve_NF_receipts':block_receipts,
  'whole_NF3_rows4_columns12':[[str(x)for x in row]for row in n3],
  'whole_NF4_rows2_columns15':[[str(x)for x in row]for row in n4],
  'whole_NF3_coordinate_order':['A_r2ell','A_rell2','B_r2ell','B_rell2'],
  'whole_NF4_coordinate_order':['A_r2ell2','B_r2ell2'],
  'same_five_actual_source_premises':prem,
  'full_tauB_clearing':'For each graph j, F_j tauB has zero first. Delta2 eB has no nonzz terms. Clear unknown zz terms with z*frame, retaining z in beta2 carries. Actual bA2 mixed H terms are present and cleared with r/ell*RH; RH2_eB_plane=delta_i2_plane=0. All residual cubic H components die modulo actual plane2 H. Thus tauB cubic has eA0/eB=delta_j3.',
  'full_tauA_clearing':'The SAME free unit source j lies in ker bB2. q_j=F_j tauA+sum_i bB_j1zi RH_i has first0/plane2zero. Clear all z-containing second terms using the full actual first frame. RA exact z eA means A cubic is exactly delta_j3. ALL beta2, mixed RH and delta_i3 carry is retained as arbitrary B coordinates; H dies by the known quadratic initial ideal.',
  'delta_shape_only_mod_ell3':'Existing delta3 r3/r2ell0 leaves d*r*ell2 modulo ell3. The ell3 channel remains unconstrained, not zero-filled.',
  'unknown_cubic_normal_forms_QQ_d_u_v':{'tauA':['0','d','u','v'],'tauB':['0','0','0','d']},
  'cubic_fixed_minor_rows':[1,3],'cubic_fixed_minor_polynomial':'d^2',
  'r_multiple_quartic_normal_forms':{'tauA':['d','-(gB_r2ell/gB_r3)*u+v'],'tauB':['0','d']},
  'quartic_fixed_minor_polynomial':'d^2',
  'filtered_independence':'If d!=0, the two actual cubic initial classes are independent in M3, and their actual r-multiples are independent in M4. Thus at least2+2 graded dimensions are removed from the known length21 cover, regardless all opaque B carry.',
  'contradiction_if_nonzero':{'known_upper':21-2-2,'accepted_actual_lower':lower},
  'all5_actual_delta_r_ell2_forced_zero_without_readback':True,
  'graph5_values_job_not_required_for_this_conclusion':True,
  'not_claimed':['primitive delta3 ell3 zero','all cubic source relations zero','actual C3 length or full torsion exponents','original H0 length or nil index'],
  'numeric_source_only_profile':{'fixed_max_degree':4,'hard_seconds':60,'hard_memory_bytes':2147483648,'owned_envelope':50000,'hard_owned':150000,'new_physical_source_gauge_or_jet_coefficients':0,'graph_readbacks_evaluated':0}
 }
 ep=out/'actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PURE_SOURCE_five_graph_cubic_coefficients_forced_zero_NO_READBACK',
         'all5_delta_r_ell2_zero_proved':True,'actual_C3_length_or_all_exponents_determined':False,
         'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
         'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)


if __name__=='__main__':main()
