#!/usr/bin/env python3
"""Independent SOURCE-only nonzero full-bracket pivot/Schur existence proof."""
import itertools,json,pathlib,re,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1];PRI='S61_COL53_NONZERO_PIVOT_CHART_SOURCE_20261105_v64';ID=PRI+'_REVIEW'
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None,parse=True):
  nonlocal raw
  p=pathlib.Path(p);k=str(p.relative_to(ROOT));s=sha(p)
  if h:assert h==s,(k,h,s)
  if k not in pins:raw+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text())if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json')
 for p,h in pr['input_sha256'].items():pin(ROOT/p,h,False)
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 def src(s):return ROOT/next(p for p in pr['input_sha256']if p.endswith(s))
 graph=pin(src('compact-face-graph.json'));owners=pin(src('exact-11-row-coefficient-operator-and-callback-owners.json'))
 common=pin(src('common-source-operands.json'));T=pin(src('actual-F10-pivot-and-full167x14-typed-operator.json'))
 vp=pin(ROOT/'reports/S61_COL53_H_FACES_VALUES_20261105_v64_REVIEW/report.json','d1708919e26ff8d484855c987d432cc4e5245cab8ae58164c5862bb44a348c26')
 assert vp['all_passed']and vp['independent_actual_signed_SOURCE_VALUE']and vp['fixed_entry_nonzero']
 vn=next(iter(vp['evidence_sha256']));v=pin(ROOT/'reports/S61_COL53_H_FACES_VALUES_20261105_v64_REVIEW'/vn,vp['evidence_sha256'][vn])
 assert v['explicit_exact_entry_nonzero_witness']=={'exponent5':[-3,0,6,-1,10],'coefficient':'241815/134217728'}
 assert graph['root']==1202 and len(graph['nodes'])==1203 and graph['nodes'][1202]['name']=='new_1013_sum'
 assert graph['nodes'][1202]['args'] # retain ALL full-DAG children, not H tied set
 assert owners['prefactor']=='P7/(U^2*F10)'and owners['stage8_vector_identity_exact']and owners['all167x14_and87plus66_RHS_P290once_historicalF9_F10_complements']
 p7=common['common_operands']['P7'];assert p7['QQ_transform_scalar']=='1'and p7['transform_monomial_exponent5']==[0]*5
 assert p7['original_literal']=={'content':'1','shift':[0]*5,'terms':[['5',0,0,0,2,0],['6',0,0,0,1,0],['1',0,0,0,0,0]]}
 assert T['actual_p10']=='F10/(353894400000*H*t0*F9)'and T['P290_applied_once_no_double_scaling']
 assert T['common_representation_ONLY']['U']=='Q^4*R^2*E'and T['common_representation_ONLY']['D10']=='U^2*F9*F10'
 full,top,sel,cols,extra,zero=[T[k]for k in['full167_row_ids','top101_row_ids','selected14_row_ids','residual14_column_ids','all87_other_top_RHS','all66_zero_constraints']]
 assert list(map(len,[full,top,sel,cols,extra,zero]))==[167,101,14,14,87,66]
 assert len(set(full))==167 and len(set(cols))==14 and 6 in sel and 53 in cols
 assert set(full)==set(top)|set(zero)and not set(top)&set(zero)
 assert set(top)==set(sel)|set(extra)and not set(sel)&set(extra)
 p=re.search(r'previous(\d+) typed units.*?(\d+) constraint units and(\d+) old units',T['inverse_substitution']);assert p
 previous,constraint,old=map(int,p.groups());assert(previous+1,constraint,old)==(10,10,49)
 assert 10+10+49+len(cols)==83 and 10+10+49+len(full)==236
 owned=count([pr,D,graph,owners,common,T,vp,v,pins])+25000;assert owned<150000 and raw+1000000<32<<20
 save(out/'positive-before-universal-identity-only.json',{'SOURCE_only_no_matrix_coefficient_update':True,'whole_numeric_upper':owned,'RAW_plus_output':raw+1000000,'no_actual_F11_full_map_recomputation':True})
 # Fresh polynomial arithmetic in12 independent variables. This checks the
 # complete antisymmetric cancellation before cancelling the OLD F9 unit.
 names=['P','f','c','d','x','y','u','v','w','z','s','t'];N=len(names)
 def sym(i):return{tuple(int(i==j)for j in range(N)):1}
 def plus(a,b,sgn=1):
  z=dict(a)
  for e,c in b.items():
   z[e]=z.get(e,0)+sgn*c
   if not z[e]:del z[e]
  return z
 def prod(a,b):
  z={}
  for e,c in a.items():
   for f,d in b.items():
    k=tuple(e[j]+f[j]for j in range(N));z[k]=z.get(k,0)+c*d
  return{k:c for k,c in z.items()if c}
 P,f,c,d,x,y,u,vv,w,z,s,tt=[sym(i)for i in range(N)]
 cf=prod(P,f)
 Nij=plus(prod(cf,x),prod(prod(c,y),u),-1)
 N58j=plus(prod(cf,vv),prod(prod(c,tt),u),-1)
 Ni56=plus(prod(cf,w),prod(prod(c,y),s),-1)
 N5856=plus(prod(cf,z),prod(prod(c,tt),s),-1)
 UF10=prod(d,N5856)
 lhs=plus(prod(UF10,Nij),prod(prod(d,Ni56),N58j),-1)
 bracket=plus(plus(prod(UF10,x),prod(prod(d,Ni56),vv),-1),prod(prod(prod(c,d),plus(prod(tt,w),prod(y,z),-1)),u))
 assert lhs==prod(cf,bracket)
 after=[len(full)-1,len(cols)-1,len(top)-1,len(sel)-1,len(extra),len(zero)]
 assert after==[166,13,100,13,87,66]
 assert D['one_prospective_pivot']['full_bracket']['all_nodes']==1203 and D['one_prospective_pivot']['full_bracket']['not_the_Hminus3_face_polynomial']
 assert D['one_prospective_pivot']['actual_entry']=='p = P7*F11_full_DAG/(U^2*F10)'
 assert D['one_prospective_pivot']['typed_ordinal']==11 and D['one_prospective_pivot']['original83_ordinal']==70
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'actual_nonzero_entry_witness_pinned_to_independent_signed_SOURCE_replay':True,'all1203_full_bracket_root1202_distinct_from411_face':True,
  'full_prefactor_antisymmetric_identity_verified_with_fresh12variable_integer_polynomials':True,
  'Schur_prefactor_derivation':'S9=N9/(U*F9). Thus p=(Nij*N5856-Ni56*N58j)/(U*F9*N5856). With U*F10=C10*N5856 and the full antisymmetric identity, this is P7*F11/(U^2*F10). F9 cancels algebraically but remains an OLD chart unit/complement.',
  'P7_old_unit':'5*q0^2+6*q0+1=(5*q0+1)*(q0+1)=Q*R. Old monomial units and Q,R,E,F9,F10,detA*detB,P290 are retained.',
  'new_chart':'D(F11_full_DAG) inside the SAME old chart. It is nonempty because all old unit factors and F11 are nonzero in the same Laurent parameter domain.',
  'explicit_pivot_inverse':'U^2*F10/(P7*F11_full_DAG)',
  'generic_unit_numbering':{'old49':49,'constraints10':10,'typed_already10':10,'prospective_typed11':11,'prospective_total70':70},
  'prospective_full166_by13_with_top100_selected13_extra87_zero66':[166,13,100,13,87,66],
  'formal_Schur_and_RHS_formula':'For all remaining rows i and columns j, Snew_ij=Sold_ij-Sold_i53*Sold_6j/p; hnew_i=hold_i-Sold_i53*hold_6/p; x53=(hold_6-sum_j Sold_6j*xj)/p. Restore prior10typed+10constraints+49old and oneP290. No coefficient or RHS update performed.',
  'same_cokernel_scope':'Unit pivot on FULL167x14 gives FULL166x13 with the same cokernel. Zero66 rows retain their consistency RHS; selected13 square is not identified with the full cokernel.',
  'all87_extra_RHS_and66_consistency_RHS_and_old_new_complements_retained':True,
  'no_new_matrix_values_or_RHS_or_fullinverse_rank_exhaustion_computed':True,'owned_upper':owned,'RAW':raw}
 ep=out/'independent-full-bracket-prefactor-unit-chart-and-all-manifest-complements.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_only':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'one_prospective_generic_pivot_exists':True,'actual_matrix_update_executed':False,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'seconds':report['seconds']}))
if __name__=='__main__':main()
