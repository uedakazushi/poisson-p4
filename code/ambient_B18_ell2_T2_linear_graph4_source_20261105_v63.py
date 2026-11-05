#!/usr/bin/env python3
"""Complete23-column linear graph and two invariant6-adjoint source gate."""
import argparse,itertools,json,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 t=read(H.PATHS['transport']);pref=read(H.PATHS['active']);res=read(H.PATHS['resolution']);g=read(H.PATHS['G2'])
 known=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json');mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json')
 read('reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_RANK_20261105_v63_REVIEW/report.json')
 adj=read('reports/AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63/actual-two-full918-source-adjoint-certificates-or-obstructions.json')
 lambdas=[{i//2:F(c) for i,c in r['lambda52']} for r in adj['results']];assert all(all(i%2 for i,c in r['lambda52']) for r in adj['results'])
 J=t['same_selected26_V_columns'];factors=list(map(H.parse,t['all56_new_column_factorizations_in_same_selectedJ_basis']));rfree=mini['remaining_free3_originalV_rows'];M=list(map(H.parse,t['actual_Mnew0_columns']))
 W=[add({r:F(1)},{J[i]:c for i,c in factors[r].items()},-1) for r in rfree];v=[F(-3375,11),F(-50),F(25,2)];w=mv(W,{i:c for i,c in enumerate(v)})
 assert w=={23:F(-3375,11),29:F(-50),30:F(25,2)} and not mv(M,w)
 selected=known['selected_bottom_columns'];bottom=list(map(H.parse,known['all642_bottom27_columns']));patterns=[];graph=[]
 for k,j in enumerate(selected):
  b=bottom[j];br={i:b.get(i*9+7,F(0)) for i in range(3) if b.get(i*9+7,F(0))};be={i:b.get(i*9+8,F(0)) for i in range(3) if b.get(i*9+8,F(0))}
  if br or be:
   cr=br.get(0,F(0))/v[0];ce=be.get(0,F(0))/v[0];assert br=={i:c*cr for i,c in enumerate(v) if c*cr} and be=={i:c*ce for i,c in enumerate(v) if c*ce};patterns.append((k,j,cr,ce))
  else:cr=ce=F(0)
  rec=mini['all642_zero_constant_combinations'][j]
  graph.append({'bottom_basis_column':k,'source642_column':j,'original768_combination':rec['original768_combination'],
   'known_b1_3_times9':srec(b),'a1_operator':'first27 rows of fixedB30inv applied to actualX1_C, allFULL9firsts alreadybound',
   'A_two_source_callbacks':'Gamma_t*(Xgamma2-Xp1*a1), onlyrell; sameactualold104/nativeC16sourceframework, mustbindnew23callbacks beforevalues',
   'four_unknown_mu_coefficients_per_channel':{'M_r_H0v':str(ce),'M_ell_H0v':str(cr)}})
 assert patterns==[(19,582,F(1),F(14,125)),(21,634,F(8,3),F(32,75))]
 # Actual G1linear in the retained physical3target and current26sourceframe.
 G2,G1a,_,_=H.full_G2(t,pref,res,g);G1=list(G1a)
 for z in pref['zero_actions']:
  for i in range(3):G1.append({i*9+z:F(1)})
 assert len(G1)==26
 mons=[(i,j) for i in range(9) for j in range(i,9) if (i,j)!=(8,8)];mi={m:i for i,m in enumerate(mons)}
 for col in G2:
  q={}
  for k,c in col.items():
   row,var=divmod(k,9)
   for p,a in G1[row].items():
    target,u=divmod(p,9);m=tuple(sorted((u,var)))
    if m!=(8,8):q=add(q,{target*44+mi[m]:a*c})
  assert not q
 # Unknown F1=M1*H0v has26*9 coefficients. These are the two genuine gaugeinvariants.
 invariants=[{},{}]
 for i,c in lambdas[0].items():invariants[0]=add(invariants[0],{i*9+7:c});invariants[1]=add(invariants[1],{i*9+8:c})
 for i,c in lambdas[1].items():invariants[1]=add(invariants[1],{i*9+7:-c})
 pairs=[(7,7),(7,8)];cols=[{} for _ in range(6)]
 for j,col in enumerate(G1):
  for k,c in col.items():
   row,u=divmod(k,9)
   for var in range(9):
    p=tuple(sorted((u,var)))
    if p in pairs:cols[2*row+pairs.index(p)]=add(cols[2*row+pairs.index(p)],{j*9+var:c})
 primitive='reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase2-all105-equations-and-primitive56-source-columns.json';hs[primitive]=H.sha(ROOT/primitive)
 slots=sum(map(len,G1+cols+invariants+W))+sum(len(r['original768_combination'])+len(r['known_b1_3_times9']) for r in graph)+3*6*6+2*234+10000;assert slots<20000
 data={'all23_linear_graph_source_records':graph,'common_bottom_direction_v':list(map(str,v)),
  'H0_free3_columns':[srec(q) for q in W],'actual_H0v_source56':srec(w),'M0_H0v_zero':True,
  'unknown_four_canonical_contractions':'mu[t,d]=lambda_t*(partial_d Mnew)*H0v, t0/1,d=rA/ell; noneassigned0',
  'exact_graph_formula':'g[t,j]=A[t,j]+ce[j]*mu[t,r]+cr[j]*mu[t,ell]',
  'gauge_source_change':'Mprime=M+G2*A0; adjustsectionbyitsunit; truecomparisonandT2unchanged, canonicalvaluesnotset0',
  'gauge_shift_order_r0_r1_ell0_ell1':'(0,u,u,v), u=b*w1+w2,v=a*w0+c*w1, w=PiY*A0*H0v',
  'gauge_constants':{'a':'1095/3652','b':'1224288/8979355','c':'1973552256/40990755575'},
  'two_gauge_invariants':'alpha=mu[0,r], beta=mu[0,ell]-mu[1,r]',
  'gauge_normal_form':'chooseu=-mu[1,r],v=-mu[1,ell] viaA0col23=(-11/3375)*k, ksupportedinfreeYrows1,38 withPiYk=(-mu[1,ell]/a,0,-mu[1,r]); thenmuprime=(alpha,0,beta,0). Thisisacoherentunitframechoice, notcanonicalzero.',
  'exact_equation_for_F1':'F1=Mnew1*H0v; G1linear*F1=(V*H0v)_2 becauseG1Mnew=V andMnew0H0v=0',
  'actual_G1linear3_by26_columns':[srec(q) for q in G1],
  'two_invariant_functionals_on_all234_F1_coefficients':[srec(q) for q in invariants],
  'six_rr_rell_source_adjoint_columns_in234_domain':[srec(q) for q in cols],
  'primitive_V2_input_provenance':{'path':primitive,'field':'all56_primitive_V_columns','source_combination':srec(w),'coefficient_monomials':'all44ordinarydegree2fullCA_F9(z7,r,ell), cachedB56; neededrhoinputsrr/rell onlyifverified'},
  'coarse_I2_closure_inspection':'keepingonlyknownG3toplinearrelationsandI²annleavesbottomfreeB³; B/I² surjectsactualraw33, so thiscoarsecoverhaslength>=99, cannotimprovetheactualupper20. Retain23graphrelationsorproveuniformclosure; noactualT2lower99claim.',
  'four_canonical_M1_values_and_A23_not_computed':True,'six_adjoint_solve_not_started':True,
  'whole_six_adjoint_slot_bound':slots,'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':150000,'bits':384}}
 ep=out/'actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'canonical_unknown_contractions':4,'gauge_invariants':2,'affected_graph_columns':[19,21],
  'six_adjoint_solve_started':False,'new_physical_or_M1_values_computed':False,'whole_six_adjoint_slot_bound':slots,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','canonical_unknown_contractions','affected_graph_columns','whole_six_adjoint_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
