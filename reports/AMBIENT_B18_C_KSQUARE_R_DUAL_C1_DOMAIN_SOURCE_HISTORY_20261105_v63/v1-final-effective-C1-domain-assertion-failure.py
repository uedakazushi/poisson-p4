#!/usr/bin/env python3
"""Rebuild actual Pi1(head0,r) on raw128 free-zero wedge inputs."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
import ambient_B18_C_Kbar_r2_initial_source_20261105_v63 as I
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KSQUARE_R_DUAL_C1_DOMAIN_SOURCE_20261105_v63'
PATHS={
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'X':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'W':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'first_peer':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63_REVIEW/report.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'minimal_peer':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW/report.json',
 'rdual':'reports/AMBIENT_B18_C_KSQUARE_R_PAIR_DUAL_SOURCE_20261105_v63/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fs={'m0':['selected_source_columns24'],'X':['all780_source_columns44_sparse_dual9','X0_pivot_columns','all780_X0_factorization_coefficients'],
 'W':['residual_original_X_columns','W_pivot_residual_columns','W_pivot_rows','W_fixed_minor_inverse'],
 'minimal':['matrices','fixed_minor_data']}
 d={k:S.selective(ROOT/p,fs[k],receipts)if k in fs else json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert d['first_peer']['all_passed']and d['minimal_peer']['all_passed']
 assert 42 not in d['m0']['selected_source_columns24']
 # Directn=E(mE)^-1 therefore hasrow42 EXACTLY0 toallorders.
 X=d['X'];a=[]
 for col in X['all780_source_columns44_sparse_dual9']:
  row=next((v for j,v in col if j==42),[]);r=dict(row);assert not Q(r.get(0,'0'));a.append(I.parse(r.get(8,'0')))
 xp=X['X0_pivot_columns'];assert len(xp)==18
 residual=d['W']['residual_original_X_columns'];wp=d['W']['W_pivot_residual_columns'];wr=d['W']['W_pivot_rows'];assert len(wp)==len(wr)==81
 iv=[[I.parse(c)for c in row]for row in d['W']['W_fixed_minor_inverse']]
 assert len(iv)==81 and all(len(row)==81 for row in iv)
 selected=[]
 for k in wp:
  j=residual[k];co=list(map(I.parse,X['all780_X0_factorization_coefficients'][j]));v=a[j]
  for c,p in zip(co,xp):v=I.addq(v,-I.mulq(c,a[p]))
  selected.append(v)
 row={}
 for k,target in enumerate(wr):
  v=Q(0)
  for j,c in enumerate(selected):v=I.addq(v,-I.mulq(c,iv[j][k]))
  if v:row[target]=v
 assert row=={27:Q(-8,3)}
 evidence={'scope':'ExactallorderPi1(head0,r)domainSOURCEonTHISraw128(0,wedge), no product paired.',
 'input_sha256':pins,'direct_n_row42_identically_zero_from_selector':True,
 'all780_rawX0_row42_zero':True,'all780_rawX1_row42_r':[str(x)for x in a],
 'all81_selected_residual_relation_head0_r':[str(x)for x in selected],
 'actual_raw128_wedge84_Pi1_head0_r_sparse':[[i,str(c)]for i,c in row.items()],
 'first24_domain_proof':'Cancelunit(n,0) usingselectedfreeErows. Itinvolvesfree44inputONLY andleaveswedge84unchanged. Onrawinputfree44=0 theoriginfree20is0; directnrow42=0 identically.',
 'second18_domain_proof':'Itsselected18unitrowcoordinatesareinthefree20block. Onfreezero/wedgew theeliminationcoefficientsareEXACTLY0, notjustorigins. Thuswedgew isunchanged; retainedfree42hasX0row0 forALL780relations.',
 'third81_domain_proof':'Retainedhead0couplingto81unitcolumns hasorigin0; itsfirst_r isthese81residualX1row42values. HencePi1(head0,r)|wedge=-a1*W0_fixed_inverse*selector81. Derivativesofunitinverse aremultipliedbycoupling0 anddrop. ThedisplayedsingleH27entryisinORIGINALraw128wedgecoordinates.',
 'r_pair_domain_premise_now_verified':'rawiota2=(free44zero,H2), so THISactualdomainidentityapplies. NocomparisonofdifferentH2quotientsorunprovedSlinearity.',
 'SOURCE_max_pre_bits':I.MAXPRE,'selective_receipts':receipts}
 owned=S.count(d)+S.count(iv)+S.count(evidence)+15000;assert owned<150000;evidence['whole_SOURCE_owned_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'literal-raw128-domain-and-all81-origin-unit-derivative-reconstruction.json';p.write_text(json.dumps(evidence,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'raw128_C1_row':[[i,str(c)]for i,c in row.items()],'input_sha256':pins,'code_sha256':sha(Path(__file__)),'evidence_sha256':{p.name:sha(p)},'whole_SOURCE_owned_upper':owned,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
