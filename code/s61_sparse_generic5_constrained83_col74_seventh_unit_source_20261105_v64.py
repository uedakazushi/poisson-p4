#!/usr/bin/env python3
"""SOURCE-only seventh typed pivot and next-column cost, no new matrix values."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,hashlib,time,signal,resource
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL74_SEVENTH_UNIT_SOURCE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def poly(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']}
def mul(A,B):
 R={}
 for e,q in A.items():
  for f,z in B.items():
   g=tuple(x+y for x,y in zip(e,f));R[g]=R.get(g,Q(0))+q*z
 return{e:q for e,q in R.items()if q}
def scale(A,q,shift=(0,0,0,0,0)):return{tuple(e[j]+shift[j]for j in range(5)):q*z for e,z in A.items()}
def main(a):
 assert a.build
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py')}
 def gate(stem,h):
  p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];inputs[p]=h;return p,r
 def ev(p,r,name):
  v=p.parent/name;assert sha(v)==r['evidence_sha256'][name];inputs[v]=sha(v);return read(v)
 p,s=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SIXTH_UNIT_COL74_SOURCE_20261105_v64','19a897e17f89aa55e68a0487d5455425f397b169e07b1cb8c6766fed54f8feb4');op=ev(p,s,'sixth-unit-and-full171x18-typed-operator.json');cost6=ev(p,s,'fixed18-source-cost-comparison.json')['all18_uniform_clearQ2_source_costs'];common=ev(p,s,'common-source-operands.json')['common_operands'];row51=ev(p,s,'rows/row-51-source-operands.json')['row_operands'];unit=ev(p,s,'old-unit2-monomial-source-prefactors.json')
 cp,cr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SIXTH_COL74_COEFFICIENTS_20261105_v64','96c24538f7f6dfce75b78a25cf755895116023b9fd46861d8a04410373d467cd');api=ev(cp,cr,'actual-cleared-col74-and-commonQ2-representation.json');values={z['row']:z for z in api['actual18_col74_numerator_values']};v=values[51];file=cp.parent/v['file'];assert sha(file)==v['sha256'];inputs[file]=sha(file);blob=file.read_bytes();assert blob[:9]==b'Q5CELL1\n\x01';literal=decode(blob[9:]);assert len(literal['terms'])==5
 Qp=poly(common['Q']['original_literal']);Q2=poly(common['Q2']['original_literal']);R={(0,0,0,1,0):Q(1),(0,0,0,0,0):Q(1)};Q3R=mul(mul(Q2,Qp),R);factor=scale(Q3R,Q(-64,3645),(2,0,-1,0,-1));assert factor==poly(literal)
 for label in('U51','V51','C51','F51'):assert not row51[label]['original_literal']['terms']
 assert not row51['A51_51']['original_literal']['terms'];raw49=poly(row51['A51_49']['original_literal']);assert raw49==scale(Qp,Q(2,6075),(1,0,0,0,0))
 # Literal independent first-two-Schur prefactors establish the simplified
 # WHOLE row/RHS alias; symbols X=A51j,Y=A11j,Z=A44j remain untouched.
 pref=[(Q(z['QQ']),tuple(z['exponent5']))for z in unit['all_four_inner_prefactors']];assert pref==[(Q(1),(0,0,0,0,0)),(Q(324000),(0,0,-1,0,-1)),(Q(36000),(0,0,-1,0,-1)),(Q(-675,8),(-1,0,0,-1,0))]
 assert Q(2,6075)*(Q(324000)-8*Q(36000))==Q(320,27)and Q(2,6075)*Q(36000)==Q(320,27)
 ep,er=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64','a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a');E={(z['row'],z['column']):z for z in ev(ep,er,'actual-fixed24-entry-source-and-exponent-plan.json')['actual_entry_sources576']}
 dp,dr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_ROWS11_44_LITERAL_SOURCE_20261105_v64','0ca442b4eb828e535bdb527357826c6f4cd7edb48f23735f08742b798ada9371');dv=ev(dp,dr,'whole-cached-rows11-44-exact-comparison.json');parse=25000+numeric(dv)+numeric(E)+numeric(op)+numeric(common)+numeric(row51)+numeric(api);assert parse<150000;nd={z['column']:len(z['difference_terms'])for z in dv['all24_complete_actual_row_comparisons']};del dv
 rows=[i for i in op['selected18_row_ids']if i!=51];cols=[j for j in op['residual18_column_ids']if j!=74];full=[i for i in op['full171_row_ids']if i!=51];top=[i for i in op['top105_row_ids']if i!=51];assert len(rows)==len(cols)==17 and len(full)==170 and len(top)==104
 nf=sum(values[i]['terms']for i in rows);assert nf==4412;costs=[]
 for z in cost6:
  j=z['column']
  if j==74:continue
  old51=3*(E[51,j]['terms']+2*(E[44,j]['terms']+nd[j]));row51j=E[51,j]['terms']+2*(E[11,j]['terms']+E[44,j]['terms']);c=3*(z['literal_tuple_majorant_clearQ2']-old51)+nf*row51j;assert c>=0;costs.append({'column':j,'sixth_source_total_excluding51':z['literal_tuple_majorant_clearQ2']-old51,'simplified_M2_row51_source_tuple_bound':row51j,'new17_fixed_source_tuple_majorant':c})
 costs.sort(key=lambda z:(z['new17_fixed_source_tuple_majorant'],z['column']));assert len(costs)==17
 evidence={'all_passed':True,'source_only':True,'pivot_row':51,'pivot_column':74,'actual5_term_numerator':literal,'literal_Qcubed_times_qplus1_factor_identity':True,'actual_S6_51_74':'-64*H^2*(5*q0+1)*(q0+1)/(3645*s0*t0)','inverse':'-3645*s0*t0/(64*H^2*(5*q0+1)*(q0+1))','additional_chart_factor':'q0+1','chart_complement_qplus1_OPEN':True,'all_previous_chart_complements_OPEN':True,'zero_accepted_U51_V51_C51_F51_proves_full_row_RHS_alias':'S6_51,j=M2_51,j, h6_51=h2_51. NOT raw row51.','source_whole_row51_formula':'M2_51,j=A51,j+(320/27)*H*(5q0+1)/(s0*t0)*(A11,j+A44,j)','full170_row_ids':full,'top104_row_ids':top,'selected17_row_ids':rows,'residual17_column_ids':cols,'all87_other_top_RHS':op['all87_other_top_RHS'],'all66_zero_constraints':op['all66_zero_constraints'],'all170x17_RF_roots':'S7_i,j=S6_i,j-(S6_i74/p7)*M2_51,j; original source callbacks retained.','RHS':'h7_i=h6_i-(S6_i74/p7)*h2_51, all87+66 consistency retained.','substitution':'x74=(h2_51-sum_j M2_51,j*xj)/p7, then restore prior six pivots, one P290 scale, ten constraints and49 units.','cokernel':'full171x18 unit-equivalent tofull170x17; full170coker=top104coker direct_sum K5^66. Selected17 square is NOT samecokernel.','future_common_representation':'Q^3*(q0+1)*S7_i,j=Q*(q0+1)*F6_i,j - J7*F6_i74*M2_51,j, J7=-3645*s0*t0/(64*H^2). Representation ONLY, no Qcancellation assumed.','all_other16_columns_unchanged_typed':True,'P290onlyonce':True,'no_new_next_column_values_rank_inverse_or_exhaustion_claim':True}
 (out/'actual-row51-seventh-unit-and-full170x17-operator.json').write_bytes(enc(evidence));(out/'remaining17-literal-source-cost-comparison.json').write_bytes(enc({'all_passed':True,'source_only':True,'all17_fixed_source_tuple_majorants':costs,'selected_minimum_in_this_frozen_source_model':costs[0]['column'],'only_metadata_no_new_QQ_matrix_cells':True,'future_support_height_custody_storage_gate_REQUIRED_before_values':True}));report={'id':ID,'all_passed':True,'source_only':True,'seventh_actual_unit_factor_verified':True,'additional_chart_factor':'q0+1','full_residual_shape':[170,17],'selected_subblock_shape':[17,17],'next_source_cheapest_column':costs[0]['column'],'next_source_tuple_majorant':costs[0]['new17_fixed_source_tuple_majorant'],'source_numeric_upper':parse,'actual_READ_RAW_bytes':sum(p.stat().st_size for p in inputs),'no_new_column_values':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',action='store_true');main(p.parse_args())
