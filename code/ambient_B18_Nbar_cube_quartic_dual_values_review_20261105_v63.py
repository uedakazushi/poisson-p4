#!/usr/bin/env python3
"""Independent fixed eight relation values, 2minor and actual quartic dual."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63_REVIEW';PINS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):PINS[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return list(map(enc,v))
 return v
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 rp='reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63/report.json';r=read(rp);assert r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 own=read('reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_SOURCE_20261105_v63_REVIEW/report.json');assert own['all_passed']and not own['rank_values_or_quartic_dual_evaluated']
 D=read('reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_SOURCE_PREFLIGHT_20261105_v63/tiny-five-relations-and-constructive-quartic-dual-source.json')
 Z=read('reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_ZERO3_SOURCE_20261105_v63/three-actual-zero-heads-and-corrected-first-source-relations.json')
 v=read('reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63/actual-eight-relations-and-complete-quartic-dual.json')
 A=[[F(c)for c in row]for row in D['source_q3_restricted_to_five_pure_r3_heads']];M=[[F(c)for c in col]for col in D['source_actual_natural_origin_columns5']];p4=[{int(i):F(c)for i,c in col}for col in D['all1650_constructive_quartic_reduction_columns3_sparse']]
 columns=[[sum(A[i][j]*col[j]for j in range(5))for i in range(3)]for col in M]
 ms=list(itertools.combinations_with_replacement(range(8),4));idx={m:i for i,m in enumerate(ms)}
 for col in Z['three_corrected_actual_M1_columns40']:
  value=[F(0)]*3
  for k,c in col.items():
   head,var=divmod(int(k),8);j=head*330+idx[tuple(sorted((var,7,7,7)))];value=[a+F(c)*p4[j].get(i,F(0))for i,a in enumerate(value)]
  columns.append(value)
 assert columns==[[F(c)for c in col]for col in v['all_eight_actual_quartic_relation_vectors3']]
 selected=v['selected_independent_source_columns'];rr=v['fixed2minor_rows'];assert selected==[4,5]and rr==[1,2]
 u,w=[columns[j]for j in selected];det=u[1]*w[2]-u[2]*w[1];assert det==F(1327744825,11226444192)==F(v['fixed2minor_determinant'])
 iv=[[w[2]/det,-w[1]/det],[-u[2]/det,u[1]/det]];mat=[[u[i],w[i]]for i in rr]
 assert all(sum(mat[i][k]*iv[k][j]for k in range(2))==int(i==j)and sum(iv[i][k]*mat[k][j]for k in range(2))==int(i==j)for i in range(2)for j in range(2))
 factors=[]
 for col in columns:
  co=[sum(row[j]*col[rr[j]]for j in range(2))for row in iv];assert [co[0]*u[i]+co[1]*w[i]for i in range(3)]==col;factors.append(co)
 assert factors==[[F(c)for c in col]for col in v['all8_fixed_factorizations']]
 assert all(col[0]==0 for col in columns);dual=[F(1),F(0),F(0)]
 assert dual==list(map(F,v['actual_complete_quartic_dual_on_three_r4_representatives']))
 phi=[[j,p.get(0,F(0))]for j,p in enumerate(p4)if p.get(0,F(0))]
 assert phi==[[int(i),F(c)]for i,c in v['actual_complete_quartic_dual1650_sparse']]
 data={'independent_eight_actual_source_columns':columns,'fixed_minor_rows':rr,'fixed_minor_columns':selected,'fixed_minor':det,'both_2inverse_products':True,'all8_span_factorizations':factors,
 'actual_complete_dual_three_r4':dual,'actual_complete_dual1650_sparse':phi,'actual_dual_source_proof':'Independent SOURCE p4 gives surjection k3→actual gr4. The eight genuine relations have exact rank2 and actual gr4 dimension1. Hence their span is the entire kernel. Dual*p4 is genuine on actual degree4 classes.',
 'newminimal_head0_r4_dual_value':1,'cube_values_computed':False,'no_new_source_gauge_jet_or_ranksearch':True}
 ep=out/'independent-eight-fixed-values-and-genuine-quartic-dual.json';ep.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'status':'PROVED_independent_actual_complete_quartic_dual','input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'actual_relation_rank':2,'actual_dual_three_r4':['1','0','0'],'cube_values_computed':False,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)
if __name__=='__main__':main()
