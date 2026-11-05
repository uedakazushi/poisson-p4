#!/usr/bin/env python3
"""SOURCE only: does the unpaired degree2 output already lie in M1*S1?"""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
import ambient_B18_C_Ksquare_leading_pair_source_20261105_v63 as L
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KSQUARE_LEADING_LINEAR_SPAN_SOURCE_20261105_v63'
PATHS={
 'leading':'reports/AMBIENT_B18_C_KSQUARE_LEADING_PAIR_SOURCE_20261105_v63/unpaired-leading-three-pair-extraction-tensor-and-rank2-sufficient-criterion.json',
 'leading_report':'reports/AMBIENT_B18_C_KSQUARE_LEADING_PAIR_SOURCE_20261105_v63/report.json',
 'quadbasis':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'point':L.PATHS['point'],'Pi':L.PATHS['Pi'],
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={'leading':['Wq44_sparse','Wbeta44_sparse','linear_B_h0','fraction_majorant_before_arithmetic'],
 'quadbasis':['exact_quotient_duals_rows_sparse'],'point':['source_m0_24_by44'],'Pi':['constant_minimal_projection5_by128_sparse']}
 d={k:S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 m=[[Q(c)for c in row]for row in d['point']['source_m0_24_by44']]
 P=[dict((i,Q(c))for i,c in row)for row in d['Pi']['constant_minimal_projection5_by128_sparse']]
 A=[dict((i,Q(c))for i,c in row)for row in d['quadbasis']['exact_quotient_duals_rows_sparse']]
 Wq=dict((i,Q(c))for i,c in d['leading']['Wq44_sparse']);Wb=dict((i,Q(c))for i,c in d['leading']['Wbeta44_sparse']);beta=dict((tuple(e),Q(c))for e,c in d['leading']['linear_B_h0'])
 columns=[];rawhash=hashlib.sha256();all_honest_pure_r=[]
 for mi,e in enumerate(L.M2):
  for qvar in range(9):
   raw={}
   if qvar<8:
    for j,c in Wq.items():L.put(raw,(j,L.shift(e,qvar)),c)
   else:
    for f,b in beta.items():
     ex=tuple(x+y for x,y in zip(e,f))
     for j,c in Wb.items():L.put(raw,(j,ex),L.mul(b,c))
   vv={}
   for (j,ex),c in raw.items():
    for i in range(24):
     if m[i][j]:L.put(vv,(i,ex),L.mul(m[i][j],c))
   hh=L.H3(vv);mini={}
   for (j,ex),c in hh.items():
    for target in range(5):
     if 44+j in P[target]:L.put(mini,36*target+L.I2[ex],L.mul(P[target][44+j],c))
   vals=[Q(0)]*6
   for i,c in mini.items():
    for a in range(6):
     if i in A[a]:vals[a]=L.add(vals[a],L.mul(A[a][i],c))
   columns.append({'monomial_index36':mi,'formal_q_or_B':qvar,'quotient6':[str(x)for x in vals]})
   rawhash.update((json.dumps([[i,str(x)]for i,x in sorted(mini.items())],separators=(',',':'))+'\n').encode())
   all_honest_pure_r.append([mi,qvar,[[j,str(c)]for (j,ex),c in sorted(hh.items())if ex==L.mon([7,7])]])
 allzero=all(all(x=='0'for x in c['quotient6'])for c in columns)
 evidence={'scope':'No LU/LS/q substitution and no signed product. Exact universal source test of membership in known174 first-linear multiples, BEFORE actual3quadratic pivots.',
 'input_sha256':pins,'all324_quotient6_columns':columns,'all324_quotient6_zero':allzero,
 'raw_minimal_degree2_stream_sha256':rawhash.hexdigest(),'all324_H_degree2_pure_r_squared_sparse':all_honest_pure_r,
 'meaning_if_zero':'All actual leading degree2 pairs lie already in M1*S1, so their coherent cubic clearing uses M2 of the35 first-relation lifts; no M3 of the3 quadratic pivots is needed solely to remove degree2.',
 'meaning_if_nonzero':'Must retain cubic coefficients of actualquadratic lifts when clearing thedegree2 output. Zero in finalquotient3 alone doesnot justify dropping thosecoefficients.',
 'next_r_pair_dual_SOURCE_only':'A whole pair liesinm3 afterclearing. Itsr-multiple ism4 and the quarticdual reads head0*r4. Its potentiallyneededpure_r cubic head0 terms come from Pi1_r*H2_rr and theactualM2R35 clearing. No output pair scalaror rank computed.',
 'bounds':d['leading']['fraction_majorant_before_arithmetic'],'source_max_stored_H':L.STORED,'source_max_pre_H':L.PRE,'selective_receipts':receipts}
 owned=S.count(d)+S.count(columns)+S.count(all_honest_pure_r)+S.count(evidence)+25000;assert owned<150000;evidence['whole_SOURCE_owned_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'unpaired-known174-membership-and-pure-r-H2-source.json';p.write_text(json.dumps(evidence,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'UNPAIRED_SOURCE_174_MEMBERSHIP_TEST','all324_quotient6_zero':allzero,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'evidence_sha256':{p.name:sha(p)},'whole_SOURCE_owned_upper':owned,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
