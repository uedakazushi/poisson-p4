#!/usr/bin/env python3
"""SOURCE: exact tenth fixed row58/col56 Schur operator, no matrix coefficients."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,hashlib,json,resource,signal,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_TENTH_UNIT_SOURCE_20261105_v64'
NINTH='S61_SPARSE_GENERIC5_CONSTRAINED83_COL82_NINTH_UNIT_SOURCE_20261105_v64'
VALUE='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64'
STRIP='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_EXISTING_QRE_STRIP_20261105_v64'
ROOTVALUE='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_ROOT_VALUE_REVIEW_20261105_v64'
QQCHECK='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_IDENTITY_REVIEW_20261105_v64'
H={NINTH:'7dbaf504004a61b9c2554b6c9c31bfffaf426bb324d15dca31221d365c009f95',VALUE:'20524cd69e98bc7f5334a07544035fc04ea4fc1486ac64ce885670c0ce70b173',STRIP:'d1f1940f49703f4c7d956b6ead3c803b25de82aa4e628d6a0a74725739462fb3',ROOTVALUE:'fff391b83d5f1ff72549982310a7d37a04a875c4d4c00936f1d6c639a79bf4f0',QQCHECK:'81c59418bd1a4409949719fdb23ba0f6d13731d4388851434f6b97a0e37d8f38'}
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(a):
 assert a.freeze
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 signal.signal(signal.SIGALRM,lambda s,f:(_ for _ in()).throw(TimeoutError('SOURCE60s')));signal.alarm(60)
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={}
 def add(p):
  p=Path(p);inputs[str(p.relative_to(ROOT))]=sha(p);return p
 def gate(stem):
  p=add(ROOT/'reports'/stem/'report.json');assert inputs[str(p.relative_to(ROOT))]==H[stem];v=json.loads(p.read_text());assert v['all_passed'];return p,v
 def evidence(p,r,name):
  f=add(p.parent/name);assert sha(f)==r['evidence_sha256'][name];return json.loads(f.read_text())
 try:
  np,nr=gate(NINTH);scope=evidence(np,nr,'actual-row35-ninth-unit-and-full168x15-operator.json')
  vp,vr=gate(VALUE);api=evidence(vp,vr,'actual15-col56-and-once-shared-row35.json')
  sp,sr=gate(STRIP);result=evidence(sp,sr,'fixed-row58-QRE-strip-and-Gauss-factorization.json')
  gp,gr=gate(ROOTVALUE);qp,qr=gate(QQCHECK)
  assert sr['only_QRE_powers']=={'Q':4,'R':2,'E':1}and qr['powers']=={'Q':4,'R':2,'E':1}
  rec=result['final_reduced_integer_polynomial'];fp=add(sp.parent/rec['file']);assert sha(fp)==rec['sha256']=='526ec2470b0ceb344ac7a837e9e45e3bfbf99e60e0a646831c5a2ae6ad917a68'
  data=fp.read_bytes();assert data[:9]==b'Q5CELL1\n\x01';p=decode(data[9:]);assert p['content']=='1'and p['shift']==[0]*5 and len(p['terms'])==317
  g=0
  for t in p['terms']:g=gcd(g,abs(int(t[0])))
  assert g==1;w=p['terms'][0];witness={'coefficient':str(w[0]),'exponent5':w[1:]}
  F10={'file':str(fp.relative_to(ROOT)),'sha256':sha(fp),'terms':317,'primitive_integer_content':1,'value_sha256':rec['value_sha256'],'nonzero_witness':witness,'no_irreducibility_or_F9_coprimality_claim':True}
  byrow={r['row']:r for r in api['actual15_ninth_col56_numerators']}
  assert scope['selected15_row_ids']==[6,79,4,46,58,77,17,124,41,78,109,5,48,26,133]and set(byrow)==set(scope['selected15_row_ids'])
  rows=scope['full168_row_ids'];top=scope['top102_row_ids'];sel=scope['selected15_row_ids'];cols=scope['residual15_column_ids'];extra=scope['all87_other_top_RHS'];zeros=scope['all66_zero_constraints']
  assert len(rows)==168 and len(top)==102 and len(sel)==len(cols)==15 and 56 in cols and 58 in top
  rows2=[i for i in rows if i!=58];top2=[i for i in top if i!=58];sel2=[i for i in sel if i!=58];cols2=[j for j in cols if j!=56]
  assert len(rows2)==167 and len(top2)==101 and len(sel2)==len(cols2)==14
  assert set(extra)==set(top2)-set(sel2)and len(extra)==87 and len(zeros)==66 and set(top2).isdisjoint(zeros)and set(top2)|set(zeros)==set(rows2)
  actual=[]
  for i in sel2:
   r=byrow[i];f=add(vp.parent/r['file']);assert sha(f)==r['sha256']
   actual.append({'row':i,'actual_N9_i56':str(f.relative_to(ROOT)),'sha256':sha(f),'terms':r['terms'],'typed_lambda':'353894400000*H*t0*N9_i56/(Q^4*R^2*E*F10)','no_RF_division_or_new_value_evaluated':True})
  pivot=byrow[58];pf=add(vp.parent/pivot['file']);assert sha(pf)==pivot['sha256'];assert pivot['terms']==951
  # A universal Laurent-polynomial verification, independent of actual entries.
  def mono(i):return{tuple(1 if j==i else 0 for j in range(6)):Q(1)}
  def summ(A,B,sign=1):
   C=dict(A)
   for e,c in B.items():
    C[e]=C.get(e,Q(0))+sign*c
    if not C[e]:del C[e]
   return C
  def times(A,B):
   C={}
   for e,c in A.items():
    for f,d in B.items():
     k=tuple(e[j]+f[j]for j in range(6));C[k]=C.get(k,Q(0))+c*d
   return{k:v for k,v in C.items()if v}
  p0,ai,pj,dij,hi,hp=map(mono,range(6));pinv={(-1,0,0,0,0,0):Q(1)}
  S=summ(dij,times(times(ai,pj),pinv),-1);hr=summ(hi,times(times(ai,hp),pinv),-1)
  assert times(p0,S)==summ(times(p0,dij),times(ai,pj),-1)
  # Substitute x56=(hp-pj*xj)/p. Constant and xj coefficients
  # of every transformed row are hr and S respectively.
  assert summ(hi,times(times(ai,hp),pinv),-1)==hr
  assert summ(dij,times(times(ai,pj),pinv),-1)==S
  op={
   'all_passed':True,'source_only':True,'pivot_row':58,'pivot_column':56,'actual_F10':F10,
   'actual_p10':'F10/(353894400000*H*t0*F9)',
   'typed_p10_inverse':'353894400000*H*t0*F9/F10',
   'additional_chart_factor':'F10 (the specified nonzero primitive317-term polynomial); no further factor asserted',
   'old_parameter_units':['old monomial units in H,L,s0,q0,t0','P290 (applied ONCE to the same four original response columns)','detA*detB','Q=5*q0+1','R=q0+1','E=32*s0+15*q0+35','F9'],
   'all_old_and_new_complements_OPEN':True,'F9_still_denominator_not_removed':True,
   'full167_row_ids':rows2,'top101_row_ids':top2,'selected14_row_ids':sel2,'residual14_column_ids':cols2,
   'all87_other_top_RHS':extra,'all66_zero_constraints':zeros,
   'full167x14_RF_roots':{'source':str((np.parent/'actual-row35-ninth-unit-and-full168x15-operator.json').relative_to(ROOT)),'schema':'S10[i,j]=S9[i,j]-lambda[i]*S9[58,j] for EVERY full167 row and remaining14 column','lambda':'S9[i,56]/p10','whole_pivot_row':'S9[58,j] is the COMPLETE current ninth-stage callback, not a raw source row, row35, or zero-filled value'},
   'actual14_selected_lambda':actual,
   'other87_lambda_UNEVALUATED':'Use the actual current S9[i,56]/p10 callback; no missing coefficient or zero is invented.',
   'zero66_lambda':'lambda=0 by the retained literal zero coefficient rows; all66 RHS consistency conditions remain.',
   'all167_residual_RHS':'h10[i]=h9[i]-lambda[i]*h9[58], including all87 extra top RHS and all66 consistency RHS',
   'incoming168_RHS_kept_for_restoration':True,
   'inverse_substitution':'x56=(h9[58]-sum_j S9[58,j]*xj)/p10, then restore previous9 typed units, oneP290,10 constraint units and49 old units',
   'full_cokernel':'full168x15 is unit-equivalent to [p10] direct_sum full167x14; sincep10unit their cokernels are isomorphic. full167coker=top101coker direct_sum K5^66. Selected14 square is NOT samecokernel.',
   'common_representation_ONLY':{'U':'Q^4*R^2*E','D9':'U*F9','D10':'U^2*F9*F10','N10':'U*F10*N9[i,j]-353894400000*H*t0*N9[i,56]*N9[58,j]','justification':'S9=N9/D9 and p10=F10/(353894400000*H*t0*F9); no additional cancellation or P290 operation.'},
   'all_other14_column_coefficients_UNEVALUATED':True,'no_new_matrix_cell_or_rank_inverse_exhaustion_claim':True,
   'P290_applied_once_no_double_scaling':True
  }
  ep=out/'actual-F10-pivot-and-full167x14-typed-operator.json';ep.write_bytes(enc(op))
  add(Path(__file__));add(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py');add(ROOT/f'research/notes/{ID}.md')
  numeric_upper=25000+8*317+8*951+12*(167+14+87+66)+2048
  assert numeric_upper==42056 and numeric_upper<150000
  raw=sum((ROOT/f).stat().st_size for f in inputs)+sum(f.stat().st_size for f in out.iterdir())+65536;assert raw<32*1024*1024
  report={'id':ID,'all_passed':True,'source_only':True,'fixed_actual_nonzero_pivot':True,'new_full_shape':[167,14],'top_shape':[101,14],'selected_square_shape':[14,14],'extra_top_RHS':87,'retained_zero_constraint_RHS':66,'all14_selected_lambda_source_binary_refs_bound':True,'all_full_RHS_and_inverse_substitution_kept':True,'no_new_matrix_coefficient_evaluated':True,'F9_denominator_and_all_complements_retained':True,'P290_once':True,'universal_Schur_and_RHS_identities_checked':True,'source_numeric_upper':numeric_upper,'fresh_READ_output_RAW_upper':raw,'input_sha256':inputs,'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
  (out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  (out/'failure.json').write_bytes(enc({'all_passed':False,'source_only':True,'error':repr(exc),'no_new_matrix_values':True}));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())
