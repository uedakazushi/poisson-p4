#!/usr/bin/env python3
"""SOURCE-only fixed35/58 shortcut test on already accepted input cells."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,hashlib,time
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_FIXED35_58_SHORTCUT_SOURCE_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def poly(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']}
def lit(z):return{tuple(e[j]+z['transform_monomial_exponent5'][j]for j in range(5)):q*Q(z['QQ_transform_scalar'])for e,q in poly(z['original_literal']).items()}
def main(a):
 assert a.freeze
 begin=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')}
 old=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64';rp=old/'report.json';assert sha(rp)=='f7f9bc898a54b9a457a59ba4146f38b0e5a8ce9f2cbda906379a71db20e5f9cb';r=read(rp);assert r['all_passed'];inputs[rp]=sha(rp);ap=old/'actual-cleared-col82-and-commonQ3RE-representation.json';assert sha(ap)==r['evidence_sha256'][ap.name];inputs[ap]=sha(ap);api={z['row']:z for z in read(ap)['actual16_col82_numerator_values']};data={}
 for row in(35,58):
  f=old/api[row]['file'];assert sha(f)==api[row]['sha256'];inputs[f]=sha(f);p=decode(f.read_bytes()[9:]);data[row]=poly(p);assert len(data[row])==111
 A=data[35];B=data[58];assert set(A)==set(B);e=max(A);k=B[e]/A[e];assert k==-1;f=next(f for f in sorted(A)if B[f]!=k*A[f]);assert A!=B
 src=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64';sr=read(src/'report.json');assert sha(src/'report.json')=='e332bea912ec182d8c6f7ed614a1cb7a25cc0dda4ab13f5d6033c85c9200b5ab';inputs[src/'report.json']=sha(src/'report.json')
 vals=[]
 for row in(35,58):
  p=src/'rows'/f'row-{row}-source-operands.json';assert sha(p)==sr['evidence_sha256'][str(p.relative_to(src))];inputs[p]=sha(p);vals.append(read(p)['row_operands'])
 prefix=[]
 for name in('U','V','C','F','N','M'):
  X=lit(vals[0][name+'35']);Y=lit(vals[1][name+'58']);prefix.append({'prefix':name,'terms35':len(X),'terms58':len(Y),'equal':X==Y,'equal_to_first_ratio_times35':Y=={e:k*q for e,q in X.items()}})
 assert all(not z['equal']for z in prefix if z['prefix']in('F','N','M'))
 current=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64';cp=current/'report.json';assert sha(cp)=='20524cd69e98bc7f5334a07544035fc04ea4fc1486ac64ce885670c0ce70b173';inputs[cp]=sha(cp);cap=current/'actual15-col56-and-once-shared-row35.json';inputs[cap]=sha(cap);candidate=next(z for z in read(cap)['actual15_ninth_col56_numerators']if z['row']==58);assert candidate['terms']==951
 evidence={'all_passed':True,'source_only':True,'fixed_pair':[35,58],'equal_old_term_counts_DO_NOT_imply_equal_values':True,'old_pivot_column82_supports_equal111':True,'old_pivot_values_equal':False,'first_constant_ratio':str(k),'first_witness':{'exponent5':list(e),'A35':str(A[e]),'B58':str(B[e])},'rejection_witness':{'exponent5':list(f),'A35':str(A[f]),'B58':str(B[f]),'first_ratio_times_A35':str(k*A[f])},'no_QQ_constant_or_parameter_monomial_proportionality':'Equal finite support forces any monomial shift to0 by coordinate minima; the two coefficient ratios disagree. This is NOT a statement about arbitrary rational K5 multipliers.','cached_prefix_comparisons':prefix,'no_raw_row58_minus_row35_shortcut_adopted':True,'new_known_col56_row58_terms951_metadata_only':candidate['terms'],'no_new_typed_pivot_or_matrix_update':True,'all168x15_other14_87plus66_RHS_oneP290_complements_kept':True}
 (out/'fixed-pair-source-counterwitness.json').write_bytes(enc(evidence));report={'id':ID,'all_passed':True,'source_only':True,'same_or_simple_proportional_shortcut_rejected':True,'new_matrix_coefficients':0,'only_fixed_pair_no_row_or_rank_search':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())
