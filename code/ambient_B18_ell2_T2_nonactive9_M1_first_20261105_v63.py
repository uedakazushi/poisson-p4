#!/usr/bin/env python3
"""Prepared known-first fold and four-invariant representation test only."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_NONACTIVE9_M1_FIRST_20261105_v63'
def main(out,gate,peer):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/gate;g=json.loads(gp.read_text());pp=ROOT/peer;ppeer=json.loads(pp.read_text());assert g['all_passed'] and ppeer['all_passed'] and H.sha(gp) in json.dumps(ppeer)
 ip=gp.parent/'pruned-known-firsts-and-positive-M1-four-invariant-test-input.json';assert H.sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text());assert d['no_signed_weighted_first_or_M1_operator_or_invariant_value_evaluated']
 for path,h in g['input_sha256'].items():assert H.sha(ROOT/path)==h
 assert count([g,ppeer])<600;H.dump(out/'phase.json',{'phase':'authorized_known_first_fold_no_M1_function_values','source_gate':H.sha(gp),'source_peer':H.sha(pp),'future_numeric_bound':g['future_whole_numeric_bound']})
 cap=g['positive_prior_bits'];assert cap<=384 and 3*cap+1<=1153;maximum=0
 def ck(x):
  nonlocal maximum
  maximum=max(maximum,abs(x.numerator).bit_length(),x.denominator.bit_length());assert maximum<=cap;return x
 def add(v,w,c=F(1)):
  q=dict(v)
  for k,x in w.items():
   y=ck(q.get(k,F(0))+ck(c*x))
   if y:q[k]=y
   else:q.pop(k,None)
  return q
 def image(A,v):
  q={}
  for k,x in v.items():q=add(q,A[k],x)
  return q
 C={j:H.parse(v) for j,v in d['selected99_original768']};Y={j:H.parse(v) for j,v in d['selected99_Y0']};R={j:H.parse(v) for j,v in d['original_r_first_columns']};E={j:H.parse(v) for j,v in d['original640_ell_first_columns']};L=list(map(H.parse,d['native_L0_columns']));Se={j:H.parse(v) for j,v in d['native_Se_columns']};S0=list(map(H.parse,d['native_S0_columns']));inv=list(map(H.parse,d['fixed_B30_inverse_columns']));Xp=list(map(H.parse,d['point_Xp27_columns']));free=d['free30_V_indices'];M=list(map(H.parse,d['M0_columns']));G2=list(map(H.parse,d['G2linear_columns']));lam=list(map(H.parse,d['lambda_two']));W=H.parse(d['common_Wv56']);p=d['fixed_extraction_lambda_row'];s=d['fixed_extraction_source_row'];assert s==29 and s!=23 and lam[0][p]*W[s]
 for j,v in Se.items():E[640+j]=image(L,v)
 Z={};sourcechecks=0
 for j,c in C.items():
  assert not image(S0,{k-640:x for k,x in c.items() if k>=640});z=[]
  for var,raw in [(7,R),(8,E)]:
   v=image(raw,c);a=image(inv,{i:v[row] for i,row in enumerate(free) if row in v});a={i:x for i,x in a.items() if i<27};q=add(v,image(Xp,a),F(-1));rhs=[{i//9:x for i,x in col.items() if i%9==var} for col in G2];assert image(M,q)==image(rhs,Y[j]);sourcechecks+=26;z.append(q)
  Z[j]=z
 weighted=[];results=[]
 for n,r in enumerate(d['formal_eighteen_weights']):
  ww=H.parse(r['weights']);v=[[{},{}],[{},{}]]
  for k,x in ww.items():
   j,ch=divmod(k,2)
   for direction in range(2):v[ch][direction]=add(v[ch][direction],Z[j][direction],x)
  weighted.append(v);support=d['positive_weighted_vector_supports'][n];assert all(set(v[ch][direction])<=set(support[ch][direction]) for ch in range(2) for direction in range(2))
  op={}
  for ch,row in enumerate(lam):
   for target,x in row.items():
    for direction,vec in [(0,v[ch][1]),(1,v[ch][0])]:
     for source,y in vec.items():op=add(op,{direction*26*56+target*56+source:y},x)
  get=lambda direction,t,i:op.get(direction*26*56+t*56+i,F(0))
  bW=ck(get(1,p,s)/(lam[0][p]*W[s]));be=ck(get(1,p,23)/lam[0][p]-ck(bW*W[23]));lp1=lam[1].get(p,F(0));rrs=ck(get(0,p,s)+ck(lp1*bW*W[s]));aW=ck(rrs/(lam[0][p]*W[s]));rr23=ck(get(0,p,23)+ck(lp1*ck(bW*W[23]+be)));ae=ck(rr23/lam[0][p]-ck(aW*W[23]));cofs=[aW,bW,ae,be]
  for c,(ds,us) in zip(cofs,d['positive_four_coefficient_bounds'][n]):
   DD,UU=int(ds),int(us);assert DD%c.denominator==0 and abs(c.numerator)*(DD//c.denominator)<=UU
  expected={}
  for target,x in lam[0].items():
   for source,w in W.items():expected=add(expected,{target*56+source:ck(x*w)},aW);expected=add(expected,{26*56+target*56+source:ck(x*w)},bW)
   expected=add(expected,{target*56+23:x},ae);expected=add(expected,{26*56+target*56+23:x},be)
  for target,x in lam[1].items():
   for source,w in W.items():expected=add(expected,{target*56+source:ck(-x*w)},bW)
   expected=add(expected,{target*56+23:-x},be)
  residue=add(op,expected,F(-1));good=not residue
  # Check every literal M1 coefficient slot, including all zero slots. Only a
  # single obstruction is serialized on failure; vectors reproduce the rest.
  assert sum(op.get(k,F(0))!=expected.get(k,F(0)) for k in range(2*26*56))==len(residue)
  obstacle=None
  if residue:
   k=min(residue);direction,rem=divmod(k,26*56);target,source=divmod(rem,56);obstacle={'direction':direction,'target26':target,'source56':source,'residual':str(residue[k])}
  results.append({'original418_source':r['original418_source'],'contraction':r['contraction'],'represented_by_four_known_invariant_functionals':good,'proposed_four_coefficients':list(map(str,cofs)),'all2912_slots_checked':True,'nonzero_residual_slot_count':len(residue),'first_exact_obstruction':obstacle,'weighted_Xeff_two_channels_r_ell':[[srec(q) for q in ch] for ch in v]})
  assert count([d,Z,weighted,op,expected,residue,results])+6000<g['future_whole_numeric_bound']
 # Full56 plane W2 verification uses only already computed first vectors;
 # coefficients also follow from the separately source-bound split theorem.
 by={(r['original418_source'],r['contraction']):v for r,v in zip(d['formal_eighteen_weights'],weighted)}
 for source in {r['original418_source'] for r in d['formal_eighteen_weights']}:
  A=by[source,'A'];B=by[source,'B'];assert not A[0][0] and not A[1][0] and not A[1][1] and not add(A[0][1],B[0][0]);assert B[1]==A[0]
 ep=out/'actual-known-first-weighted-M1-operators-and-four-invariant-tests.json';H.dump(ep,{'all18_results':results,'all34_full_source_chain_equations_checked':sourcechecks,'full56_W2_plane_zero_for_all9':True,'four_basis_functionals':d['four_basis_definitions'],'no_M1_Phi_native_S2_pivot_G2quadratic_or_actual_invariant_values_queried':True,'scope':'Known-first coefficient identities only. A failedspankeepsanexactnewsourceoperator; a successpermitsusingexistingfourinvariantvaluesonlyinaseparatecompletevaluegate. Actualpsi/quadraticclasses/HF2/T2lengthnotcomputed.'})
 successes=sum(r['represented_by_four_known_invariant_functionals'] for r in results);assert sum(q.stat().st_size for q in out.iterdir())<g['future_RAW_bound']
 report={'id':STEM,'all_passed':True,'known_M1_operator_representations_succeeded':successes,'known_M1_operator_representations_failed':18-successes,'all2912_times18_operator_slots_checked':True,'full56_W2_plane_zero_all9':True,'actual_M1_psi_or_quadratic_values_not_evaluated':True,'input_sha256':{gate:H.sha(gp),peer:H.sha(pp),str(ip.relative_to(ROOT)):H.sha(ip)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:H.sha(ep)},'source_committed_bound':cap,'source_temporary_bound':3*cap+1,'actual_maximum_bits':maximum,'future_whole_numeric_bound':g['future_whole_numeric_bound'],'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','known_M1_operator_representations_succeeded','known_M1_operator_representations_failed','actual_maximum_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);p.add_argument('--source-peer',required=True);a=p.parse_args();main(a.output.resolve(),a.gate,a.source_peer)
