#!/usr/bin/env python3
"""SOURCE only: fixed row35 numerator stripping by Q, q+1 and E only."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,json,hashlib,time,signal,resource
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL82_ROW35_EXISTING_UNIT_STRIP_SOURCE_V2_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def serial(A):return[[list(e),str(q)]for e,q in sorted(A.items(),key=lambda z:(sum(z[0]),z[0]),reverse=True)]
def primitive(A):
 d=1
 for q in A.values():d=d*q.denominator//gcd(d,q.denominator)
 Z={e:(q*d).numerator for e,q in A.items()};g=0
 for z in Z.values():g=gcd(g,abs(z))
 return{'QQ_content':str(Q(g,d)),'integer_terms':[[list(e),str(z//g)]for e,z in sorted(Z.items(),key=lambda z:(sum(z[0]),z[0]),reverse=True)]}
def mul(A,B):
 R={}
 for e,q in A.items():
  for f,z in B.items():
   g=tuple(x+y for x,y in zip(e,f));R[g]=R.get(g,Q(0))+q*z
 return{e:q for e,q in R.items()if q}
def stats(A):
 d=1
 for q in A.values():d=d*q.denominator//gcd(d,q.denominator)
 return d,max([abs((q*d).numerator)for q in A.values()]+[0])
def main(a):
 assert a.build
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);cp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64/report.json';assert sha(cp)=='f7f9bc898a54b9a457a59ba4146f38b0e5a8ce9f2cbda906379a71db20e5f9cb';cr=read(cp);assert cr['all_passed'];api=cp.parent/'actual-cleared-col82-and-commonQ3RE-representation.json';assert sha(api)==cr['evidence_sha256'][api.name];v=next(z for z in read(api)['actual16_col82_numerator_values']if z['row']==35);file=cp.parent/v['file'];assert sha(file)==v['sha256'];p=decode(file.read_bytes()[9:]);assert len(p['terms'])==111;shift=[min(p['shift'][j]+t[j+1]for t in p['terms'])for j in range(5)];assert shift==[0,0,-1,-1,-2];P={tuple(p['shift'][j]+t[j+1]-shift[j]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']};original=dict(P)
 gatep=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL82_ROW35_EXISTING_UNIT_STRIP_SOURCE_20261105_v64/changed-fibre-support-SOURCE-gate.json';assert sha(gatep)=='45208513c7b602e3c26130830812e65667192d1acac3e4c9617a84229932c353';g=read(gatep);assert g['all_passed']
 groups={tuple(k):d for k,d in g['fixed_H_L_t_fibres']};universe={(h,l,sq,qq,t)for (h,l,t),D in groups.items()for sq in range(D+1)for qq in range(D-sq+1)};assert len(universe)==200 and set(P)<=universe;DM=max(groups.values());tries=DM+3;d0,z0=stats(P);Zsum=sum(abs((q*d0).numerator)for q in P.values());Dglobal=d0*160**(DM*tries);Mglobal=Zsum*(DM+1)**tries*51**((DM+1)*tries)*160**(DM*tries)*82**(DM+1);coeffchars=(Dglobal.bit_length()+Mglobal.bit_length())//3+32;RAWbound=3*200*tries*(5*16+coeffchars+128)+200000;assert RAWbound<32*1024*1024
 factors=[('Q',3,{(0,0,0,1,0):Q(5),(0,0,0,0,0):Q(1)}),('R=q0+1',3,{(0,0,0,1,0):Q(1),(0,0,0,0,0):Q(1)}),('E',2,{(0,0,1,0,0):Q(32),(0,0,0,1,0):Q(15),(0,0,0,0,0):Q(35)})];powers={};traces=[];peak=44625;rawupper=RAWbound;bits=max(Dglobal.bit_length(),Mglobal.bit_length());completed=0
 def bound(A):
  assert set(A)<=universe and all(Dglobal%q.denominator==0 and abs((q*Dglobal).numerator)<=Mglobal for q in A.values())
 for name,j,F in factors:
  lead=tuple(1 if h==j else 0 for h in range(5));lc=F[lead];power=0
  while True:
   degree=max([e[j]for e in P]+[0]);n=len(P);support_upper=len(universe);owned=44625;assert owned<150000;bitupper=bits;assert completed<tries
   R=dict(P);U={};steps=0
   while True:
    eligible=[e for e in R if e[j]>=1]
    if not eligible:break
    e=max(eligible,key=lambda e:(e[j],e));q=R[e]/lc;u=tuple(e[h]-lead[h]for h in range(5));U[u]=U.get(u,Q(0))+q
    for f,c in F.items():
     g=tuple(u[h]+f[h]for h in range(5));R[g]=R.get(g,Q(0))-q*c
     if not R[g]:del R[g]
    steps+=1;assert len(R)+len(U)<=4*support_upper
   bound(R);bound(U);back=mul(F,U)
   for e,q in R.items():back[e]=back.get(e,Q(0))+q
   back={e:q for e,q in back.items()if q};bound(back);assert back==P and all(e[j]==0 for e in R)
   trace={'factor':name,'successful_division':not R,'input':serial(P),'quotient':serial(U),'remainder':serial(R),'fixed_linear_leading_variable':j,'leading_coefficient':str(lc),'division_steps':steps,'prearithmetic_support_upper':support_upper,'prearithmetic_owned_upper':owned,'prearithmetic_bit_upper':bitupper,'literal_remultiply_plus_remainder':True};tf=out/f'linear-division-{completed:02d}.json';tf.write_bytes(enc(trace));traces.append({'factor':name,'successful_division':not R,'file':tf.name,'sha256':sha(tf)});completed+=1;del trace,back
   if R:break
   P=U;power+=1;assert P
  powers[name]=power
 assert P and all(e>=0 for x in P for e in x);rebuilt=dict(P)
 for name,j,F in reversed(factors):
  for k in range(powers[name]):rebuilt=mul(rebuilt,F)
 assert rebuilt==original;pr=primitive(P);assert len(pr['integer_terms'])==len(P)
 evidence={'all_passed':True,'source_only':True,'fixed_existing_cell':[35,82],'actual_input111_terms':111,'actual_Laurent_monomial_shift5':shift,'only_existing_unit_factors':powers,'known_linear_division_traces':traces,'reduced_nonzero_polynomial':pr,'exact_rebuild_of_original111_terms':True,'actual_entry_formula':'content*H^m0*L^m1*s0^m2*q0^m3*t0^m4 * Q^(powerQ-3)*(q0+1)^(powerR-1)*E^(powerE-1)*primitive_reduced_polynomial','new_chart_if_future_pivot':'nonzero reduced polynomial; all previous chart factors and complements retained. No pivot performed.','no_general_polynomial_gcd_factor_rank_or_row_search':True,'changed_fibre_gate_sha256':sha(gatep),'global_denominator_bound':str(Dglobal),'global_scaled_integer_bound':str(Mglobal),'no_new_matrix_cell_or_coefficient_job':True}
 (out/'fixed-row35-existing-unit-divisions.json').write_bytes(enc(evidence));inputs={cp:sha(cp),api:sha(api),file:sha(file),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py'),gatep:sha(gatep)};report={'id':ID,'all_passed':True,'source_only':True,'fixed_cell':[35,82],'only_existing_unit_powers':powers,'reduced_polynomial_terms':len(P),'reduced_polynomial_nonzero':True,'source_owned_proven_upper':peak,'all_linear_division_remultiplications_pass':True,'source_RAW_proven_upper':sum(p.stat().st_size for p in inputs)+rawupper+100000,'prearithmetic_integer_bit_upper':bits,'no_pivot_matrix_or_new_coefficient_executed':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};assert report['source_RAW_proven_upper']<32*1024*1024;(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',action='store_true');main(p.parse_args())
