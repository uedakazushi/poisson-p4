#!/usr/bin/env python3
"""Prepared fixed951-term SOURCE strip: integer pseudo-division by Q/R/E only."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,json,hashlib,time,signal,resource,re,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
from s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64 import write_cell,validate_cell,descending_alias_keys
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_EXISTING_QRE_STRIP_20261105_v64';PREFLIGHT='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_STRIP_PREFLIGHT_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def scalars(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',v)else 0
 if isinstance(v,list):return sum(map(scalars,v))
 if isinstance(v,dict):return sum(map(scalars,v.values()))
 return 0
def main(a):
 assert a.execute and a.root_review_file and a.root_review_sha
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 def timeout(s,f):raise TimeoutError('hard60s; preserve all completed SOURCE traces')
 signal.signal(signal.SIGALRM,timeout);signal.alarm(60);rootreview=ROOT/a.root_review_file;assert sha(rootreview)==a.root_review_sha and rootreview.stat().st_size<=65536;rootdoc=read(rootreview);assert rootdoc['all_passed']and scalars(rootdoc)<=2000
 pp=ROOT/f'reports/{PREFLIGHT}/report.json';pr=read(pp);assert pr['all_passed']and pr['prepared_fixed_QRE_worker_safe'];bp=pp.parent/'integer-pseudo-division-budget-and-proof.json';assert sha(bp)==pr['evidence_sha256'][bp.name];b=read(bp);assert sha(Path(__file__))==b['worker_sha256']
 for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);groups={tuple(k):D for k,D in b['fixed_H_L_t_fibres']};U={(h,l,s,q,t)for(h,l,t),D in groups.items()for s in range(D+1)for q in range(D-s+1)};assert len(U)==1739;Dmax=b['max_s_q_total_degree'];Mdiv=int(b['all_pseudo_maps_integer_l1_bound']);Mrebuild=int(b['all_rebuild_maps_integer_l1_bound']);Mglobal=max(Mdiv,Mrebuild);Midentity=83*Mdiv;META=25000;LIMIT=b['whole_owned_upper'];assert LIMIT<=68784;pool={};peak=META;traceindex=[];completed=0;successes=0;identities=0;maxbits=0
 def own(extra=128):
  nonlocal peak
  n=META+7*sum(len(z)for z in pool.values())+extra;peak=max(peak,n);assert n<=LIMIT,(n,LIMIT);return n
 def meta(decoded=None):
  n=5*len(U)+scalars(rootdoc)+scalars(pr)+scalars(b)+scalars(traceindex)+(scalars(decoded)if decoded is not None else 0)+3000;assert n<=META,n
 def put(k,e,c,bound):
  nonlocal maxbits
  assert e in U and isinstance(c,int)and abs(c)<=bound;A=pool[k]
  if c:
   if e not in A:own(7+128)
   A[e]=c
  else:A.pop(e,None)
  maxbits=max(maxbits,abs(c).bit_length())
 def output(k,name):
  A=pool[k];own(3*len(A)+2048);path=out/name;rec=write_cell(path,Q(1),A,[0]*5);assert validate_cell(path,Q(1),A,rec)and path.stat().st_size<=b['one_integer_cell_RAW_upper'];return{'file':path.name,**rec,'all_coefficients_stream_readback':True}
 def logical_hash(A):
  own(3*len(A)+2048);h=hashlib.sha256();h.update((str(len(A))+'\n').encode());keys=descending_alias_keys(A)
  for e in keys:h.update((str(A[e])+' '+' '.join(map(str,e))+'\n').encode())
  return h.hexdigest()
 try:
  ip=ROOT/b['original_binary_path'];blob=ip.read_bytes();p=decode(blob[9:]);meta(p);qOrigin=Q(p['content']);shift=b['actual_Laurent_monomial_valuation5'];assert str(qOrigin)==b['original_QQ_content'];pool['P']={};G=0
  for t in p['terms']:
   e=tuple(p['shift'][j]+t[j+1]-shift[j]for j in range(5));c=int(t[0]);assert c and e not in pool['P'];put('P',e,c,Mglobal);G=gcd(G,abs(c))
  assert len(pool['P'])==951 and G==1;del p,blob;meta();initial=output('P','initial-integer-polynomial.bin');current_ref={'file':initial['file'],'sha256':initial['sha256'],'encoded_value_integer_divisor':'1','logical_integer_value_sha256':initial['value_sha256']};powers={};factors=[('Q',3,5,{(0,0,0,1,0):5,(0,0,0,0,0):1}),('R',3,1,{(0,0,0,1,0):1,(0,0,0,0,0):1}),('E',2,32,{(0,0,1,0,0):32,(0,0,0,1,0):15,(0,0,0,0,0):35})];factor_records=[]
  for label,j,aLead,F in factors:
   count=0
   while True:
    assert completed<b['max_trials'];meta();P=pool['P'];oldtotal=max(e[2]+e[3]for e in P);assert oldtotal<=Dmax-successes and sum(abs(c)for c in P.values())<=int(b['all_normalized_quotient_l1_bound']);d=max(e[j]for e in P);A=aLead**d;pool['R']={};pool['V']={}
    for e,c in P.items():put('R',e,A*c,Mdiv)
    steps=0;lead=tuple(1 if h==j else 0 for h in range(5))
    while True:
     eligible=(e for e in pool['R']if e[j]>0);e=max(eligible,key=lambda z:(z[j],z),default=None)
     if e is None:break
     c=pool['R'][e];assert c%aLead==0;v=c//aLead;u=tuple(e[h]-lead[h]for h in range(5));assert u not in pool['V'];put('V',u,v,Mdiv)
     for f,coef in F.items():
      g=tuple(u[h]+f[h]for h in range(5));term=coef*v;assert abs(term)<=Midentity;new=pool['R'].get(g,0)-term;assert new+term==pool['R'].get(g,0);put('R',g,new,Mdiv)
     steps+=1;assert steps<=b['support_universe_terms']
    assert all(e[j]==0 for e in pool['R'])
    def check(e):
     nonlocal identities,maxbits
     left=A*P.get(e,0);right=pool['R'].get(e,0)
     for f,coef in F.items():
      u=tuple(e[h]-f[h]for h in range(5));term=coef*pool['V'].get(u,0);right+=term;assert abs(term)<=Midentity and abs(right)<=Midentity
     assert left==right;identities+=1;maxbits=max(maxbits,abs(right).bit_length(),abs(left).bit_length())
    for e in P:check(e)
    for e in pool['R']:check(e)
    for u in pool['V']:
     for f in F:check(tuple(u[h]+f[h]for h in range(5)))
    qrec=output('V',f'trial-{completed:02d}-{label}-scaled-quotient.bin');rrec=output('R',f'trial-{completed:02d}-{label}-scaled-remainder.bin');success=not pool['R'];trace={'factor':label,'trial':completed,'leading_variable':j,'leading_coefficient':aLead,'input_integer_polynomial':current_ref,'degree':d,'pseudo_multiplier':str(A),'identity':'a^d*P=F*V+R','quotient':qrec,'remainder':rrec,'zero_remainder':success,'coefficientwise_identity_checks_total':identities,'long_division_steps':steps,'only_QRE_no_F9_division':True};tf=out/f'trial-{completed:02d}-{label}.json';tf.write_bytes(enc(trace));traceindex.append({'factor':label,'trial':completed,'file':tf.name,'sha256':sha(tf),'zero_remainder':success});completed+=1;del trace
    if not success:
     del pool['R'],pool['V'];break
    # Primitive F and P imply primitive integer quotient by Gauss. This
    # converts the successful pseudo-quotient by a^d and preserves qOrigin.
    assert qrec['G']==A
    for e,c in pool['V'].items():assert c%A==0;pool['V'][e]=c//A
    del pool['P'],pool['R'];pool['P']=pool.pop('V');P=pool['P'];assert P and max(e[2]+e[3]for e in P)==oldtotal-1;g=0
    for c in P.values():g=gcd(g,abs(c))
    assert g==1 and all(min(e[h]for e in P)==0 for h in range(5));successes+=1;count+=1;assert successes<=Dmax;current_ref={'file':qrec['file'],'sha256':qrec['sha256'],'encoded_value_integer_divisor':str(A),'logical_integer_value_sha256':logical_hash(P)};assert time.monotonic()-start<55
   powers[label]=count;factor_records.append({'factor':label,'power':count})
  final=output('P','final-reduced-integer-polynomial.bin');finalhash=final['value_sha256'];del P;pool['B']=dict(pool['P']);own();del pool['P']
  # Exact reverse multiplication, never a decomposition or a gcd search.
  for label,j,aLead,F in reversed(factors):
   for _ in range(powers[label]):
    pool['T']={}
    for e,c in pool['B'].items():
     for f,v in F.items():
      g=tuple(e[h]+f[h]for h in range(5));new=pool['T'].get(g,0)+c*v;put('T',g,new,Mrebuild)
    del pool['B'];pool['B']=pool.pop('T')
  assert logical_hash(pool['B'])==initial['value_sha256'];del pool['B'];assert not pool
  result={'all_passed':True,'only_existing_unit_powers':powers,'original_QQ_content_unchanged_by_Gauss':str(qOrigin),'actual_Laurent_monomial_valuation5':shift,'final_reduced_integer_polynomial':final,'input_and_final_factorization_chain':traceindex,'original_rebuilt_value_sha256':initial['value_sha256'],'final_reduced_nonzero':final['terms']>0,'no_F9_division_or_new_matrix_value':True,'full168x15_87plus66_RHS_other14_oneP290_complements_untouched':True};(out/'fixed-row58-QRE-strip-and-Gauss-factorization.json').write_bytes(enc(result));raw=b['actual_READ_input_RAW_bytes']+pp.stat().st_size+bp.stat().st_size+rootreview.stat().st_size+sum(p.stat().st_size for p in out.iterdir());assert raw+65536<=b['fresh_RAW_upper'];report={'id':ID,'all_passed':True,'source_only':True,'existing_input_terms':951,'only_QRE_powers':powers,'final_integer_terms':final['terms'],'source_trials':completed,'all_coefficientwise_identities':identities,'tracked_map_and_identity_integer_bits_actual':maxbits,'owned_numeric_peak':peak,'owned_limit':LIMIT,'actual_RAW_before_report':raw,'Gauss_primitive_content_checks':True,'original_reconstruction_passed':True,'no_F9_or_general_GCD_or_row_search_or_matrix_update':True,'input_sha256':{**pr['input_sha256'],str(pp.relative_to(ROOT)):sha(pp),str(bp.relative_to(ROOT)):sha(bp),str(rootreview.relative_to(ROOT)):sha(rootreview)},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));assert raw+(out/'report.json').stat().st_size<=b['fresh_RAW_upper'];signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  (out/'failure.json').write_bytes(enc({'all_passed':False,'error':repr(exc),'completed_trials':traceindex,'no_automatic_retry':True,'all_completed_trace_files_preserved':True,'seconds':time.monotonic()-start}));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true');p.add_argument('--root-review-file');p.add_argument('--root-review-sha');main(p.parse_args())
