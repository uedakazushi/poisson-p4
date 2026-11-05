#!/usr/bin/env python3
"""Emit an already-computed pure-ell restriction; never solve source or gauge."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63'
SOURCE='reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
HEIGHT_STOP_HANDLER=None
def halt(reason):
 if HEIGHT_STOP_HANDLER is not None:HEIGHT_STOP_HANDLER(reason)
 raise RuntimeError('OPERATIONAL_HEIGHT_STOP_before_Fraction_operation:'+str(reason))
def bits(x):return max(abs(x.numerator).bit_length(),x.denominator.bit_length())
def qconst(n,d=1):
 if d==0 or max(abs(n).bit_length(),abs(d).bit_length())>2048:halt(['input_ratio',abs(n).bit_length(),abs(d).bit_length()])
 z=Q(n,d)
 if bits(z)>2048:halt(['input_reduced_ratio',bits(z)])
 return z
def qparse(x):
 if isinstance(x,Q):
  if bits(x)>2048:halt(['provided_fraction',bits(x)])
  return x
 z=str(x).split('/')
 return qconst(int(z[0]),int(z[1])if len(z)==2 else 1)
def qadd(x,y):
 if not x:return y
 if not y:return x
 nx,dx=abs(x.numerator).bit_length(),x.denominator.bit_length();ny,dy=abs(y.numerator).bit_length(),y.denominator.bit_length()
 bound=max(max(nx+dy,ny+dx)+1,dx+dy)
 if bound>2048:halt(['add_bound',bound])
 z=x+y;assert bits(z)<=2048;return z
def qmul(x,y):
 if not x or not y:return qconst(0)
 nx,dx=abs(x.numerator).bit_length(),x.denominator.bit_length();ny,dy=abs(y.numerator).bit_length(),y.denominator.bit_length()
 bound=max(nx+ny,dx+dy)
 if bound>2048:halt(['multiply_bound',bound])
 z=x*y;assert bits(z)<=2048;return z
def qneg(x):assert bits(x)<=2048;return -x
def plus(a,b):return tuple(qadd(x,y)for x,y in zip(a,b))
def neg(a):return tuple(qneg(x)for x in a)
def times(a,b):
 ans=[]
 for i in range(3):
  z=qconst(0)
  for j in range(i+1):z=qadd(z,qmul(a[j],b[i-j]))
  ans.append(z)
 return tuple(ans)
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return list(map(enc,v))
 return v
def count(v):
 if isinstance(v,(int,Q)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(map(count,v.values()))
 if isinstance(v,(tuple,list)):return sum(map(count,v))
 return 0
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source-peer-report',required=True);a=ap.parse_args()
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 paths=[SOURCE+'/report.json',SOURCE+'/compiled-original164-DAG-and-provided-pure-ell-inputs.json',SOURCE+'/pure-ell-readback-custody-height-and-whole-budget.json',a.source_peer_report]
 pins={p:sha(ROOT/p)for p in paths};report=json.loads((ROOT/paths[0]).read_text());D=json.loads((ROOT/paths[1]).read_text());gate=json.loads((ROOT/paths[2]).read_text());peer=json.loads((ROOT/a.source_peer_report).read_text())
 assert report['all_passed']and peer['all_passed']and peer['input_sha256'][paths[0]]==pins[paths[0]]
 assert report['evidence_sha256'][Path(paths[1]).name]==pins[paths[1]]and report['evidence_sha256'][Path(paths[2]).name]==pins[paths[2]]
 assert gate['future_static_cap']==130000 and gate['source_bit_bound']==2048
 global HEIGHT_STOP_HANDLER
 records=D['operator_records_originalID_op_args'];vals=[None]*(max(r[0]for r in records)+1);inversechecks=[];maxbits=0
 def operational_checkpoint(reason):
  partial=[[j,vals[j]]for j in D['requested_root_nodes164']if vals[j]is not None]
  (out/'root-operational-height-stop.json').write_text(json.dumps(enc({'reason':reason,'computed_requested_roots_only':partial,'stop_before_offending_Fraction_operation':True,'no_cap_raise':True}),indent=2)+'\n')
 HEIGHT_STOP_HANDLER=operational_checkpoint
 phi={i:tuple(map(qparse,v))for i,v in D['provided105_originalID_pure_ell_coefficients']}
 for r in records:
  i,op,*args=r
  if op=='input':v=phi[i]
  elif op=='c':v=(qconst(args[0],args[1]),qconst(0),qconst(0))
  elif op=='a':v=plus(vals[args[0]],vals[args[1]])
  elif op=='n':v=neg(vals[args[0]])
  elif op=='m':v=times(vals[args[0]],vals[args[1]])
  else:
   assert op=='i';x=vals[args[0]];assert x[0]==1
   v=(qconst(1),qneg(x[1]),qadd(qmul(x[1],x[1]),qneg(x[2])));assert times(x,v)==(qconst(1),qconst(0),qconst(0));inversechecks.append([i,args[0],'1'])
  assert all(max(abs(x.numerator).bit_length(),x.denominator.bit_length())<=2048 for x in v)
  maxbits=max(maxbits,*(max(abs(x.numerator).bit_length(),x.denominator.bit_length())for x in v));vals[i]=v
 assert inversechecks==D['unit_inverse_centers6']
 targets={p:(Q(0),Q(1,5),Q(0))if j==8 else(Q(0),Q(0),Q(0))for j,p in enumerate(D['CA_target_positions'])}
 for j,i in enumerate(D['equation105_order']):assert vals[i]==targets.get(j,(Q(0),Q(0),Q(0))),(j,i,vals[i])
 for i,v in D['expected164_constant_and_ell_coefficients_from_C2']:assert vals[i][:2]==tuple(map(qparse,v)),i
 nus=[vals[i]for i in D['nu_roots_CA_CB']];assert nus==[tuple(map(qparse,v))for v in D['expected_nu_pure_ell3_from_Q20']]
 if 2*maxbits+2>2048:
  checkpoint={'all164_root_values3':[[i,vals[i]]for i in D['requested_root_nodes164']],'all105_equations_and_prior_overlaps_passed':True,'stop_before_physical_delta_subtractions':True,'actual_root_bit_bound':maxbits,'required_subtraction_bound':2*maxbits+2}
  (out/'root-readback-height-stop.json').write_text(json.dumps(enc(checkpoint),indent=2)+'\n')
  raise RuntimeError('SOURCE_HEIGHT_STOP_before_delta_subtractions_no_cap_raise')
 physical=[]
 for r in D['physical12_coefficient_roots']:
  delta=plus(vals[r['aA_node']],neg(vals[r['aB_node']]))
  assert all(max(abs(x.numerator).bit_length(),x.denominator.bit_length())<=2048 for x in delta)
  physical.append({'physical_row':r['physical_row'],'source_roots':r,'delta_jet3':delta,'bA_jet3':vals[r['bA_node']],'bB_jet3':vals[r['bB_node']]})
 roots=[[i,vals[i]]for i in D['requested_root_nodes164']]
 data={'parameter':'zA0..6=0,rA=0,sA=ell/5; ordinarycoeffs inQQ[ell]/ell³','all164_actual_pure_ell_root_values3':roots,'all12_actual_physical_coefficients_pure_ell3':physical,'actual_nuA_B_pure_ell3':nus,'all105_original_equations_verified':True,'all164_constant_and_first_C2_overlaps_verified':True,'all6_actual_inverse_products_verified':True,'both_Q20_nu_overlaps_verified':True,'source_phi105_never_changed_or_solved':True,'no_IFT_gauge_higherjet_or_C3_rank_exponent_result':True,'maximum_actual_coefficient_bits':maxbits}
 owned=count([report,D,gate,peer,vals,phi,data]);assert owned<=130000,owned
 ev=out/'all164-cached-pure-ell-readbacks-and-physical12-second-values.json';blob=json.dumps(enc(data),indent=2)+'\n';assert sum((ROOT/p).stat().st_size for p in paths)+len(blob.encode())<2<<20;ev.write_text(blob)
 result={'id':STEM,'all_passed':True,'status':'PROVED_existing164_pure_ell_readback_values_only','physical_source_gauge_or_higherjet_solves':0,'C3_rank_torsion_exponents_decided':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'owned_numeric_slots':owned,'max_actual_bits':maxbits,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));signal.alarm(0)
if __name__=='__main__':main()
