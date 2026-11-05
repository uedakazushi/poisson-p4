#!/usr/bin/env python3
"""Prepare fixed j3 source cut and same-Haxis ring worker; NO tail values."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,importlib.util,json,resource,signal,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def numeric(x):
 if isinstance(x,(int,float)):return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except ValueError:return 0
 if isinstance(x,list):return sum(numeric(y) for y in x)
 if isinstance(x,dict):return len(x)+sum(numeric(y) for y in x.values())
 return 0
def main(out):
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60);start=time.monotonic()
 out.mkdir(parents=True,exist_ok=False);pins={};runtime={}
 def pin(s,h=None,run=False):
  p=ROOT/s;v=sha(p)
  if h:assert v==h
  pins[s]=v
  if run:runtime[s]=v
  return p
 helper=pin('code/s61_nonsquare_prePOWER_tail1_worker_custody_20261105_v63.py')
 spec=importlib.util.spec_from_file_location('lexical_only_reader',helper);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);Spans=mod.Spans
 base='reports/S61_NONSQUARE_PREPOWER_TAIL8_COMPLEMENT_SOURCE_20261105_v63/'
 pin(base+'report.json','5e72b7f654a2a375edf1e48219343a9b453571b14ab53f72d692a08a9b4e5c48')
 family=Spans(pin(base+'typed-family-and-source-plans.json','d5f447ba15fcacba48ae193c9dc395d596894d553a1ca795f9aac32bcb3e3e98'))
 plan=family.value(family.index(family.obj('source_plans'),2));family.close()
 assert plan['target']==3 and len(plan['nodes'])==162 and plan['source_current15_rows']==[2,4,5,8]
 N=plan['nodes'];cols=plan['column_roots'];ks=plan['cached_kappa_nodes'];assert len(ks)==4
 needed=set();todo=cols[:]
 while todo:
  i=todo.pop()
  if i not in needed:needed.add(i);todo+=N[i]['args']
 assert len(needed)==149
 cp=Spans(pin('reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json'))
 pd=Spans(pin('reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json'))
 q=Spans(pin('reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json'))
 cps=cp.obj('matrix25x32');pds=pd.obj('Lframe_dual6x2');qd=q.obj('scalar_DAG');qn=q.obj('nodes',qd);qa=q.obj('new_unit_atoms_num_den',qd)
 targets=q.value(q.obj('new22_target_original32_columns'));j=targets.index(3)
 remaining=q.obj('remaining15x22_scalar_roots')
 assert [q.value(q.index(remaining,r,j)) for r in [2,4,5,8]]==[2551,2761,2866,3181]
 B=[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
 active=needed|set(ks)
 # The plan's literal one node is unused after simplification. Copy exactly
 # the old149-column closure plus cached4, matching the frozen source count.
 atoms=[i for i,n in enumerate(N) if i in active and n['op']=='external_atom'];witness=[];amin=128;amax=-128;alln=colterms=kterms=0;maxatom=0;prepbits=rombits=0
 native=out/'native-input.txt';literal=out/'literal-source-records.jsonl'
 with native.open('w') as f,literal.open('w') as g:
  f.write('J3HAX001\n'+str(len(atoms))+'\n')
  for i in atoms:
   n=N[i];ref=n['reference'];den=[[[0,0,0],'1']]
   if 'kind' in ref:
    definition=ref['definition'];assert q.value(q.index(qn,ref['old_node']))==definition
    if ref['kind']=='checkpoint':raw,den=cp.value(cp.index(cps,*definition[1:]))
    elif ref['kind']=='accepted_pair_second':raw,den=pd.value(pd.index(pds,definition[1],1))
    else:assert ref['kind']=='new_unit_atom';raw,den=q.value(q.index(qa,definition[1]))
   elif 'cached_kappa_root' in ref:raw=json.loads(pin(ref['path']).read_text())['polynomial']
   elif 'named_FULLU_unit' in ref:assert ref['named_FULLU_unit']=='B';raw=B
   else:raw=[] if ref['QQ']=='0' else [[[0,0,0],ref['QQ']]]
   maxatom=max(maxatom,numeric(raw)+numeric(den));D=int(n['QQ_D']);assert D%32003
   u=Q(n['QQ_divisor']);assert u and u.numerator%32003 and u.denominator%32003
   if 'kind' in ref:
    inv=next(x for x in plan['atom_denominator_inventory'] if x['old_node']==ref['old_node'])
    assert Q(inv['QQ_unit'])==u and inv['a_shift']==n['a_shift']
    if inv['B_power']==0:assert den==[[[-n['a_shift'],0,0],str(u)]]
    else:
     assert inv['B_power']==1
     assert {tuple(m):Q(v) for m,v in den}=={tuple(m):u*Q(v) for m,v in B}
   f.write(f'{i} {D} {len(raw)}\n');g.write(json.dumps(dict(node=i,reference=ref,original_num_den=[raw,den],QQ_divisor=str(u),a_shift=n['a_shift']),separators=(',',':'))+'\n')
   L=0
   for m,v in raw:
    c=Q(v);h=max(abs(c.numerator).bit_length()+u.denominator.bit_length()+D.bit_length(),c.denominator.bit_length()+abs(u.numerator).bit_length())+2
    assert h<=2048;num=c.numerator*u.denominator*D;dd=c.denominator*u.numerator;assert dd and num%dd==0
    z=num//dd;assert z*dd==num;prepbits=max(prepbits,h,abs(z).bit_length()+abs(dd).bit_length());assert prepbits<=2048
    a,e,d=m[0]+n['a_shift'],m[1],m[2];assert -128<=a<=128 and 0<=e<=128 and 0<=d<=128
    amin=min(amin,a);amax=max(amax,a);rombits=max(rombits,abs(z).bit_length());L+=abs(z);f.write(f'{a} {e} {d} {z}\n')
   assert L==int(n['integer_l1_upper']);alln+=len(raw)
   if i in ks:kterms+=len(raw)
   else:colterms+=len(raw)
   witness.append(dict(node=i,terms=len(raw),D=str(D),D_punit=True,QQ_unit_punit=True,integer_normalization_reverse=True))
 assert colterms==9832 and kterms==4920 and alln==14752 and rombits<=384
 spans=2*sum(len(v) for ss in (cp,pd,q) for v in ss.cache.values())
 for ss in (cp,pd,q):ss.close()
 prepared_owned=numeric(plan)+spans+maxatom+10000;assert prepared_owned<150000
 instructions=[[{'external_atom':0,'add_bound_only':1,'sub_bound_only':2,'mul_bound_only':3}[n['op']],*(n['args'] if n['args'] else [-1,-1])] for n in N]
 pin('reports/S61_TWO_TAIL_H_D12_AXIS_QUOTIENT_ENTRY_20261105_v63/report.json','23646500355809a05d19af1774fa33ce75a8ec56f0e31151f52fa56ed76063c8',True)
 pin('reports/S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63/report.json','49990187b3b0b8edaf8939e4d13ac297f4f2ecba7f8b9f5326639089fa043ce2',True)
 pin('reports/S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63_REVIEW/report.json','6640c4f2be0e409a8f4854fd458fbbd700b4e28e7d6a360c153b4f2f67bcd25b',True)
 pin('reports/S61_NONSQUARE_PREPOWER_TAIL1_FULL57_UNION_SOURCE_20261105_v63/report.json','3dbea3f3095fff06ffde0a809dd870d7370940d3567b149eccf2d650934b194c',True)
 pin('reports/S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63_REVIEW/report.json','f946ca6312c0b5167aa5cc257434d9338727a0ff3a6142496bc2ffa55d3890c8',True)
 # All actual new j3 literal values remain UNEVALUATED as a polynomial or ring element.
 worker=ROOT/'code/s61_two_tail_H_tail3_axis_20261105_v63.py';pins[str(worker.relative_to(ROOT))]=sha(worker)
 note=ROOT/'research/notes/S61_TWO_TAIL_H_D12_GAP_THIRD_SOURCE_20261105_v63.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 x=subprocess.run([sys.executable,str(worker),'--help'],capture_output=True,text=True,timeout=10);assert x.returncode==0
 (out/'worker-help.txt').write_text(x.stdout);(out/'worker-help.stderr').write_text(x.stderr)
 powerN=max(amax,1)-min(amin,-1)+1;ringM=sum(op==3 for op,a,b in instructions)
 productupper=4*alln+len(atoms)*35+28*(powerN+ringM+3)+100;assert productupper<100000
 proof=dict(preparation_only=True,no_tail_or_ring_values_executed=True,source_rows=[2,4,5,8],original32_target=3,old_roots=[2551,2761,2866,3181],
  column149=True,total_typed_nodes=162,cached_kappa4_no_inverse=True,literal_column_terms=colterms,cached_kappa_terms=kterms,
  all_literal_denominators_punit=True,atom_normalization=witness,preparation_owned=prepared_owned,preparation_temp_bits=prepbits,existing_ROM_integer_bits=rombits,
  source_global_identity='A_anchor*D15_3=-T3*D12 over COMPLETE U6 only; genuine original15 minor/normalizations/all tails retained',
  no_oldOmega_or_kappa_or_anchor_chart=True,conditional_gap='after independent D12 nondivision: squarefreeH and gcd(H,U6)=1 imply dim2 nonempty V(H)∩D(U6*D12)',
  combined_test_REQUIRED='D12_remainder*T3_remainder must be nonzero in SAME ring; separate nonzero tests invalid in zero-divisor quotient',
  combined_positive_scope='after independent peer some SAME Hcomponent has nonempty dim2 open with actual full57; no wholeH/global coverage',
  zero_combined_inconclusive=True,ring='Fp32003[a]/[17430,11013,30453,27809,1]',not_field_assumed=True,
  scalar_products_upper=productupper,hard_scalar_products=100000,hard_field_bits=30,hard_owned=12000,
  source_owned_runtime=6000+162*4+powerN*5,hard_runtime_RAW=8388608,hard_prepare_RAW=16777216,
  soft_seconds=45,hard_seconds=60,address_space=2000000000,D12_independent_value_peer_REQUIRED_before_adoption=True)
 (out/'source-and-whole-guard-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 runtime[str(native.relative_to(ROOT))]=sha(native)
 rep=dict(id=out.name,all_passed=True,frozen=True,prepared_only=True,no_actual_j3_or_combined_values=True,
  worker_sha256=sha(worker),runtime_input_sha256=runtime,input_sha256=pins,native_input=str(native.relative_to(ROOT)),
  literal_atom_count=len(atoms),instructions=instructions,active_column_and_cached_nodes=sorted(active),output_root=plan['output_root'],minimum_a_exponent=amin,maximum_a_exponent=amax,
  monic_H_axis=[17430,11013,30453,27809,1],D12_report='reports/S61_TWO_TAIL_H_D12_AXIS_QUOTIENT_ENTRY_20261105_v63/report.json',
  independent_D12_value_peer_pending=True,evidence_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()})
 p=out/'report.json';p.write_text(json.dumps(rep,indent=2)+'\n');raw=sum(x.stat().st_size for x in out.iterdir())+sum((ROOT/s).stat().st_size for s in pins);assert raw<16777216
 print(json.dumps(dict(SHA256=sha(p),RAW_prepare=raw,scalar_upper=productupper,owned_runtime=proof['source_owned_runtime'],atoms=len(atoms),a_box=[amin,amax],ROM_bits=rombits),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
