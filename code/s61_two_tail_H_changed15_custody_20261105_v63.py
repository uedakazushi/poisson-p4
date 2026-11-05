#!/usr/bin/env python3
"""SOURCE/custody for changed15 pair4; copy SAVED axis coefficients, no ring values."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,importlib.util,json,math,resource,signal,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main(out):
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 pins={};runtime={}
 def pin(s,h=None,run=False):
  p=ROOT/s;v=sha(p)
  if h:assert v==h
  pins[s]=v
  if run:runtime[s]=v
  return p
 def rep(st,h=None,run=False):return json.loads(pin('reports/'+st+'/report.json',h,run).read_text())
 direct=rep('S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63_REVIEW','f946ca6312c0b5167aa5cc257434d9338727a0ff3a6142496bc2ffa55d3890c8',True)
 name='reports/S61_NONSQUARE_PREPOWER_DIRECT_M12_FULL15_SOURCE_20261105_v63_REVIEW/independent-source-normalization-and-A-D15-minusT-D12-proof.json'
 dp=pin(name,direct['evidence_sha256'][Path(name).name],True);d=json.loads(dp.read_text())
 assert d['all450_row_inverse_products'] and d['all12_pair_H3_entries_literal_zero'] and d['pair_rank1_on_exact9_pure_columns']==[22,11,13,25,15,19,14,18,17]
 rep('S61_NONSQUARE_PREPOWER_GLOBAL_FULL57_SOURCE_20261105_v63_REVIEW','6640c4f2be0e409a8f4854fd458fbbd700b4e28e7d6a360c153b4f2f67bcd25b',True)
 rep('S61_TWO_TAIL_H_D12_ROOT_NEWTON_REVIEW_20261105_v63','041d9234ccf704286eacd0da6e0f6161d8b268b27dfaa678e2080657f5f32193',True)
 rep('S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63','49990187b3b0b8edaf8939e4d13ac297f4f2ecba7f8b9f5326639089fa043ce2',True)
 k=rep('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63',run=True)
 assert k['all_passed']
 rep('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63_REVIEW','fb6343721086e85f9ee44c7cd361d585249bcbbb35fc9a7bc35c3337a5576726',True)
 t=rep('S61_NONSQUARE_PREPOWER_TAIL1_AXIS_REVIEW_20261105_v63','fee8b8a767cba862a47ceecf0f4e3c3bdf79772c9420dea42ca3912839f2332b',True)
 assert t['four_column_axes_and_complete_T1_axis_match']
 data=pin('reports/S61_NONSQUARE_PREPOWER_TAIL1_CACHED_KAPPA_ENTRY_20261105_v63/integer-column4-and-result.txt','a40a01d3bac906dd10023d086090029756f9918a2219a89713faa51c717070ec',True)
 tail1=[];native_counts=[]
 with data.open() as f:
  assert f.readline().strip()=='J1KV001'
  for expect in (66,94,122,150,155):
   node,D=map(int,f.readline().split());assert node==expect and D>0;n=int(f.readline());axis={}
   for _ in range(n):
    key,v=map(int,f.readline().split())
    if node!=155 and key%1048576==0:
     a=key//1048576;assert a not in axis;axis[a]=Q(v,D)
   native_counts.append(dict(node=node,full_terms=n,axis_terms=len(axis)))
   if node!=155:tail1.append(axis)
  assert not f.read().strip()
 helper=pin('code/s61_nonsquare_prePOWER_tail1_worker_custody_20261105_v63.py')
 spec=importlib.util.spec_from_file_location('lexical_reader_only',helper);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);Spans=mod.Spans
 hcname='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json'
 hc=Spans(pin(hcname,k['input_sha256'][hcname],True));hs=hc.obj('raw_compact14x5');rs=hc.obj('known_row_units_rho14')
 bn='reports/S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63/native.json'
 bot=Spans(pin(bn,k['input_sha256'][bn],True));bs=bot.obj('bottom4x6')
 P=[[None]*4 for _ in range(4)];refs=[]
 def axis(poly,scale=Q(1)):
  z={}
  for (a,e,d),v in poly:
   if e==d==0:
    assert a not in z;qv=Q(v)
    assert abs(qv.numerator).bit_length()+abs(scale.numerator).bit_length()<=2048
    assert qv.denominator.bit_length()+scale.denominator.bit_length()<=2048
    z[a]=qv*scale
  return z
 for i,r in enumerate([1,3,4,7]):
  rho=Q(hc.value(hc.obj('constant',hc.index(rs,r))));P[i][1]=axis(hc.value(hc.index(hs,r,0)),rho)
  P[i][0]=axis(bot.value(bot.index(bs,i,0)))
  P[i][2]=tail1[i]
  n=[66,94,122,150][i];p=pin('reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63/root'+str(n)+'.json',k['evidence_sha256']['root'+str(n)+'.json'],True)
  P[i][3]=axis(json.loads(p.read_text())['polynomial'])
  refs.append(dict(row=i,current15_row=[2,4,5,8][i],anchor=[bn,'bottom4x6',i,0],H27=[hcname,'raw_compact14x5',r,0],H27_rho_constant=str(rho),tail1=[str(data.relative_to(ROOT)),n],tail2=[str(p.relative_to(ROOT)),'polynomial'],no_extra_a_or_L_factor=True))
 hc.close();bot.close();inputs=[];amin=128;amax=-128;rom=denbits=prepbits=0;terms=0
 for i in range(4):
  for j in range(4):
   q=P[i][j];D=math.lcm(*(v.denominator for v in q.values())) if q else 1;assert D%32003
   rows=[]
   for a,v in sorted(q.items()):
    bound=abs(v.numerator).bit_length()+D.bit_length()+2;assert bound<=2048
    z=v.numerator*(D//v.denominator);assert Q(z,D)==v
    prepbits=max(prepbits,bound);rom=max(rom,abs(z).bit_length());amin=min(amin,a);amax=max(amax,a);rows.append([a,str(z)])
   denbits=max(denbits,D.bit_length());terms+=len(rows);inputs.append(dict(row=i,column=j,D=str(D),coefficients=rows))
 assert rom<=384 and denbits<=384
 perm=list(range(11))+[14,11,12,13];inverse_perm={old:new for new,old in enumerate(perm)}
 assert sum(perm[i]>perm[j] for i in range(15) for j in range(i+1,15))==3
 rowop=[d['row_operator_QQ_Laurent_a'][i] for i in perm]
 rowinv=[{str(inverse_perm[int(j)]):p for j,p in row.items()} for row in d['row_operator_inverse_QQ_Laurent_a']]
 order=[0,1,3,6,2,11,12,9,7,13,10,14,4,5,8]
 cols=[22,11,13,25,15,19,14,18,28,30,29,17,27,1,2]
 mathematical=dict(genuine_current15_rows=order,target_original32_columns=cols,pair_rows=[2,4,5,8],pair_columns=['anchor22','H27/L','B*tail1','B*tail2'],
  old_M12_rows=order[:12],old_M12_columns=cols[:12],all15_row_operator=rowop,all15_row_operator_inverse=rowinv,
  relative_row_parity_to_old15=-1,row_determinant='minus old accepted row determinant; QQstar*a^-1',
  target_normalization={'28':'1/L','30':'1/L','29':'1/L','27':'1/L','1':'B','2':'B','others':'1','L':'a*F*C^2'},column_unit_product='B^2/L^4',
  clear_lower_rows_positions=[12,13,14],extra_columns_positions=[12,13,14],Laplace_index_sum=78,Laplace_sign=1,
  universal_block_identity='det(A*V-c*u)=A^2*detPair4 (accepted old universal identity, no new determinant values)',
  domain_identity='A*D15_new=+detPair4*D12',anchor_nonzero_source_terms=156,anchor_NOT_inverted=True,
  GLOBAL_full57_bridge='accepted complete U6 rank(old25)=10+rank(current15) and original63/raw91 product41/Cdelta lift; all tails retained',
  all_inversions_only='QQstar,a,d,F,B,C,J',no_private_kappa_Omega_assumption=True,
  combined_required='D12*detPair4 mod SAME Haxis, separate nonzero factors insufficient in zero-divisor ring',
  positive_after_independent_peer='some Hcomponent has nonempty dimension2 parameter-open with actual full57',
  zero_inconclusive_no_fallback=True,source_entry_references=refs,tail1_source_block_counts=native_counts)
 (out/'genuine15-source-identity-and-frame.json').write_text(json.dumps(mathematical,indent=2)+'\n')
 powerN=max(amax,1)-min(amin,-1)+1;prod=4*terms+16*35+28*(powerN+32+3)+100;assert prod<10000
 own=2048+64+64+5*powerN;assert own<4096
 budget=dict(source_only=True,no_ring_or_determinant_values=True,axis_input_terms=terms,a_box=[amin,amax],existing_ROM_integer_bits=rom,existing_D_bits=denbits,
  preparation_scalar_temp_bits=prepbits,preparation_owned_upper=50000,preparation_owned_proof='one original1257-term root JSON plus16short axis polynomials,accepted24term universal proof/rowGL and lexical spans/control <50000',field_temporary_bits=30,scalar_product_upper=prod,hard_products=10000,owned_prior=own,hard_owned=4096,
  owned_proof='2048 input/source/control/failure reserve +64DP+64matrix+5 per Laurent power',
  runtime_RAW_cap=4194304,prepare_RAW_cap=8388608,soft_seconds=45,hard_seconds=60,address_space=2000000000,no_value_worker_run=True)
 (out/'source-and-whole-budget.json').write_text(json.dumps(budget,indent=2)+'\n')
 worker=ROOT/'code/s61_two_tail_H_changed15_axis_20261105_v63.py';pins[str(worker.relative_to(ROOT))]=sha(worker)
 note=ROOT/'research/notes/S61_TWO_TAIL_H_CHANGED15_SOURCE_20261105_v63.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 h=subprocess.run([sys.executable,str(worker),'--help'],capture_output=True,text=True,timeout=10);assert h.returncode==0
 (out/'worker-help.txt').write_text(h.stdout);(out/'worker-help.stderr').write_text(h.stderr)
 repout=dict(id=out.name,all_passed=True,frozen=True,prepared_only=True,no_determinant_or_ring_values_executed=True,
  worker_sha256=sha(worker),runtime_input_sha256=runtime,input_sha256=pins,axis_input_polynomials=inputs,minimum_a_exponent=amin,maximum_a_exponent=amax,
  D12_peer='reports/S61_TWO_TAIL_H_D12_ROOT_NEWTON_REVIEW_20261105_v63/report.json',
  evidence_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()})
 p=out/'report.json';p.write_text(json.dumps(repout,indent=2)+'\n')
 raw=sum(x.stat().st_size for x in out.iterdir())+sum((ROOT/s).stat().st_size for s in pins);assert raw<8388608
 rr=sum((ROOT/s).stat().st_size for s in runtime)+p.stat().st_size+worker.stat().st_size+131072;assert rr<4194304
 print(json.dumps(dict(SHA256=sha(p),source_RAW=raw,runtime_RAW_prior=rr,products_upper=prod,owned=own,axis_terms=terms,a_box=[amin,amax],ROM_bits=rom,D_bits=denbits),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output.resolve())
