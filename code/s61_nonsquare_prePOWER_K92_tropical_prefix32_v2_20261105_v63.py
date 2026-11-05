#!/usr/bin/env python3
"""ONE fixed QQ coefficient job at w=(1,2244,68), K32; no fallback."""
import argparse,gc,hashlib,importlib.util,json,math,resource,shlex,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_TROPICAL_PREFIX32_V2_20261105_v63'
W=(1,2244,68);K=32
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 gate='reports/S61_NONSQUARE_PREPOWER_K92_TROPICAL_SOURCE_PREFLIGHT_V2_20261105_v63';gr=src(gate+'/report.json');assert gr['all_passed']and pins[gate+'/report.json']=='2904fcec4f1563e62c8484ab7b26e44bb33f7d5a272d21f2e636266a16d9b1f6'
 profiles=src(gate+'/six-support-orders.json')['profiles'];profile=next(x for x in profiles if tuple(x['weight'])==W);assert profile['tropical_degree']==232363;rowpot=profile['row_potentials'];colpot=profile['column_potentials'];cost=next(x for x in profile['prefix_costs']if x['K']==K);assert cost['selected_existing_source_terms']==78 and cost['QQinteger_coefficient_bit_upper']==592
 helper=ROOT/'code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py';pins[str(helper.relative_to(ROOT))]=sha(helper);sp=importlib.util.spec_from_file_location('sourcegate',helper);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
 binding=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02'
 M,inventory,kept,skipped=h.selected_inputs(state,binding['actual_pure_matrix_roots']);assert kept==21690
 rows=[];cleared=[[{}for j in range(11)]for i in range(11)];records=[];rowL1=[]
 for i,row in enumerate(M):
  selected=[];D=1
  for j,p in enumerate(row):
   for idx,(mon,v)in enumerate(p['co']):
    deficit=rowpot[i]+colpot[j]-sum(x*y for x,y in zip(mon,W));assert deficit>=0
    if deficit<=K:selected.append((j,p,idx,mon,v,deficit));D=math.lcm(D,p['c'].denominator)
  L=0
  for j,p,idx,mon,v,gap in selected:
   value=v*p['c'].numerator*(D//p['c'].denominator);assert value;cleared[i][j][gap]=cleared[i][j].get(gap,0)+value;L+=abs(value);records.append(dict(row=i,column=j,state_polynomial_id=p['state_id'],source_index=idx,original_monomial=mon,original_integer=str(v),common_QQfactor=str(p['c']),deficit=gap,row_LCM=str(D),cleared_integer=str(value)))
  rows.append(D);rowL1.append(L);assert L.bit_length()==cost['rows'][i]['row_integer_L1_bits']
 assert len(records)==78 and all(v for row in cleared for entry in row for v in entry.values());bound=sum(x.bit_length()for x in rowL1);assert bound==592;Dprod=math.prod(rows)
 ip=out/'native-input.txt'
 with ip.open('w')as f:
  f.write('11 32 592\n');f.write(' '.join(map(str,rows))+'\n')
  for row in cleared:
   for entry in row:
    f.write(str(len(entry))+'\n')
    for gap,v in sorted(entry.items()):f.write(f'{gap} {v}\n')
 dump(out/'selected78-source.json',dict(weight=W,K=K,tropical_degree=232363,row_potentials=rowpot,column_potentials=colpot,matching=profile['matching'],row_LCMs=list(map(str,rows)),row_L1_bits=[x.bit_length()for x in rowL1],row_LCM_product=str(Dprod),source_terms=records,
  identity='det(cleared shifted M)(t)=(product rowLCM)*t^232363*det(M(a=t^-1,e=t^-2244,d=t^-68)); all deficits nonnegative; discard deficits>32 only',sign='ordinary determinant roworder/columnorder unchanged; subsettransition sign=(-1)^popcount(mask>>(j+1))',NO_source_point=True,NO_new_prime=True,all_source572_NOT_replayed=True))
 dump(out/'resource.json',dict(K=K,source_parse_terms=21690,selected78_terms=78,native_DP_initialized_slots_upper=924*33,native_total_plus_scratch_upper=924*33+78+15000,lastuse_source_released_BEFORE_native=True,integer_bit_bound=592,one_weight_only=W,NO_Fp_fallback=True))
 cpp=Path(__file__).with_suffix('.cpp');note=ROOT/'research/notes'/f'{ID}.md';exe=out/'native-worker'
 for p in[Path(__file__),cpp,note]:pins[str(p.relative_to(ROOT))]=sha(p)
 del profiles,profile,M,inventory,binding,cleared,records,gr,cost;gc.collect();dump(out/'phase.json',dict(phase='source78 frozen/source21690released beforecompile'))
 cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+shlex.quote(str(cpp))+' -o '+shlex.quote(str(exe))+' $(pkg-config --cflags --libs flint)';r=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,55-(time.monotonic()-start)));(out/'compile.log').write_text(r.stdout)
 if r.returncode:
  dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,arithmetic_STARTED=False,compile_failure_preserved=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()}));signal.alarm(0);return
 dump(out/'phase.json',dict(phase='ONEfixedK32 integerQQ coefficient job begins'))
 cmd='source /workspace/cloud-setup/activate.sh && exec '+shlex.quote(str(exe))+' '+shlex.quote(str(ip))+' '+shlex.quote(str(out));r=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,58-(time.monotonic()-start)));(out/'native.log').write_text(r.stdout)
 if r.returncode:
  dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,exact_prefix_complete=False,partial_state_preserved=True,NO_fallback=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()}));signal.alarm(0);return
 nr=read(out/'native.json');assert nr['all_passed']and int(nr['row_denominator_product'])==Dprod;integers=list(map(int,nr['cleared_integer_coefficients']));assert len(integers)==33;coefficients=[Q(c,Dprod)for c in integers];nz=next((i for i,c in enumerate(integers)if c),None);mon=None
 if nz is not None:
  target=232363-nz;n=target-5;a=5+n%68;n//=68;d=n%33;e=n//33;mon=(a,e,d);assert 5<=a<=72 and 0<=e<=107 and 0<=d<=32 and sum(x*y for x,y in zip(mon,W))==target
 dump(out/'exact33-coefficients.json',dict(weight=W,K=K,tropical_degree=232363,row_LCM_product=str(Dprod),cleared_integer_coefficients=list(map(str,integers)),QQ_coefficients=list(map(str,coefficients)),first_nonzero_deficit=nz,unique_original_monomial=mon,
  original_det11_nonzero_primary=nz is not None,allzero_is_ONLY_prefix_zero=True,NO_second_weight_widerK_or_Fp=True))
 dump(out/'phase.json',dict(phase='FINAL ALL33 exactQQ coefficients computed; NO fallback'))
 total=sum(p.stat().st_size for p in out.iterdir());assert total<32*1024*1024 and nr['peak_plus15000scratch']<=150000
 report=dict(id=ID,all_passed=True,frozen=True,ONE_fixedweight_QQ_job=True,weight=W,K=K,all33_exact_coefficients_saved=True,all33_zero=nz is None,first_nonzero_deficit=nz,unique_original_monomial=mon,actual_det11_nonzero_primary=nz is not None,independent_peer_pending=True,
  source21690_released_BEFORE_DP=True,all572source_NOT_replayed=True,NO_numeric_parameter_point=True,NO_Fp_or_otherweight_fallback=True,no_global57_claim=True,coefficient_updates=nr['coefficient_updates'],peak_integer_slots_plus_scratch=nr['peak_plus15000scratch'],total_RAW_before_report=total,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
