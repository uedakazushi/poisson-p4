#!/usr/bin/env python3
"""SOURCE ONLY: all72 remaining A/B/C exact linear entries and known-zero filters."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_EXACT_ABC_SOURCE_PREFLIGHT_20261105_v63';W=(1,2244,68);INF=10**12
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 profbase='reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63';rank=src(profbase+'/report.json');assert rank['all_passed']and pins[profbase+'/report.json']=='adc8ab36ddb71baae8aa885841287e865b63b35ff8a5701b32a015535fac3450';known=src(profbase+'/all37-and32-factor-identities.json')['all32_integer_linear_identities'];assert len(known)==32
 base='reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63';rt=src(base+'/retropicalized-source.json');row=src(base+'/actual-row11.json');R=[[Q(v)for v in rr]for rr in rt['row_transform']];Qmat=[[Q(v)for v in rr]for rr in rt['column_transform']];DR=math.lcm(*(v.denominator for rr in R for v in rr));DQ=math.lcm(*(v.denominator for rr in Qmat for v in rr));DD=DR*DQ;RZ=[[int(v*DR)for v in rr]for rr in R];QZ=[[int(v*DQ)for v in rr]for rr in Qmat];Ds=list(map(int,rt['constant_clearing_rowLCMs']));rp=rt['row_potentials'];cp=rt['column_potentials']
 hp=ROOT/'code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py';pins[str(hp.relative_to(ROOT))]=sha(hp);sp=importlib.util.spec_from_file_location('sourcegate',hp);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h);bind=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02';M,inventory,kept,skip=h.selected_inputs(state,bind['actual_pure_matrix_roots'])
 for p in row['columns']:M[1][p['column']]=dict(state_id=-1001-p['column'],c=Q(p['QQ_common_factor']),co=[(tuple(m),int(v))for m,v in p['primitive_integer_polynomial']])
 known_map={(q['global_row'],q['global_column'],q['gap']):Q(q['QQ_value'])for q in known};entries=[];maxsupport=totalsupport=ops=raw=maxbits=0;updated_min=[[INF]*11 for i in range(11)];knownzeros=0
 for u in range(11):
  for v in range(11):
   if u>=4 and v>=4:continue
   kind='A'if u<4 and v<4 else'B'if u<4 else'C';br=u if u<4 else u-4;bc=v if v<4 else v-4;keys=set();L=0;eo=0;contributors=[]
   for i in range(11):
    if not RZ[u][i]:continue
    for j in range(11):
     if not QZ[j][v]:continue
     p=M[i][j]
     if not p['co']:continue
     scale=Q(RZ[u][i]*QZ[j][v]*Ds[i])*p['c'];assert scale.denominator==1;gaps=[rp[i]+cp[j]-sum(a*b for a,b in zip(mon,W))for mon,z in p['co']];assert all(g>=0 for g in gaps)and len(gaps)==len(set(gaps));keys.update(gaps);L+=abs(scale.numerator)*sum(abs(z)for mon,z in p['co']);eo+=len(p['co']);contributors.append(dict(source_row=i,source_column=j,state_polynomial_id=p['state_id'],integer_scale=str(scale.numerator),source_terms=len(p['co']),source_QQfactor=str(p['c']),gap_shift=rp[i]+cp[j]))
   expected0=Q(u==v and u<4);checks=[dict(gap=0,exact_QQ_value=str(expected0),integer_numerator=str(int(expected0*DD)))];zeros={0}if not expected0 else set()
   for (i,j,g),value in known_map.items():
    if(i,j)==(u,v):assert(value*DD).denominator==1;checks.append(dict(gap=g,exact_QQ_value=str(value),integer_numerator=str(int(value*DD))));zeros.update([g]if not value else[])
   knownzeros+=len(zeros-{0});updated_min[u][v]=min((g for g in keys if g>0 and g not in zeros),default=INF);bits=L.bit_length();size=len(keys)*(35+(bits*30103+99999)//100000)+2000;maxsupport=max(maxsupport,len(keys));totalsupport+=len(keys);ops+=eo;raw+=size;maxbits=max(maxbits,bits)
   entries.append(dict(entry_index=len(entries),block=kind,block_row=br,block_column=bc,global_row=u,global_column=v,source_contributors=contributors,source_scalar_additions=eo,union_exponent_support_upper=len(keys),max_exponent=max(keys,default=0),integer_coefficient_L1_bits=bits,exact_integer_L1_upper=str(L),raw_serialization_upper=size,known_coefficients_to_recheck=checks,known_zero_filtered_positive_gap_lower=updated_min[u][v],NO_new_coefficient_sums=True))
 assert len(entries)==72 and {b:sum(e['block']==b for e in entries)for b in['A','B','C']}=={'A':16,'B':28,'C':28};live=kept+row['actual_total_terms']+2*maxsupport+15000;newraw=raw+500000;external=sum((ROOT/p).stat().st_size for p in pins);assert live<150000 and external+newraw<32*1024*1024
 dump(out/'all72-linear-operator.json',dict(fixed_weight=W,formal_source='Bwhole(t)=diag(rowLCM)*diag(t^rowpot)*actualM(t^-w)*diag(t^colpot); transformed=R*Bwhole*Q',blocks={'A':[4,4],'B':[4,7],'C':[7,4]},QQ_common_denominator=str(DD),row_LCMs=list(map(str,Ds)),row_potentials=rp,column_potentials=cp,R_integer=RZ,R_denominator=str(DR),Q_integer=QZ,Q_denominator=str(DQ),entries=entries,no_polynomial_convolution_or_inverse=True,no_source572_replay=True))
 dump(out/'knownzero-filter-source.json',dict(positive_C_B_coefficients_already_known=32,known_zero_positive_slots=knownzeros,updated_boolean_min_positive_gaps=updated_min,
  these_are_STILL_lower_bounds_until_full_exactABC=True,next_goal='Exact72coefficients removes ALL Boolean cancellation ambiguity in A/B/C supports. Then updated Ainv-minplus routes and typed residual6 SOURCE valuation; no nextseriescoefficient in this budget.'))
 dump(out/'resource-gate.json',dict(source_input_terms=kept,normalized_row_terms=row['actual_total_terms'],max_one_entry_union_support=maxsupport,all72_union_support=totalsupport,total_scalar_additions=ops,integer_coefficient_bit_upper=maxbits,whole_live_upper=live,
  all72_coefficient_RAW_upper=raw,new_artifacts_plus500k_upper=newraw,external_read_RAW=external,total_read_plus_output_upper=external+newraw,hard_seconds=60,hard_cpu_guard=45,hard_address_bytes=2000000000,hard_scalar_terms=150000,hard_total_RAW=33554432,
  entry_stream_and_release=True,integer_and_reverseQQ_wholeidentity_replay=True,allA0I4_B0C0zeros_and32existing_coefficients_rechecked=True,NO_coefficient_job_executed=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);total=sum(p.stat().st_size for p in out.iterdir());assert total<500000
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,blocks_A16_B28_C28_TOTAL72=True,NO_coefficient_sums_or_series_inverse_det=True,max_single_entry_union_support=maxsupport,total_union_support=totalsupport,total_scalar_additions=ops,
  max_integer_bits=maxbits,whole_live_upper=live,new_RAW_upper=newraw,external_read_RAW=external,total_read_plus_new_upper=external+newraw,positive_known_zero_slots=knownzeros,no_newweight_K_point_prime_or_source572=True,
  no_det11_zero_nonzero_or_global57_claim=True,total_RAW_before_report=total,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
