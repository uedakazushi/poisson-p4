#!/usr/bin/env python3
"""SOURCE ONLY: stream exact formal D7x7 linear polynomials; no coefficient sums."""
import argparse,hashlib,importlib.util,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_EXACT_D_SOURCE_PREFLIGHT_20261105_v63';W=(1,2244,68)
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 f7=src('reports/S61_NONSQUARE_PREPOWER_K92_FINITE_LEADING7_20261105_v63/report.json');assert f7['all_passed']and f7['exact_leading7_rank']==0;f7entries=src('reports/S61_NONSQUARE_PREPOWER_K92_FINITE_LEADING7_20261105_v63/all41-entry-identities.json')
 base='reports/S61_NONSQUARE_PREPOWER_K92_ONE_ROW_PRECONDITION_20261105_v63';r=src(base+'/report.json');assert r['all_passed'];x=src(base+'/retropicalized-source.json');ro=src(base+'/actual-row11.json');R=[[Q(v)for v in rr]for rr in x['row_transform']];C=[[Q(v)for v in rr]for rr in x['column_transform']];assert x['new_leading_constant_rank']==4
 hp=ROOT/'code/s61_nonsquare_prePOWER_K92_tropical_source_preflight_v2_20261105_v63.py';pins[str(hp.relative_to(ROOT))]=sha(hp);sp=importlib.util.spec_from_file_location('sourcegate',hp);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
 bind=src('reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/source-binding.json');state=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_DET11_GUARDED_ATTEMPT_20261104_v63/state.bin';pins[str(state.relative_to(ROOT))]=sha(state);assert pins[str(state.relative_to(ROOT))]=='cd5fb2431c2f65649153d82890d52b52ce026f9d0e478bc30592687f17b20f02';M,inventory,kept,skip=h.selected_inputs(state,bind['actual_pure_matrix_roots'])
 for p in ro['columns']:M[1][p['column']]=dict(state_id=-1001-p['column'],c=Q(p['QQ_common_factor']),co=[(tuple(m),int(v))for m,v in p['primitive_integer_polynomial']])
 DR=math.lcm(*(v.denominator for rr in R for v in rr));DQ=math.lcm(*(v.denominator for rr in C for v in rr));DD=DR*DQ;RZ=[[int(v*DR)for v in rr]for rr in R];QZ=[[int(v*DQ)for v in rr]for rr in C];Ds=list(map(int,x['constant_clearing_rowLCMs']));rp=x['row_potentials'];cp=x['column_potentials'];profiles=[];max_support=total_support=ops=raw=0;maxbits=0
 for u in range(4,11):
  for v in range(4,11):
   keys=set();L=0;entryops=0;contributors=[]
   for i in range(11):
    if not RZ[u][i]:continue
    for j in range(11):
     if not QZ[j][v]:continue
     p=M[i][j];scale=Q(RZ[u][i]*QZ[j][v]*Ds[i])*p['c'];assert scale.denominator==1
     if not p['co']:continue
     gaps=[rp[i]+cp[j]-sum(a*b for a,b in zip(mon,W))for mon,z in p['co']];assert all(g>=0 for g in gaps)and len(gaps)==len(set(gaps));keys.update(gaps);l1=sum(abs(z)for mon,z in p['co']);L+=abs(scale.numerator)*l1;entryops+=len(p['co']);contributors.append(dict(source_row=i,source_column=j,state_polynomial_id=p['state_id'],integer_scale=str(scale.numerator),source_terms=len(p['co']),source_QQfactor=str(p['c']),gap_shift=rp[i]+cp[j]))
   bits=L.bit_length();size=len(keys)*(35+(bits*30103+99999)//100000)+2000;knownzeros=[0]+[q['gap']for q in f7entries['entries']if q['residual_row']==u-4 and q['residual_column']==v-4];max_support=max(max_support,len(keys));total_support+=len(keys);ops+=entryops;raw+=size;maxbits=max(maxbits,bits)
   profiles.append(dict(entry_index=len(profiles),residual_row=u-4,residual_column=v-4,global_row=u,global_column=v,integer_numerator_common_denominator=str(DD),source_contributors=contributors,
    source_scalar_additions=entryops,union_exponent_support_upper=len(keys),max_exponent=max(keys,default=0),integer_coefficient_L1_bits=bits,exact_integer_L1_upper=str(L),raw_serialization_upper=size,certified_zero_coefficients_to_recheck=sorted(set(knownzeros)),NO_new_coefficient_sums=True))
   dump(out/'phase.json',dict(phase='SOURCE ONLY exactD linear support unions',entries_profiled=len(profiles)))
 assert len(profiles)==49
 live=kept+ro['actual_total_terms']+2*max_support+15000;newraw=raw+500000;external=sum((ROOT/p).stat().st_size for p in pins);assert live<150000 and newraw+external<32*1024*1024
 dump(out/'all49-linear-operator.json',dict(formal_source='B(t)=diag(rowLCM)*diag(t^rowpot)*actualM(t^-w)*diag(t^colpot)',D_definition='lower7x7 block of R*B(t)*Q; exactfull sparseunivariate polynomials, NOT determinant/Schurinverse',
  fixed_weight=W,QQ_common_denominator=str(DD),row_LCMs=list(map(str,Ds)),row_potentials=rp,column_potentials=cp,R_integer=RZ,R_denominator=str(DR),Q_integer=QZ,Q_denominator=str(DQ),entries=profiles,
  no_polynomial_convolution=True,NO_unknown_parameter_inverse=True,NO_numeric_point_or_prime=True,no_source572_replay=True))
 dump(out/'resource-gate.json',dict(source_input_terms=kept,normalized_one_row_terms=ro['actual_total_terms'],max_ONEentry_union_support=max_support,all49_union_support=total_support,total_scalar_additions=ops,integer_coefficient_bit_upper=maxbits,
  whole_live_upper=live,max_entry_RAW_upper=max(p['raw_serialization_upper']for p in profiles),all49_coefficients_RAW_upper=raw,new_artifacts_plus500k_upper=newraw,actual_external_read_RAW=external,total_read_plus_output_bound=external+newraw,
  hard_seconds=60,hard_native_clock=45,hard_address_bytes=2000000000,hard_owned_scalar_terms=150000,hard_total_RAW_bytes=33554432,
  lastuse='one entry dictionary at a time; stream final sparse coefficient list directly to JSONL then release BEFORE nextentry; never retain all49',
  stop='before any new coefficient/bit/filegrowth guard, preserve currententry exact dictionary+sourcecursor and previouscompleted49entry stream; no auto continuation',
  measured_gate_passes=True,NO_coefficient_job_executed=True))
 dump(out/'completion-and-next-source.json',dict(next_allowed_if_ROOT_GO='ONE all49 exactlinear coefficients, allsource integer scale identities, normalized QQcontents, allD0 and41previousleadingzero checks',
  completion_criterion='49 complete entry records with exact source sum identities and hashes, no determinant or inverse series',next_source_only='Use EXACT D entry valuations, retaining original A/B/C support and ALL C*AinvB correction routes; re-evaluate sparsegap Schur minimum then derive finite next coefficient request.',
  no_actual_solution_or_global57_claim=True,formal_D_zero_does_NOT_mean_physical_parameterpolynomial_entry_zero_without_separate_injectivity=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);total=sum(p.stat().st_size for p in out.iterdir());assert total<500000
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,all49_exactD_linear_operator_bound=True,NO_coefficient_sum_or_det_or_inverse_series=True,max_single_entry_union_support=max_support,all49_union_support=total_support,total_integer_scalar_additions=ops,
  integer_bit_upper=maxbits,whole_live_upper=live,new_RAW_upper=newraw,external_read_RAW=external,read_plus_new_upper=external+newraw,no_otherweight_K_point_prime_or_source572=True,no_det11_zero_nonzero_or_global57_claim=True,
  total_RAW_before_report=total,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
